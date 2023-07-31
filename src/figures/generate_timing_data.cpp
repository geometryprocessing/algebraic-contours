// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "apply_transformation.h"
#include "compute_boundaries.h"
#include "contour_network.h"
#include "generate_transformation.h"
#include "twelve_split_spline.h"
#include <igl/Timer.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <sys/resource.h>
#include <CLI/CLI.hpp>

int main(int argc, char *argv[]) {
  // Set stack size
  const rlim_t kStackSize = 512 * 1024 * 1024;   // min stack size = 512 MB
  struct rlimit rl;
  int result;
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
  {
      if (rl.rlim_cur < kStackSize)
      {
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0)
          {
              fprintf(stderr, "setrlimit returned result = %d\n", result);
          }
      }
  }

  // Get command line arguments
  CLI::App app{"Generate perspective distortion figure images for a given mesh and camera."};
  std::string input_filename = "";
  std::string output_dir = "./";
  int num_tests = 0;
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  app.add_option("--num_tests", num_tests, "Number of tests to run")
    ->check(CLI::NonNegativeNumber);
  CLI11_PARSE(app, argc, argv);

  // Set logging level
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd initial_V, V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, initial_V, uv, N, F, FT, FN);

  // Start timer for quadratic spline generation
  igl::Timer timer;
  timer.start();

  // Generate quadratic spline
  OptimizationParameters optimization_params;
  std::vector<std::vector<int>> face_to_patch_indices;
  std::vector<int> patch_to_face_indices;
  Eigen::SparseMatrix<double> fit_matrix;
  Eigen::SparseMatrix<double> energy_hessian;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>
      energy_hessian_inverse;
  AffineManifold affine_manifold(F, uv, FT);
  TwelveSplitSplineSurface spline_surface(
      initial_V, affine_manifold,
      optimization_params, face_to_patch_indices, patch_to_face_indices,
      fit_matrix, energy_hessian, energy_hessian_inverse);

  // Get spline construction time
  double spline_surface_time = timer.getElapsedTime();
  timer.stop();

  // Time and get the boundary edges
  timer.start();
  std::vector<std::pair<int, int>> patch_boundary_edges(0);
  compute_twelve_split_spline_patch_boundary_edges(F, face_to_patch_indices, patch_boundary_edges);
  double compute_patch_boundary_time = timer.getElapsedTime();
  timer.stop();

  // Start the timer for the contour generation
  timer.start();

  // Build the contours
  IntersectionParameters intersect_params;
  InvisibilityParameters invisibility_params;
  invisibility_params.write_contour_soup = false;
  invisibility_params.invisibility_method = InvisibilityMethod::propagation;
  invisibility_params.check_propagation = false;
  ContourNetwork contour_network(spline_surface,
                                 intersect_params, invisibility_params,
                                 patch_boundary_edges);

  // Get contour network construction time
  double initial_contour_network_time = timer.getElapsedTime();
  timer.stop();

  // Write view independent timing data
  std::ofstream out_view_independent(join_path(output_dir, "view_independent.csv"),
                                     std::ios::app);
  out_view_independent << input_filename << "," << F.rows() << ","
                       << spline_surface_time << ","
                       << compute_patch_boundary_time << ",\n";
                    
  // Write view dependent timing data
  std::ofstream out_per_view(join_path(output_dir, "per_view.csv"), std::ios::app);
  out_per_view << input_filename << ","
               << "1 0 0 0 1 0 0 0 1"
               << "," << initial_contour_network_time << "," << 0.0 << ","
               << contour_network.compute_contour_time << ","
               << contour_network.compute_cusp_time << ","
               << contour_network.compute_intersection_time << ","
               << contour_network.compute_visibility_time << ","
               << contour_network.compute_projected_time << ","
               << contour_network.segment_number << ","
               << contour_network.interior_cusp_number << ","
               << contour_network.boundary_cusp_number << ","
               << contour_network.intersection_call << ","
               << contour_network.ray_intersection_call << ","
               << spline_surface.num_patches() << ",\n";

  // Generate a random number generator for angles
  std::mt19937 e2(0);
  std::uniform_real_distribution<double> angle_distribution(0.0, 2.0 * M_PI);

  // Run tests
  for (int i = 0; i < num_tests; ++i) {
    // Generate random rotation matrix
    double theta_x = angle_distribution(e2);
    double theta_y = angle_distribution(e2);
    double theta_z = angle_distribution(e2);
    MatrixXr x_rotation_matrix = x_axis_rotation_projective_matrix(theta_x);
    MatrixXr y_rotation_matrix = y_axis_rotation_projective_matrix(theta_y);
    MatrixXr z_rotation_matrix = z_axis_rotation_projective_matrix(theta_z);
    MatrixXr rotation_matrix =
        x_rotation_matrix * y_rotation_matrix * z_rotation_matrix;
    Matrix3x3r frame = rotation_matrix.block(0, 0, 3, 3);

    // Apply transformation
    apply_camera_frame_transformation_to_vertices(initial_V, frame, V);

    // Time update
    timer.start();
    spline_surface.update_positions(V, fit_matrix, energy_hessian_inverse);
    ContourNetwork contour_network(spline_surface,
                                    intersect_params, invisibility_params,
                                    patch_boundary_edges);
    double total_time = timer.getElapsedTime();

    // Write timing data
    out_per_view << input_filename << "," << rotation_matrix.row(0) << " "
                 << rotation_matrix.row(1) << " " << rotation_matrix.row(2)
                 << "," << total_time << ","
                 << contour_network.surface_update_position_time << ","
                 << contour_network.compute_contour_time << ","
                 << contour_network.compute_cusp_time << ","
                 << contour_network.compute_intersection_time << ","
                 << contour_network.compute_visibility_time << ","
                 << contour_network.compute_projected_time << ","
                 << contour_network.segment_number << ","
                 << contour_network.interior_cusp_number << ","
                 << contour_network.boundary_cusp_number << ","
                 << contour_network.intersection_call << ","
                 << contour_network.ray_intersection_call << ","
                 << spline_surface.num_patches() << ",\n";
  }

  // Close output files
  out_per_view.close();
  out_view_independent.close();
}
