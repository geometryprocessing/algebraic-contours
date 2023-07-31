// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "apply_transformation.h"
#include "compute_boundaries.h"
#include "contour_network.h"
#include "twelve_split_spline.h"
#include "generate_transformation.h"
#include "globals.cpp"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

int main(int argc, char *argv[]) {
  CLI::App app{"Generate rotation animation for a given mesh."};
  std::string input_filename = "";
  std::string output_dir = "./";
  bool view_surface = false;
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  app.add_option("--view_surface", view_surface, "View mesh and surface rotation");
  CLI11_PARSE(app, argc, argv);

  // Set logging and discretization level
  DISCRETIZATION_LEVEL = 2;
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd initial_V, V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, initial_V, uv, N, F, FT, FN);

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

  // Get the boundary edges
  std::vector<std::pair<int, int>> patch_boundary_edges(0);
  compute_twelve_split_spline_patch_boundary_edges(F, face_to_patch_indices, patch_boundary_edges);

  MatrixXr pert_frame(3, 3);
  pert_frame << 0.9218127938, -0.3785297887, 0.08352467839, 0.3815369887,
       0.9240681946, -0.0229673378, -0.06848867707, 0.0530393400,
       0.9962410000;
  double pert = 1e-2;
  double theta = 0;
  double dx = 1;
  std::string output_path;
  for (size_t i = 0; i < 360; ++i)
  {
    // Rotate the mesh
    theta = dx * i + pert;
    MatrixXr rotation_matrix = y_axis_rotation_projective_matrix(theta);
    Matrix3x3r frame = rotation_matrix.block(0, 0, 3, 3) * pert_frame;
    apply_camera_frame_transformation_to_vertices(initial_V, frame, V);

    // Optionally screenshot the mesh
    if (view_surface)
    {
      output_path = join_path(output_dir, "mesh_animation_frame_" + std::to_string(i) + ".png");
      screenshot_mesh(V, F, output_path, SpatialVector(0, 0, -1), SpatialVector(0, 0, 0), true);
    }

    // Update the spline surface vertex positions
    spline_surface.update_positions(
        V, fit_matrix, energy_hessian_inverse);

    // Optionally screenshot the surface
    if (view_surface)
    {
      output_path = join_path(output_dir, "surface_animation_frame_" + std::to_string(i) + ".png");
      spline_surface.screenshot(output_path, SpatialVector(0, 0, -1), SpatialVector(0, 0, 0), true);
    }

    // Build the contours
    IntersectionParameters intersect_params;
    InvisibilityParameters invisibility_params;
    ContourNetwork contour_network(spline_surface,
                                   intersect_params, invisibility_params,
                                   patch_boundary_edges);

    // Optionally write the raster contours if viewing the surface
    if (view_surface)
    {
      output_path = join_path(output_dir, "animation_frame_" +
        std::to_string(i) + ".png"); contour_network.view_contours();
      contour_network.write_rasterized_contours(output_path);
    }

    // Write the contours
    output_path = join_path(output_dir, "animation_frame_" + std::to_string(i) + ".svg");
    contour_network.write(output_path);
  }
}
