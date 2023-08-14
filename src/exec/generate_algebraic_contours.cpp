// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "globals.cpp"
#include "apply_transformation.h"
#include "generate_transformation.h"
#include "compute_boundaries.h"
#include "contour_network.h"
#include "twelve_split_spline.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>


int main(int argc, char *argv[])
{
  // Build maps from strings to enums
  std::map<std::string, spdlog::level::level_enum> log_level_map {
    {"trace",    spdlog::level::trace},
    {"debug",    spdlog::level::debug},
    {"info",     spdlog::level::info},
    {"warn",     spdlog::level::warn},
    {"critical", spdlog::level::critical},
    {"off",      spdlog::level::off},
  };
  std::map<std::string, InvisibilityMethod> invisibility_method_map{
    {"none",        InvisibilityMethod::none},
    {"direct",      InvisibilityMethod::direct},
    {"chaining",    InvisibilityMethod::chaining},
    {"propagation", InvisibilityMethod::propagation},
  };
  std::map<std::string, SVGOutputMode> svg_output_mode_map{
    {"visible",    SVGOutputMode::uniform_visible_curves},
    {"simplified", SVGOutputMode::uniform_simplified_visible_curves},
    {"closed",     SVGOutputMode::uniform_closed_curves},
    {"contrast",   SVGOutputMode::contrast_invisible_segments},
    {"chains",     SVGOutputMode::random_chains},
  };

  // Get command line arguments
  CLI::App app{"Generate smooth occluding contours for a mesh."};
  std::string input_filename = "";
  std::string output_dir = "./";
  std::string camera_filename = "";
  spdlog::level::level_enum log_level = spdlog::level::off;
  SVGOutputMode svg_output_mode = SVGOutputMode::uniform_visible_curves;
  OptimizationParameters optimization_params;
  IntersectionParameters intersect_params;
  InvisibilityParameters invisibility_params;
  double weight = optimization_params.position_difference_factor;
  double trim = intersect_params.trim_amount;
  double pad = invisibility_params.pad_amount;
  InvisibilityMethod invisibility_method = invisibility_params.invisibility_method;
  bool show_nodes = false;
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  app.add_option("-c,--camera", camera_filename, "Camera filepath")
    ->check(CLI::ExistingFile);
  app.add_option("--log_level", log_level, "Level of logging")
    ->transform(CLI::CheckedTransformer(log_level_map, CLI::ignore_case));
  app.add_option("--invisibility_method", invisibility_method, "Method for invisibility tests of contours")
    ->transform(CLI::CheckedTransformer(invisibility_method_map, CLI::ignore_case));
  app.add_option("--svg_mode", svg_output_mode, "Output mode for SVG contours")
    ->transform(CLI::CheckedTransformer(svg_output_mode_map, CLI::ignore_case));
  app.add_option("-w,--weight", weight, "Fitting weight for the quadratic surface approximation")
    ->check(CLI::PositiveNumber);
  app.add_option("--trim", trim, "Trimming for contour intersection checks")
    ->check(CLI::NonNegativeNumber);
  app.add_option("--pad", pad, "Padding for contour chaining checks")
    ->check(CLI::NonNegativeNumber);
  app.add_flag("--show_nodes", show_nodes, "Show important nodes in the contours");
  CLI11_PARSE(app, argc, argv);

  // Set logger level
  spdlog::set_level(log_level);

  // Set optimization parameters
  optimization_params.position_difference_factor = weight;

  // Set intersection parameters
  intersect_params.trim_amount = trim;

  // Set invisibility parameters
  invisibility_params.pad_amount = pad;
  invisibility_params.invisibility_method = invisibility_method;

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);

  // Set up the camera
  if (camera_filename == "")
  {
    Matrix3x3r frame(3, 3);
    frame <<
      1, 0, 0,
      0, 1, 0,
      0, 0, 1;
    spdlog::info("Projecting onto frame:\n{}", frame);
    MatrixXr V_copy = V;
    apply_camera_frame_transformation_to_vertices(V_copy, frame, V);
  } else
  {
    Eigen::Matrix<double, 4, 4> camera_matrix, projection_matrix;
    double camera_to_plane_distance = 1.0;
    read_camera_matrix(camera_filename, camera_matrix);
    spdlog::info("Using camera matrix:\n{}", camera_matrix);

    projection_matrix = origin_to_infinity_projective_matrix(camera_to_plane_distance);
    projection_matrix = projection_matrix * camera_matrix;
    apply_transformation_to_vertices_in_place(V, projection_matrix);
  }

  // Generate quadratic spline
  spdlog::info("Computing spline surface");
  std::vector<std::vector<int>> face_to_patch_indices;
  std::vector<int> patch_to_face_indices;
  Eigen::SparseMatrix<double> fit_matrix;
  Eigen::SparseMatrix<double> energy_hessian;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>
      energy_hessian_inverse;
  AffineManifold affine_manifold(F, uv, FT);
  TwelveSplitSplineSurface spline_surface(
      V, affine_manifold,
      optimization_params, face_to_patch_indices, patch_to_face_indices,
      fit_matrix, energy_hessian, energy_hessian_inverse);

  // Get the boundary edges
	std::vector<std::pair<int, int>> patch_boundary_edges(0);
  compute_twelve_split_spline_patch_boundary_edges(F, face_to_patch_indices, patch_boundary_edges);

  // Build the contours
  spdlog::info("Computing contours");
  ContourNetwork contour_network(
    spline_surface,
    intersect_params,
    invisibility_params,
    patch_boundary_edges
  );

	// Save the contours to file
  spdlog::info("Saving contours");
  std::string contour_network_file ("contours.svg");
  std::string contour_network_path = join_path(output_dir, contour_network_file);
  contour_network.write(contour_network_path, svg_output_mode, show_nodes);
}
