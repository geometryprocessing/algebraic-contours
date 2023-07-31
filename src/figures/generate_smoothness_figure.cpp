// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "apply_transformation.h"
#include "generate_transformation.h"
#include "quadratic_spline_surface.h"
#include "twelve_split_spline.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

int main(int argc, char *argv[]) {
  // Get command line arguments
  CLI::App app{"Generate smoothness figure images for a given mesh and camera."};
  std::string input_filename = "";
  std::string output_dir = "./";
  std::string camera_filename = "";
  double w_f = 1.0;
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-c,--camera", camera_filename, "Camera filepath")
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  app.add_option("-w, --weight", w_f, "Fitting weight")
    ->check(CLI::NonNegativeNumber);
  CLI11_PARSE(app, argc, argv);

  // Set logger level
  spdlog::set_level(spdlog::level::off);

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
    double theta_x = 0;
    double theta_y = 160;
    double theta_z = 0;
    MatrixXr rotation_matrix = axis_rotation_projective_matrix(theta_x, theta_y, theta_z);
    spdlog::info("Projecting onto frame:\n{}", frame);
    frame = rotation_matrix.block(0, 0, 3, 3);
    MatrixXr V_copy = V;
    bool orthographic = false;
    apply_camera_frame_transformation_to_vertices(V_copy, frame, V, orthographic);
  } else
  {
    Eigen::Matrix<double, 4, 4> camera_matrix, projection_matrix;
    read_camera_matrix(camera_filename, camera_matrix);
    spdlog::info("Using camera matrix:\n{}", camera_matrix);
    apply_transformation_to_vertices_in_place(V, camera_matrix);
  }

  // View the initial mesh
  screenshot_mesh(
    V,
    F,
    join_path(output_dir, "mesh.png"),
    SpatialVector(0, 0, 0),
    SpatialVector(0, 0, 1)
  );

  // Generate quadratic spline
  OptimizationParameters optimization_params;
  optimization_params.position_difference_factor = w_f;
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

  // View the quadratic surface
  polyscope::init();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  spline_surface.screenshot(
    join_path(output_dir, "spline_surface.png"),
    SpatialVector(0, 0, 0),
    SpatialVector(0, 0, 1),
    false
  );
}
