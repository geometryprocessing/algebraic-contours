// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "apply_transformation.h"
#include "common.h"
#include "generate_transformation.h"
#include "twelve_split_spline.h"
#include "polyscope/polyscope.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

int main(int argc, char *argv[]) {
  // Get command line arguments
  CLI::App app{"Generate perspective distortion figure images for a given mesh and camera."};
  std::string input_filename = "";
  std::string output_dir = "./";
  std::string camera_filename = "";
  double translation = 0.0;
  double perspective_fov = 45;
  double orthographic_fov = 45;
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-c,--camera", camera_filename, "Camera filepath")
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  app.add_option("--translation", translation, "translation amount");
  app.add_option("--perspective_fov", perspective_fov, "Field of view for perspective rendering");
  app.add_option("--orthographic_fov", orthographic_fov, "Field of view for orthographic rendering");
  CLI11_PARSE(app, argc, argv);

  // Set logging and discretization level
  DISCRETIZATION_LEVEL = 2;
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);

  // Set up polyscope
  polyscope::init();

  // Set up the camera
  if (camera_filename == "")
  {
    Matrix3x3r frame(3, 3);
    frame <<
      1, 0, 0,
      0, 1, 0,
      0, 0, 1;
    MatrixXr rotation_matrix = axis_rotation_projective_matrix(100, 10, 0);
    Matrix3x3r rotation_frame = rotation_matrix.block(0, 0, 3, 3);
    frame = rotation_frame;
    spdlog::info("Projecting onto frame:\n{}", frame);
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

  // Apply additional translation
  int num_vertices = V.rows();
  for (int i = 0; i < num_vertices; ++i) {
    V.row(i) = V.row(i) + SpatialVector(0, 0, translation);
  }

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
      V, affine_manifold,
      optimization_params, face_to_patch_indices, patch_to_face_indices,
      fit_matrix, energy_hessian, energy_hessian_inverse);
      
  // View the perspective spline surface
  polyscope::view::fov = perspective_fov;
  spline_surface.screenshot(
    join_path(output_dir, "perspective_spline_surface.png"),
    SpatialVector(0, 0, 0),
    SpatialVector(0, 0, 1),
    false
  );

  // Send the camera to infinity and update the vertex positions
  MatrixXr projection_matrix = origin_to_infinity_projective_matrix(1.0);
  MatrixXr orthographic_V;
  apply_transformation_to_vertices(V, projection_matrix, orthographic_V);
  spline_surface.update_positions(orthographic_V, fit_matrix,
                                  energy_hessian_inverse);
  
  // View the orthographic spline surface
  polyscope::view::fov = orthographic_fov;
  spline_surface.screenshot(
    join_path(output_dir, "orthographic_spline_surface.png"),
    SpatialVector(0, 0, -5),
    SpatialVector(0, 0, 0),
    true
  );
}
