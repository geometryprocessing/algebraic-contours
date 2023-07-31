// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "apply_transformation.h"
#include "compute_boundaries.h"
#include "contour_network.h"
#include "generate_transformation.h"
#include "twelve_split_spline.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

/// \file generate_example_figure.cpp
///
/// Executable to generate a pipeline example from a mesh and camera specification.
/// Generates an image of the mesh, the quadratic surface (with front camera view),
/// the full contours with cusps and intersections, and the final occluding contours.

// The following rotation angles were originally used for the below meshes, in
// order x, y, z.
//
// Spot 1 -30 1
// Bigguy: 0 150 0
// Blub: 0 150 0
// Fertility: 0 150 0
// Bob: 15 -30 15
// Killaroo: 90 30 0
// Monsterfrog: 0 150 0
// Pawn: 90 0 0
// Pipes: 0 170 0
// Ogre 0 150 0
// Toad 75 130 0


int main(int argc, char *argv[]) {
  // Get command line arguments
  CLI::App app{"Generate example figure images for a given mesh and camera."};
  std::string input_filename = "";
  std::string output_dir = "./";
  std::string camera_filename = "";
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-c,--camera", camera_filename, "Camera filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  CLI11_PARSE(app, argc, argv);

  // Set logger level
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);

  // Set up the camera
  Eigen::Matrix<double, 4, 4> camera_matrix, projection_matrix;
  double camera_to_plane_distance = 1.0;
  read_camera_matrix(camera_filename, camera_matrix);
  spdlog::info("Using camera matrix:\n{}", camera_matrix);

  // Apply camera and perspective projection transformations
  projection_matrix = origin_to_infinity_projective_matrix(camera_to_plane_distance);
  projection_matrix = projection_matrix * camera_matrix;
  apply_transformation_to_vertices_in_place(V, projection_matrix);

  // Generate quadratic spline
  spdlog::info("Computing spline surface");
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


  // Get the boundary edges
  std::vector<std::pair<int, int>> patch_boundary_edges(0);
  compute_twelve_split_spline_patch_boundary_edges(F, face_to_patch_indices, patch_boundary_edges);

  // Build the contours
  spdlog::info("Computing contours");
  IntersectionParameters intersect_params;
  InvisibilityParameters invisibility_params;
  ContourNetwork contour_network(spline_surface,
                                 intersect_params, invisibility_params,
                                 patch_boundary_edges);

  // View the contours
  affine_manifold.screenshot(
    join_path(output_dir, "mesh.png"),
    V,
    SpatialVector(0, 0, -1),
    SpatialVector(0, 0, 0),
    true
  );
  contour_network.screenshot(
    join_path(output_dir, "spline_surface.png"),
    spline_surface,
    SpatialVector(0, 0, -1),
    SpatialVector(0, 0, 0),
    true
  );
  contour_network.write(
    join_path(output_dir, "full_contours.svg"),
    SVGOutputMode::contrast_invisible_segments,
    true
  );
  contour_network.write(
    join_path(output_dir, "contours.svg"),
    SVGOutputMode::uniform_visible_chains,
    false
  );
}
