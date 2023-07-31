// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "apply_transformation.h"
#include "common.h"
#include "generate_transformation.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include <igl/readOBJ.h>
#include <CLI/CLI.hpp>

// Generate a camera frustum curve network
void generate_frustum(double distance, double height, double width,
                      MatrixXr &points,
                      std::vector<std::array<int, 2>> &edges) {
  points.resize(8, 3);
  std::vector<std::vector<int>> polylines(0);
  double x0 = -width / 2.0;
  double x1 = width / 2.0;
  double y0 = -height / 2.0;
  double y1 = height / 2.0;
  double d = distance;

  // Add points for bounding box
  points.row(0) = Eigen::Vector3d(x0 / d, y0 / d, 1.0);
  points.row(1) = Eigen::Vector3d(x1 / d, y0 / d, 1.0);
  points.row(2) = Eigen::Vector3d(x1 / d, y1 / d, 1.0);
  points.row(3) = Eigen::Vector3d(x0 / d, y1 / d, 1.0);
  points.row(4) = Eigen::Vector3d(x0, y0, d);
  points.row(5) = Eigen::Vector3d(x1, y0, d);
  points.row(6) = Eigen::Vector3d(x1, y1, d);
  points.row(7) = Eigen::Vector3d(x0, y1, d);

  // Add polylines for bounding box
  std::vector<int> polyline_bottom(5);
  polyline_bottom[0] = 0;
  polyline_bottom[1] = 1;
  polyline_bottom[2] = 2;
  polyline_bottom[3] = 3;
  polyline_bottom[4] = 0;
  std::vector<int> polyline_top(5);
  polyline_top[0] = 4;
  polyline_top[1] = 5;
  polyline_top[2] = 6;
  polyline_top[3] = 7;
  polyline_top[4] = 4;
  std::vector<int> polyline_edge_0(2);
  polyline_edge_0[0] = 0;
  polyline_edge_0[1] = 4;
  std::vector<int> polyline_edge_1(2);
  polyline_edge_1[0] = 1;
  polyline_edge_1[1] = 5;
  std::vector<int> polyline_edge_2(2);
  polyline_edge_2[0] = 2;
  polyline_edge_2[1] = 6;
  std::vector<int> polyline_edge_3(2);
  polyline_edge_3[0] = 3;
  polyline_edge_3[1] = 7;
  polylines.push_back(polyline_bottom);
  polylines.push_back(polyline_top);
  polylines.push_back(polyline_edge_0);
  polylines.push_back(polyline_edge_1);
  polylines.push_back(polyline_edge_2);
  polylines.push_back(polyline_edge_3);
  edges = convert_polylines_to_edges(polylines);
}

// Initialize vertices for frustum view
void initialize_vertices(
  const MatrixXr &input_V,
  const MatrixXr &frame,
  Eigen::MatrixXd &output_V
) {
  // Get bounding box of mesh
  VectorXr min_point;
  VectorXr max_point;
  compute_point_cloud_bounding_box(input_V, min_point, max_point);
  VectorXr mesh_midpoint = 0.5 * (max_point + min_point);
  VectorXr bounding_box_diagonal = max_point - min_point;
  spdlog::info("Initial mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Initial mesh midpoint: {}", mesh_midpoint);

  // Normalize the vertices
  double scale_factor =  bounding_box_diagonal.norm();
  size_t num_vertices = input_V.rows();
  output_V.resize(num_vertices, 3);
  for (size_t i = 0; i < num_vertices; ++i)
  {
    output_V.row(i) = 4.0 * (input_V.row(i).transpose() - mesh_midpoint) / scale_factor;
  }
  compute_point_cloud_bounding_box(input_V, min_point, max_point);
  mesh_midpoint = 0.5 * (max_point + min_point);
  bounding_box_diagonal = max_point - min_point;
  spdlog::info("Normalized mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Normalized mesh midpoint: {}", mesh_midpoint);

  // Generate rotation matrix
  spdlog::info("Projecting onto frame:\n{}", frame);
  MatrixXr frame_rotation_matrix = rotate_frame_projective_matrix(frame);
  
  // Generate translation matrix
  double z_distance = 3.0;
  VectorXr translation(3);
  translation << 0.0,
                 0.0,
                 z_distance;
  MatrixXr translation_matrix = translation_projective_matrix(translation);

  // Apply the transformations
  MatrixXr projective_transformation = translation_matrix * frame_rotation_matrix;
  apply_transformation_to_vertices_in_place(output_V, projective_transformation);
}

int main(int argc, char *argv[]) {
  // Get command line arguments
  CLI::App app{"Generate perspective figure example for a given mesh."};
  std::string input_filename = "";
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  CLI11_PARSE(app, argc, argv);

  // Set logger level
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd input_V, V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, input_V, uv, N, F, FT, FN);

  // Initialize vertices
  MatrixXr rotation_matrix = axis_rotation_projective_matrix(0, -90, 0);
  Matrix3x3r rotation_frame = rotation_matrix.block(0, 0, 3, 3);
  initialize_vertices(input_V, rotation_frame, V);

  // Build camera frustum
  double distance = 4.0;
  double height = 5.0;
  double width = 5.0;
  MatrixXr perspective_points;
  std::vector<std::array<int, 2>> edges;
  generate_frustum(distance, height, width, perspective_points, edges);

  // Build perspective viewer
  polyscope::init();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::registerCurveNetwork("perspective_frustum", perspective_points,
                                  edges);
  polyscope::getCurveNetwork("perspective_frustum")
      ->setColor(glm::vec3(0.0, 0.0, 0.0))
      ->setRadius(0.002);
  polyscope::registerSurfaceMesh("perspective_mesh", V, F);
  polyscope::getSurfaceMesh("perspective_mesh")
      ->setSurfaceColor(glm::vec3(0.670, 0.673, 0.292));

  // Look at viewer
  glm::vec3 glm_camera_position = {0.0, 0.0, 0.0};
  glm::vec3 glm_camera_target = {0.0, 0.0, 4.0};
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::show();
  polyscope::removeAllStructures();

  // Project mesh and frustum to send the camera to infinity while fixing the
  // plane z = 1
  MatrixXr projection_matrix = origin_to_infinity_projective_matrix(1.0);
  MatrixXr orthographic_points, orthographic_V;
  apply_transformation_to_vertices(perspective_points, projection_matrix,
                                   orthographic_points);
  apply_transformation_to_vertices(V, projection_matrix, orthographic_V);

  // Build orthographic viewer
  polyscope::registerCurveNetwork("orthographic_frustum", orthographic_points,
                                  edges);
  polyscope::registerSurfaceMesh("orthographic_mesh", orthographic_V, F);
  polyscope::getCurveNetwork("orthographic_frustum")
      ->setColor(glm::vec3(0.0, 0.0, 0.0))
      ->setRadius(0.002);
  polyscope::getSurfaceMesh("orthographic_mesh")
      ->setSurfaceColor(glm::vec3(0.670, 0.673, 0.292));

  // Look at viewer
  glm_camera_position = { 0.0, 0.0, -1.0 };
  glm_camera_target = { 0.0, 0.0, 1.0 };
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::show();
  polyscope::removeAllStructures();
}
