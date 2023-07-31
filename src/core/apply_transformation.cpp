// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "apply_transformation.h"

// Convert point to homogeneous coordinates
void
convert_point_to_homogeneous_coords(
  const SpatialVector& point,
  Eigen::Matrix<double, 4, 1>& homogeneous_coords)
{
  // Simply append a homogeneous coordinate of value 1 to the point
  homogeneous_coords.head(3) = point;
  homogeneous_coords(3) = 1.0;
}

// Convert homogeneous coordinates to a point
void
convert_homogeneous_coords_to_point(
  const Eigen::Matrix<double, 4, 1>& homogeneous_coords,
  SpatialVector& point)
{
  // Extract homogenous coordinates
  double x = homogeneous_coords(0);
  double y = homogeneous_coords(1);
  double z = homogeneous_coords(2);
  double w = homogeneous_coords(3);

  // If w is zero, return
  if (float_equal(w, 0.0))
    return;

  // Divide by the homogeneous coordinate
  point(0) = x / w;
  point(1) = y / w;
  point(2) = z / w;
}

void
apply_transformation_to_point(
  const SpatialVector& point,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  SpatialVector& transformed_point)
{
  // Get homogenous coordinates for the point
  Eigen::Matrix<double, 4, 1> homogeneous_coords;
  convert_point_to_homogeneous_coords(point, homogeneous_coords);

  // Transform the homogenous coordinates
  Eigen::Matrix<double, 4, 1> transformed_coords =
    projective_transformation * homogeneous_coords;

  // Convert transformed homogenous coordinates to a point
  convert_homogeneous_coords_to_point(transformed_coords, transformed_point);
}

void
apply_transformation_to_points(
  const std::vector<SpatialVector>& points,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  std::vector<SpatialVector>& transformed_points)
{
  transformed_points.resize(points.size());

  // Apply transformation to each point individually
  for (size_t i = 0; i < points.size(); ++i) {
    apply_transformation_to_point(
      points[i], projective_transformation, transformed_points[i]);
  }
}

void
apply_transformation_to_points_in_place(
  std::vector<SpatialVector>& points,
  const Eigen::Matrix<double, 4, 4>& projective_transformation)
{
  // Apply transformation to each point individually
  for (size_t i = 0; i < points.size(); ++i) {
    SpatialVector point = points[i];
    apply_transformation_to_point(point, projective_transformation, points[i]);
  }
}

void
apply_transformation_to_control_points(
  const std::vector<std::vector<SpatialVector>>& input_control_point_grid,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  std::vector<std::vector<SpatialVector>>& output_control_point_grid)
{
  output_control_point_grid.resize(input_control_point_grid.size());

  // Apply transformation to each control point individually
  for (size_t i = 0; i < input_control_point_grid.size(); ++i) {
    apply_transformation_to_points(input_control_point_grid[i],
                                   projective_transformation,
                                   output_control_point_grid[i]);
  }
}

void
apply_transformation_to_control_points_in_place(
  std::vector<std::vector<SpatialVector>>& control_point_grid,
  const Eigen::Matrix<double, 4, 4>& projective_transformation)
{
  // Apply transformation to each control point individually
  for (size_t i = 0; i < control_point_grid.size(); ++i) {
    apply_transformation_to_points_in_place(control_point_grid[i],
                                            projective_transformation);
  }
}

void
apply_transformation_to_vertices(
  const MatrixXr& input_V,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  MatrixXr& output_V)
{
  // Apply transformation to each vertex point individually
  output_V.resize(input_V.rows(), input_V.cols());
  for (Eigen::Index i = 0; i < input_V.rows(); ++i) {
    SpatialVector v = input_V.row(i);
    SpatialVector Tv;
    apply_transformation_to_point(v, projective_transformation, Tv);
    output_V.row(i) = Tv;
  }
}

void
apply_transformation_to_vertices_in_place(
  MatrixXr& V,
  const Eigen::Matrix<double, 4, 4>& projective_transformation)
{
  // Apply transformation to each vertex point individually
  for (Eigen::Index i = 0; i < V.rows(); ++i) {
    SpatialVector v = V.row(i);
    SpatialVector Tv;
    apply_transformation_to_point(v, projective_transformation, Tv);
    V.row(i) = Tv;
  }
}

void
generate_projective_transformation(
  const Matrix3x3r& input_frame,
  const TransformationParameters& transformation_params,
  Eigen::Matrix<double, 4, 4>& projective_transformation)
{
  // Generate the matrix to align the frame with the standard frame
  spdlog::trace("Aligning frame {}", input_frame);
  Eigen::Matrix<double, 4, 4> frame_rotation_matrix =
    rotate_frame_projective_matrix(input_frame);
  spdlog::trace("Frame rotation matrix:\n{}", frame_rotation_matrix);
  projective_transformation = frame_rotation_matrix;

  // Generate the rotation matrix for rotation along the axes
  spdlog::trace("Applying rotation ({}, {}, {})",
                transformation_params.x_rotation,
                transformation_params.y_rotation,
                transformation_params.z_rotation);
  Eigen::Matrix<double, 4, 4> axis_rotation_matrix =
    axis_rotation_projective_matrix(transformation_params.x_rotation,
                                    transformation_params.y_rotation,
                                    transformation_params.z_rotation);
  spdlog::trace("Axis rotation matrix:\n{}", axis_rotation_matrix);
  projective_transformation = axis_rotation_matrix * projective_transformation;

  // Get translation of control point grid from the origin
  spdlog::trace("Applying translation of {}",
                transformation_params.translation);
  double z_distance = transformation_params.translation;

  // Generate the translation matrix
  SpatialVector translation(3);
  translation << 0.0, 0.0, z_distance;
  Eigen::Matrix<double, 4, 4> translation_matrix =
    translation_projective_matrix(translation);
  spdlog::trace("Translation matrix:\n{}", translation_matrix);
  projective_transformation = translation_matrix * projective_transformation;

  // Optionally apply a projective transformation to send the camera to infinity
  if (transformation_params.send_camera_to_infinity) {
    spdlog::trace("Sending the camera to infinity");
    Eigen::Matrix<double, 4, 4> projection_matrix =
      origin_to_infinity_projective_matrix(
        transformation_params.camera_to_plane_distance);
    spdlog::trace("Projection matrix:\n{}", projection_matrix);
    projective_transformation = projection_matrix * projective_transformation;
  }
}

void
initialize_control_points(
  const std::vector<std::vector<SpatialVector>>& input_control_point_grid,
  const Matrix3x3r& input_frame,
  const TransformationParameters& transformation_params,
  std::vector<std::vector<SpatialVector>>& output_control_point_grid,
  Matrix3x3r& output_frame)
{
  // Get the standard frame
  output_frame.resize(3, 3);
  output_frame.setIdentity();

  // Generate the transformation for the control points
  Eigen::Matrix<double, 4, 4> projective_transformation;
  generate_projective_transformation(
    input_frame, transformation_params, projective_transformation);

  // Apply transformations to the control point grid
  apply_transformation_to_control_points(input_control_point_grid,
                                         projective_transformation,
                                         output_control_point_grid);
}

void
initialize_vertices(const MatrixXr& input_vertices,
                    const Matrix3x3r& input_frame,
                    const TransformationParameters& transformation_params,
                    MatrixXr& output_vertices,
                    Matrix3x3r& output_frame)
{
  // Get the standard frame
  output_frame.resize(3, 3);
  output_frame.setIdentity();

  // Generate the transformation for the vertices
  Eigen::Matrix<double, 4, 4> projective_transformation;
  generate_projective_transformation(
    input_frame, transformation_params, projective_transformation);

  // Apply transformations to the control point grid
  apply_transformation_to_vertices(
    input_vertices, projective_transformation, output_vertices);
}

void
apply_camera_frame_transformation_to_vertices(const MatrixXr& input_V,
                                              const Matrix3x3r& frame,
                                              Eigen::MatrixXd& output_V,
                                              bool orthographic,
                                              bool recenter_mesh)
{
  // Compute mesh midpoint and bounding box diagonal
  SpatialVector min_point;
  SpatialVector max_point;
  compute_point_cloud_bounding_box(input_V, min_point, max_point);
  SpatialVector mesh_midpoint = 0.5 * (max_point + min_point);
  SpatialVector bounding_box_diagonal = max_point - min_point;
  spdlog::info("Initial mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Initial mesh midpoint: {}", mesh_midpoint);

  // Normalize the vertices
  double scale_factor = bounding_box_diagonal.maxCoeff();
  size_t num_vertices = input_V.rows();
  output_V.resize(num_vertices, 3);
  for (size_t i = 0; i < num_vertices; ++i) {
    output_V.row(i) = 2.0 * (input_V.row(i) - mesh_midpoint) / scale_factor;
  }
  compute_point_cloud_bounding_box(input_V, min_point, max_point);
  mesh_midpoint = 0.5 * (max_point + min_point);
  bounding_box_diagonal = max_point - min_point;
  spdlog::info("Normalized mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Normalized mesh midpoint: {}", mesh_midpoint);

  // Generate rotation matrix
  spdlog::info("Projecting onto frame:\n{}", frame);
  Eigen::Matrix<double, 4, 4> frame_rotation_matrix =
    rotate_frame_projective_matrix(frame);

  // Generate translation matrix
  double z_distance = 4.0;
  SpatialVector translation(3);
  translation << 0.0, 0.0, z_distance;
  Eigen::Matrix<double, 4, 4> translation_matrix =
    translation_projective_matrix(translation);
  
  // Optionally generate matrix to send the origin to infinity
  Eigen::Matrix<double, 4, 4> projective_transformation;
  if (orthographic)
  {
    double camera_to_plane_distance = 1.0;
    Eigen::Matrix<double, 4, 4> projection_matrix =
      origin_to_infinity_projective_matrix(camera_to_plane_distance);
    projective_transformation = projection_matrix * translation_matrix * frame_rotation_matrix;
  }
  else
  {
    projective_transformation = translation_matrix * frame_rotation_matrix;
  }

  // Apply the transformations
  spdlog::info("Apply transformation:\n{}", projective_transformation);
  apply_transformation_to_vertices_in_place(output_V,
                                            projective_transformation);

  if (recenter_mesh) {
    // Renormalize the projected vertices
    compute_point_cloud_bounding_box(output_V, min_point, max_point);
    mesh_midpoint = 0.5 * (max_point + min_point);
    spdlog::info("Projected mesh bounding box: {}, {}", min_point, max_point);
    spdlog::info("Projected mesh midpoint: {}", mesh_midpoint);
    bounding_box_diagonal = max_point - min_point;
    scale_factor = bounding_box_diagonal.norm();
    for (size_t i = 0; i < num_vertices; ++i) {
      output_V.row(i) = (output_V.row(i) - mesh_midpoint) / scale_factor;
    }
  }

  // Check final midpoint location
  compute_point_cloud_bounding_box(output_V, min_point, max_point);
  mesh_midpoint = 0.5 * (max_point + min_point);
  spdlog::info("Final mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Final mesh midpoint: {}", mesh_midpoint);
}