// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>

#include "common.h"
#include "generate_transformation.h"

/// \file apply_transformation.h
///
/// Methods to apply projective transformation matrices to various point data
/// types.

struct TransformationParameters
{
  // Distance from the camera to the view plane
  double camera_to_plane_distance = 1;

  // Translation to apply to the surface
  double translation = 1;

  // Rotation parameters to apply to the surface
  double x_rotation = 0;
  double y_rotation = 0;
  double z_rotation = 0;

  // If true, apply a projective transformation to send the camera at the origin
  // to infinity
  bool send_camera_to_infinity = true;
};

void
convert_point_to_homogeneous_coords(
  const SpatialVector& point,
  Eigen::Matrix<double, 4, 1>& homogeneous_coords);

// Convert homogeneous coordinates to a point
void
convert_homogeneous_coords_to_point(
  const Eigen::Matrix<double, 4, 1>& homogeneous_coords,
  SpatialVector& point);

/// Apply projective transformation to a single point
///
/// @param[in] point: point to transform
/// @param[in] projective_transformation: transformation to apply to the point
/// @param[out] transformed_point: point after transformation
void
apply_transformation_to_point(
  const SpatialVector& point,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  SpatialVector& transformed_point);

/// Apply projective transformation to a vector of points.
///
/// @param[in] points: points to transform
/// @param[in] projective_transformation: transformation to apply to the vector
/// @param[out] transformed_points: transformed points
void
apply_transformation_to_points(
  const std::vector<SpatialVector>& points,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  std::vector<SpatialVector>& transformed_points);

/// Apply projective transformation to a vector of points in place.
///
/// @param[in, out] points: points to transform
/// @param[in] projective_transformation: transformation to apply to the vector
void
apply_transformation_to_points_in_place(
  std::vector<SpatialVector>& points,
  const Eigen::Matrix<double, 4, 4>& projective_transformation);

/// Apply projective transformation to a grid of control points.
///
/// @param[in] input_control_point_grid: control point grid to transform
/// @param[in] projective_transformation: transformation to apply to the vector
/// @param[out] output_control_point_grid: transformed control point grid
void
apply_transformation_to_control_points(
  const std::vector<std::vector<SpatialVector>>& input_control_point_grid,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  std::vector<std::vector<SpatialVector>>& output_control_point_grid);

/// Apply projective transformation to a grid of control points in place.
///
/// @param[in, out] control_point_grid: control point grid to transform
/// @param[in] projective_transformation: transformation to apply to the vector
void
apply_transformation_to_control_points_in_place(
  std::vector<std::vector<SpatialVector>>& control_point_grid,
  const Eigen::Matrix<double, 4, 4>& projective_transformation);

/// Apply projective transformation to a matrix of vertices.
///
/// @param[in] input_V: vertex matrix to transform
/// @param[in] projective_transformation: transformation to apply to the vector
/// @param[out] output_V: transformed vertex matrix
void
apply_transformation_to_vertices(
  const MatrixXr& input_V,
  const Eigen::Matrix<double, 4, 4>& projective_transformation,
  MatrixXr& output_V);

/// Apply projective transformation to a matrix of vertices in place.
///
/// @param[in, out] V: vertex matrix to transform
/// @param[in] projective_transformation: transformation to apply to the vector
void
apply_transformation_to_vertices_in_place(
  MatrixXr& V,
  const Eigen::Matrix<double, 4, 4>& projective_transformation);

/// Apply transformations to the control point grid and frame to align the view
/// direction with the z axis, set viewing distances and angles, and send the
/// camera to infinity.
///
/// @param[in] input_control_point_grid: control points defining surface before
/// transformation
/// @param[in] input_frame: frame defining the contour view direction before
/// transformation
/// @param[in] transformation_params: parameters for the transformation to apply
/// @param[out] output_control_point_grid: control points defining surface after
/// transformation
/// @param[out] output_frame: frame defining the contour view direction after
/// transformation
void
initialize_control_points(
  const std::vector<std::vector<SpatialVector>>& input_control_point_grid,
  const Matrix3x3r& input_frame,
  const TransformationParameters& transformation_params,
  std::vector<std::vector<SpatialVector>>& output_control_point_grid,
  Matrix3x3r& output_frame);

/// Apply transformations to a mesh to set the initial camera view direction
/// with optional conversion to an orthographic perspective
///
/// The vertices are optimally first normalized to be in a bounding box of
/// diagonal 1 at the origin.
///
/// @param[in] input_V: mesh vertices before the transformation
/// @param[in] camera_direction: view direction for the mesh
/// @param[out] output_V: mesh vertices after transformation
/// @param[in] orthographic: project camera to infinity if true
/// @param[in] normalize_initial_positions: normalize the initial vertices to a
/// bounding box if true
void
apply_camera_frame_transformation_to_vertices(const MatrixXr& input_V,
                                              const Matrix3x3r& frame,
                                              Eigen::MatrixXd& output_V,
                                              bool orthographic = true,
                                              bool recenter_mesh = false);