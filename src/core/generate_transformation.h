// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// \file generate_transformation.h
///
/// Methods to generate projective transformation matrices.

/// Get a 3x3 matrix representing an orientation preserving frame aligned with
/// tau.
///
/// @param[in] tau: vector determining the frame (up to rotation)
/// @return 3x3 matrix with frame vectors as the rows, including the direction
///         determined by tau as the third vector
Matrix3x3r
get_frame(const Eigen::Matrix<double, 3, 1>& tau);

/// Generate the projective matrix that sends the origin to infinity while
/// fixing the plane z = plane_distance
///
/// @param[in] plane_distance: distance from the origin to the plane
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
origin_to_infinity_projective_matrix(double plane_distance);

/// Generate the projective matrix that sends a point at infinity to the origin
/// while fixing the plane z = plane_distance.
///
/// This is the inverse of the map sending the origin to infinity.
///
/// @param[in] plane_distance: distance from the origin to the plane
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
infinity_to_origin_projective_matrix(double plane_distance);

/// Generate the rotation matrix that sends the given frame to the standard
/// frame.
///
/// @param[in] frame: 3x3 frame matrix to align with the standard frame
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
rotate_frame_projective_matrix(const Matrix3x3r& frame);

/// Generate the projective matrix representing translation by the given
/// translation vector.
///
/// @param[in] translation: 1x3 translation vector
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
translation_projective_matrix(const SpatialVector& translation);

/// Generate the projective matrix for a given uniform scaling.
///
/// @param[in] scale: uniform scale factor
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
scaling_projective_matrix(double scale);

/// Generate the projective matrix for rotation around the x axis.
///
/// @param[in] degree: degree of rotation around the x axis
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
x_axis_rotation_projective_matrix(double degree);

/// Generate the projective matrix for rotation around the y axis.
///
/// @param[in] degree: degree of rotation around the y axis
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
y_axis_rotation_projective_matrix(double degree);

/// Generate the projective matrix for rotation around the z axis.
///
/// @param[in] degree: degree of rotation around the z axis
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
z_axis_rotation_projective_matrix(double degree);

/// Generate the projective matrix for chained rotation around the standard
/// axes.
///
/// The order of rotation is z axis -> y axis -> x axis
///
/// @param[in] x_degree: degree of rotation around the x axis
/// @param[in] y_degree: degree of rotation around the y axis
/// @param[in] z_degree: degree of rotation around the z axis
/// @return 4x4 projective matrix for the transformation
Eigen::Matrix<double, 4, 4>
axis_rotation_projective_matrix(double x_degree,
                                double y_degree,
                                double z_degree);
