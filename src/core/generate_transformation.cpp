// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "generate_transformation.h"

Matrix3x3r
get_frame(const Eigen::Matrix<double, 3, 1>& tau)
{
  Matrix3x3r frame;
  frame.setZero();

  // Get standard frame
  Eigen::Matrix<double, 3, 1> e1(1, 0, 0);
  Eigen::Matrix<double, 3, 1> e2(0, 1, 0);
  Eigen::Matrix<double, 3, 1> e3(0, 0, 1);

  // Set the last frame vector to normalized tau
  frame.col(2) = tau / tau.norm();

  // Get candidate frame vectors by crossing tau with the standard basis
  SpatialVector candidate_frame_vector_1 =
    cross_product<double>(frame.col(2), e1);
  SpatialVector candidate_frame_vector_2 =
    cross_product<double>(frame.col(2), e2);
  SpatialVector candidate_frame_vector_3 =
    cross_product<double>(frame.col(2), e3);

  // Set the frame to the candidate with the largest magnitude
  frame.col(0) = candidate_frame_vector_1;
  if (candidate_frame_vector_2.norm() > frame.row(0).norm()) {
    frame.col(0) = candidate_frame_vector_2;
  }
  if (candidate_frame_vector_3.norm() > frame.row(0).norm()) {
    frame.col(0) = candidate_frame_vector_3;
  }

  // Normalize the first frame vector
  frame.col(0) /= frame.col(0).norm();

  // Get last frame vector from the cross product
  frame.col(1) = cross_product<double>(frame.col(2), frame.col(0));

  // TODO Ensure orientation is preserved

  return frame;
}

Eigen::Matrix<double, 4, 4>
origin_to_infinity_projective_matrix(double plane_distance)
{
  Eigen::Matrix<double, 4, 4> projection_matrix;
  projection_matrix.setZero();

  // Scale by plane distance on x, y coordinates in the fixed plane
  projection_matrix(0, 0) = plane_distance;
  projection_matrix(1, 1) = plane_distance;

  // Ensure plane points remain in the plane (and invert z coordinate)
  projection_matrix(2, 3) = -plane_distance * plane_distance;

  // The homogeneous coordinate is the original z coordinate
  projection_matrix(3, 2) = 1.0;

  return projection_matrix;
}

Eigen::Matrix<double, 4, 4>
infinity_to_origin_projective_matrix(double plane_distance)
{
  Eigen::Matrix<double, 4, 4> inverse_projection_matrix =
    origin_to_infinity_projective_matrix(plane_distance);
  return inverse_projection_matrix.inverse();
}

Eigen::Matrix<double, 4, 4>
rotate_frame_projective_matrix(const Matrix3x3r& frame)
{
  Eigen::Matrix<double, 4, 4> rotation_matrix;
  rotation_matrix.setZero();

  // The desired rotation is the transpose of the frame
  rotation_matrix.topLeftCorner(3, 3) = frame.transpose();

  // No homogeneous scaling for the rotation
  rotation_matrix(3, 3) = 1;

  return rotation_matrix;
}

Eigen::Matrix<double, 4, 4>
translation_projective_matrix(const SpatialVector& translation)
{
  assert(translation.size() == 3);

  // Initialize matrix to the identity
  Eigen::Matrix<double, 4, 4> translation_matrix;
  translation_matrix.setIdentity();

  // Add translation using homogenous coordinates
  translation_matrix.topRightCorner(3, 1) = translation.transpose();

  return translation_matrix;
}

Eigen::Matrix<double, 4, 4>
scaling_projective_matrix(double scale)
{
  Eigen::Matrix<double, 4, 4> scaling_matrix;
  scaling_matrix.setZero();

  // Scale coordinates (except homogenous coordinate)
  scaling_matrix(0, 0) = scale;
  scaling_matrix(1, 1) = scale;
  scaling_matrix(2, 2) = scale;
  scaling_matrix(3, 3) = 1.0;

  return scaling_matrix;
}

Eigen::Matrix<double, 4, 4>
x_axis_rotation_projective_matrix(double degree)
{
  Eigen::Matrix<double, 4, 4> rotation_matrix;
  rotation_matrix.setZero();
  double angle = (2 * M_PI / 360.0) * degree;

  // Leave x coordinate fixed
  rotation_matrix(0, 0) = 1.0;

  // Rotate y and z coordinates
  rotation_matrix(1, 1) = cos(angle);
  rotation_matrix(1, 2) = -sin(angle);
  rotation_matrix(2, 1) = sin(angle);
  rotation_matrix(2, 2) = cos(angle);

  // No homogeneous coordinate change
  rotation_matrix(3, 3) = 1.0;

  return rotation_matrix;
}

Eigen::Matrix<double, 4, 4>
y_axis_rotation_projective_matrix(double degree)
{
  Eigen::Matrix<double, 4, 4> rotation_matrix;
  rotation_matrix.setZero();
  double angle = (2 * M_PI / 360.0) * degree;

  // Leave y coordinate fixed
  rotation_matrix(1, 1) = 1.0;

  // Rotate x and z coordinates
  rotation_matrix(0, 0) = cos(angle);
  rotation_matrix(0, 2) = -sin(angle);
  rotation_matrix(2, 0) = sin(angle);
  rotation_matrix(2, 2) = cos(angle);

  // No homogeneous coordinate change
  rotation_matrix(3, 3) = 1.0;

  return rotation_matrix;
}

Eigen::Matrix<double, 4, 4>
z_axis_rotation_projective_matrix(double degree)
{
  Eigen::Matrix<double, 4, 4> rotation_matrix;
  rotation_matrix.setZero();
  double angle = (2 * M_PI / 360.0) * degree;

  // Leave y coordinate fixed
  rotation_matrix(2, 2) = 1.0;

  // Rotate x and z coordinates
  rotation_matrix(0, 0) = cos(angle);
  rotation_matrix(0, 1) = -sin(angle);
  rotation_matrix(1, 0) = sin(angle);
  rotation_matrix(1, 1) = cos(angle);

  // No homogeneous coordinate change
  rotation_matrix(3, 3) = 1.0;

  return rotation_matrix;
}

Eigen::Matrix<double, 4, 4>
axis_rotation_projective_matrix(double x_degree,
                                double y_degree,
                                double z_degree)
{
  Eigen::Matrix<double, 4, 4> rotation_matrix =
    z_axis_rotation_projective_matrix(z_degree);
  rotation_matrix =
    y_axis_rotation_projective_matrix(y_degree) * rotation_matrix;
  rotation_matrix =
    x_axis_rotation_projective_matrix(x_degree) * rotation_matrix;

  return rotation_matrix;
}
