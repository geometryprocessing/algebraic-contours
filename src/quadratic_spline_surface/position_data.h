// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "affine_manifold.h"
#include "common.h"
// #include "conformal_ideal_delaunay/ConformalInterface.hh"
#include "line_segment.h"
#include "polynomial_function.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "rational_function.h"
#include "vertex_circulator.h"

/// \file evaluate_surface.h
///
/// Methods to generate a per triangle local degrees of freedom from the global
/// degrees of freedom

/// Position and derivative data at the corner of a triangle.
template<typename VectorX>
struct TriangleCornerData
{
  VectorX function_value; // position value vector at the corner
  VectorX first_edge_derivative; // derivative in the ccw edge direction
  VectorX second_edge_derivative; // derivate in the clockwise edge direction

  /// Default constructor
  ///
  /// @param[in] input_function_value: position value vector at the corner
  /// @param[in] input_first_edge_derivative: derivative in the ccw edge direction
  /// @param[in] input_second_edge_derivative: derivate in the clockwise edge direction
  TriangleCornerData(const VectorX& input_function_value,
                     const VectorX& input_first_edge_derivative,
                     const VectorX& input_second_edge_derivative)
  {
    function_value = input_function_value;
    first_edge_derivative = input_first_edge_derivative;
    second_edge_derivative = input_second_edge_derivative;
  }

  /// Trivial constructor
  TriangleCornerData()
  {
    function_value = VectorX();
    first_edge_derivative = VectorX();
    second_edge_derivative = VectorX();
  }
};

/// Derivative data at the midpoint of a triangle
template<typename VectorX>
struct TriangleMidpointData
{
  VectorX normal_derivative; // derivative in the direction of the opposite corner

  /// Default constructor
  ///
  /// @param[in] input_normal_derivative: derivative in the direction of the opposite corner
  TriangleMidpointData(const VectorX& input_normal_derivative)
  {
    normal_derivative = input_normal_derivative;
  }

  /// Trivial constructor
  TriangleMidpointData() { normal_derivative = VectorX(); }
};

// Typedef for spatial vector case
typedef TriangleCornerData<SpatialVector> TriangleCornerFunctionData;
typedef TriangleMidpointData<SpatialVector> TriangleMidpointFunctionData;

/// Generate corner position data for a mesh with an affine manifold structure and
/// per vertex position and gradients.
///
/// @param[in] V: mesh vertex embedding
/// @param[in] affine_manifold: mesh topology and affine manifold structure
/// @param[in] gradients: per vertex uv gradients in the local charts
/// @param[out] corner_data: quadratic vertex position and derivative data
void
generate_affine_manifold_corner_data(
  const Eigen::MatrixXd& V,
  const AffineManifold& affine_manifold,
  const std::vector<Matrix2x3r>& gradients,
  std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data);

/// Generate midpoint position data for a mesh with an affine manifold structure and
/// per edge gradients.
///
/// @param[in] affine_manifold: mesh topology and affine manifold structure
/// @param[in] edge_gradients: per edge gradients in the local charts
/// @param[out] midpoint_data: quadratic edge midpoint derivative data
void
generate_affine_manifold_midpoint_data(
  const AffineManifold& affine_manifold,
  const std::vector<std::array<Matrix2x3r, 3>>& edge_gradients,
  std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data);

/// Given corner data for the endpoints of the edge, compute the midpoint
/// and the edge aligned midpoint gradient of the corresponding Powell-Sabin
/// quadratic spline patch
///
/// @param[in] edge_origin_corner_data: corner data at the origin corner of the
/// oriented edge
/// @param[in] edge_dest_corner_data: corner data at the destination corner of
/// the oriented edge
/// @param[out] midpoint: quadratic edge function midpoint
/// @param[out] midpoint_edge_gradient: edge aligned gradient
void
compute_edge_midpoint_with_gradient(
  const TriangleCornerFunctionData& edge_origin_corner_data,
  const TriangleCornerFunctionData& edge_dest_corner_data,
  SpatialVector& midpoint,
  SpatialVector& midpoint_edge_gradient);

/// Given per corner position data, rearrange it into matrices with row i
/// corresponding to the data for vertex i.
///
/// @param[in] corner_data: quadratic vertex position and derivative data
/// @param[out] position_matrix: matrix with position data as rows
/// @param[out] first_derivative_matrix: matrix with first derivative data as rows
/// @param[out] second_derivative_matrix: matrix with second derivative data as rows
void
generate_corner_data_matrices(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  MatrixXr& position_matrix,
  MatrixXr& first_derivative_matrix,
  MatrixXr& second_derivative_matrix);

/// Given position data, rearrange the edge midpoint data into matrices with row i
/// corresponding to the data for edge i.
///
/// Note that the normal derivative matrix is extracted directly from the midpoint data
/// and the position and tangent derivative matrices are inferred from the corner data.
///
/// @param[in] corner_data: quadratic vertex position and derivative data
/// @param[in] midpoint_data: quadratic edge midpoint derivative data
/// @param[out] position_matrix: matrix with midpoint position data as rows
/// @param[out] tangent_derivative_matrix: matrix with edge tangent derivative data as rows
/// @param[out] normal_derivative_matrix: matrix with normal derivative data as rows
void
generate_midpoint_data_matrices(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  const std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data,
  MatrixXr& position_matrix,
  MatrixXr& tangent_derivative_matrix,
  MatrixXr& normal_derivative_matrix);
