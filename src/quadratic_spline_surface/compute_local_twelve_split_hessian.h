// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "PS12_patch_coeffs.h"


/// Compute the matrix to convert from global degrees of freedom
/// [f0, gu_0, gv_0, f1, gu_1, gv_1, f2, gu_2, gv_2, gm_12, gm_20, gm_01]
/// to the local triangle degrees of freedom
/// [f0, f1, f2, d01, d02, d10, d12, d20, d21, h01, h12, h20]
///
/// Here, dij is the derivative of f in the direction of edge eij and
/// hij is the derivative of f from the edge midpoint to the opposite vertex
///
/// @param[in] uv: global uv vertex positions
/// @param[in] corner_to_corner_uv_positions: per vertex matrices with the
/// local vertex chart edge directions as rows
/// @param[in] reverse_edge_orientations: per edge booleans indicating if the
/// edge midpoint needs to be flipped from frame orientation consistency
/// @return matrix mapping global to local degrees of freedom
Eigen::Matrix<double, 12, 12>
get_C_gl(const Eigen::Matrix<double, 3, 2>& uv,
         const std::array<Matrix2x2r, 3>& corner_to_corner_uv_positions,
         const std::array<bool, 3>& reverse_edge_orientations);

/// Compute the matrix to go from local triangle degrees of freedom
/// [f0, f1, f2, d01, d02, d10, d12, d20, d21, h01, h12, h20]
/// to the second derivatives
/// [paa, pbb, pab]_{1,...,12}
/// of the twelve subtriangle quadratic functions defined by the Powell-
/// Sabin interpolant with respect to the barycentric coordinates of
/// the domain triangle
///
/// The full chain of theoretical operations is:
///     Csub: local triangle dof -> Bezier points
///     Cder: Bezier points -> second derivatives in subtriangle coords
///     Tquad: second derivatives in subtriangle coords
///            -> second derivatives in triangle coords
/// However, this matrix is constant and is thus hard coded
/// 
/// @return matrix mapping local dof to second derivatives 
Eigen::Matrix<double, 36, 12>
get_Tquad_Cder_Csub();

/// Compute the matrix to go from second derivatives 
/// [paa, pbb, pab]_{1,...,12}
/// of the twelve subtriangle quadratic functions defined by the Powell-
/// Sabin interpolant with respect to the barycentric coordinates of
/// the domain triangle to the second derivatives with respect to a 
/// parametric domain triangle defined by uv
///
/// @param[in] uv: uv coordinates of the domain triangle
/// @return matrix mapping second derivatives wrt barycentric coordinates
///     to second derivatives wrt uv coordinates
Eigen::Matrix<double, 36, 36>
get_R_quad(const Eigen::Matrix<double, 3, 2>& uv);

/// Compute the diagonal uv triangle area weighting matrix for the Hessian
///
/// @param[in] uv: uv coordinates of the domain triangle
/// @return diagonal Hessian weight matrix
Eigen::Matrix<double, 36, 36>
get_S(const Eigen::Matrix<double, 3, 2>& uv);

/// Compute the thin plate energy Hessian matrix for a single 12-split 
/// Powell-Sabin element
///
/// @param[in] uv: uv coordinates of the domain triangle
/// @param[in] corner_to_corner_uv_positions: per vertex matrices with the
/// local vertex chart edge directions as rows
/// @param[in] reverse_edge_orientations: per edge booleans indicating if the
/// edge midpoint needs to be flipped from frame orientation consistency
/// @return Hessian matrix in terms of Powell-Sabin degrees of freedom
Eigen::Matrix<double, 12, 12>
build_local_smoothness_hessian(
  const Eigen::Matrix<double, 3, 2>& uv,
  const std::array<Matrix2x2r, 3>& corner_to_corner_uv_positions,
  const std::array<bool, 3>& reverse_edge_orientations);