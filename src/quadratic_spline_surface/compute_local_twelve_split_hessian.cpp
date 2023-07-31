// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_local_twelve_split_hessian.h"

Eigen::Matrix<double, 12, 12>
get_C_gl(const Eigen::Matrix<double, 3, 2>& uv,
         const std::array<Matrix2x2r, 3>& corner_to_corner_uv_positions,
         const std::array<bool, 3>& reverse_edge_orientations)
{
  // Get global uv vertex positions
  Eigen::Vector2d q0 = uv.row(0);
  Eigen::Vector2d q1 = uv.row(1);
  Eigen::Vector2d q2 = uv.row(2);

  // Compute global edge directions and squared lengths
  Eigen::Vector2d e01_global = q1 - q0;
  Eigen::Vector2d e12_global = q2 - q1;
  Eigen::Vector2d e20_global = q0 - q2;
  double l01sq = e01_global.dot(e01_global);
  double l12sq = e12_global.dot(e12_global);
  double l20sq = e20_global.dot(e20_global);

  // Compute global midpoint to corner directions
  Eigen::Vector2d e01_m = q2 - (q0 + q1) / 2;
  Eigen::Vector2d e12_m = q0 - (q2 + q1) / 2;
  Eigen::Vector2d e20_m = q1 - (q0 + q2) / 2;

  // Compute global edge perpendicular directions
  Eigen::Vector2d perpe01(-e01_global[1], e01_global[0]);
  Eigen::Vector2d perpe12(-e12_global[1], e12_global[0]);
  Eigen::Vector2d perpe20(-e20_global[1], e20_global[0]);

  // Compute local edge midpoint directions in frames defined by eij, perpeij
  Eigen::Vector2d e01_m_loc(e01_m.dot(e01_global) / l01sq,
                            e01_m.dot(perpe01) / l01sq);
  Eigen::Vector2d e12_m_loc(e12_m.dot(e12_global) / l12sq,
                            e12_m.dot(perpe12) / l12sq);
  Eigen::Vector2d e20_m_loc(e20_m.dot(e20_global) / l20sq,
                            e20_m.dot(perpe20) / l20sq);

  // Extract vertex chart rotated corner to corner uv directions
  // Note that eij = R(qj - qi), where R is some rigid transformation mapping
  // the global uv triangle to some local layout containing the given corner
  Eigen::Vector2d e01 = corner_to_corner_uv_positions[0].row(0); // q1 - q0
  Eigen::Vector2d e02 = corner_to_corner_uv_positions[0].row(1); // q2 - q0
  Eigen::Vector2d e12 = corner_to_corner_uv_positions[1].row(0); // q2 - q1
  Eigen::Vector2d e10 = corner_to_corner_uv_positions[1].row(1); // q0 - q1
  Eigen::Vector2d e20 = corner_to_corner_uv_positions[2].row(0); // q0 - q2
  Eigen::Vector2d e21 = corner_to_corner_uv_positions[2].row(1); // q1 - q2

  // Assign labels to the elementary row basis vectors in R^12 corresponding
  // to global degree of freedom
  Eigen::Matrix<double, 12, 1> f0, f1, f2, gu_0, gv_0, gu_1, gv_1, gu_2, gv_2,
    gm_12, gm_20, gm_01, d01, d02, d10, d12, d20, d21, ge_01, ge_12, ge_20, h01,
    h12, h20;
  f0    << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  gu_0  << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  gv_0  << 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  f1    << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
  gu_1  << 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
  gv_1  << 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
  f2    << 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
  gu_2  << 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
  gv_2  << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
  gm_12 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
  gm_20 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
  gm_01 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  // Compute the derivative of f in the directions of the parametric edges eij
  // Note that these (and subsequent quantities) are actually row vectors that extract
  // the desired values from the column vector of global degrees of freedom by multiplication
  d01 = e01(0) * gu_0 + e01(1) * gv_0;
  d02 = e02(0) * gu_0 + e02(1) * gv_0;
  d10 = e10(0) * gu_1 + e10(1) * gv_1;
  d12 = e12(0) * gu_1 + e12(1) * gv_1;
  d20 = e20(0) * gu_2 + e20(1) * gv_2;
  d21 = e21(0) * gu_2 + e21(1) * gv_2;

  // Compute the derivative along eij at the edge midpoint
  ge_01 = -2 * f0 - 0.5 * d01 + 2 * f1 + 0.5 * d10;
  ge_12 = -2 * f1 - 0.5 * d12 + 2 * f2 + 0.5 * d21;
  ge_20 = -2 * f2 - 0.5 * d20 + 2 * f0 + 0.5 * d02;

  // Reverse (global) edge orientations as needed for consistency of the 
  // midpoint reference frames
  if (reverse_edge_orientations[0])
    gm_12 = -gm_12;
  if (reverse_edge_orientations[1])
    gm_20 = -gm_20;
  if (reverse_edge_orientations[2])
    gm_01 = -gm_01;

  // Compute the derivative of f in the direction of the edge midpoint to opposite
  // vertex
  h01 = e01_m_loc[0] * ge_01 + e01_m_loc[1] * gm_01;
  h12 = e12_m_loc[0] * ge_12 + e12_m_loc[1] * gm_12;
  h20 = e20_m_loc[0] * ge_20 + e20_m_loc[1] * gm_20;

  // Assemble matrix
  Eigen::Matrix<double, 12, 12> C_gl;
  C_gl.row(0) = f0;
  C_gl.row(1) = f1;
  C_gl.row(2) = f2;
  C_gl.row(3) = d01;
  C_gl.row(4) = d02;
  C_gl.row(5) = d10;
  C_gl.row(6) = d12;
  C_gl.row(7) = d20;
  C_gl.row(8) = d21;
  C_gl.row(9) = h01;
  C_gl.row(10) = h12;
  C_gl.row(11) = h20;

  return C_gl;
}

Eigen::Matrix<double, 36, 12>
get_Tquad_Cder_Csub()
{
  // Autogenerate matrix entries as 12 subtriangle matrices
  // Note that unused values (first derivatives and constants) are also generated
  // but are unused
  double patch_coeffs[12][6][12];
  PS12_patch_coeffs(patch_coeffs);

  // Combine 12 subtriangle matrices, flattening the second derivatives per subtriangle
  Eigen::Matrix<double, 36, 12> Tquad_Cder_Csub;
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 12; k++) {
        Tquad_Cder_Csub(i * 3 + j, k) = patch_coeffs[i][j + 3][k];
      }
    }
  }

  return Tquad_Cder_Csub;
}

Eigen::Matrix<double, 36, 36>
get_R_quad(const Eigen::Matrix<double, 3, 2>& uv)
{
  // Compute 2x2 matrix mapping parametric to barycentric coordinates
  Eigen::Matrix2d R_inv;
  R_inv << uv(0, 0) - uv(2, 0), uv(1, 0) - uv(2, 0), //
    uv(0, 1) - uv(2, 1), uv(1, 1) - uv(2, 1);

  // Compute inverse 2x2 matrix mapping barycentric to parametric coordinates
  Eigen::Matrix2d R;
  R << R_inv(1, 1), -R_inv(0, 1), //
    -R_inv(1, 0), R_inv(0, 0);
  R /= (R_inv(0, 0) * R_inv(1, 1) - R_inv(0, 1) * R_inv(1, 0));

  // Compute 3x3 matrix q_l mapping second derivatives wrt barycentric coordinates
  // to second derivatives wrt parametric coordinates
  Eigen::Matrix3d q_l;
  q_l << R(0, 0) * R(1, 1) + R(0, 1) * R(1, 0), 2 * R(0, 0) * R(0, 1),
    2 * R(1, 0) * R(1, 1),                                               //
    2 * R(0, 0) * R(1, 0), 2 * R(0, 0) * R(0, 0), 2 * R(1, 0) * R(1, 0), //
    2 * R(0, 1) * R(1, 1), 2 * R(0, 1) * R(0, 1), 2 * R(1, 1) * R(1, 1);

  // Build R_quad as 12 block copies of q_l
  Eigen::Matrix<double, 36, 36> R_quad;
  R_quad.fill(0);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        R_quad(i * 3 + j, i * 3 + k) = q_l(j, k);
      }
    }
  }
  return R_quad;
}

// Get barycentric patch triangle area matrix weighted by some arbitrary area A
Eigen::Matrix<double, 36, 36>
get_S_weighted(double A)
{
  Eigen::Matrix<double, 36, 36> S;
  S.fill(0);
  for (int i = 0; i < 36; i++) {
    if (i % 3 == 0) {
      S(i, i) = 2 * A;
    } else {
      S(i, i) = A;
    }
  }

  for (int i = 0; i < 18; i++) {
    S(i, i) /= 24;
  }
  for (int i = 18; i < 36; i++) {
    S(i, i) /= 8;
  }

  return S;
}

Eigen::Matrix<double, 36, 36>
get_S(const Eigen::Matrix<double, 3, 2>& uv)
{
  // Compute the area of the uv triangle
  Eigen::Matrix3d tri;
  tri << uv(0, 0), uv(0, 1), 1, //
    uv(1, 0), uv(1, 1), 1,      //
    uv(2, 0), uv(2, 1), 1;
  double A = tri.determinant() / 2;

  // Compute the corresponding Hessian weighting matrix
  return get_S_weighted(A * A);
}

Eigen::Matrix<double, 12, 12>
build_local_smoothness_hessian(
  const Eigen::Matrix<double, 3, 2>& uv,
  const std::array<Matrix2x2r, 3>& corner_to_corner_uv_positions,
  const std::array<bool, 3>& reverse_edge_orientations)
{
  // Build all elementary matrices composing the Hessian
  Eigen::Matrix<double, 12, 12> C_gl =
    get_C_gl(uv, corner_to_corner_uv_positions, reverse_edge_orientations);
  Eigen::Matrix<double, 36, 12> Tquad_Cder_Csub = get_Tquad_Cder_Csub();
  Eigen::Matrix<double, 36, 36> R_quad = get_R_quad(uv);
  Eigen::Matrix<double, 36, 36> S = get_S(uv);

  // Assemble matrices into the Hessian
  Eigen::Matrix<double, 36, 12> G = R_quad * Tquad_Cder_Csub * C_gl;
  Eigen::Matrix<double, 12, 12> hessian = G.transpose() * S * G;

  return hessian;
}