// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "evaluate_surface_normal.h"

#include "bivariate_quadratic_function.h"
#include "polynomial_function.h"

void
generate_quadratic_surface_normal_coeffs(
  const Matrix6x3r& surface_mapping_coeffs,
  Matrix6x3r& normal_mapping_coeffs)
{
  // Get directional derivatives
  Eigen::Matrix<double, 3, 6> D_u = u_derivative_matrix();
  Eigen::Matrix<double, 3, 6> D_v = v_derivative_matrix();
  Eigen::Matrix<double, 3, 3> u_derivative_coeffs =
    D_u * surface_mapping_coeffs;
  Eigen::Matrix<double, 3, 3> v_derivative_coeffs =
    D_v * surface_mapping_coeffs;

  // Compute normal from the cross product
  normal_mapping_coeffs =
    compute_quadratic_cross_product(u_derivative_coeffs, v_derivative_coeffs);
}
