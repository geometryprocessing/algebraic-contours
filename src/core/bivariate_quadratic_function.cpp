// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "bivariate_quadratic_function.h"

#include <unsupported/Eigen/Polynomials>

Eigen::Matrix<double, 1, 3>
generate_linear_monomials(const PlanarPoint& domain_point)
{
  double u = domain_point[0];
  double v = domain_point[1];

  Eigen::Matrix<double, 1, 3> w;
  w << 1, u, v;
  return w;
}

Eigen::Matrix<double, 1, 6>
generate_quadratic_monomials(const PlanarPoint& domain_point)
{
  double u = domain_point[0];
  double v = domain_point[1];

  Eigen::Matrix<double, 1, 6> w;
  w << 1, u, v, u * v, u * u, v * v;
  return w;
}

double
evaluate_line(const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& line_coeffs,
              const PlanarPoint& domain_point)
{
  OneFormXr w = generate_linear_monomials(domain_point);
  return w * line_coeffs;
}

double
evaluate_quadratic(
  const Eigen::Ref<const Eigen::Matrix<double, 6, 1>>& quadratic_coeffs,
  const PlanarPoint& domain_point)
{
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(domain_point);
  return w * quadratic_coeffs;
}

Eigen::Matrix<double, 6, 1>
compute_linear_product(
  const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& L1_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& L2_coeffs)
{
  Eigen::Matrix<double, 6, 1> product_coeffs;

  product_coeffs(0) = L1_coeffs(0) * L2_coeffs(0);
  product_coeffs(1) = L1_coeffs(1) * L2_coeffs(0) + L1_coeffs(0) * L2_coeffs(1);
  product_coeffs(2) = L1_coeffs(2) * L2_coeffs(0) + L1_coeffs(0) * L2_coeffs(2);
  product_coeffs(3) = L1_coeffs(1) * L2_coeffs(2) + L1_coeffs(2) * L2_coeffs(1);
  product_coeffs(4) = L1_coeffs(1) * L2_coeffs(1);
  product_coeffs(5) = L1_coeffs(2) * L2_coeffs(2);

  return product_coeffs;
}

Matrix6x3r
compute_quadratic_cross_product(const Eigen::Matrix<double, 3, 3>& V_coeffs,
                                const Eigen::Matrix<double, 3, 3>& W_coeffs)
{
  Matrix6x3r N_coeffs;

  // 1 coefficient
  N_coeffs.row(0) = cross_product<double>(V_coeffs.row(0), W_coeffs.row(0));

  // u coefficient
  N_coeffs.row(1) = cross_product<double>(V_coeffs.row(0), W_coeffs.row(1)) +
                    cross_product<double>(V_coeffs.row(1), W_coeffs.row(0));

  // v coefficient
  N_coeffs.row(2) = cross_product<double>(V_coeffs.row(0), W_coeffs.row(2)) +
                    cross_product<double>(V_coeffs.row(2), W_coeffs.row(0));

  // uv coefficient
  N_coeffs.row(3) = cross_product<double>(V_coeffs.row(1), W_coeffs.row(2)) +
                    cross_product<double>(V_coeffs.row(2), W_coeffs.row(1));

  // u^2 coefficient
  N_coeffs.row(4) = cross_product<double>(V_coeffs.row(1), W_coeffs.row(1));

  // v^2 coefficient
  N_coeffs.row(5) = cross_product<double>(V_coeffs.row(2), W_coeffs.row(2));

  return N_coeffs;
}

// FIXME Make static
Eigen::Matrix<double, 3, 6>
u_derivative_matrix()
{
  Eigen::Matrix<double, 3, 6> D_u;
  D_u.setZero();

  // Set nonzero elements explicitly
  D_u(0, 1) = 1;
  D_u(1, 4) = 2;
  D_u(2, 3) = 1;

  return D_u;
}

// FIXME Make static
Eigen::Matrix<double, 3, 6>
v_derivative_matrix()
{
  Eigen::Matrix<double, 3, 6> D_v;
  D_v.setZero();

  // Set nonzero elements explicitly
  D_v(0, 2) = 1;
  D_v(1, 3) = 1;
  D_v(2, 5) = 2;

  return D_v;
}

void
generate_bezier_to_monomial_matrix(
  Eigen::Matrix<double, 6, 6>& change_of_basis_matrix)
{
  change_of_basis_matrix << 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, -2, 0, 2, 0, 0, 0,
    -2, 2, -2, -2, 0, 0, 2, 0, 0, -2, 1, 0, 1, 0, -2, 0, 0, 1, 1;
};

void
generate_monomial_to_bezier_matrix(
  Eigen::Matrix<double, 6, 6>& change_of_basis_matrix)
{
  change_of_basis_matrix << 1, 0.5, 0.5, 0.5, 0, 0, 1, 0, 0.5, 0, 0, 0, 1, 0.5,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0;
};

/// FIXME This is not fully general
bool
is_conic_standard_form(const VectorXr& C_coeffs)
{
  // Mixed term must be zero
  if (!float_equal(C_coeffs(3), 0.0))
    return false;

  return true;
}
