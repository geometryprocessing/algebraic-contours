// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "polynomial_function.h"

template<int dimension>
double&
const_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
            int col)
{
  assert(col > dimension);
  return quadratic_coeffs(0, col);
}

template<int dimension>
double&
u_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
        int col)
{
  assert(col >= 0);
  assert(col < dimension);
  return quadratic_coeffs(1, col);
}

template<int dimension>
double&
v_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
        int col)
{
  assert(col >= 0);
  assert(col < dimension);
  return quadratic_coeffs(2, col);
}

template<int dimension>
double&
uv_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
         int col)
{
  assert(col >= 0);
  assert(col < dimension);
  return quadratic_coeffs(3, col);
}

template<int dimension>
double&
uu_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
         int col)
{
  assert(col >= 0);
  assert(col < dimension);
  return quadratic_coeffs(4, col);
}

template<int dimension>
double&
vv_coeff(Eigen::Ref<Eigen::Matrix<double, 6, dimension>> quadratic_coeffs,
         int col)
{
  assert(col >= 0);
  assert(col < dimension);
  return quadratic_coeffs(5, col);
}

/// Generate monomial variable terms for a quadratic in order [1, u, v, uv, uu,
/// vv]
///
/// @param[in] domain_point: uv coordinates to generate the monomials for
/// @return quadratic monomials
Eigen::Matrix<double, 1, 6>
generate_quadratic_monomials(const PlanarPoint& domain_point);

/// Generate monomial variable terms for a line in order [1, u, v]
///
/// @param[in] domain_point: uv coordinates to generate the monomials for
/// @return linear monomials
Eigen::Matrix<double, 1, 3>
generate_linear_monomials(const PlanarPoint& domain_point);

/// Evaluate a quadratic bivariate equation with scalar coefficients.
///
/// @param[in] quadratic_coeffs: quadratic coefficients in order [1, u, v, uv,
/// uu, vv]
/// @param[in] domain_point: uv coordinates to evaluate the quadratic at
/// @param[out] quadratic_evaluation: quadratic function evaluation
template<int dimension>
void
evaluate_quadratic_mapping(
  const Eigen::Ref<const Eigen::Matrix<double, 6, dimension>>& quadratic_coeffs,
  const PlanarPoint& domain_point,
  Eigen::Ref<Eigen::Matrix<double, 1, dimension>> quadratic_evaluation)
{
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(domain_point);
  quadratic_evaluation = w * quadratic_coeffs;
}

/// Evaluate a quadratic bivariate equation with scalar coefficients.
///
/// @param[in] quadratic_coeffs: quadratic coefficients in order [1, u, v, uv,
/// uu, vv]
/// @param[in] domain_point: uv coordinates to evaluate the quadratic at
/// @return quadratic function evaluation
double
evaluate_quadratic(
  const Eigen::Ref<const Eigen::Matrix<double, 6, 1>>& quadratic_coeffs,
  const PlanarPoint& domain_point);

/// Evaluate a linear bivariate equation with scalar coefficients.
///
/// @param[in] line_coeffs: line coefficients in order [1, u, v]
/// @param[in] domain_point: uv coordinates to evaluate the line at
/// @return linear function evaluation
double
evaluate_line(const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& line_coeffs,
              const PlanarPoint& domain_point);

/// Compute the quadratic coefficients for the scalar product of two linear
/// scalar functions V(u,v) and W(u,v) with coefficients in order [1, u, v]
///
/// @param[in] V_coeffs: coefficients for the first linear vector function
/// @param[in] W_coeffs: coefficients for the second linear vector function
/// @return coefficients for the quadratic product function
Eigen::Matrix<double, 6, 1>
compute_linear_product(
  const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& L1_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, 3, 1>>& L2_coeffs);

/// Compute the quadratic coefficients for the cross product of two linear row
/// vector functions V(u,v) and W(u,v) with coefficients in order [1, u, v]
///
/// @param[in] V_coeffs: coefficients for the first linear vector function
/// @param[in] W_coeffs: coefficients for the second linear vector function
/// @return coefficients for the quadratic cross function
Matrix6x3r
compute_quadratic_cross_product(const Eigen::Matrix<double, 3, 3>& V_coeffs,
                                const Eigen::Matrix<double, 3, 3>& W_coeffs);

/// Build matrix from quadratic coefficients to linear coefficients representing
/// the derivative in the u direction
///
/// @return u derivative matrix
Eigen::Matrix<double, 3, 6>
u_derivative_matrix();

/// Build matrix from quadratic coefficients to linear coefficients representing
/// the derivative in the v direction
///
/// @return v derivative matrix
Eigen::Matrix<double, 3, 6>
v_derivative_matrix();

/// Generate the matrix to go from Bezier control points to quadratic
/// coefficients over the standard u + v <= 1 triangle in the positive quadrant.
///
/// @param[out] change_of_basis_matrix: matrix going from bezier points to
/// monomial coefficients
void
generate_bezier_to_monomial_matrix(
  Eigen::Matrix<double, 6, 6>& change_of_basis_matrix);

/// Generate the matrix to go from quadratic coefficients over the standard u +
/// v <= 1 triangle in the positive quadrant to Bezier control points.
///
/// @param[out] change_of_basis_matrix: matrix going from monomial coefficients
/// to bezier points
void
generate_monomial_to_bezier_matrix(
  Eigen::Matrix<double, 6, 6>& change_of_basis_matrix);

/// Return true iff the conic with quadratic coefficients C_coeffs is in
/// standard form with no mixed terms
///
/// @param[in] C_coeffs: quadratic coefficients for the conic
/// @return: true iff C_coeffs is a conic in standard form
bool
is_conic_standard_form(const VectorXr& C_coeffs);

/// Generate a human readable format of a quadratic mapping
///
/// @param[in] quadratic_coeffs: quadratic coefficients in order [1, u, v, uv,
/// uu, vv]
/// @return formatted quadratic mapping
template<int dimension>
std::string
formatted_bivariate_quadratic_mapping(
  const Eigen::Ref<const Eigen::Matrix<double, 6, dimension>>& quadratic_coeffs,
  size_t precision = 16)
{
  std::stringstream quadratic_string;
  for (int i = 0; i < quadratic_coeffs.cols(); ++i) {
    quadratic_string << std::fixed << std::setprecision(precision)
                     << quadratic_coeffs(0, i);
    quadratic_string << formatted_term(quadratic_coeffs(1, i), "u", precision);
    quadratic_string << formatted_term(quadratic_coeffs(2, i), "v", precision);
    quadratic_string << formatted_term(quadratic_coeffs(3, i), "uv", precision);
    quadratic_string << formatted_term(
      quadratic_coeffs(4, i), "u^2", precision);
    quadratic_string << formatted_term(
      quadratic_coeffs(5, i), "v^2", precision);
    quadratic_string << "\n";
  }

  return quadratic_string.str();
}

/// Generate a human readable format of a linear mapping
///
/// @param[in] line_coeffs: linear coefficients in order [1, u, v]
/// @return formatted linear mapping
template<int dimension>
std::string
formatted_bivariate_linear_mapping(
  const Eigen::Ref<const Eigen::Matrix<double, 3, dimension>>& line_coeffs,
  size_t precision = 16)
{
  std::stringstream line_string;
  for (int i = 0; i < dimension; ++i) {
    line_string << std::fixed << std::setprecision(precision)
                << line_coeffs(0, i);
    line_string << formatted_term(line_coeffs(1, i), "u", precision);
    line_string << formatted_term(line_coeffs(2, i), "v", precision);
    line_string << "\n";
  }

  return line_string.str();
}

/// Given an affine transformation [u, v]^T = A*[u', v']^T + b of R^2, generate
/// the change of basis matrix C for the bivariate quadratic monomial
/// coefficients vector Q with respect to u, v so that Q' = C * Q is the
/// coefficient vector for the bivariate quadratic monomials with respect to u',
/// v'.
///
/// @param[in] linear_transformation: linear part of the affine transformation
/// @param[in] translation: translation part of the affine transformation
/// @param[out] change_of_basis_matrix: change of coefficient basis matrix
template<typename Scalar>
void
generate_quadratic_coordinate_affine_transformation_matrix(
  const Eigen::Matrix<Scalar, 2, 2>& linear_transformation,
  const Eigen::Matrix<Scalar, 1, 2>& translation,
  Eigen::Matrix<Scalar, 6, 6>& change_of_basis_matrix)
{
  // Get matrix information
  Scalar b11 = linear_transformation(0, 0);
  Scalar b12 = linear_transformation(0, 1);
  Scalar b21 = linear_transformation(1, 0);
  Scalar b22 = linear_transformation(1, 1);
  Scalar b1 = translation[0];
  Scalar b2 = translation[1];
  Scalar one(1);
  Scalar zero(0);

  // Set matrix
  change_of_basis_matrix << one, b1, b2, b1 * b2, b1 * b1, b2 * b2, zero, b11,
    b12, b11 * b2 + b12 * b1, 2 * b1 * b11, 2 * b2 * b12, zero, b21, b22,
    b22 * b1 + b21 * b2, 2 * b1 * b21, 2 * b2 * b22, zero, zero, zero,
    b11 * b22 + b21 * b12, 2 * b11 * b21, 2 * b12 * b22, zero, zero, zero,
    b11 * b12, b11 * b11, b12 * b12, zero, zero, zero, b21 * b22, b21 * b21,
    b22 * b22;
}

/// Generate the matrix to transform quadratic monomial coefficients for
/// coordinates (u, v) to coefficients for translated coordinates (u', v') = (u
/// + du, v + dv).
///
/// @param[in] du: change in u coordinate
/// @param[in] dv: change in v coordinate
/// @param[out] change_of_basis_matrix: change of coefficient basis matrix
template<typename Scalar>
void
generate_quadratic_coordinate_translation_matrix(
  double du,
  double dv,
  Eigen::Matrix<Scalar, 6, 6>& change_of_basis_matrix)
{
  // Build (inverse) translation and linear identity matrix
  Eigen::Matrix<Scalar, 2, 1> translation = { -du, -dv };
  Eigen::Matrix<Scalar, 2, 2> identity;
  identity << Scalar(1), Scalar(0), Scalar(0), Scalar(1);

  // Generate the translation transformation as a special case of an affine
  // transformation
  generate_quadratic_coordinate_affine_transformation_matrix<Scalar>(
    identity, translation, change_of_basis_matrix);
}

/// Given an barycentric transformation [u, v, w]^T = A*[u', v', w']^T of RP^2
/// with coordinates normalized so that u + v + w = u' + v' + w' = 1, generate
/// the change of basis matrix C for the bivariate quadratic monomial
/// coefficients vector Q with respect to u, v so that Q' = C * Q is the
/// coefficient vector for the bivariate quadratic monomials with respect to u',
/// v'.
///
/// @param[in] barycentric_transformation: transformation of the coordinates
/// @param[out] change_of_basis_matrix: change of coefficient basis matrix
template<typename Scalar>
void
generate_quadratic_coordinate_barycentric_transformation_matrix(
  const Eigen::Matrix<Scalar, 3, 3>& barycentric_transformation,
  Eigen::Matrix<Scalar, 6, 6>& change_of_basis_matrix)
{
  Eigen::Matrix<Scalar, 2, 2> linear_transformation;
  Eigen::Matrix<Scalar, 2, 1> translation;

  // Get relevant barycentric matrix components. Since the w coordinate is
  // discarded, we only use the first two rows of the barycentric transformation
  Scalar a11 = barycentric_transformation(0, 0);
  Scalar a12 = barycentric_transformation(0, 1);
  Scalar a13 = barycentric_transformation(0, 2);
  Scalar a21 = barycentric_transformation(1, 0);
  Scalar a22 = barycentric_transformation(1, 1);
  Scalar a23 = barycentric_transformation(1, 2);

  // Build affine transformation from barycentric transformation with w = 1 - u
  // - v
  linear_transformation << a11 - a13, a21 - a23, a12 - a13, a22 - a23;
  translation << a13, a23;

  // Get change of basis matrix from the affine transformation
  generate_quadratic_coordinate_affine_transformation_matrix<Scalar>(
    linear_transformation, translation, change_of_basis_matrix);
}

/// Given vertex positions in R^2 for a domain triangle, generate the change
/// of basis matrix C for the bivariate quadratic monomial coefficients vector
/// Q with respect to u, v so that Q' = C * Q is the coefficient vector for the
/// surface mapping over the triangle in the positive quadrant with u + v <= 1
/// that has the same image as the surface mapping over the input domain.
///
/// @param[in] v0: first domain triangle vertex position
/// @param[in] v1: second domain triangle vertex position
/// @param[in] v2: third domain triangle vertex position
/// @param[out] change_of_basis_matrix: change of coefficient basis matrix
template<typename Scalar>
void
generate_quadratic_coordinate_domain_triangle_normalization_matrix(
  const Eigen::Matrix<Scalar, 2, 1>& v0,
  const Eigen::Matrix<Scalar, 2, 1>& v1,
  const Eigen::Matrix<Scalar, 2, 1>& v2,
  Eigen::Matrix<Scalar, 6, 6>& change_of_basis_matrix)
{
  // Generate affine transformation mapping the standard triangle to the new
  // triangle
  Eigen::Matrix<Scalar, 2, 2> linear_transformation;
  Eigen::Matrix<Scalar, 2, 1> translation;
  linear_transformation.row(0) = v1 - v0;
  linear_transformation.row(1) = v2 - v0;
  translation = v0;

  // Get change of basis matrix from the affine transformation
  generate_quadratic_coordinate_affine_transformation_matrix<Scalar>(
    linear_transformation, translation, change_of_basis_matrix);
}

/// Deprecated
///
/// TODO: This is incomplete and is only correct up to translation. This is fine
/// for the quadratic surface energy, but not for the general mapping needed to
/// build boundary boxes. It can be replaced by
/// generate_quadratic_coordinate_domain_triangle_normalization_matrix.
template<typename Scalar>
void
generate_reparametrization(
  const Eigen::Matrix<Scalar, 2, 1>& v0,
  const Eigen::Matrix<Scalar, 2, 1>& v1,
  const Eigen::Matrix<Scalar, 2, 1>& v2,
  Eigen::Matrix<Scalar, 2, 2>& general_to_canonical_reparametrization,
  Eigen::Matrix<Scalar, 2, 2>& canonical_to_general_reparametrization)
{
  Scalar l0 = (v1 - v0).norm();
  ;
  Scalar l1 = (v2 - v1).norm();
  ;
  Scalar l2 = (v0 - v2).norm();
  ;
  Scalar cos_angle = (l0 * l0 + l2 * l2 - l1 * l1) / (2 * l0 * l2);
  Scalar sin_angle = sqrt(1 - cos_angle * cos_angle);
  general_to_canonical_reparametrization << Scalar(1) / l0, Scalar(0),
    -cos_angle / (l0 * sin_angle), Scalar(1) / (l2 * sin_angle);
  canonical_to_general_reparametrization << l0, l2 * cos_angle, Scalar(0),
    l2 * sin_angle;
}
