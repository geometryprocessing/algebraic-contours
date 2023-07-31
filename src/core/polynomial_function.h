// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include <unsupported/Eigen/Polynomials>

/// @file polynomial_function.h
///
/// Methods to construct and evaluate polynomials in one variable. We represent
/// polynomials as column vectors of coefficients with basis order 1, t,
/// t^2,.... For vector valued polynomial functions (i.e. functions f: R -> R^n
/// where each coordinate projection is a polynomial), we use matrices with the
/// column vectors as columns. We refer to these vector valued polynomial
/// functions as polynomial mappings to distinguish them from scalar valued
/// polynomials. Evaluation of the polynomials and polynomial mappings is then
/// simply the left product of the row vector [1 t t^2 ... t^n] with the column
/// vector or matrix.

/// @brief Generate the row vector T of n + 1 monomials 1, t, ... , t^n.
///
/// @tparam degree: maximum monomial degree
/// @param[in] t: evaluation point for the monomials
/// @param[out] T: row vector of monomials
template<int degree>
void
generate_monomials(double t, Eigen::Ref<Eigen::Matrix<double, 1, degree + 1>> T)
{
  T[0] = 1.0;
  for (int i = 1; i <= degree; ++i) {
    T[i] = t * T[i - 1];
  }
}

/// @brief Evaluate the polynomial with given coefficients at t.
///
/// @tparam degree: maximum monomial degree
/// @param[in] polynomial_coeffs: coefficients of the polynomial
/// @param[in] t: evaluation point for the polynomial
/// @return evaluation of the polynomial
template<int degree>
double
evaluate_polynomial(
  const Eigen::Ref<const Eigen::Matrix<double, degree + 1, 1>>&
    polynomial_coeffs,
  double t)
{
  Eigen::Matrix<double, 1, degree + 1> T;
  generate_monomials<degree>(t, T);
  return T * polynomial_coeffs;
}

/// @brief Evaluate the polynomial with given coefficients at t.
///
/// @tparam dimension: polynomial dimension
/// @tparam degree: maximum monomial degree
/// @param[in] polynomial_coeffs: coefficients of the polynomial
/// @param[in] t: evaluation point for the polynomial
/// @param[out] polynomial_evaluation: evaluation of the polynomial
template<int degree, int dimension>
void
evaluate_polynomial_mapping(
  const Eigen::Ref<const Eigen::Matrix<double, degree + 1, dimension>>&
    polynomial_coeffs,
  double t,
  Eigen::Ref<Eigen::Matrix<double, 1, dimension>> polynomial_evaluation)
{
  Eigen::Matrix<double, 1, degree + 1> T;
  generate_monomials<degree>(t, T);
  polynomial_evaluation = T * polynomial_coeffs;
}

/// @brief Generate the polynomial coefficients for the kronecker product of two
/// polynomials of the same dimension.
///
/// @tparam dimension: polynomial dimension
/// @tparam first_degree: maximum monomial degree of the first polynomial
/// @tparam second_degree: maximum monomial degree of the second polynomial
/// @param[in] first_polynomial_coeffs: coefficients of the first polynomial
/// @param[in] second_polynomial_coeffs: coefficients of the second polynomial
/// @param[out] product_polynomial_coeffs: product polynomial coefficients
template<int first_degree, int second_degree, int dimension>
void
compute_polynomial_mapping_product(
  const Eigen::Ref<const Eigen::Matrix<double, first_degree + 1, dimension>>&
    first_polynomial_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, second_degree + 1, dimension>>&
    second_polynomial_coeffs,
  Eigen::Ref<Eigen::Matrix<double, first_degree + second_degree + 1, dimension>>
    product_polynomial_coeffs)
{
  // Compute the new polynomial coefficients by convolution
  product_polynomial_coeffs.setZero();
  for (int i = 0; i <= first_degree; ++i) {
    for (int j = 0; j <= second_degree; ++j) {
      for (int k = 0; k < dimension; ++k) {
        product_polynomial_coeffs(i + j, k) +=
          first_polynomial_coeffs(i, k) * second_polynomial_coeffs(j, k);
      }
    }
  }
}

/// @brief Generate the polynomial coefficients for the product of a
/// scalar polynomial and a vector valued polynomial mapping.
///
/// @tparam dimension: polynomial mapping dimension
/// @tparam first_degree: maximum monomial degree of the first polynomial
/// @tparam second_degree: maximum monomial degree of the second polynomial
/// @param[in] scalar_polynomial_coeffs: coefficients of the scalar polynomial
/// @param[in] polynomial_coeffs: coefficients of the vector valued polynomial
/// @param[out] product_polynomial_coeffs: product vector valued polynomial
/// mapping coefficients
template<int first_degree, int second_degree, int dimension>
void
compute_polynomial_mapping_scalar_product(
  const Eigen::Ref<const Eigen::Matrix<double, first_degree + 1, 1>>&
    scalar_polynomial_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, second_degree + 1, dimension>>&
    polynomial_coeffs,
  Eigen::Ref<Eigen::Matrix<double, first_degree + second_degree + 1, dimension>>
    product_polynomial_coeffs)
{
  // Compute the new polynomial mapping coefficients by convolution
  product_polynomial_coeffs.setZero();
  for (int i = 0; i <= first_degree; ++i) {
    for (int j = 0; j <= second_degree; ++j) {
      for (int k = 0; k < dimension; ++k) {
        product_polynomial_coeffs(i + j, k) +=
          scalar_polynomial_coeffs(i) * polynomial_coeffs(j, k);
      }
    }
  }
}

/// @brief Generate the polynomial coefficients for the cross product of two
/// vector valued polynomial mappings with range R^3.
///
/// @tparam first_degree: maximum monomial degree of the first polynomial
/// @tparam second_degree: maximum monomial degree of the second polynomial
/// @param[in] first_polynomial_coeffs: coefficients of the first polynomial
/// @param[in] second_polynomial_coeffs: coefficients of the second polynomial
/// @param[out] product_polynomial_coeffs: product vector valued polynomial
/// mapping coefficients
template<int first_degree, int second_degree>
void
compute_polynomial_mapping_cross_product(
  const Eigen::Ref<const Eigen::Matrix<double, first_degree + 1, 3>>&
    first_polynomial_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, second_degree + 1, 3>>&
    second_polynomial_coeffs,
  Eigen::Ref<Eigen::Matrix<double, first_degree + second_degree + 1, 3>>
    product_polynomial_coeffs)
{
  // Compute cross product terms
  Eigen::Matrix<double, first_degree + second_degree + 1, 1> A0B1, A0B2, A1B0,
    A1B2, A2B0, A2B1;
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(0), second_polynomial_coeffs.col(1), A0B1);
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(0), second_polynomial_coeffs.col(2), A0B2);
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(1), second_polynomial_coeffs.col(0), A1B0);
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(1), second_polynomial_coeffs.col(2), A1B2);
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(2), second_polynomial_coeffs.col(0), A2B0);
  compute_polynomial_mapping_product<first_degree, second_degree, 1>(
    first_polynomial_coeffs.col(2), second_polynomial_coeffs.col(1), A2B1);

  // Assemble the cross product from the terms
  product_polynomial_coeffs.col(0) = A1B2 - A2B1;
  product_polynomial_coeffs.col(1) = A2B0 - A0B2;
  product_polynomial_coeffs.col(2) = A0B1 - A1B0;
}

/// @brief Generate the polynomial coefficients for the cross product of two
/// vector valued polynomial mappings with the same range.
///
/// @tparam dimension: polynomial mapping dimension
/// @tparam first_degree: maximum monomial degree of the first polynomial
/// @tparam second_degree: maximum monomial degree of the second polynomial
/// @param[in] first_polynomial_coeffs: coefficients of the first polynomial
/// @param[in] second_polynomial_coeffs: coefficients of the second polynomial
/// @param[out] product_polynomial_coeffs: product vector valued polynomial
/// mapping coefficients
template<int first_degree, int second_degree, int dimension>
void
compute_polynomial_mapping_dot_product(
  const Eigen::Ref<const Eigen::Matrix<double, first_degree + 1, dimension>>&
    first_polynomial_coeffs,
  const Eigen::Ref<const Eigen::Matrix<double, second_degree + 1, dimension>>&
    second_polynomial_coeffs,
  Eigen::Ref<Eigen::Matrix<double, first_degree + second_degree + 1, 1>>
    product_polynomial_coeffs)
{
  product_polynomial_coeffs.setZero();
  for (int i = 0; i < dimension; ++i) {
    Eigen::Matrix<double, first_degree + second_degree + 1, 1>
      product_polynomial_term_coeffs;
    compute_polynomial_mapping_product<first_degree, second_degree, dimension>(
      first_polynomial_coeffs.col(i),
      second_polynomial_coeffs.col(i),
      product_polynomial_term_coeffs);
    product_polynomial_coeffs += product_polynomial_term_coeffs;
  }
}

/// @brief Generate the polynomial coefficients for the derivative of a
/// polynomial mapping.
///
/// @param[in] polynomial_coeffs: coefficients of the polynomial mapping
/// @param[out] derivative-polynomial_coeffs: derivative polynomial mapping
/// coefficients
template<int degree, int dimension>
void
compute_polynomial_mapping_derivative(
  const Eigen::Ref<const Eigen::Matrix<double, degree + 1, dimension>>&
    polynomial_coeffs,
  Eigen::Ref<Eigen::Matrix<double, degree, dimension>>
    derivative_polynomial_coeffs)
{
  for (int i = 1; i <= degree; ++i) {
    for (int j = 0; j < dimension; ++j) {
      derivative_polynomial_coeffs(i - 1, j) = i * polynomial_coeffs(i, j);
    }
  }
}

/// @brief Compute the real roots of a quadratic polynomial.
///
/// @param[in] quadratic_coeffs: coefficients of the polynomial
/// @param[out] solutions: real roots of the polynomial
/// @param[out] num_solutions: solution count
/// @param[in] eps: threshold for zero comparisons
void
quadratic_real_roots(const Eigen::Matrix<double, 3, 1>& quadratic_coeffs,
                     std::array<double, 2>& solutions,
                     int& num_solutions,
                     double eps = 1e-10);

/// @brief Compute the real roots of a polynomial.
///
/// @param[in] A_coeffs: coefficients of the polynomial
/// @return real roots of the polynomial
std::vector<double>
polynomial_real_roots(const VectorXr& A_coeffs);

/// @brief Construct a formatted string for a variable raised to some power
///
/// @param[in] variable: variable
/// @param[in] degree: power to raise the variable to
/// @return formatted monomial string
std::string
formatted_monomial(const std::string& variable, int degree);

/// @brief Construct a formatted string for a term with given coefficient
/// and variable.
///
/// @param[in] coefficient: coefficient of the term
/// @param[in] variable: variable of the term
/// @param[in] precision: floating point precision
/// @return formatted term string
std::string
formatted_term(double coefficient, std::string variable, int precision = 16);

/// @brief Construct a formatted string for a polynomial with given coefficients
/// TODO Implement separate method for polynomial mappings
///
/// @param[in] P_coeffs: coefficients of the polynomial
/// @param[in] precision: floating point precision
/// @return formatted polynomial string
template<int degree, int dimension>
std::string
formatted_polynomial(
  const Eigen::Ref<const Eigen::Matrix<double, degree + 1, dimension>>&
    polynomial_coeffs,
  int precision = 16)
{
  std::stringstream polynomial_string;

  // Handle trivial case
  if (polynomial_coeffs.cols() == 0)
    return "";

  for (int i = 0; i < polynomial_coeffs.cols(); ++i) {
    polynomial_string << std::fixed << std::setprecision(precision)
                      << polynomial_coeffs(0, i);
    for (int j = 1; j < polynomial_coeffs.rows(); ++j) {
      std::string monomial_string = formatted_monomial("t", j);
      polynomial_string << formatted_term(
        polynomial_coeffs(j, i), monomial_string, precision);
    }
    polynomial_string << "\n";
  }

  return polynomial_string.str();
}

/// @brief Substitute some variable value that supports addition,
/// multiplication, and double multiplication into a polynomial.
///
/// Abstractly, this represents the polynomial evaluation homomorphism between
/// R[X] and some R-Algebra S.
///
/// @param[in] A_coeffs: coefficients of the polynomial function
/// @param[in] t: evaluation point
/// @return evaluated polynomial value
template<typename VariableType>
VariableType
substitute_polynomial(const VectorXr& A_coeffs, VariableType t)
{
  if (A_coeffs.size() == 0)
    return VariableType();

  // Initialize polynomial evaluation with the constant term
  VariableType A = VariableType(A_coeffs[0]);

  // Add the higher order terms
  int degree = A_coeffs.size() - 1;
  VariableType power_of_t = t;
  for (int i = 1; i <= degree; ++i) {
    A = A + VariableType(A_coeffs[i]) * power_of_t;
    power_of_t = power_of_t * t;
  }

  return A;
}
