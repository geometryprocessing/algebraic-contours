// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "bivariate_quadratic_function.h"
#include "common.h"
#include "conic.h"
#include "polynomial_function.h"
#include "rational_function.h"

/// @brief Representation of a line segment, which is a rational function
/// with degree 1 numerator and degree 0 denominator.
class LineSegment : public Conic
{
public:
  LineSegment() {}

  LineSegment(const Matrix2x2r& numerator_coeffs)
  {
    init_conic_coefficients(numerator_coeffs);
  }

  LineSegment(const Matrix2x2r& numerator_coeffs, const Interval& input_domain)
  {
    init_conic_coefficients(numerator_coeffs);
    domain() = input_domain;
  }

  template<size_t dimension>
  void pullback_linear_function(
    const Eigen::Matrix<double, 3, dimension>& F_coeffs,
    RationalFunction<1, dimension>& pullback_function) const
  {
    SPDLOG_TRACE("Pulling back line segment by linear function {}",
                 formatted_bivariate_linear_mapping(F_coeffs));

    // Separate the individual polynomial coefficients from the rational
    // function
    Matrix3x2r P_coeffs = get_numerators();
    Eigen::Matrix<double, 3, 1> u_coeffs = P_coeffs.col(0);
    Eigen::Matrix<double, 3, 1> v_coeffs = P_coeffs.col(1);
    Eigen::Matrix<double, 2, 1> Q_coeffs(1.0, 0.0);
    SPDLOG_TRACE("u function before pullback: {}",
                 formatted_polynomial(u_coeffs));
    SPDLOG_TRACE("v function before pullback: {}",
                 formatted_polynomial(v_coeffs));

    // Combine quadratic monomial functions into a matrix
    Matrix2x3r monomial_coeffs;
    monomial_coeffs.setZero();
    monomial_coeffs(0, 0) = 1.0;
    monomial_coeffs(0, 1) = u_coeffs[0];
    monomial_coeffs(1, 1) = u_coeffs[1];
    monomial_coeffs(0, 2) = v_coeffs[0];
    monomial_coeffs(1, 2) = v_coeffs[1];
    spdlog::trace("Monomial coefficient matrix:\n{}", monomial_coeffs);

    // Compute the pulled back rational function numerator
    spdlog::trace("Linear coefficient matrix:\n{}", F_coeffs);
    Eigen::Matrix<double, 2, dimension> pullback_coeffs =
      monomial_coeffs * F_coeffs;
    SPDLOG_TRACE("Pullback function: {}",
                 formatted_polynomial(pullback_coeffs));

    pullback_function =
      RationalFunction<1, dimension>(pullback_coeffs, Q_coeffs, domain());
  }

private:
  // Given line segment coefficients, expand them to conic coefficients
  void init_conic_coefficients(const Matrix2x2r& numerator_coeffs)
  {
    // Build conic numerator with trivial quadratic term
    Matrix3x2r conic_numerator_coeffs;
    conic_numerator_coeffs.block(0, 0, 2, 2) = numerator_coeffs;
    conic_numerator_coeffs.block(2, 0, 1, 2).setZero();

    // Build constant 1 denominator
    Eigen::Matrix<double, 3, 1> conic_denominator_coeffs;
    conic_denominator_coeffs << 1.0, 0.0, 0.0;

    // Set conic coefficients
    set_numerators(conic_numerator_coeffs);
    set_denominator(conic_denominator_coeffs);
  }
};