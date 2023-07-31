// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <Eigen/Core>
#include <iostream>

#include "bivariate_quadratic_function.h"
#include "common.h"
#include "polynomial_function.h"
#include "rational_function.h"

enum class ConicType
{
  ellipse,
  hyperbola,
  parabola,
  parallel_lines,
  intersecting_lines,
  line,
  point,
  empty,
  plane,
  error,
  unknown
};

/// @brief Explicit representation of a conic segment
// class Conic : private RationalFunction {
class Conic : public RationalFunction<2, 2>
{
public:
  /// Build a trivial conic;
  Conic() { m_type = ConicType::unknown; }

  Conic(const Matrix3x2r& numerator_coeffs,                   // 3x2
        const Eigen::Matrix<double, 3, 1>& denominator_coeffs // 3
        )
    : RationalFunction<2, 2>(numerator_coeffs, denominator_coeffs)
  {
    m_type = ConicType::unknown;
    assert(is_valid());
  }

  Conic(const Matrix3x2r& numerator_coeffs,                    // 3x2
        const Eigen::Matrix<double, 3, 1>& denominator_coeffs, // 3
        const ConicType& type)
    : RationalFunction<2, 2>(numerator_coeffs, denominator_coeffs)
    , m_type(type)
  {
    assert(is_valid());
  }

  Conic(const Matrix3x2r& numerator_coeffs,
        const Eigen::Matrix<double, 3, 1>& denominator_coeffs,
        const Interval domain)
    : RationalFunction<2, 2>(numerator_coeffs, denominator_coeffs, domain)
  {
    m_type = ConicType::unknown;
    assert(is_valid());
  }

  Conic(const Matrix3x2r& numerator_coeffs,
        const Eigen::Matrix<double, 3, 1>& denominator_coeffs,
        const Interval domain,
        const ConicType& type)
    : RationalFunction<2, 2>(numerator_coeffs, denominator_coeffs, domain)
    , m_type(type)
  {
    assert(is_valid());
  }

  /// @brief Get the type (e.g. hyperbola, line, etc.) of the conic.
  ///
  /// @return type identifier
  ConicType get_type() const;

  void transform(const Matrix2x2r& rotation, const PlanarPoint& translation);

  template<size_t dimension>
  void pullback_quadratic_function(
    const Eigen::Matrix<double, 6, dimension>& F_coeffs,
    RationalFunction<4, dimension>& pullback_function) const
  {
    SPDLOG_TRACE("Pulling back conic by quadratic function {}",
                 formatted_bivariate_quadratic_mapping(F_coeffs));

    // Separate the individual polynomial coefficients from the rational
    // function
    Matrix3x2r P_coeffs = get_numerators();
    Eigen::Matrix<double, 3, 1> u_coeffs = P_coeffs.col(0);
    Eigen::Matrix<double, 3, 1> v_coeffs = P_coeffs.col(1);
    Eigen::Matrix<double, 3, 1> Q_coeffs = get_denominator();
    SPDLOG_TRACE("u function before pullback: ({})/({})",
                 formatted_polynomial(u_coeffs),
                 formatted_polynomial(Q_coeffs));
    SPDLOG_TRACE("v function before pullback: ({})/({})",
                 formatted_polynomial(v_coeffs),
                 formatted_polynomial(Q_coeffs));

    // Compute (homogenized) polynomial coefficients for the quadratic terms
    Eigen::Matrix<double, 5, 1> QQ_coeffs, Qu_coeffs, Qv_coeffs, uv_coeffs,
      uu_coeffs, vv_coeffs;
    compute_polynomial_mapping_product<2, 2, 1>(Q_coeffs, Q_coeffs, QQ_coeffs);
    compute_polynomial_mapping_product<2, 2, 1>(Q_coeffs, u_coeffs, Qu_coeffs);
    compute_polynomial_mapping_product<2, 2, 1>(Q_coeffs, v_coeffs, Qv_coeffs);
    compute_polynomial_mapping_product<2, 2, 1>(u_coeffs, v_coeffs, uv_coeffs);
    compute_polynomial_mapping_product<2, 2, 1>(u_coeffs, u_coeffs, uu_coeffs);
    compute_polynomial_mapping_product<2, 2, 1>(v_coeffs, v_coeffs, vv_coeffs);

    // Combine quadratic monomial functions into a matrix
    Eigen::Matrix<double, 5, 6> monomial_coeffs;
    monomial_coeffs.setZero();
    monomial_coeffs.col(0).head(QQ_coeffs.size()) = QQ_coeffs;
    monomial_coeffs.col(1).head(Qu_coeffs.size()) = Qu_coeffs;
    monomial_coeffs.col(2).head(Qv_coeffs.size()) = Qv_coeffs;
    monomial_coeffs.col(3).head(uv_coeffs.size()) = uv_coeffs;
    monomial_coeffs.col(4).head(uu_coeffs.size()) = uu_coeffs;
    monomial_coeffs.col(5).head(vv_coeffs.size()) = vv_coeffs;
    spdlog::trace("Monomial coefficient matrix:\n{}", monomial_coeffs);

    // Compute the pulled back rational function numerator
    spdlog::trace("Quadratic coefficient matrix:\n{}", F_coeffs);
    Eigen::Matrix<double, 5, dimension> pullback_coeffs =
      monomial_coeffs * F_coeffs;
    SPDLOG_TRACE("Pullback numerator: {}",
                 formatted_polynomial(pullback_coeffs));
    SPDLOG_TRACE("Pullback denominator: {}", formatted_polynomial(QQ_coeffs));

    pullback_function =
      RationalFunction<4, dimension>(pullback_coeffs, QQ_coeffs, domain());
  }

  friend std::ostream& operator<<(std::ostream& out, const Conic& F);

private:
  ConicType m_type;

  bool is_valid() const;

  std::string formatted_conic() const;
};

std::ostream&
operator<<(std::ostream& out, const Conic& F);
