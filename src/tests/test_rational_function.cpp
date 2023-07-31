#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "rational_function.h"
#include "generate_shapes.h"


TEST_CASE ( "The derivative of a rational function can be found", "[rational_function]")
{

  SECTION ( "Zero function")
  {
    Eigen::Matrix<double, 2, 1> P_coeffs;
    Eigen::Matrix<double, 2, 1> Q_coeffs;
    P_coeffs << 0, 0;
    Q_coeffs << 1, 0;
    RationalFunction<1, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<2, 1> F_derivative;
    F.compute_derivative(F_derivative);

    REQUIRE( float_equal(F_derivative(-1.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), 0.0) );
  }

  SECTION ( "Constant function")
  {
    Eigen::Matrix<double, 2, 1> P_coeffs;
    Eigen::Matrix<double, 2, 1> Q_coeffs;
    P_coeffs << 1, 0;
    Q_coeffs << 1, 0;
    RationalFunction<1, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<2, 1> F_derivative;
    F.compute_derivative(F_derivative);

    REQUIRE( float_equal(F_derivative(-1.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), 0.0) );
  }

  SECTION ( "Linear function")
  {
    Eigen::Matrix<double, 2, 1> P_coeffs;
    Eigen::Matrix<double, 2, 1> Q_coeffs;
    P_coeffs << -1, 2;
    Q_coeffs << 1, 0;
    RationalFunction<1, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<2, 1> F_derivative;
    F.compute_derivative(F_derivative);

    REQUIRE( float_equal(F_derivative(-1.0)(0), 2.0) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 2.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), 2.0) );
  }

  SECTION ( "Quadratic function")
  {
    Eigen::Matrix<double, 3, 1> P_coeffs;
    Eigen::Matrix<double, 3, 1> Q_coeffs;
    P_coeffs << 1, -2, 1;
    Q_coeffs << 1, 0;
    RationalFunction<2, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> F_derivative;
    F.compute_derivative(F_derivative);

    REQUIRE( float_equal(F_derivative(-1.0)(0), -4.0) );
    REQUIRE( float_equal(F_derivative(0.0)(0), -2.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), 0.0) );
  }

  SECTION ( "Inverse monomial function")
  {
    Eigen::Matrix<double, 3, 1> P_coeffs;
    Eigen::Matrix<double, 3, 1> Q_coeffs;
    P_coeffs << 1, 0, 0;
    Q_coeffs << 0, 0, 1;
    RationalFunction<2, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> F_derivative;
    F.compute_derivative(F_derivative);

    REQUIRE( float_equal(F_derivative(-1.0)(0), 2.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), -2.0) );
    REQUIRE( float_equal(F_derivative(2.0)(0), -0.25) );
  }

  SECTION ( "Inverse quadratic function")
  {
    Eigen::Matrix<double, 3, 1> P_coeffs;
    Eigen::Matrix<double, 3, 1> Q_coeffs;
    P_coeffs << 1, 0, 0;
    Q_coeffs << 1, 0, 1;
    RationalFunction<2, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> F_derivative;
    F.compute_derivative(F_derivative);

    // -2t / (1 + t^2)^2
    REQUIRE( float_equal(F_derivative(-1.0)(0), 0.5) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), -0.5) );
    REQUIRE( float_equal(F_derivative(2.0)(0), -0.16) );
  }

  SECTION ( "Rational function")
  {
    Eigen::Matrix<double, 3, 1> P_coeffs;
    Eigen::Matrix<double, 3, 1> Q_coeffs;
    P_coeffs << 1, 1, 0;
    Q_coeffs << 1, 0, 1;
    RationalFunction<2, 1> F(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> F_derivative;
    F.compute_derivative(F_derivative);

    // (1 + t^2) - (1 + t)(2t)
    // (1 - 2t - t^2) / (1 + t^2)^2
    REQUIRE( float_equal(F_derivative(-1.0)(0), 0.5) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 1.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), -0.5) );
    REQUIRE( float_equal(F_derivative(2.0)(0), -0.28) );
  }

  SECTION ( "Planar rational function")
  {
    Eigen::Matrix<double, 3, 2> P_coeffs;
    Eigen::Matrix<double, 3, 1> Q_coeffs;
    P_coeffs << 1, 1, 0, 1;
    Q_coeffs << 1, 0, 1;
    RationalFunction<2, 2> F(P_coeffs, Q_coeffs);
    RationalFunction<4, 2> F_derivative;
    F.compute_derivative(F_derivative);
    spdlog::debug("Function {}", F);
    spdlog::debug("Derivative {}", F_derivative);

    // -2t / (1 + t^2)^2 
    // (1 - 2t - t^2) / (1 + t^2)^2
    REQUIRE( float_equal(F_derivative(-1.0)(0), 0.5) );
    REQUIRE( float_equal(F_derivative(0.0)(0), 0.0) );
    REQUIRE( float_equal(F_derivative(1.0)(0), -0.5) );
    REQUIRE( float_equal(F_derivative(2.0)(0), -0.16) );

    REQUIRE( float_equal(F_derivative(-1.0)(1), 0.5) );
    REQUIRE( float_equal(F_derivative(0.0)(1), 1.0) );
    REQUIRE( float_equal(F_derivative(1.0)(1), -0.5) );
    REQUIRE( float_equal(F_derivative(2.0)(1), -0.28) );
  }

}