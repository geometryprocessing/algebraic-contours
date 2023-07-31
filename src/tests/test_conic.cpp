#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "conic.h"


TEST_CASE( "A conic can be pulled back by a quadratic", "[conic]" ) {
  Eigen::Matrix<double, 3, 2> P_coeffs;
  Eigen::Matrix<double, 3, 1> Q_coeffs;
  Eigen::Matrix<double, 6, 1> F_coeffs;

  SECTION ( "Zero case" )
  {
    P_coeffs << 0, 2,
                0, 3,
                0, 1;
    Q_coeffs << 1, 0, 0;
    F_coeffs << 0, 0, 0, 0, 0, 0;
    Conic conic(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> pullback;
    conic.pullback_quadratic_function<1>(F_coeffs, pullback);
    REQUIRE( float_equal(pullback(-1.0)(0), 0.0) );
    REQUIRE( float_equal(pullback(0.0)(0), 0.0) );
    REQUIRE( float_equal(pullback(1.0)(0), 0.0) );
  }

  SECTION ( "Unit pullback case" )
  {
    P_coeffs << 0, 2,
                0, 3,
                0, 1;
    Q_coeffs << 1, 0, 0;
    F_coeffs << 1, 0, 0, 0, 0, 0;
    Conic conic(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> pullback;
    conic.pullback_quadratic_function<1>(F_coeffs, pullback);
    REQUIRE( float_equal(pullback(-1.0)(0), 1.0) );
    REQUIRE( float_equal(pullback(0.0)(0), 1.0) );
    REQUIRE( float_equal(pullback(1.0)(0), 1.0) );
  }

  SECTION ( "u projection case" )
  {
    P_coeffs << 1, 1,
                2, -2,
                1, 1;
    Q_coeffs << 1, 0, 1;
    F_coeffs << 0, 1, 0, 0, 0, 0;
    Conic conic(P_coeffs, Q_coeffs);
    RationalFunction<4, 1> pullback;
    conic.pullback_quadratic_function<1>(F_coeffs, pullback);
    REQUIRE( float_equal(pullback(-1.0)(0), 0.0) );
    REQUIRE( float_equal(pullback(0.0)(0), 1.0) );
    REQUIRE( float_equal(pullback(1.0)(0), 2.0) );
  }

}