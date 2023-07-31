#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "polynomial_function.h"
#include "generate_shapes.h"
#include "quadratic_spline_surface.h"


TEST_CASE ( "The cross product of polynomial mappings can be found", "[polynomial_function]")
{

  SECTION ( "Elementary constant functions")
  {
    Eigen::Matrix<double, 1, 3> A_coeffs;
    Eigen::Matrix<double, 1, 3> B_coeffs;
    A_coeffs << 1, 0, 0;
    B_coeffs << 0, 1, 0;
    Eigen::Matrix<double, 1, 3> cross_product_coeffs;
    compute_polynomial_mapping_cross_product<0, 0>(A_coeffs, B_coeffs, cross_product_coeffs);

    REQUIRE( float_equal(cross_product_coeffs(0, 0), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 1), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 2), 1.0) );
  }

  SECTION ( "Elementary linear functions")
  {
    Eigen::Matrix<double, 2, 3> A_coeffs;
    Eigen::Matrix<double, 2, 3> B_coeffs;
    A_coeffs << 2, 0, 0,
                1, 0, 0;
    B_coeffs << 0, 1, 0,
                0, 1, 0;
    Eigen::Matrix<double, 3, 3> cross_product_coeffs;
    compute_polynomial_mapping_cross_product<1, 1>(A_coeffs, B_coeffs, cross_product_coeffs);

    REQUIRE( float_equal(cross_product_coeffs(0, 0), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 1), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 2), 2.0) );
    REQUIRE( float_equal(cross_product_coeffs(1, 2), 3.0) );
    REQUIRE( float_equal(cross_product_coeffs(2, 2), 1.0) );
  }

  SECTION ( "General constant functions")
  {
    Eigen::Matrix<double, 1, 3> A_coeffs;
    Eigen::Matrix<double, 1, 3> B_coeffs;
    A_coeffs << 1, 2, 3;
    B_coeffs << 4, 5, 6;
    Eigen::Matrix<double, 1, 3> cross_product_coeffs;
    compute_polynomial_mapping_cross_product<0, 0>(A_coeffs, B_coeffs, cross_product_coeffs);

    REQUIRE( float_equal(cross_product_coeffs(0, 0), -3.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 1), 6.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 2), -3.0) );
  }

  SECTION ( "Cancelling linear functions")
  {
    Eigen::Matrix<double, 2, 3> A_coeffs;
    Eigen::Matrix<double, 2, 3> B_coeffs;
    A_coeffs << 1, 2, 3,
                1, 1, 1;
    B_coeffs << 4, 5, 6,
                1, 1, 1;
    Eigen::Matrix<double, 3, 3> cross_product_coeffs;
    compute_polynomial_mapping_cross_product<1, 1>(A_coeffs, B_coeffs, cross_product_coeffs);

    REQUIRE( float_equal(cross_product_coeffs(0, 0), -3.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 1), 6.0) );
    REQUIRE( float_equal(cross_product_coeffs(0, 2),-3.0) );
    
    REQUIRE( float_equal(cross_product_coeffs(1, 0), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(1, 1), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(1, 2), 0.0) );

    REQUIRE( float_equal(cross_product_coeffs(2, 0), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(2, 1), 0.0) );
    REQUIRE( float_equal(cross_product_coeffs(2, 2), 0.0) );
  }
}


TEST_CASE ( "The roots of a polynomial can be found", "[polynomial_function]")
{

  SECTION ( "Linear function")
  {
    VectorXr A_coeffs(2);
    A_coeffs << 1, 1;
    std::vector<double> roots = polynomial_real_roots(A_coeffs);
    
    REQUIRE( roots.size() == 1 );
    REQUIRE( float_equal(roots[0], -1.0) );
  }

  SECTION ( "Quadratic function with roots")
  {
    VectorXr A_coeffs(3);
    A_coeffs << -1, 0, 1;
    std::vector<double> roots = polynomial_real_roots(A_coeffs);
    
    REQUIRE( roots.size() == 2 );
    REQUIRE( ((float_equal(roots[0], -1.0)) || (float_equal(roots[0], 1.0))) );
    REQUIRE( ((float_equal(roots[1], -1.0)) || (float_equal(roots[1], 1.0))) );
  }

  SECTION ( "Quadratic function without roots")
  {
    VectorXr A_coeffs(3);
    A_coeffs << 1, 0, 1;
    std::vector<double> roots = polynomial_real_roots(A_coeffs);
    
    REQUIRE( roots.size() == 0 );
  }
}
