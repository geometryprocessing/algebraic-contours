#include "common.h"
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "convert_conic.h"

TEST_CASE( "The eigen decomposition of the identity is trivial", "[parametrize_contour]" ) {
  Matrix2x2r A(2,2);
  A << 1, 0,
       0, 1;
  std::array<double, 2> eigenvalues;
  Matrix2x2r rotation;
  compute_symmetric_matrix_eigen_decomposition(A, eigenvalues, rotation);

  REQUIRE(float_equal(eigenvalues[0], 1.0));
  REQUIRE(float_equal(eigenvalues[1], 1.0));
  REQUIRE(float_equal(rotation(0, 0), 1.0));
  REQUIRE(float_equal(rotation(0, 1), 0.0));
  REQUIRE(float_equal(rotation(1, 0), 0.0));
  REQUIRE(float_equal(rotation(1, 1), 1.0));
}

TEST_CASE( "The eigen decomposition of a diagonal matrix is trivial", "[parametrize_contour]" ) {
  Matrix2x2r A(2,2);
  A << 5, 0,
       0, 0.5;
  std::array<double, 2> eigenvalues;
  Matrix2x2r rotation;
  compute_symmetric_matrix_eigen_decomposition(A, eigenvalues, rotation);

  REQUIRE(float_equal(eigenvalues[0], 5.0));
  REQUIRE(float_equal(eigenvalues[1], 0.5));
  // Not generally unique
  //REQUIRE(float_equal(rotation(0, 0), 1.0)); 
  //REQUIRE(float_equal(rotation(0, 1), 0.0));
  //REQUIRE(float_equal(rotation(1, 0), 0.0));
  //REQUIRE(float_equal(rotation(1, 1), 1.0));
}

TEST_CASE( "The eigen decomposition of U^T D U is D, U", "[parametrize_contour]" ) {
  Matrix2x2r D(2,2);
  Matrix2x2r U(2,2);
  Matrix2x2r A(2,2);

  SECTION("Singular values (5, 0.5), rotation angle 0.1") {
    D << 5, 0, 0, 0.5;
    double theta = 0.1;
    U << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
    A = U.transpose() * D * U;
    std::array<double, 2> eigenvalues;
    Matrix2x2r rotation;
    compute_symmetric_matrix_eigen_decomposition(A, eigenvalues, rotation);

    spdlog::set_level(spdlog::level::debug);
    spdlog::debug("Expected rotation is {}", U);
    spdlog::debug("Computed rotation is {}", rotation);
    REQUIRE(float_equal(eigenvalues[0], 5.0));
    REQUIRE(float_equal(eigenvalues[1], 0.5));
    // Not generally unique
    //REQUIRE(float_equal(rotation(0, 0), U(0, 0)));
    //REQUIRE(float_equal(rotation(0, 1), U(0, 1)));
    //REQUIRE(float_equal(rotation(1, 0), U(1, 0)));
    //REQUIRE(float_equal(rotation(1, 1), U(1, 1)));
  }

  SECTION("Singular values (-5, -0.5), rotation angle 1") {
    D << -5, 0, 0, -0.5;
    double theta = 1;
    U << std::cos(theta), -std::sin(theta), std::sin(theta), std::cos(theta);
    A = U.transpose() * D * U;
    std::array<double, 2> eigenvalues;
    Matrix2x2r rotation;
    compute_symmetric_matrix_eigen_decomposition(A, eigenvalues, rotation);

    SPDLOG_DEBUG("Expected rotation is {}", U);
    SPDLOG_DEBUG("Computed rotation is {}", rotation);
    REQUIRE(float_equal(eigenvalues[0], -0.5));
    REQUIRE(float_equal(eigenvalues[1], -5.0));
    // Not generally unique
    //REQUIRE(float_equal(rotation(0, 0), U(0, 1)));
    //REQUIRE(float_equal(rotation(1, 0), U(1, 1)));
    //REQUIRE(float_equal(rotation(0, 1), -U(1, 1)));
    //REQUIRE(float_equal(rotation(1, 1), U(0, 1)));
  }
}
