#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include "bivariate_quadratic_function.h"

TEST_CASE( "Patch monomials are computed at (0,0)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(0, 0));

  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 0.0 );
  REQUIRE( w[2] == 0.0 );
  REQUIRE( w[3] == 0.0 );
  REQUIRE( w[4] == 0.0 );
  REQUIRE( w[5] == 0.0 );
}

TEST_CASE( "Patch monomials are computed at (1,0)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(1, 0));

  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 1.0 );
  REQUIRE( w[2] == 0.0 );
  REQUIRE( w[3] == 0.0 );
  REQUIRE( w[4] == 1.0 );
  REQUIRE( w[5] == 0.0 );
}

TEST_CASE( "Patch monomials are computed at (0,1)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(0, 1));

  REQUIRE( w.size() == 6 );
  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 0.0 );
  REQUIRE( w[2] == 1.0 );
  REQUIRE( w[3] == 0.0 );
  REQUIRE( w[4] == 0.0 );
  REQUIRE( w[5] == 1.0 );
}

TEST_CASE( "Patch monomials are computed at (1,1)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(1, 1));

  REQUIRE( w.size() == 6 );
  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 1.0 );
  REQUIRE( w[2] == 1.0 );
  REQUIRE( w[3] == 1.0 );
  REQUIRE( w[4] == 1.0 );
  REQUIRE( w[5] == 1.0 );
}

TEST_CASE( "Patch monomials are computed at (0.5,1)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(0.5, 1));

  REQUIRE( w.size() == 6 );
  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 0.5 );
  REQUIRE( w[2] == 1.0 );
  REQUIRE( w[3] == 0.5 );
  REQUIRE( w[4] == 0.25 );
  REQUIRE( w[5] == 1.0 );
}

TEST_CASE( "Patch monomials are computed at (0.5,0.5)", "[patch]" ) {
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(PlanarPoint(0.5, 0.5));

  REQUIRE( w[0] == 1.0 );
  REQUIRE( w[1] == 0.5 );
  REQUIRE( w[2] == 0.5 );
  REQUIRE( w[3] == 0.25 );
  REQUIRE( w[4] == 0.25 );
  REQUIRE( w[5] == 0.25 );
}

//TEST_CASE( "Patches are constant for constant control points", "[patch]" ) {
//  std::vector<std::vector<SpatialVector>> P(3);
//  for (int i = 0; i < 3; ++i) {
//    P[i].resize(3);
//    for (int j = 0; j < 3; ++j) {
//      P[i][j].setZero();
//    }
//  }
//
//  SECTION( "Evaluation at the center gives 0 ") {
//    SpatialVector p;
//    PatchIndices patch_indices = { 1, 1, 0 };
//    PlanarPoint v;
//    v << 0, 0;
//    evaluate_spline_surface_patch(P, patch_indices, v, p);
//
//    REQUIRE( abs(p(0)) < 1e-12 );
//    REQUIRE( abs(p(1)) < 1e-12 );
//    REQUIRE( abs(p(2)) < 1e-12 );
//  }
//
//  SECTION( "Evaluation at a corner gives 0 ") {
//    SpatialVector p;
//    PatchIndices patch_indices = { 0, 0, 0 };
//    PlanarPoint v;
//    v << 0, 0;
//    evaluate_spline_surface_patch(P, patch_indices, v, p);
//
//    REQUIRE( abs(p(0)) < 1e-12 );
//    REQUIRE( abs(p(1)) < 1e-12 );
//    REQUIRE( abs(p(2)) < 1e-12 );
//  }
//
//  SECTION( "Evaluation at an arbitrary point gives 0 ") {
//    SpatialVector p;
//    PatchIndices patch_indices = { 1, 2, 1 };
//    PlanarPoint v;
//    v << 0.85, 0.45;
//    evaluate_spline_surface_patch(P, patch_indices, v, p);
//
//    REQUIRE( abs(p(0)) < 1e-12 );
//    REQUIRE( abs(p(1)) < 1e-12 );
//    REQUIRE( abs(p(2)) < 1e-12 );
//  }
//}
