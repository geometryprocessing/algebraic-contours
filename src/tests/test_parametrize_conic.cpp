#include "common.h"
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "parametrize_conic.h"
#include "polynomial_function.h"

bool _test_parametrization(const Eigen::Matrix<double, 6, 1> &C_coeffs) {
  spdlog::info("Testing conic with equation {}",
               formatted_bivariate_quadratic_mapping<1>(C_coeffs));

  std::vector<Conic> conics;
  parametrize_conic(C_coeffs, conics);
  for (int i = 0; i < conics.size(); ++i) {
    std::vector<PlanarPoint> points;
    conics[i].sample_points(10, points);
    for (int j = 0; j < points.size(); ++j) {
      PlanarPoint p = points[j];
      if (!float_equal(evaluate_quadratic(C_coeffs, p), 0.0)) {
        return false;
      }
    }
  }

  return true;
}

TEST_CASE("An ellipse is recognized", "[parametrize_contour]") {
  Eigen::Matrix<double, 6, 1> C_coeffs(6);
  C_coeffs << 1, 0, 0, 0, -1, -1;

  SECTION("Circle") { REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse); }
  SECTION("Vertically stretched ellipse") {
    C_coeffs(4) = -2;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("Horizontally stretched ellipse") {
    C_coeffs(5) = -2;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("Large circle") {
    C_coeffs(4) = -0.1;
    C_coeffs(5) = -0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("General standard form") {
    C_coeffs(0) = 5;
    C_coeffs(4) = -10;
    C_coeffs(5) = -0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("Negative standard form") {
    C_coeffs(0) = -5;
    C_coeffs(4) = 10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("Linear terms") {
    C_coeffs(0) = -5;
    C_coeffs(1) = -5;
    C_coeffs(2) = 5;
    C_coeffs(4) = 10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
  SECTION("General form") {
    C_coeffs(0) = -5;
    C_coeffs(1) = -5;
    C_coeffs(2) = 5;
    C_coeffs(3) = 1.9;
    C_coeffs(4) = 10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::ellipse);
  }
}

TEST_CASE("A hyperbola is recognized", "[parametrize_contour]") {
  Eigen::Matrix<double, 6, 1> C_coeffs(6);
  C_coeffs << 1, 0, 0, 0, 1, -1;

  SECTION("Vertically symmetric") {
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("Horizontally symmetric") {
    C_coeffs(4) = -1;
    C_coeffs(5) = 1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("Vertically stretched hyperbola") {
    C_coeffs(4) = 2;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("Horizontally stretched hyperbola") {
    C_coeffs(5) = -2;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("General standard form") {
    C_coeffs(0) = 5;
    C_coeffs(4) = -10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("Negative standard form") {
    C_coeffs(0) = -5;
    C_coeffs(4) = 10;
    C_coeffs(5) = -0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("Linear terms") {
    C_coeffs(0) = -5;
    C_coeffs(1) = -5;
    C_coeffs(2) = 5;
    C_coeffs(4) = -10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
  SECTION("General form") {
    C_coeffs(0) = -5;
    C_coeffs(1) = -5;
    C_coeffs(2) = 5;
    C_coeffs(3) = 2.1;
    C_coeffs(4) = 10;
    C_coeffs(5) = 0.1;
    REQUIRE(identify_conic(C_coeffs) == ConicType::hyperbola);
  }
}

TEST_CASE("An ellipse is parametrized") {
  Eigen::Matrix<double, 6, 1> C_coeffs(6);
  C_coeffs.setZero(6);
  C_coeffs(0) = -1.0;
  C_coeffs(4) = 1.0;
  C_coeffs(5) = 1.0;

  SECTION("Unit circle") { REQUIRE(_test_parametrization(C_coeffs)); }

  SECTION("Small circle") {
    C_coeffs(4) = 2;
    C_coeffs(5) = 2;
    REQUIRE(_test_parametrization(C_coeffs));
  }

  SECTION("General circle") {
    C_coeffs(0) = -0.5;
    C_coeffs(1) = -2;
    C_coeffs(2) = 3;
    C_coeffs(4) = 2;
    C_coeffs(5) = 2;
    REQUIRE(_test_parametrization(C_coeffs));
  }

  SECTION("General ellipse") {
    C_coeffs(0) = -0.5;
    C_coeffs(1) = -2;
    C_coeffs(2) = 3;
    C_coeffs(4) = 0.1;
    C_coeffs(5) = 2;
    REQUIRE(_test_parametrization(C_coeffs));
  }
}

TEST_CASE("A hyperbola is parametrized") {
  Eigen::Matrix<double, 6, 1> C_coeffs(6);
  C_coeffs.setZero(6);
  C_coeffs(0) = -1.0;
  C_coeffs(4) = 1.0;
  C_coeffs(5) = -1.0;

  SECTION("Unit hyperbola") { REQUIRE(_test_parametrization(C_coeffs)); }

  SECTION("Translated hyperbola") {
    C_coeffs(1) = -2;
    C_coeffs(2) = 3;

    REQUIRE(_test_parametrization(C_coeffs));
  }

  SECTION("Stretched hyperbola") {
    C_coeffs(4) = -2.0;
    C_coeffs(5) = 0.5;
    REQUIRE(_test_parametrization(C_coeffs));
  }

  SECTION("General hyperbola") {
    C_coeffs(0) = 0.5;
    C_coeffs(1) = -2;
    C_coeffs(2) = 3;
    C_coeffs(4) = 0.1;
    C_coeffs(5) = -2.0;
    REQUIRE(_test_parametrization(C_coeffs));
  }
}