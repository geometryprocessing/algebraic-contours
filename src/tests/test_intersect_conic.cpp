#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "intersect_conic.h"
#include "zwart_powell_spline.h"
#include "generate_shapes.h"


TEST_CASE ( "Points in the triangle patches are recognized", "[parametrize_contour]")
{
  std::array<Eigen::Matrix<double, 3, 1>, 3> triangle_boundary_coeffs;
  double r = 0.5;
  double eps = 1e-5;


  SECTION ( "First patch")
  {
    generate_zwart_powell_spline_patch_boundaries(0, triangle_boundary_coeffs);
    ConvexPolygon polygon(triangle_boundary_coeffs);

    REQUIRE ( polygon.contains(Eigen::Vector2d(-r/2.0, 0.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(-r+eps, 0.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(-r+eps, -r+eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(-r+eps, r-eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0-eps, 0.0)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r-eps, 0.0)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r-eps, -r-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r-eps, r+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0+eps, 0.0)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r+eps, -r)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r+eps, r)) );
  }
  SECTION ( "Second patch")
  {
    generate_zwart_powell_spline_patch_boundaries(1, triangle_boundary_coeffs);
    ConvexPolygon polygon(triangle_boundary_coeffs);

    REQUIRE ( polygon.contains(Eigen::Vector2d(r/2.0, 0.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(r-eps, 0.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(r-eps, -r+eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(r-eps, r-eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0+eps, 0.0)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(r+eps, 0.0)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(r+eps, -r-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(r+eps, r+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0-eps, 0.0)) );
  }
  SECTION ( "Third patch")
  {
    generate_zwart_powell_spline_patch_boundaries(2, triangle_boundary_coeffs);
    ConvexPolygon polygon(triangle_boundary_coeffs);

    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, -r/2.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, -r+eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(-r+eps, -r+eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(r-eps, -r+eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, 0.0-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0, -r-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r-eps, -r-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(r+eps, -r-eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0, 0.0+eps)) );
  }
  SECTION ( "Fourth patch")
  {
    generate_zwart_powell_spline_patch_boundaries(3, triangle_boundary_coeffs);
    ConvexPolygon polygon(triangle_boundary_coeffs);

    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, r/2.0)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, r-eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(-r+eps, r-eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(r-eps, r-eps)) );
    REQUIRE ( polygon.contains(Eigen::Vector2d(0.0, 0.0+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0, r+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(-r-eps, r+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(r+eps, r+eps)) );
    REQUIRE ( !polygon.contains(Eigen::Vector2d(0.0, 0.0-eps)) );
  }
}
