#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "compute_intersections.h"
#include "generate_shapes.h"



TEST_CASE ( "A potentially self intersecting curve is split", "[compute_intersections]")
{
  std::vector<std::vector<double>> split_points;
  Eigen::Matrix<double, 5, 2> P_coeffs;
  Eigen::Matrix<double, 5, 1> Q_coeffs;
  SECTION ( "Circle" )
  {
    P_coeffs <<  1, 0,
                       0, 1,
                       -1, 0,
                       0, 0,
                       0, 0;
    Q_coeffs << 1, 0, 1, 0, 0;
    Interval domain(-2, 2);

    RationalFunction<4, 2> planar_segment(P_coeffs, Q_coeffs, domain);
    std::vector<RationalFunction<4, 2>> planar_curves = { planar_segment };
    split_planar_curves_no_self_intersection(planar_curves, split_points);
    REQUIRE(split_points.size() == 1);
    REQUIRE(split_points[0].size() == 1);
  }


  SECTION ( "Spot bug" )
  {
    P_coeffs << -0.73120161816068385, -2.25202965564466018,
                -0.40059603916677894,  0.12008492684645028,
                -1.39905702192626835,  4.47876303596620318,
                -4.50593422038637392, -5.30152851610740328,
                -0.57491127825807053, -3.18194386122352757;
    Interval domain(-0.09280896795608366, -0.07905005296651063);
    // Observed to split at 32272 knots due to potential self intersections
    Q_coeffs << 1, 0, -2, 0, 1;
    RationalFunction<4, 2> planar_segment(P_coeffs, Q_coeffs, domain);
    std::vector<RationalFunction<4, 2>> planar_curves = { planar_segment };
    split_planar_curves_no_self_intersection(planar_curves, split_points);
    REQUIRE(split_points.size() == 1);
    REQUIRE(split_points[0].size() == -1);
  }
}