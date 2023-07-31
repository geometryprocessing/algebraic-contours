#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "compute_intersections.h"
#include "generate_shapes.h"

TEST_CASE ( "The intersections of two curves can be found", "[compute_intersections]")
{
  Eigen::Matrix<double, 5, 2> first_P_coeffs;
  Eigen::Matrix<double, 5, 1> first_Q_coeffs;
  Eigen::Matrix<double, 5, 2> second_P_coeffs;
  Eigen::Matrix<double, 5, 1> second_Q_coeffs;
  std::vector<double> intersections;
  std::vector<double> second_curve_intersections;
  IntersectionStats intersection_stats;
  IntersectionParameters intersect_params;
  intersect_params.use_heuristics=false;

  SECTION ( "Simple linear functions" )
  {
    first_P_coeffs << 0, 0,
                      1, 0,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << 0, 0,
                       0, 1,
                       0, 0,
                       0, 0,
                       0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 0.0) );
  }
  SECTION ( "Integer linear functions" )
  {
    first_P_coeffs << 1, 0,
                      2, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  4, 0,
                       -1, 1,
                        0, 0,
                        0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }
  SECTION ( "General linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                        0, 0,
                        0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Bounded linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                       0, 0,
                       0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;
    Interval bounds(-2, 2);

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, bounds);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, bounds);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Perturbed linear functions" )
  {
    first_P_coeffs << 0, 0,
                      0, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << 0, 1e-6,
                       1,    0,
                       0, 0,
                       0, 0,
                       0,    0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1e-6) );
  }

  SECTION ( "General linear functions" )
  {
    first_P_coeffs <<  1, 0,
                       0, 1,
                       0, 0,
                       0, 0,
                      -1, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << -1, 0,
                        0, 1,
                        0, 0,
                        0, 0,
                        1, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 2 );
  }
}

TEST_CASE ( "The intersections of two curves can be found with bezier clipping", "[compute_intersections]")
{
  Eigen::Matrix<double, 5, 2> first_P_coeffs;
  Eigen::Matrix<double, 5, 1> first_Q_coeffs;
  Eigen::Matrix<double, 5, 2> second_P_coeffs;
  Eigen::Matrix<double, 5, 1> second_Q_coeffs;
  std::vector<double> intersections;
  std::vector<double> second_curve_intersections;
  IntersectionStats intersection_stats;
  IntersectionParameters intersect_params;
  intersect_params.use_heuristics=false;
  Interval domain(-2, 2);

  SECTION ( "Simple linear functions" )
  {
    first_P_coeffs << 0, 0,
                      1, 0,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << 0, 0,
                       0, 1,
                       0, 0,
                       0, 0,
                       0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 0.0) );
  }
  SECTION ( "Integer linear functions" )
  {
    first_P_coeffs << 1, 0,
                      2, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  4, 0,
                       -1, 1,
                        0, 0,
                        0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }
  SECTION ( "General linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                        0, 0,
                        0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Bounded linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                        0, 0,
                        0, 0,
                        0, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;
    Interval bounds(-2, 2);

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Perturbed linear functions" )
  {
    first_P_coeffs << 0, 0,
                      0, 1,
                      0, 0,
                      0, 0,
                      0, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << 0, 1e-6,
                       1,    0,
                       0,    0,
                       0,    0,
                       0,    0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1e-6) );
  }

  SECTION ( "General linear functions" )
  {
    first_P_coeffs <<  1, 0,
                       0, 1,
                       0, 0,
                       0, 0,
                      -1, 0;
    first_Q_coeffs << 1, 0, 0, 0, 0;
    second_P_coeffs << -1, 0,
                        0, 1,
                        0, 0,
                        0, 0,
                        1, 0;
    second_Q_coeffs << 1, 0, 0, 0, 0;

    RationalFunction<4, 2> first_image_segment(first_P_coeffs, first_Q_coeffs, domain);
    RationalFunction<4, 2> second_image_segment(second_P_coeffs, second_Q_coeffs, domain);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersect_params, intersections, second_curve_intersections, intersection_stats);

    REQUIRE ( intersections.size() == 2 );
  }
}

TEST_CASE ( "A potentially self intersecting curve is split", "[compute_intersections]")
{
  Eigen::Matrix<double, 5, 2> first_P_coeffs;
  Eigen::Matrix<double, 5, 1> first_Q_coeffs;
  std::vector<std::vector<double>> split_points;
  SECTION ( "Circle" )
  {
    first_P_coeffs <<  1, 0,
                       0, 1,
                       -1, 0,
                       0, 0,
                       0, 0;
    first_Q_coeffs << 1, 0, 1, 0, 0;
    Interval domain(-2, 2);

    RationalFunction<4, 2> planar_segment(first_P_coeffs, first_Q_coeffs, domain);
    std::vector<RationalFunction<4, 2>> planar_curves = { planar_segment };
    split_planar_curves_no_self_intersection(planar_curves, split_points);
  }
}