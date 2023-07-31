#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "compute_intersections.h"
#include "generate_shapes.h"


TEST_CASE ( "The intersections of two curves can be found", "[compute_intersections]")
{
  MatrixXr first_P_coeffs(3,2);
  VectorXr first_Q_coeffs(1);
  MatrixXr second_P_coeffs(3,2);
  VectorXr second_Q_coeffs(1);
  std::vector<double> intersections;
  IntersectionStats intersection_stats;

  SECTION ( "Simple linear functions" )
  {
    first_P_coeffs << 0, 0,
                      1, 0,
                      0, 0;
    first_Q_coeffs << 1;
    second_P_coeffs << 0, 0,
                       0, 1,
                       0, 0;
    second_Q_coeffs << 1;

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 0.0) );
  }
  SECTION ( "Integer linear functions" )
  {
    first_P_coeffs << 1, 0,
                      2, 1,
                      0, 0;
    first_Q_coeffs << 1;
    second_P_coeffs <<  4, 0,
                       -1, 1,
                        0, 0;
    second_Q_coeffs << 1;

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }
  SECTION ( "General linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0;
    first_Q_coeffs << 1;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                        0, 0;
    second_Q_coeffs << 1;

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Bounded linear functions" )
  {
    first_P_coeffs << 1, 0,
                      0.5, 1,
                      0, 0;
    first_Q_coeffs << 1;
    second_P_coeffs <<  2.5, 0,
                       -1, 1,
                        0, 0;
    second_Q_coeffs << 1;
    Interval bounds(-2, 2);

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs, bounds);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs, bounds);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1.0) );
  }

  SECTION ( "Perturbed linear functions" )
  {
    first_P_coeffs << 0, 0,
                      0, 1,
                      0, 0;
    first_Q_coeffs << 1;
    second_P_coeffs << 0, 1e-6,
                       1,    0,
                       0,    0;
    second_Q_coeffs << 1;

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 1 );
    REQUIRE ( float_equal(intersections[0], 1e-6) );
  }

  SECTION ( "General linear functions" )
  {
    first_P_coeffs <<  1, 0,
                       0, 1,
                      -1, 0;
    first_Q_coeffs << 1;
    second_P_coeffs << -1, 0,
                        0, 1,
                        1, 0;
    second_Q_coeffs << 1;

    RationalFunction first_image_segment(first_P_coeffs, first_Q_coeffs);
    RationalFunction second_image_segment(second_P_coeffs, second_Q_coeffs);
    compute_planar_curve_intersections(first_image_segment, second_image_segment, intersections, intersection_stats);

    REQUIRE ( intersections.size() == 2 );
  }
}