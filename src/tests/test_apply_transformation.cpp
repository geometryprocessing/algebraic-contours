#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "apply_transformation.h"


TEST_CASE ( "Vectors can be converted to homogeneous coordinates", "[apply_transformation]")
{

  SECTION ( "Zero case" )
  {
    SpatialVector vector(3);
    vector << 0, 0, 0;
    Eigen::Matrix<double, 4, 1> homogeneous_coords;
    convert_point_to_homogeneous_coords(vector, homogeneous_coords);

    REQUIRE( float_equal(homogeneous_coords(0), 0.0) );
    REQUIRE( float_equal(homogeneous_coords(1), 0.0) );
    REQUIRE( float_equal(homogeneous_coords(2), 0.0) );
    REQUIRE( float_equal(homogeneous_coords(3), 1.0) );
  }

  SECTION ( "General case" )
  {
    SpatialVector vector(3);
    vector << 1, 2, 3;
    Eigen::Matrix<double, 4, 1> homogeneous_coords;
    convert_point_to_homogeneous_coords(vector, homogeneous_coords);

    REQUIRE( float_equal(homogeneous_coords(0), 1.0) );
    REQUIRE( float_equal(homogeneous_coords(1), 2.0) );
    REQUIRE( float_equal(homogeneous_coords(2), 3.0) );
    REQUIRE( float_equal(homogeneous_coords(3), 1.0) );
  }
}

TEST_CASE ( "Homogeneous coordinates can be translated", "[apply_transformation]")
{
  SpatialVector translation(3);
  Eigen::Matrix<double, 4, 1>  e1, e2, e3;
  convert_point_to_homogeneous_coords(elementary_basis_vector<3>(0), e1);
  convert_point_to_homogeneous_coords(elementary_basis_vector<3>(1), e2);
  convert_point_to_homogeneous_coords(elementary_basis_vector<3>(2), e3);

  SECTION ( "Zero case" )
  {
    translation << 0, 0, 0;
    Eigen::Matrix<double, 4, 4> translation_matrix = translation_projective_matrix(translation);
    Eigen::Matrix<double, 4, 1> v1 = translation_matrix * e1;
    Eigen::Matrix<double, 4, 1> v2 = translation_matrix * e2;
    Eigen::Matrix<double, 4, 1> v3 = translation_matrix * e3;

    REQUIRE( float_equal(v1(0), 1.0) );
    REQUIRE( float_equal(v1(1), 0.0) );
    REQUIRE( float_equal(v1(2), 0.0) );
    REQUIRE( float_equal(v2(0), 0.0) );
    REQUIRE( float_equal(v3(0), 0.0) );
  }
}