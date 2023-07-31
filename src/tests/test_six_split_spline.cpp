#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "generate_shapes.h"
#include "six_split_spline.h"


TEST_CASE ( "Surface mappings can be found", "[six_split_spline]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  SECTION ( "Constant" )
  {
    SpatialVector p(1.0, 2.0, 3.0);
    SpatialVector zero(0.0, 0.0, 0.0);
    std::array<TriangleCornerFunctionData, 3> position_data;
    position_data[0] = TriangleCornerFunctionData(p, zero, zero);
    position_data[1] = TriangleCornerFunctionData(p, zero, zero);
    position_data[2] = TriangleCornerFunctionData(p, zero, zero);
    std::array<Eigen::Matrix<double, 6, 3>, 6> surface_mappings;
    generate_six_split_spline_patch_surface_mapping<double>(
      position_data,
      surface_mappings
    );
    SpatialVector q;
    evaluate_quadratic_mapping<3>(surface_mappings[0], PlanarPoint(0.25, 0.25), q);

    REQUIRE ( surface_mappings.size() == 6 );
    REQUIRE ( 
      vector_equal<3>(
        q,
        p
      )
    );
  }

  SECTION ( "Linear" )
  {
    SpatialVector p0(1.0, 0.0, 1.0);
    SpatialVector p1(0.0, 1.0, 1.0);
    SpatialVector p2(0.0, 0.0, 1.0);
    SpatialVector zero(0.0, 0.0, 0.0);
    std::array<TriangleCornerFunctionData, 3> position_data;
    position_data[0] = TriangleCornerFunctionData(p0, p1 - p0, p2 - p0);
    position_data[1] = TriangleCornerFunctionData(p1, p2 - p1, p0 - p1);
    position_data[2] = TriangleCornerFunctionData(p2, p0 - p2, p1 - p2);
    std::array<Eigen::Matrix<double, 6, 3>, 6> control_points;
    generate_six_split_spline_patch_control_points<double>(
      position_data,
      control_points
    );

    std::array<Eigen::Matrix<double, 6, 3>, 6> surface_mappings;
    generate_six_split_spline_patch_surface_mapping<double>(
      position_data,
      surface_mappings
    );

    REQUIRE ( control_points.size() == 6 );
    REQUIRE ( 
      vector_equal<3>(
        control_points[0].row(5),
        (p0 + p1 + p2) / 3.0
      )
    );
  }

  SECTION ( "Quadratic" )
  {
    spdlog::set_level(spdlog::level::debug);
    // z = 3xy + x^2 + 2y^2
    // Note: Order is important here. Swapping the order changes which two coordinates
    // of u, v, and w are used
    std::array<TriangleCornerFunctionData, 3> position_data;
    position_data[0] = TriangleCornerFunctionData(
      SpatialVector( 1, 0,  1),
      SpatialVector(-1, 1,  1),
      SpatialVector(-1, 0, -2)
    );
    position_data[1] = TriangleCornerFunctionData(
      SpatialVector(0,  1,  2),
      SpatialVector(0, -1, -4),
      SpatialVector(1, -1, -1)
    );
    position_data[2] = TriangleCornerFunctionData(
      SpatialVector(0, 0, 0),
      SpatialVector(1, 0, 0),
      SpatialVector(0, 1, 0)
    );
    std::array<Eigen::Matrix<double, 6, 3>, 6> control_points;
    generate_six_split_spline_patch_control_points<double>(
      position_data,
      control_points
    );

    std::array<Eigen::Matrix<double, 6, 3>, 6> surface_mappings;
    generate_six_split_spline_patch_surface_mapping<double>(
      position_data,
      surface_mappings
    );
    MatrixXr expected_surface_mapping(6, 3);
    expected_surface_mapping <<
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 3,
      0, 0, 1,
      0, 0, 2;
    
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[0]);
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[1]);
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[2]);
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[3]);
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[4]);
    spdlog::debug("Computed surface mapping:\n{}", surface_mappings[5]);

    REQUIRE ( control_points.size() == 6 );
    REQUIRE ( matrix_equal(surface_mappings[0], expected_surface_mapping) );
  }
}
