#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <algorithm>

#include "generate_shapes.h"
#include "generate_position_data.h"
#include "quadratic_spline_surface.h"
#include "twelve_split_spline.h"


// Test that a quadratic surface can be reproduced from analytic corner and midpoint data
bool test_twelve_split_quadratic_reproduction(
    double uv_coeff,
    double uu_coeff,
    double vv_coeff
) {
  Eigen::MatrixXd V(3, 2);
  Eigen::MatrixXi F(1, 3);
  V <<
    1.0,  0.0,
    0.0,  1.0,
    0.0,  0.0;
  F <<
    0, 1, 2;

  ParametricAffineManifold parametric_affine_manifold(F, V);
  QuadraticPositionFunction position_func(uv_coeff, uu_coeff, vv_coeff);
  QuadraticGradientFunction gradient_func(uv_coeff, uu_coeff, vv_coeff);

  // Generate function data
  std::vector<std::array<TriangleCornerFunctionData, 3>> corner_data;
  std::vector<std::array<TriangleMidpointFunctionData, 3>> midpoint_data;
  generate_parametric_affine_manifold_corner_data(
    position_func,
    gradient_func,
    parametric_affine_manifold,
    corner_data
  );
  generate_parametric_affine_manifold_midpoint_data(
    gradient_func,
    parametric_affine_manifold,
    midpoint_data
  );

  std::array<Eigen::Matrix<double, 6, 3>, 12> surface_mappings;
  generate_twelve_split_spline_patch_surface_mapping<double>(
    corner_data[0],
    midpoint_data[0],
    surface_mappings
  );
  SpatialVector q;
  evaluate_quadratic_mapping<3>(surface_mappings[0], PlanarPoint(0.2, 0.3), q);

  if ( surface_mappings.size() != 12 ) return false;
  if ( 
    !vector_equal(
      q,
      position_func(0.2, 0.3)
    )
  ) {
    return false;
  }

  return true;
}


TEST_CASE ( "Twelve split surface mappings can be found", "[twelve_split_spline]")
{


  SECTION ( "Constant" )
  {
    // Build constant function triangle data
    SpatialVector p(1.0, 2.0, 3.0);
    SpatialVector zero(0.0, 0.0, 0.0);
    std::array<TriangleCornerFunctionData, 3> corner_data;
    corner_data[0] = TriangleCornerFunctionData(p, zero, zero);
    corner_data[1] = TriangleCornerFunctionData(p, zero, zero);
    corner_data[2] = TriangleCornerFunctionData(p, zero, zero);
    std::array<TriangleMidpointFunctionData, 3> midpoint_data;
    midpoint_data[0] = TriangleMidpointFunctionData(zero);
    midpoint_data[1] = TriangleMidpointFunctionData(zero);
    midpoint_data[2] = TriangleMidpointFunctionData(zero);
    std::array<Eigen::Matrix<double, 6, 3>, 12> surface_mappings;
    generate_twelve_split_spline_patch_surface_mapping<double>(
      corner_data,
      midpoint_data,
      surface_mappings
    );
    SpatialVector q;
    evaluate_quadratic_mapping<3>(surface_mappings[0], PlanarPoint(0.25, 0.25), q);

    REQUIRE ( surface_mappings.size() == 12 );
    REQUIRE ( vector_equal(q, p) );
  }

  SECTION ( "Constant" )
  {
    // Build constant function triangle data
    SpatialVector p(1, 2, 3);
    SpatialVector zero(0, 0, 0);
    std::array<TriangleCornerData<SpatialVector>, 3> corner_data;
    corner_data[0] = TriangleCornerData<SpatialVector>(p, zero, zero);
    corner_data[1] = TriangleCornerData<SpatialVector>(p, zero, zero);
    corner_data[2] = TriangleCornerData<SpatialVector>(p, zero, zero);
    std::array<TriangleMidpointData<SpatialVector>, 3> midpoint_data;
    midpoint_data[0] = TriangleMidpointData<SpatialVector>(zero);
    midpoint_data[1] = TriangleMidpointData<SpatialVector>(zero);
    midpoint_data[2] = TriangleMidpointData<SpatialVector>(zero);
    std::array<Eigen::Matrix<double, 6, 3>, 12> surface_mappings;
    generate_twelve_split_spline_patch_surface_mapping<double>(
      corner_data,
      midpoint_data,
      surface_mappings
    );

    REQUIRE ( surface_mappings.size() == 12 );
  }

  SECTION ( "Linear" )
  {
    // Build linear "quadratic" functionals
    REQUIRE (test_twelve_split_quadratic_reproduction(0.0, 0.0, 0.0));
  }

  SECTION ( "Quadratic" )
  {
    // Test linear "quadratic" functionals
    REQUIRE (test_twelve_split_quadratic_reproduction(1.0, 0.0, 0.0));
    REQUIRE (test_twelve_split_quadratic_reproduction(0.0, 1.0, 0.0));
    REQUIRE (test_twelve_split_quadratic_reproduction(0.0, 0.0, 1.0));
    REQUIRE (test_twelve_split_quadratic_reproduction(1.0, 2.0, -1.0));
  }
}

