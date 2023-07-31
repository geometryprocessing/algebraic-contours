#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "generate_shapes.h"
#include "quadratic_spline_surface.h"


//bool test_quadratic_reproduction(
//  double u_curvature = 0.0,
//  double v_curvature = 0.0
//) {
//  int resolution = 10;
//  std::vector< std::vector<SpatialVector> > param_grid;
//  Eigen::MatrixXd param_V;
//  Eigen::MatrixXi param_F;
//  std::vector<std::vector<double>> param_l;
//  std::vector< std::vector<SpatialVector> > control_point_grid;
//  Eigen::MatrixXd control_V;
//  Eigen::MatrixXi control_F;
//  std::vector<std::vector<double>> l;
//  std::vector<std::vector<int>> face_to_patch_indices;
//  std::vector<int> patch_to_face_indices;
//
//  std::vector< std::vector<PlanarPoint> > layout_point_grid;
//  generate_global_layout_grid(resolution, layout_point_grid);
//  generate_quadratic_grid(layout_point_grid, u_curvature, v_curvature, 0, control_point_grid);
//  generate_mesh_from_grid(control_point_grid, control_V, control_F, l, true);
//  generate_plane_grid(resolution, 0.0, 0.0, param_grid);
//  generate_mesh_from_grid(param_grid, param_V, param_F, param_l, true);
//  QuadraticSplineParameters quadratic_spline_params;
//  quadratic_spline_params.spline_type = SplineType::powell_sabin_six_split;
//  OptimizationParameters optimization_params;
//  optimization_params.method == OptimizationParameters::least_squares;
//  QuadraticSplineSurface spline_surface(
//    control_V,
//    control_F,
//    l,
//    quadratic_spline_params,
//    optimization_params,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//  QuadraticSplineSurface param_surface(
//    param_V,
//    param_F,
//    param_l,
//    quadratic_spline_params,
//    optimization_params,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//
//  for (int j = 0; j < face_to_patch_indices[50].size(); ++j)
//  {
//    int patch_index = face_to_patch_indices[50][j];
//    std::vector<PlanarPoint> domain_points;
//    spline_surface.get_patch(patch_index).get_domain().sample(
//      5,
//      domain_points
//    );
//    for (size_t i = 0; i < domain_points.size(); ++i)
//    {
//      SpatialVector param_point;
//      SpatialVector surface_point;
//      param_surface.evaluate_patch(patch_index, domain_points[i], param_point);
//      spline_surface.evaluate_patch(patch_index, domain_points[i], surface_point);
//      REQUIRE ( float_equal(param_point[0], surface_point[0]) );
//      REQUIRE ( float_equal(param_point[1], surface_point[1]) );
//      double u = param_point[0];
//      double v = param_point[1];
//      double f = 0.5 * u_curvature * u * u + 0.5 * v_curvature * v * v;
//      spdlog::info(
//        "Evaluated surface at ({}, {}) with actual value {} and expected value {}",
//        u,
//        v,
//        surface_point[2],
//        f
//      );
//      if ( !float_equal(surface_point[2], f) ) return false;
//    }
//  }
//
//  return true;
//}

TEST_CASE ( "Powell Sabin can reproduce quadratic surfaces", "[powell_sabin_spline]")
{
  std::vector< std::vector<SpatialVector> > param_grid;
  Eigen::MatrixXd param_V;
  Eigen::MatrixXi param_F;
  std::vector<std::vector<double>> param_l;
  std::vector< std::vector<SpatialVector> > control_point_grid;
  Eigen::MatrixXd control_V;
  Eigen::MatrixXi control_F;
  std::vector<std::vector<double>> l;

  SECTION ( "Linear" )
  {
//    int resolution = 10;
//    double u_slope = 1.0;
//    double v_slope = 2.0;
//    std::vector<std::vector<int>> face_to_patch_indices;
//    std::vector<int> patch_to_face_indices;
//    generate_plane_grid(resolution, u_slope, v_slope, control_point_grid);
//    generate_mesh_from_grid(control_point_grid, control_V, control_F, l, true);
//    generate_plane_grid(resolution, 0.0, 0.0, param_grid);
//    generate_mesh_from_grid(param_grid, param_V, param_F, param_l, true);
//    OptimizationParameters optimization_params;
//    optimization_params.method == OptimizationParameters::least_squares;
//    QuadraticSplineSurface spline_surface(
//      control_V,
//      control_F,
//      l,
//      quadratic_spline_params,
//      optimization_params,
//      face_to_patch_indices,
//      patch_to_face_indices
//    );
//    QuadraticSplineSurface param_surface(
//      param_V,
//      param_F,
//      param_l,
//      quadratic_spline_params,
//      optimization_params,
//      face_to_patch_indices,
//      patch_to_face_indices
//    );
//
//    int patch_index = face_to_patch_indices[30][0];
//    std::vector<PlanarPoint> domain_points;
//    spline_surface.get_patch(patch_index).get_domain().sample(
//      5,
//      domain_points
//    );
//    for (size_t i = 0; i < domain_points.size(); ++i)
//    {
//      SpatialVector param_point;
//      SpatialVector surface_point;
//      param_surface.evaluate_patch(patch_index, domain_points[i], param_point);
//      spline_surface.evaluate_patch(patch_index, domain_points[i], surface_point);
//      REQUIRE ( float_equal(param_point[0], surface_point[0]) );
//      REQUIRE ( float_equal(param_point[1], surface_point[1]) );
//      double u = param_point[0];
//      double v = param_point[1];
//      double f = u_slope * u + v_slope * v;
//      spdlog::info(
//        "Evaluated surface at ({}, {}) with actual value {} and expected value {}",
//        u,
//        v,
//        surface_point[2],
//        f
//      );
//      REQUIRE ( float_equal(surface_point[2], f) );
//    }
  }


  SECTION ( "Quadratic" )
  {
//    REQUIRE( test_quadratic_reproduction(0.0, 0.0) );
//    REQUIRE( test_quadratic_reproduction(1.0, 0.0) );
//    REQUIRE( test_quadratic_reproduction(0.0, 1.0) );
//    REQUIRE( test_quadratic_reproduction(1.0, 1.0) );
//    REQUIRE( test_quadratic_reproduction(1.0, 2.0) );
//    REQUIRE( test_quadratic_reproduction(-2.0, 1.0) );
//    REQUIRE( test_quadratic_reproduction(-2.0, -1.0) );
  }
}
