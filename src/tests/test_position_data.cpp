#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "generate_shapes.h"
#include "position_data.h"


TEST_CASE ( "Gradients can be found", "[position_data]")
{
  Eigen::MatrixXd V(4, 3);
  Eigen::MatrixXi F(3, 3);
  int vertex_index = 0;
  std::vector<int> vertex_one_ring = { 1, 2, 3};
  std::vector<int> face_one_ring = { 0, 1, 2};
  MatrixXr one_ring_uv_positions(3, 2);
  Matrix2x3r gradient;
  F <<
    0, 1, 2,
    0, 2, 3,
    0, 3, 1;

  one_ring_uv_positions <<
    1.0, 0.0,
    -std::sqrt(3) / 2.0, 0.5,
    -std::sqrt(3) / 2.0, -0.5;

  SECTION ( "Constant" )
  {
    V <<
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0;

    compute_least_squares_vertex_gradient(
      V,
      vertex_index,
      vertex_one_ring,
      one_ring_uv_positions,
      gradient
    );
    spdlog::info(gradient);

  }

  SECTION ( "Linear" )
  {
    V <<
      0.0, 0.0, 0.0,
      1.0, 0.0, 1.0,
      -std::sqrt(3) / 2.0, 0.5, -std::sqrt(3) / 2.0,
      -std::sqrt(3) / 2.0, -0.5, -std::sqrt(3) / 2.0;

    compute_least_squares_vertex_gradient(
      V,
      vertex_index,
      vertex_one_ring,
      one_ring_uv_positions,
      gradient
    );
    spdlog::info(gradient);

  }
}
