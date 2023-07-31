#include <catch2/catch_test_macros.hpp>
#include "common.h"

#include "affine_manifold.h"
#include "generate_shapes.h"


TEST_CASE ( "An affine manifold can be built for a triangle" )
{
  Eigen::MatrixXd V(3, 2);
  Eigen::MatrixXi F(1, 3);
  V <<
    0.0, 0.0,
    3.0, 0.0,
    0.0, 4.0;
  F <<
    0, 1, 2;
  double width = 3;
  double height = 4;
  std::vector<std::vector<double>> l(1);
  l[0] = { 5, 4, 3 };

// REQUIRE GLOBAL UV
//  SECTION ( "Build from lengths" )
//  {
//    AffineManifold affine_manifold(F, l);
//
//    // Check basic manifold information
//    REQUIRE( affine_manifold.num_faces() == 1 );
//    REQUIRE( affine_manifold.num_vertices() == 3 );
//
//    // Check vertex chart at vertex 0
//    VertexManifoldChart chart = affine_manifold.get_vertex_chart(0);
//    REQUIRE( chart.vertex_index == 0 );
//    REQUIRE( chart.vertex_one_ring.size() == 2 );
//    REQUIRE( chart.vertex_one_ring[0] == 1 );
//    REQUIRE( chart.vertex_one_ring[1] == 2 );
//    REQUIRE( chart.face_one_ring.size() == 1 );
//    REQUIRE( chart.face_one_ring[0] == 0 );
//    REQUIRE( float_equal(chart.one_ring_uv_positions(0, 0), 3.0) );
//    REQUIRE( float_equal(chart.one_ring_uv_positions(0, 1), 0.0) );
//    REQUIRE( float_equal(chart.one_ring_uv_positions(1, 0), 0.0) );
//    REQUIRE( float_equal(chart.one_ring_uv_positions(1, 1), 4.0) );
//
//    // Check face corner charts for face 0
//    std::array<Matrix2x2r, 3> corner_uv_positions;
//    affine_manifold.get_face_corner_charts(0, corner_uv_positions);
//    REQUIRE( float_equal(corner_uv_positions[0](0, 0), 3.0) );
//    REQUIRE( float_equal(corner_uv_positions[0](0, 1), 0.0) );
//    REQUIRE( float_equal(corner_uv_positions[0](1, 0), 0.0) );
//    REQUIRE( float_equal(corner_uv_positions[0](1, 1), 4.0) );
//    REQUIRE( float_equal(corner_uv_positions[1](0, 0), 5.0) );
//    REQUIRE( float_equal(corner_uv_positions[1](0, 1), 0.0) );
//    REQUIRE( float_equal(corner_uv_positions[2](0, 0), 4.0) );
//    REQUIRE( float_equal(corner_uv_positions[2](0, 1), 0.0) );
//
//    // Check edge chart for vertex 2 in face 0
//    EdgeManifoldChart edge_chart = affine_manifold.get_edge_chart(0, 2);
//    assert( edge_chart.top_face_index == 0 );
//    assert( edge_chart.left_vertex_index == 0 );
//    assert( edge_chart.right_vertex_index == 1 );
//    assert( edge_chart.top_vertex_index == 2 );
//    assert( float_equal(edge_chart.left_vertex_uv_position[0], -0.5) );
//    assert( float_equal(edge_chart.left_vertex_uv_position[1], 0.0) );
//    assert( float_equal(edge_chart.right_vertex_uv_position[0], 0.5) );
//    assert( float_equal(edge_chart.right_vertex_uv_position[1], 0.0) );
//    assert( float_equal(edge_chart.top_vertex_uv_position[0], -0.5) );
//    assert( float_equal(edge_chart.top_vertex_uv_position[1], 4.0 / 3.0) );
//
//  }

  SECTION ( "Build from global uvs" )
  {
    ParametricAffineManifold affine_manifold(F, V);


    // Check basic manifold information
    REQUIRE( affine_manifold.num_faces() == 1 );
    REQUIRE( affine_manifold.num_vertices() == 3 );

    // Check vertex chart at vertex 0
    VertexManifoldChart chart = affine_manifold.get_vertex_chart(0);
    REQUIRE( chart.vertex_index == 0 );
    REQUIRE( chart.vertex_one_ring.size() == 2 );
    REQUIRE( chart.vertex_one_ring[0] == 1 );
    REQUIRE( chart.vertex_one_ring[1] == 2 );
    REQUIRE( chart.face_one_ring.size() == 1 );
    REQUIRE( chart.face_one_ring[0] == 0 );
    REQUIRE( float_equal(chart.one_ring_uv_positions(0, 0), 3.0) );
    REQUIRE( float_equal(chart.one_ring_uv_positions(0, 1), 0.0) );
    REQUIRE( float_equal(chart.one_ring_uv_positions(1, 0), 0.0) );
    REQUIRE( float_equal(chart.one_ring_uv_positions(1, 1), 4.0) );

    // Check face corner charts for face 0
    std::array<Matrix2x2r, 3> corner_uv_positions;
    affine_manifold.get_face_corner_charts(0, corner_uv_positions);
    REQUIRE( float_equal(corner_uv_positions[0](0, 0),  3.0) );
    REQUIRE( float_equal(corner_uv_positions[0](0, 1),  0.0) );
    REQUIRE( float_equal(corner_uv_positions[0](1, 0),  0.0) );
    REQUIRE( float_equal(corner_uv_positions[0](1, 1),  4.0) );
    REQUIRE( float_equal(corner_uv_positions[1](0, 0), -3.0) );
    REQUIRE( float_equal(corner_uv_positions[1](0, 1),  4.0) );
    REQUIRE( float_equal(corner_uv_positions[1](1, 0), -3.0) );
    REQUIRE( float_equal(corner_uv_positions[1](1, 1),  0.0) );
    REQUIRE( float_equal(corner_uv_positions[2](0, 0),  0.0) );
    REQUIRE( float_equal(corner_uv_positions[2](0, 1), -4.0) );
    REQUIRE( float_equal(corner_uv_positions[2](1, 0),  3.0) );
    REQUIRE( float_equal(corner_uv_positions[2](1, 1), -4.0) );

    // Check global uv
    PlanarPoint uv0, uv1, uv2;
    affine_manifold.get_vertex_global_uv(0, uv0);
    affine_manifold.get_vertex_global_uv(1, uv1);
    affine_manifold.get_vertex_global_uv(2, uv2);
    REQUIRE( matrix_equal(affine_manifold.get_global_uv(), V) );
    REQUIRE( vector_equal<2>(uv0, V.row(0)) );
    REQUIRE( vector_equal<2>(uv1, V.row(1)) );
    REQUIRE( vector_equal<2>(uv2, V.row(2)) );
  }
}