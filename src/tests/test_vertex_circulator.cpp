#include <catch2/catch_test_macros.hpp>
#include "common.h"
#include <vector>
#include <iostream>
#include <Eigen/Core>

#include "generate_shapes.h"
#include "vertex_circulator.h"


TEST_CASE ( "One rings of vertices can be found", "[vertex_circulator]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  SECTION ( "Tetrahedron" )
  {
    generate_tetrahedron_VF(V, F);
    VertexCirculator vertex_circulator(F);
    std::vector<int> vertex_one_ring;
    std::vector<int> face_one_ring;
    vertex_circulator.get_one_ring(0, vertex_one_ring, face_one_ring);


    REQUIRE ( vertex_one_ring.size() == 4 );
    REQUIRE ( face_one_ring.size() == 3 );
    REQUIRE ( vector_contains(vertex_one_ring, 1) );
    REQUIRE ( vector_contains(vertex_one_ring, 2) );
    REQUIRE ( vector_contains(vertex_one_ring, 3) );
    REQUIRE ( vector_contains(face_one_ring, 0) );
    REQUIRE ( vector_contains(face_one_ring, 1) );
    REQUIRE ( vector_contains(face_one_ring, 2) );
  }

  SECTION ( "Torus" )
  {
    generate_minimal_torus_VF(V, F);
    VertexCirculator vertex_circulator(F);
    std::vector<int> vertex_one_ring;
    std::vector<int> face_one_ring;
    vertex_circulator.get_one_ring(4, vertex_one_ring, face_one_ring);


    REQUIRE ( vertex_one_ring.size() == 7 );
    REQUIRE ( face_one_ring.size() == 6 );
    REQUIRE ( vector_contains(vertex_one_ring, 1) );
    REQUIRE ( vector_contains(vertex_one_ring, 2) );
    REQUIRE ( vector_contains(vertex_one_ring, 3) );
    REQUIRE ( vector_contains(vertex_one_ring, 5) );
    REQUIRE ( vector_contains(vertex_one_ring, 6) );
    REQUIRE ( vector_contains(vertex_one_ring, 7) );
    REQUIRE ( vector_contains(face_one_ring, 1) );
    REQUIRE ( vector_contains(face_one_ring, 2) );
    REQUIRE ( vector_contains(face_one_ring, 3) );
    REQUIRE ( vector_contains(face_one_ring, 6) );
    REQUIRE ( vector_contains(face_one_ring, 7) );
    REQUIRE ( vector_contains(face_one_ring, 8) );
  }

  SECTION ( "Triangle" )
  {
    F.resize(1, 3);
    F <<
      0, 1, 2;

    VertexCirculator vertex_circulator(F);
    std::vector<int> vertex_one_ring;
    std::vector<int> face_one_ring;

    // First vertex
    vertex_circulator.get_one_ring(0, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 2 );
    REQUIRE ( face_one_ring.size() == 1 );
    REQUIRE ( vertex_one_ring[0] == 1 );
    REQUIRE ( vertex_one_ring[1] == 2 );
    REQUIRE ( face_one_ring[0] == 0 );

    // Second vertex
    vertex_circulator.get_one_ring(1, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 2 );
    REQUIRE ( face_one_ring.size() == 1 );
    REQUIRE ( vertex_one_ring[0] == 2 );
    REQUIRE ( vertex_one_ring[1] == 0 );
    REQUIRE ( face_one_ring[0] == 0 );

    // Third vertex
    vertex_circulator.get_one_ring(2, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 2 );
    REQUIRE ( face_one_ring.size() == 1 );
    REQUIRE ( vertex_one_ring[0] == 0 );
    REQUIRE ( vertex_one_ring[1] == 1 );
    REQUIRE ( face_one_ring[0] == 0 );
  }

  SECTION ( "Square" )
  {
    F.resize(2, 3);
    F <<
      0, 1, 2,
      1, 3, 2;

    VertexCirculator vertex_circulator(F);
    std::vector<int> vertex_one_ring;
    std::vector<int> face_one_ring;

    // First vertex
    vertex_circulator.get_one_ring(0, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 2 );
    REQUIRE ( face_one_ring.size() == 1 );
    REQUIRE ( vertex_one_ring[0] == 1 );
    REQUIRE ( vertex_one_ring[1] == 2 );
    REQUIRE ( face_one_ring[0] == 0 );

    // Second vertex
    vertex_circulator.get_one_ring(1, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 3 );
    REQUIRE ( face_one_ring.size() == 2 );
    REQUIRE ( vertex_one_ring[0] == 3 );
    REQUIRE ( vertex_one_ring[1] == 2 );
    REQUIRE ( vertex_one_ring[2] == 0 );
    REQUIRE ( face_one_ring[0] == 1 );
    REQUIRE ( face_one_ring[1] == 0 );

    // Third vertex
    vertex_circulator.get_one_ring(2, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 3 );
    REQUIRE ( face_one_ring.size() == 2 );
    REQUIRE ( vertex_one_ring[0] == 0 );
    REQUIRE ( vertex_one_ring[1] == 1 );
    REQUIRE ( vertex_one_ring[2] == 3 );
    REQUIRE ( face_one_ring[0] == 0 );
    REQUIRE ( face_one_ring[1] == 1 );

    // Third vertex
    vertex_circulator.get_one_ring(3, vertex_one_ring, face_one_ring);
    REQUIRE ( vertex_one_ring.size() == 2 );
    REQUIRE ( face_one_ring.size() == 1 );
    REQUIRE ( vertex_one_ring[0] == 2 );
    REQUIRE ( vertex_one_ring[1] == 1 );
    REQUIRE ( face_one_ring[0] == 1 );
  }
}