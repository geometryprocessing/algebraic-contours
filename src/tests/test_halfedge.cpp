#include <catch2/catch_test_macros.hpp>
#include "common.h"

#include "generate_shapes.h"
#include "halfedge.h"


TEST_CASE ( "Halfedge can be built from VF", "[halfedge]")
{
  Eigen::MatrixXi F;

  SECTION ( "One triangle" )
  {
    F.resize(1, 3);
    F << 0, 1, 2;
    std::vector<std::vector<Halfedge::Index>> corner_to_he;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> he_to_corner;
    Halfedge mesh(F, corner_to_he, he_to_corner);

    // Check size information
    REQUIRE ( corner_to_he.size() == 1 );
    REQUIRE ( corner_to_he[0].size() == 3 );
    REQUIRE ( he_to_corner.size() == 3 );
    REQUIRE ( mesh.num_halfedges() == 3 );
    REQUIRE ( mesh.num_faces() == 1 );
    REQUIRE ( mesh.num_vertices() == 3 );
    REQUIRE ( mesh.num_edges() == 3 );

    // Get mesh elements
    size_t he0 = corner_to_he[0][0];
    size_t he1 = corner_to_he[0][1];
    size_t he2 = corner_to_he[0][2];
    size_t f0 = 0;
    size_t v0 = 0;
    size_t v1 = 1;
    size_t v2 = 2;

    // Check next 
    REQUIRE ( mesh.next_halfedge(he0) == he1 );
    REQUIRE ( mesh.next_halfedge(he1) == he2 );
    REQUIRE ( mesh.next_halfedge(he2) == he0 );

    // Check face
    REQUIRE ( mesh.halfedge_to_face(he0) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he1) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he2) == f0 );

    // Check vertex
    REQUIRE ( mesh.halfedge_to_tail_vertex(he0) == v1 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he0) == v2 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he1) == v2 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he1) == v0 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he2) == v0 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he2) == v1 );
  }

  SECTION ( "Two closed triangles" )
  {
    F.resize(2, 3);
    F << 0, 1, 2,
         0, 2, 1;
    std::vector<std::vector<Halfedge::Index>> corner_to_he;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> he_to_corner;
    Halfedge mesh(F, corner_to_he, he_to_corner);

    // Check size information
    REQUIRE ( corner_to_he.size() == 2 );
    REQUIRE ( corner_to_he[0].size() == 3 );
    REQUIRE ( corner_to_he[1].size() == 3 );
    REQUIRE ( he_to_corner.size() == 6 );
    REQUIRE ( mesh.num_halfedges() == 6 );
    REQUIRE ( mesh.num_faces() == 2 );
    REQUIRE ( mesh.num_vertices() == 3 );
    REQUIRE ( mesh.num_edges() == 3 );

    // Get mesh elements
    // Halfedges are indexed by face and global vertex index
    size_t he00 = corner_to_he[0][0];
    size_t he01 = corner_to_he[0][1];
    size_t he02 = corner_to_he[0][2];
    size_t he10 = corner_to_he[1][0];
    size_t he11 = corner_to_he[1][2];
    size_t he12 = corner_to_he[1][1];
    size_t f0 = 0;
    size_t f1 = 1;
    size_t v0 = 0;
    size_t v1 = 1;
    size_t v2 = 2;

    // Check next 
    REQUIRE ( mesh.next_halfedge(he00) == he01 );
    REQUIRE ( mesh.next_halfedge(he01) == he02 );
    REQUIRE ( mesh.next_halfedge(he02) == he00 );
    REQUIRE ( mesh.next_halfedge(he10) == he12 );
    REQUIRE ( mesh.next_halfedge(he11) == he10 );
    REQUIRE ( mesh.next_halfedge(he12) == he11 );

    // Check opposite 
    REQUIRE ( mesh.opposite_halfedge(he00) == he10 );
    REQUIRE ( mesh.opposite_halfedge(he01) == he11 );
    REQUIRE ( mesh.opposite_halfedge(he02) == he12 );
    REQUIRE ( mesh.opposite_halfedge(he10) == he00 );
    REQUIRE ( mesh.opposite_halfedge(he11) == he01 );
    REQUIRE ( mesh.opposite_halfedge(he12) == he02 );

    // Check face
    REQUIRE ( mesh.halfedge_to_face(he00) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he01) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he02) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he10) == f1 );
    REQUIRE ( mesh.halfedge_to_face(he11) == f1 );
    REQUIRE ( mesh.halfedge_to_face(he12) == f1 );

    // Check vertex
    REQUIRE ( mesh.halfedge_to_tail_vertex(he00) == v1 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he00) == v2 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he01) == v2 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he01) == v0 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he02) == v0 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he02) == v1 );

    REQUIRE ( mesh.halfedge_to_tail_vertex(he10) == v2 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he10) == v1 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he11) == v0 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he11) == v2 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he12) == v1 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he12) == v0 );
  }

  SECTION ( "Two open triangles" )
  {
    F.resize(2, 3);
    F << 0, 1, 2,
         0, 3, 1;
    std::vector<std::vector<Halfedge::Index>> corner_to_he;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> he_to_corner;
    Halfedge mesh(F, corner_to_he, he_to_corner);

    // Check size information
    REQUIRE ( corner_to_he.size() == 2 );
    REQUIRE ( corner_to_he[0].size() == 3 );
    REQUIRE ( corner_to_he[1].size() == 3 );
    REQUIRE ( he_to_corner.size() == 6 );
    REQUIRE ( mesh.num_halfedges() == 6 );
    REQUIRE ( mesh.num_faces() == 2 );
    REQUIRE ( mesh.num_vertices() == 4 );
    REQUIRE ( mesh.num_edges() == 5 );

    // Get mesh elements
    // Halfedges are indexed by face and global vertex index
    size_t he00 = corner_to_he[0][0];
    size_t he01 = corner_to_he[0][1];
    size_t he02 = corner_to_he[0][2];
    size_t he10 = corner_to_he[1][0];
    size_t he11 = corner_to_he[1][2];
    size_t he13 = corner_to_he[1][1];
    size_t f0 = 0;
    size_t f1 = 1;
    size_t v0 = 0;
    size_t v1 = 1;
    size_t v2 = 2;
    size_t v3 = 3;

    // Check next 
    REQUIRE ( mesh.next_halfedge(he00) == he01 );
    REQUIRE ( mesh.next_halfedge(he01) == he02 );
    REQUIRE ( mesh.next_halfedge(he02) == he00 );
    REQUIRE ( mesh.next_halfedge(he10) == he13 );
    REQUIRE ( mesh.next_halfedge(he11) == he10 );
    REQUIRE ( mesh.next_halfedge(he13) == he11 );

    // Check opposite 
    REQUIRE ( mesh.opposite_halfedge(he02) == he13 );
    REQUIRE ( mesh.opposite_halfedge(he13) == he02 );

    // Check face
    REQUIRE ( mesh.halfedge_to_face(he00) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he01) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he02) == f0 );
    REQUIRE ( mesh.halfedge_to_face(he10) == f1 );
    REQUIRE ( mesh.halfedge_to_face(he11) == f1 );
    REQUIRE ( mesh.halfedge_to_face(he13) == f1 );

    // Check vertex
    REQUIRE ( mesh.halfedge_to_tail_vertex(he00) == v1 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he00) == v2 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he01) == v2 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he01) == v0 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he02) == v0 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he02) == v1 );

    REQUIRE ( mesh.halfedge_to_tail_vertex(he10) == v3 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he10) == v1 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he11) == v0 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he11) == v3 );
    REQUIRE ( mesh.halfedge_to_tail_vertex(he13) == v1 );
    REQUIRE ( mesh.halfedge_to_head_vertex(he13) == v0 );
  }
}

