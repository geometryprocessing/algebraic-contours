// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "halfedge.h"

#include "vertex_circulator.h"
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>

void
build_halfedge_to_edge_maps(
  const std::vector<Halfedge::Index>& opp,
  std::vector<Halfedge::Index>& he2e,
  std::vector<std::pair<Halfedge::Index, Halfedge::Index>>& e2he)
{
  Halfedge::Index num_he = opp.size();
  he2e.resize(num_he);
  e2he.reserve(num_he);

  // Iterate over halfedges to build maps between halfedges and edges
  Halfedge::Index e_index = 0;
  for (Halfedge::Index he = 0; he < num_he; ++he) {
    // Check if the halfedge is on the boundary
    bool is_boundary = ((opp[he] < 0) || (opp[he] >= num_he));

    // Skip interior halfedges with lower index, but always process a boundary
    // halfedge
    if ((he >= opp[he]) || (is_boundary)) {
      // Build maps for the current edge
      e2he.push_back(std::make_pair(he, opp[he]));
      he2e[he] = e_index;

      // Only valid for interior edges
      if (!is_boundary) {
        he2e[opp[he]] = e_index;
      }

      // Update current edge index
      e_index++;
    }
  }
}

/// Build halfedge mesh from mesh faces F with
Halfedge::Halfedge(
  const Eigen::MatrixXi& F,
  std::vector<std::vector<Index>>& corner_to_he,
  std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner)
  : m_F(F)
{
  Eigen::Index num_faces = F.rows();
  int num_vertices = F.maxCoeff() + 1;
  int num_halfedges = 3 * num_faces;

  // Check input
#if CHECK_VALIDITY
  if (!is_manifold(F)) {
    spdlog::error("Input mesh is not manifold");
    clear();
    return;
  }
#endif

  // Build maps between corners and halfedges
  build_corner_to_he_maps(num_faces, corner_to_he, he_to_corner);

  // Iterate over faces to build next, face, and to arrays
  m_next.resize(num_halfedges, invalid_halfedge_index());
  m_face.resize(num_halfedges, invalid_face_index());
  m_to.resize(num_halfedges, invalid_vertex_index());
  m_from.resize(num_halfedges, invalid_vertex_index());
  for (Eigen::Index face_index = 0; face_index < num_faces; ++face_index) {
    for (Eigen::Index i = 0; i < 3; ++i) {
      Index current_he = corner_to_he[face_index][i];
      Index next_he = corner_to_he[face_index][(i + 1) % 3];
      m_next[current_he] = next_he;
      m_face[current_he] = face_index;
      m_to[current_he] = F(face_index, (i + 2) % 3);
      m_from[current_he] = F(face_index, (i + 1) % 3);
    }
  }

  // Build out and f2he arrays
  m_out.resize(num_vertices, -1);
  m_f2he.resize(num_faces, -1);
  for (Index he_index = 0; he_index < num_halfedges; ++he_index) {
    m_out[m_to[he_index]] = m_next[he_index];
    m_f2he[m_face[he_index]] = he_index;
  }

  // Iterate over vertices to build opp using a vertex circulator
  // Note that this is the main difficulty in constructing halfedge from VF
  VertexCirculator vertex_circulator(F);
  m_opp.resize(3 * num_faces, -1);
  for (Index vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
    // Get vertex one ring
    std::vector<int> vertex_one_ring;
    std::vector<int> face_one_ring;
    vertex_circulator.get_one_ring(
      vertex_index, vertex_one_ring, face_one_ring);

    // Determine if we are in a boundary case
    bool is_boundary = (vertex_one_ring.front() != vertex_one_ring.back());
    Index num_adjacent_faces = face_one_ring.size();
    Index num_interior_edges =
      is_boundary ? num_adjacent_faces - 1 : num_adjacent_faces;

    // Build opposite arrays
    for (Index i = 0; i < num_interior_edges; ++i) {
      // Get current face prev (cw) halfedge from the vertex
      Index fi = face_one_ring[i];
      Index fi_vertex_index = find_face_vertex_index(F.row(fi), vertex_index);
      Index current_he = corner_to_he[fi][(fi_vertex_index + 1) % 3];

      // Get next (ccw) face next (ccw) halfedge from the vertex
      Index fj = face_one_ring[(i + 1) % num_adjacent_faces];
      Index fj_vertex_index = find_face_vertex_index(F.row(fj), vertex_index);
      Index opposite_he = corner_to_he[fj][(fj_vertex_index + 2) % 3];

      // Assign opposite halfedge
      m_opp[current_he] = opposite_he;
    }
  }

  // Build maps between edges and halfedges
  build_halfedge_to_edge_maps(m_opp, m_he2e, m_e2he);

  // Set sizes
  m_num_halfedges = num_halfedges;
  m_num_faces = num_faces;
  m_num_vertices = num_vertices;
  m_num_edges = m_e2he.size();

  // Check validity
#if CHECK_VALIDITY
  if (!is_valid()) {
    spdlog::error("Could not build halfedge");
    clear();
    return;
  }
#endif
}