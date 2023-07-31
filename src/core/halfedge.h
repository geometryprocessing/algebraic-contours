// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// \file halfedge.h
///
/// Class to build halfedge from VF

// TODO The cleanest way to handle this is to fill boundaries with faces
// and handle boundary cases with these to avoid invalid operations. The
// interface should be chosen with care to balance elegance and versatility.

/// Mesh halfedge representation. Supports meshes with boundary and basic
/// topological information. Can be initialized from face topology information.
class Halfedge
{
public:
  typedef int Index;

  /// Default trivial halfedge
  Halfedge() { clear(); }

  /// Build halfedge mesh from mesh faces F with
  Halfedge(const Eigen::MatrixXi& F,
           std::vector<std::vector<Index>>& corner_to_he,
           std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner);

  // **************
  // Element counts
  // **************

  Index num_halfedges() const { return m_num_halfedges; }

  Index num_faces() const { return m_num_faces; }

  Index num_vertices() const { return m_num_vertices; };

  Index num_edges() const { return m_num_edges; };

  // *********
  // Adjacency
  // *********

  Index next_halfedge(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_halfedge_index();
    return m_next[he];
  }

  Index opposite_halfedge(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_halfedge_index();
    return m_opp[he];
  }

  Index halfedge_to_face(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_face_index();
    return m_face[he];
  }

  Index halfedge_to_head_vertex(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_vertex_index();
    return m_to[he];
  }

  Index halfedge_to_tail_vertex(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_vertex_index();
    return m_from[he];
  }

  // ********************
  // Edge Representations
  // ********************

  Index halfedge_to_edge(Index he) const
  {
    if (!is_valid_halfedge_index(he))
      return invalid_edge_index();
    return m_he2e[he];
  }

  std::pair<Index, Index> edge_to_halfedge(Index e) const { return m_e2he[e]; }

  Index edge_to_first_halfedge(Index e) const { return m_e2he[e].first; }

  Index edge_to_second_halfedge(Index e) const { return m_e2he[e].second; }

  std::vector<Index> const& get_halfedge_to_edge_map() const { return m_he2e; }

  std::vector<std::pair<Index, Index>> const& get_edge_to_halfedge_map() const
  {
    return m_e2he;
  }

  // ******************
  // Element predicates
  // ******************

  bool is_boundary_edge(Index e) const
  {
    if (!is_valid_halfedge_index(edge_to_first_halfedge(e)))
      return true;
    if (!is_valid_halfedge_index(edge_to_second_halfedge(e)))
      return true;
    return false;
  }

  bool is_boundary_halfedge(Index he) const
  {
    return is_boundary_edge(halfedge_to_edge(he));
  }

  void build_boundary_edge_list(std::vector<Index>& boundary_edges) const
  {
    boundary_edges.reserve(num_edges());
    for (Index ei = 0; ei < num_edges(); ++ei) {
      if (is_boundary_edge(ei)) {
        boundary_edges.push_back(ei);
      }
    }
  }

  void build_boundary_halfedge_list(
    std::vector<Index>& boundary_halfedges) const
  {
    boundary_halfedges.reserve(num_halfedges());
    for (Index hi = 0; hi < num_halfedges(); ++hi) {
      if (is_boundary_halfedge(hi)) {
        boundary_halfedges.push_back(hi);
      }
    }
  }

  void clear()
  {
    m_next.clear();
    m_opp.clear();
    m_he2e.clear();
    m_e2he.clear();
    m_to.clear();
    m_from.clear();
    m_face.clear();
    m_out.clear();
    m_f2he.clear();
    m_F.resize(0, 0);
  }

private:
  std::vector<Index> m_next;
  std::vector<Index> m_opp;
  std::vector<Index> m_he2e;
  std::vector<std::pair<Index, Index>> m_e2he;
  std::vector<Index> m_to;
  std::vector<Index> m_from;
  std::vector<Index> m_face;
  std::vector<Index> m_out;
  std::vector<Index> m_f2he;
  Eigen::MatrixXi m_F;
  Index m_num_vertices;
  Index m_num_faces;
  Index m_num_halfedges;
  Index m_num_edges;

  void build_corner_to_he_maps(
    Eigen::Index num_faces,
    std::vector<std::vector<Index>>& corner_to_he,
    std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner) const
  {
    corner_to_he.resize(num_faces);
    he_to_corner.resize(3 * num_faces);

    // Iterate over faces to build corner to he maps
    Index he_index = 0;
    for (Eigen::Index face_index = 0; face_index < num_faces; ++face_index) {
      corner_to_he[face_index].resize(3);
      for (Eigen::Index i = 0; i < 3; ++i) {
        // Assign indices
        corner_to_he[face_index][i] = he_index;
        he_to_corner[he_index] = std::make_pair(face_index, i);

        // Update current face index
        he_index++;
      }
    }
  }

  // *********************
  // Index validity checks
  // *********************

  bool is_valid_halfedge_index(Index he) const
  {
    if (he < 0)
      return false;
    if (he >= num_halfedges())
      return false;
    return true;
  }

  bool is_valid_vertex_index(Index vertex_index) const
  {
    if (vertex_index < 0)
      return false;
    if (vertex_index >= num_vertices())
      return false;
    return true;
  }

  bool is_valid_face_index(Index face_index) const
  {
    if (face_index < 0)
      return false;
    if (face_index >= num_faces())
      return false;
    return true;
  }

  bool is_valid_edge_index(Index edge_index) const
  {
    if (edge_index < 0)
      return false;
    if (edge_index >= num_edges())
      return false;
    return true;
  }

  // **************************
  // Invalid index constructors
  // **************************

  Index invalid_halfedge_index() const { return -1; }

  Index invalid_vertex_index() const { return -1; }

  Index invalid_face_index() const { return -1; }

  Index invalid_edge_index() const { return -1; }

  // **********************
  // Global validity checks
  // **********************

  bool is_valid() const
  {
    // Check sizes
    if (static_cast<Index>(m_next.size()) != num_halfedges()) {
      spdlog::error("next domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_opp.size()) != num_halfedges()) {
      spdlog::error("opp domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_he2e.size()) != num_halfedges()) {
      spdlog::error("he2e domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_to.size()) != num_halfedges()) {
      spdlog::error("to domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_from.size()) != num_halfedges()) {
      spdlog::error("from domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_face.size()) != num_halfedges()) {
      spdlog::error("face domain not in bijection with halfedges");
      return false;
    }
    if (static_cast<Index>(m_e2he.size()) != num_edges()) {
      spdlog::error("e2he domain not in bijection with edges");
      return false;
    }
    if (static_cast<Index>(m_out.size()) != num_vertices()) {
      spdlog::error("out domain not in bijection with vertices");
      return false;
    }
    if (static_cast<Index>(m_f2he.size()) != num_faces()) {
      spdlog::error("f2he domain not in bijection with faces");
      return false;
    }
    if (static_cast<Index>(m_F.rows()) != num_faces()) {
      spdlog::error("F rows not in bijection with faces");
      return false;
    }

    // Check opp
    // TODO: This requires a more robust halfedge that correctly assigns
    // faces to fill holes so, e.g., opp[opp[he]] == he in all cases
    for (Index he = 0; he < num_halfedges(); ++he) {
      // if (m_opp[m_opp[he]] != he)
      //{
      //   spdlog::error("Broken edge opp pair {} and {}", he, m_opp[he]);
      //   return false;
      // }
    }

    return true;
  }
};

/// Build a map from halfedges to edges and the inverse from opposite halfedge
/// information.
///
/// @param[in] opp: array mapping halfedges to their opposite halfedge pair
/// @param[out] he2e: array mapping halfedges to edge indices
/// @param[out] e2he: array mapping edges to their pair of corresponding
/// halfedge indices
void
build_halfedge_to_edge_maps(
  const std::vector<Halfedge::Index>& opp,
  std::vector<Halfedge::Index>& he2e,
  std::vector<std::pair<Halfedge::Index, Halfedge::Index>>& e2he);