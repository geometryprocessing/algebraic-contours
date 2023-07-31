// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "halfedge.h"

void
compute_face_boundary_edges(
  const Eigen::MatrixXi& F,
  std::vector<std::pair<int, int>>& face_boundary_edges)
{
  // Build halfedge for the mesh
  std::vector<std::vector<Halfedge::Index>> corner_to_he;
  std::vector<std::pair<Eigen::Index, Eigen::Index>> he_to_corner;
  Halfedge halfedge(F, corner_to_he, he_to_corner);

  // Get boundary halfedges
  std::vector<Halfedge::Index> boundary_halfedges;
  halfedge.build_boundary_halfedge_list(boundary_halfedges);

  // Get boundary face corners opposite halfedges
  face_boundary_edges.reserve(boundary_halfedges.size());
  for (size_t i = 0; i < boundary_halfedges.size(); ++i) {
    face_boundary_edges.push_back(he_to_corner[boundary_halfedges[i]]);
  }
}

void
compute_boundary_vertices(const Eigen::MatrixXi& F,
                          std::vector<int>& boundary_vertices)
{
  // Get face boundary edges
  std::vector<std::pair<int, int>> face_boundary_edges;
  compute_face_boundary_edges(F, face_boundary_edges);

  // Get boolean array of boundary indices
  size_t num_vertices = F.maxCoeff() + 1;
  std::vector<bool> is_boundary_vertex(num_vertices, false);
  for (size_t i = 0; i < face_boundary_edges.size(); ++i) {
    // Mark boundary edge endpoints as boundary vertices
    int face_index = face_boundary_edges[i].first;
    int face_vertex_index = face_boundary_edges[i].second;
    is_boundary_vertex[F(face_index, (face_vertex_index + 1) % 3)] = true;
    is_boundary_vertex[F(face_index, (face_vertex_index + 2) % 3)] = true;
  }

  // Convert boolean array to index vector
  std::vector<size_t> unsigned_boundary_vertices;
  convert_boolean_array_to_index_vector(is_boundary_vertex,
                                        unsigned_boundary_vertices);
  convert_unsigned_vector_to_signed(unsigned_boundary_vertices,
                                    boundary_vertices);
}

