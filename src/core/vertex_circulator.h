// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// \file vertex_circulator.h
///
/// Class to build circulators around vertices in VF representation

class VertexCirculator
{
public:
  /// @brief Constructor for the vertex circulator from the faces of the mesh.
  ///
  /// @param[in] F: input mesh faces
  VertexCirculator(const Eigen::MatrixXi& F);

  /// @brief Get the one ring of a vertex.
  ///
  /// The one ring of both faces and vertices counter clockwise around the
  /// vertex are returned. For boundary vertices, the faces and vertices start
  /// at the right boundary and traverse the faces in order to the left
  /// boundary. For interior vertices, an arbitrary start face is chosen, and
  /// the vertex one ring is closed so that v_0 = v_n.
  ///
  /// @param[in] vertex_index: index of the vertex to get the one ring for
  /// @param[out] vertex_one_ring: vertices ccw around the one ring
  /// @param[out] face_one_ring: faces ccw around the one ring
  void get_one_ring(int vertex_index,
                    std::vector<int>& vertex_one_ring,
                    std::vector<int>& face_one_ring) const;

private:
  Eigen::MatrixXi m_F;
  std::vector<std::vector<int>> m_all_adjacent_faces;
  std::vector<std::vector<int>> m_all_vertex_one_rings;
  std::vector<std::vector<int>> m_all_face_one_rings;
};
