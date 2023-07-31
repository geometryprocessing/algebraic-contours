// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "vertex_circulator.h"

// Return true iff the face contains the edge { vertex_0, vertex_1 }
bool
contains_edge(Eigen::VectorXi face, int vertex_0, int vertex_1)
{
  return ((contains_vertex(face, vertex_0)) &&
          (contains_vertex(face, vertex_1)));
}

// Return true iff the face is to the left of the given edge
bool
is_left_face(Eigen::VectorXi face, int vertex_0, int vertex_1)
{
  if ((face[0] == vertex_0) && (face[1] == vertex_1))
    return true;
  if ((face[1] == vertex_0) && (face[2] == vertex_1))
    return true;
  if ((face[2] == vertex_0) && (face[0] == vertex_1))
    return true;

  return false;
}

// Return true iff the face is to the right of the given edge
bool
is_right_face(Eigen::VectorXi face, int vertex_0, int vertex_1)
{
  if ((face[1] == vertex_0) && (face[0] == vertex_1))
    return true;
  if ((face[2] == vertex_0) && (face[1] == vertex_1))
    return true;
  if ((face[0] == vertex_0) && (face[2] == vertex_1))
    return true;

  return false;
}

// Get the index of the vertex in the face ccw from the given vertex
int
find_next_vertex(Eigen::VectorXi face, int vertex)
{
  if (face[0] == vertex)
    return face[1];
  if (face[1] == vertex)
    return face[2];
  if (face[2] == vertex)
    return face[0];

  return -1;
}

// Get the index of the vertex in the face clockwise from the given vertex
int
find_prev_vertex(Eigen::VectorXi face, int vertex)
{
  if (face[0] == vertex)
    return face[2];
  if (face[1] == vertex)
    return face[0];
  if (face[2] == vertex)
    return face[1];

  return -1;
}

// Return true iff the two faces are adjacent
bool
are_adjacent(Eigen::VectorXi face_0, Eigen::VectorXi face_1)
{
  if (contains_edge(face_0, face_1[0], face_1[1]))
    return true;
  if (contains_edge(face_0, face_1[1], face_1[2]))
    return true;
  if (contains_edge(face_0, face_1[2], face_1[0]))
    return true;
  return false;
}

// Get list of all faces adjacent to each vertex
void
compute_adjacent_faces(const Eigen::MatrixXi& F,
                       std::vector<std::vector<int>>& all_adjacent_faces)
{
  // Initialize adjacent faces list
  size_t num_vertices = F.maxCoeff() + 1;
  all_adjacent_faces.resize(num_vertices);
  for (size_t i = 0; i < all_adjacent_faces.size(); ++i) {
    all_adjacent_faces[i].clear();
  }

  for (Eigen::Index i = 0; i < F.rows(); ++i) {
    for (Eigen::Index j = 0; j < F.cols(); ++j) {
      all_adjacent_faces[F(i, j)].push_back(i);
    }
  }
}

// Compute the first face of the vertex one ring, which should be right
// boundary face for a boundary vertex.
int
compute_vertex_one_ring_first_face(const Eigen::MatrixXi& F,
                                   int vertex_index,
                                   const std::vector<int>& adjacent_faces)
{
  if (adjacent_faces.empty())
    return -1;

  // Get arbitrary adjacent face to start and vertex on the face
  int current_face = adjacent_faces[0];
  int current_vertex = find_next_vertex(F.row(current_face), vertex_index);
  SPDLOG_TRACE("Starting search for first face from vertex {} on face {}",
               current_vertex,
               F.row(current_face));

  // Cycle clockwise to a starting face
  for (size_t i = 1; i < adjacent_faces.size(); ++i) {
    // Get previous face or return if none exists

    int prev_face = -1;
    for (size_t j = 0; j < adjacent_faces.size(); ++j) {
      int f = adjacent_faces[j];
      if (is_right_face(F.row(f), vertex_index, current_vertex)) {
        prev_face = f;
        break;
      }
    }

    // Return current face if no previous face found
    if (prev_face == -1) {
      return current_face;
    }

    // Get previous face and vertex
    current_face = prev_face;
    current_vertex = find_prev_vertex(F.row(current_face), current_vertex);
  }

  // If we have not returned yet, this is an interior vertex, and we return
  // the current face as an arbitrary choice
  return current_face;
}

// Compute the vertex one ring for a vertex index using adjacent faces.
void
compute_vertex_one_ring(const Eigen::MatrixXi& F,
                        int vertex_index,
                        const std::vector<int>& adjacent_faces,
                        std::vector<int>& vertex_one_ring,
                        std::vector<int>& face_one_ring)
{
  size_t num_faces = adjacent_faces.size();
  vertex_one_ring.resize(num_faces + 1);
  face_one_ring.resize(num_faces);
  if (adjacent_faces.empty())
    return;

  // Get first face and vertex
  face_one_ring[0] =
    compute_vertex_one_ring_first_face(F, vertex_index, adjacent_faces);
  vertex_one_ring[0] = find_next_vertex(F.row(face_one_ring[0]), vertex_index);

  // Get remaining one ring faces and vertices
  for (size_t i = 1; i < num_faces; ++i) {
    // Get next vertex
    vertex_one_ring[i] =
      find_next_vertex(F.row(face_one_ring[i - 1]), vertex_one_ring[i - 1]);

    // Get next face
    for (size_t j = 0; j < num_faces; ++j) {
      int f = adjacent_faces[j];
      if (is_left_face(F.row(f), vertex_index, vertex_one_ring[i])) {
        face_one_ring[i] = f;
      }
    }
  }

  // Get final vertex (same as first for closed loop)
  SPDLOG_TRACE("Adding last vertex for face {} from vertex {}",
               F.row(face_one_ring[num_faces - 1]),
               vertex_one_ring[num_faces - 1]);
  vertex_one_ring[num_faces] = find_next_vertex(
    F.row(face_one_ring[num_faces - 1]), vertex_one_ring[num_faces - 1]);
  spdlog::trace("Last vertex: {}", vertex_one_ring[num_faces]);
}

VertexCirculator::VertexCirculator(const Eigen::MatrixXi& F)
  : m_F(F)
{
  // Initialize adjacent faces list
  size_t num_vertices = F.maxCoeff() + 1;
  compute_adjacent_faces(F, m_all_adjacent_faces);

  // Compute face and vertex one rings
  m_all_vertex_one_rings.resize(num_vertices);
  m_all_face_one_rings.resize(num_vertices);
  for (size_t i = 0; i < num_vertices; ++i) {
    compute_vertex_one_ring(F,
                            i,
                            m_all_adjacent_faces[i],
                            m_all_vertex_one_rings[i],
                            m_all_face_one_rings[i]);
  }
}

void
VertexCirculator::get_one_ring(int vertex_index,
                               std::vector<int>& vertex_one_ring,
                               std::vector<int>& face_one_ring) const
{
  vertex_one_ring = m_all_vertex_one_rings[vertex_index];
  face_one_ring = m_all_face_one_rings[vertex_index];
}