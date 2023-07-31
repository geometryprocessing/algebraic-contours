// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "powell_sabin_local_to_global_indexing.h"

// The variables we are optimizing need to be linearized and initialized with
// some values. Similarly, the final optimized results need to be extracted.
// Furthermore, for the optimization it is useful to combine the variable values
// into a single vector so that, e.g., gradient descent or Newton's method can
// be applied.

// *******************
// Local block indices
// *******************

// Get the local start index of the block of position variable indices for a
// given vertex
size_t
generate_local_vertex_position_variables_start_index(size_t vertex_index,
                                                     size_t dimension = 3)
{
  size_t relative_vertex_index = 3 * dimension * vertex_index;

  return relative_vertex_index;
}

// Get the local start index of the block of gradient variable indices for a
// given vertex
size_t
generate_local_vertex_gradient_variables_start_index(size_t vertex_index,
                                                     size_t dimension = 3)
{
  size_t relative_vertex_index = 3 * dimension * vertex_index;
  size_t position_block_size = dimension;

  return relative_vertex_index + position_block_size;
}

// Get the local start index of the block of gradient variable indices for a
// given edge
size_t
generate_local_edge_gradient_variables_start_index(size_t edge_index,
                                                   size_t dimension = 3)
{
  size_t vertex_block_size = 9 * dimension;
  size_t relative_edge_index = dimension * edge_index;

  return vertex_block_size + relative_edge_index;
}

// **********************
// Local variable indices
// **********************

// Get the local index of the position variable indices for a given coordinate
// and vertex index
size_t
generate_local_vertex_position_variable_index(size_t face_vertex_index,
                                              size_t coord,
                                              size_t dimension)
{
  size_t start_index = generate_local_vertex_position_variables_start_index(
    face_vertex_index, dimension);

  return start_index + coord;
}

// Get the local index of the gradient variable indices for a given matrix index
// pair and vertex index
size_t
generate_local_vertex_gradient_variable_index(size_t face_vertex_index,
                                              size_t row,
                                              size_t col,
                                              size_t dimension)
{
  size_t start_index = generate_local_vertex_gradient_variables_start_index(
    face_vertex_index, dimension);
  size_t matrix_index =
    generate_local_variable_matrix_index(row, col, dimension);

  return start_index + matrix_index;
}

// Get the local index of the gradient variable indices for a given coordinate
// and edge index pair
size_t
generate_local_edge_gradient_variable_index(size_t face_edge_index,
                                            size_t coord,
                                            size_t dimension)
{
  size_t start_index = generate_local_edge_gradient_variables_start_index(
    face_edge_index, dimension);

  return start_index + coord;
}


// ********************
// Global block indices
// ********************

// Get the start index of the block of vertex position variable indices
int
generate_global_vertex_position_variables_block_start_index()
{
  return 0;
}

// Get the start index of the block of vertex gradient variable indices
int
generate_global_vertex_gradient_variables_block_start_index(
  int num_variable_vertices,
  int dimension)
{
  // There are dimension many position variables per variable vertex
  return dimension * num_variable_vertices;
}

// Get the start index of the block of edge gradient variable indices
int
generate_global_edge_gradient_variables_block_start_index(
  int num_variable_vertices,
  int dimension)
{
  // There are dimension many position variables and 2 * dimension many vector
  // gradient variables per variable vertex
  return 3 * dimension * num_variable_vertices;
}

// Get the start index of the block of position variable indices for a given
// vertex
int
generate_global_vertex_position_variables_start_index(int vertex_index,
                                                      int dimension)
{
  int start_index =
    generate_global_vertex_position_variables_block_start_index();
  int relative_vertex_index = dimension * vertex_index;

  return start_index + relative_vertex_index;
}

// Get the start index of the block of gradient variable indices for a given
// vertex
int
generate_global_vertex_gradient_variables_start_index(int num_variable_vertices,
                                                      int vertex_index,
                                                      int dimension)
{
  int start_index = generate_global_vertex_gradient_variables_block_start_index(
    num_variable_vertices, dimension);
  int relative_vertex_index = 2 * dimension * vertex_index;

  return start_index + relative_vertex_index;
}

// Get the start index of the block of gradient variable indices for a given
// edge
int
generate_global_edge_gradient_variables_start_index(int num_variable_vertices,
                                                    int edge_index,
                                                    int dimension)
{
  int start_index = generate_global_edge_gradient_variables_block_start_index(
    num_variable_vertices, dimension);
  int relative_edge_index = dimension * edge_index;

  return start_index + relative_edge_index;
}

// ***********************
// Global variable indices
// ***********************

// Get the global index of the position variable indices for a given coordinate
// and vertex index
int
generate_global_vertex_position_variable_index(int vertex_index,
                                               int coord,
                                               int dimension)
{
  int start_index = generate_global_vertex_position_variables_start_index(
    vertex_index, dimension);

  return start_index + coord;
}

// Get the index of the gradient variable indices for a given matrix index pair
// and vertex index
int
generate_global_vertex_gradient_variable_index(int num_variable_vertices,
                                               int vertex_index,
                                               int row,
                                               int col,
                                               int dimension)
{
  int start_index = generate_global_vertex_gradient_variables_start_index(
    num_variable_vertices, vertex_index, dimension);
  int matrix_index = generate_local_variable_matrix_index(row, col, dimension);

  return start_index + matrix_index;
}

// Get the index of the gradient variable indices for a given coordinate and
// edge index pair
int
generate_global_edge_gradient_variable_index(int num_variable_vertices,
                                             int edge_index,
                                             int coord,
                                             int dimension)
{
  int start_index = generate_global_edge_gradient_variables_start_index(
    num_variable_vertices, edge_index, dimension);

  return start_index + coord;
}

// *******************
// Variable flattening
// *******************

// Get flat vector of all current variable values for the six-split
// NOTE: Also used as a subroutine to generate the twelve split maps
void
generate_six_split_variable_value_vector(
  const std::vector<SpatialVector>& vertex_positions,
  const std::vector<Matrix2x3r>& vertex_gradients,
  const std::vector<int>& variable_vertices,
  VectorXr& variable_values)
{
  int num_variable_vertices = variable_vertices.size();
  if (num_variable_vertices == 0) {
    variable_values.resize(0);
    spdlog::warn("Building value vector for zero variable vertices");
    return;
  }
  int dimension = 3;
  variable_values.resize(3 * dimension * num_variable_vertices);

  // Get position values
  for (int vertex_index = 0; vertex_index < num_variable_vertices;
       ++vertex_index) {
    int start_index = generate_global_vertex_position_variables_start_index(
      vertex_index, dimension);
    vertex_positions[variable_vertices[vertex_index]].size();
    for (int i = 0; i < dimension; ++i) {
      variable_values[start_index + i] =
        vertex_positions[variable_vertices[vertex_index]][i];
    }
  }

  // Get gradient values
  for (int vertex_index = 0; vertex_index < num_variable_vertices;
       ++vertex_index) {
    int start_index = generate_global_vertex_gradient_variables_start_index(
      num_variable_vertices, vertex_index, dimension);
    Matrix2x3r const& variable_matrix =
      vertex_gradients[variable_vertices[vertex_index]];
    for (Eigen::Index i = 0; i < variable_matrix.rows(); ++i) {
      for (Eigen::Index j = 0; j < variable_matrix.cols(); ++j) {
        int local_index = generate_local_variable_matrix_index(i, j, dimension);
        variable_values[start_index + local_index] = variable_matrix(i, j);
      }
    }
  }
}

// Get flat vector of all current variable values for the twelve-split
void
generate_twelve_split_variable_value_vector(
  const std::vector<SpatialVector>& vertex_positions,
  const std::vector<Matrix2x3r>& vertex_gradients,
  const std::vector<std::array<SpatialVector, 3>>& edge_gradients,
  const std::vector<int>& variable_vertices,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  VectorXr& variable_values)
{
  // Get the variable values shared with the six-split
  VectorXr six_split_variable_values;
  generate_six_split_variable_value_vector(vertex_positions,
                                           vertex_gradients,
                                           variable_vertices,
                                           six_split_variable_values);

  // Add six split to the variable value vector
  // Build a halfedge representation to get unique edge values
  int num_variable_vertices = variable_vertices.size();
  int num_variable_edges = variable_edges.size();
  int dimension = 3;
  variable_values.resize(3 * dimension * num_variable_vertices +
                         dimension * num_variable_edges);
  variable_values.head(six_split_variable_values.size()) =
    six_split_variable_values;

  // Get flat values for edge gradients
  for (int variable_edge_index = 0; variable_edge_index < num_variable_edges;
       ++variable_edge_index) {
    // Get one corner for the given edge
    Halfedge::Index edge_index = variable_edges[variable_edge_index];
    Halfedge::Index halfedge_index =
      halfedge.edge_to_first_halfedge(edge_index);
    Eigen::Index face_index = he_to_corner[halfedge_index].first;
    Eigen::Index face_vertex_index = he_to_corner[halfedge_index].second;

    // Extract each coordinate for the corner
    int num_coordinates = 3;
    for (int coord = 0; coord < num_coordinates; ++coord) {
      // Get edge variable index
      int variable_index = generate_global_edge_gradient_variable_index(
        num_variable_vertices, variable_edge_index, coord);

      // Extract variable value of the edge to the corner
      variable_values[variable_index] =
        edge_gradients[face_index][face_vertex_index][coord];
    }
  }
}

// Map local triangle vertex indices to their global variable indices
// NOTE: Also used as a subroutine to generate the twelve split maps
void
generate_six_split_local_to_global_map(
  const std::array<int, 3>& global_vertex_indices,
  int num_variable_vertices,
  std::vector<int>& local_to_global_map)
{
  int dimension = 3;
  local_to_global_map.resize(27);
  for (int local_vertex_index = 0; local_vertex_index < 3;
       ++local_vertex_index) {
    int global_vertex_index = global_vertex_indices[local_vertex_index];

    // Add vertex position index values
    for (int coord = 0; coord < dimension; ++coord) {
      int local_index = generate_local_vertex_position_variable_index(
        local_vertex_index, coord, dimension);
      int global_index;
      if (global_vertex_index < 0) {
        global_index = -1;
      } else {
        global_index = generate_global_vertex_position_variable_index(
          global_vertex_index, coord, dimension);
      }
      local_to_global_map[local_index] = global_index;
    }

    // Add vertex gradient index values
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < dimension; ++col) {
        int local_index = generate_local_vertex_gradient_variable_index(
          local_vertex_index, row, col, dimension);
        int global_index;
        if (global_vertex_index < 0) {
          global_index = -1;
        } else {
          global_index = generate_global_vertex_gradient_variable_index(
            num_variable_vertices, global_vertex_index, row, col, dimension);
        }
        local_to_global_map[local_index] = global_index;
      }
    }
  }
}

// Map local triangle vertex and edge indices to their global variable indices
void
generate_twelve_split_local_to_global_map(
  const std::array<int, 3>& global_vertex_indices,
  const std::array<int, 3>& global_edge_indices,
  int num_variable_vertices,
  std::vector<int>& local_to_global_map)
{
  // Get index map for the Powell-Sabin shared variables
  int dimension = 3;
  local_to_global_map.resize(36);
  std::vector<int> six_split_local_to_global_map;
  generate_six_split_local_to_global_map(global_vertex_indices,
                                         num_variable_vertices,
                                         six_split_local_to_global_map);
  std::copy(six_split_local_to_global_map.begin(),
            six_split_local_to_global_map.end(),
            local_to_global_map.begin());

  for (int local_edge_index = 0; local_edge_index < 3; ++local_edge_index) {
    int global_edge_index = global_edge_indices[local_edge_index];

    // Add edge gradient index values
    for (int coord = 0; coord < dimension; ++coord) {
      int local_index = generate_local_edge_gradient_variable_index(
        local_edge_index, coord, dimension);
      int global_index;
      if (global_edge_index < 0) {
        global_index = -1;
      } else {
        global_index = generate_global_edge_gradient_variable_index(
          num_variable_vertices, global_edge_index, coord, dimension);
      }
      local_to_global_map[local_index] = global_index;
    }
  }
}

// Update variables in a vector from the vector of all variable values from some
// start index
void
update_independent_variable_vector(const VectorXr& variable_values,
                                   SpatialVector& variable_vector,
                                   int start_index)
{
  for (Eigen::Index i = 0; i < variable_vector.size(); ++i) {
    int variable_index = start_index + i;
    variable_vector[i] = variable_values[variable_index];
  }
}

// Update variables in a matrix from the vector of all variable values from some
// start index The flattening of the matrix is assumed to be row major.
void
update_independent_variable_matrix(const VectorXr& variable_values,
                                   Matrix2x3r& variable_matrix,
                                   int start_index)
{
  int dimension = variable_matrix.cols();
  for (Eigen::Index i = 0; i < variable_matrix.rows(); ++i) {
    for (Eigen::Index j = 0; j < variable_matrix.cols(); ++j) {
      int local_index = generate_local_variable_matrix_index(i, j, dimension);
      int variable_index = start_index + local_index;
      variable_matrix(i, j) = variable_values[variable_index];
    }
  }
}

// Update all position variables
void
update_position_variables(const VectorXr& variable_values,
                          const std::vector<int>& variable_vertices,
                          std::vector<SpatialVector>& vertex_positions)
{
  if (vertex_positions.empty())
    return;
  int num_variable_vertices = variable_vertices.size();
  int dimension = 3;
  for (int vertex_index = 0; vertex_index < num_variable_vertices;
       ++vertex_index) {
    int start_index = generate_global_vertex_position_variables_start_index(
      vertex_index, dimension);
    update_independent_variable_vector(
      variable_values, vertex_positions[variable_vertices[vertex_index]], start_index);
  }
}

// Update all vertex gradient variables
void
update_vertex_gradient_variables(const VectorXr& variable_values,
                                 const std::vector<int>& variable_vertices,
                                 std::vector<Matrix2x3r>& vertex_gradients)
{
  int num_variable_vertices = variable_vertices.size();
  int dimension = 3;
  for (int vertex_index = 0; vertex_index < num_variable_vertices;
       ++vertex_index) {
    int start_index = generate_global_vertex_gradient_variables_start_index(
      num_variable_vertices, vertex_index, dimension);
    update_independent_variable_matrix(
      variable_values, vertex_gradients[variable_vertices[vertex_index]], start_index);
  }
}

// Update all edge gradient variables
void
update_edge_gradient_variables(
  const VectorXr& variable_values,
  const std::vector<int>& variable_vertices,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  std::vector<std::array<SpatialVector, 3>>& edge_gradients)
{
  int dimension = 3;
  int num_variable_vertices = variable_vertices.size();
  int num_variable_edges = variable_edges.size();
  // Get flat values for edge gradients
  for (int variable_edge_index = 0; variable_edge_index < num_variable_edges;
       ++variable_edge_index) {
    // Get corner for the given edge
    int edge_index = variable_edges[variable_edge_index];
    int first_halfedge_index = halfedge.edge_to_first_halfedge(edge_index);
    int first_face_index = he_to_corner[first_halfedge_index].first;
    int first_face_vertex_index = he_to_corner[first_halfedge_index].second;

    // Get index in the flattened variable vector
    int start_index = generate_global_edge_gradient_variables_start_index(
      num_variable_vertices, variable_edge_index, dimension);

    // Update the gradients for the first corner
    update_independent_variable_vector(
      variable_values,
      edge_gradients[first_face_index][first_face_vertex_index],
      start_index);

    // Update the gradients for the second corner if it exists
    if (!halfedge.is_boundary_edge(edge_index)) {
      int second_halfedge_index = halfedge.edge_to_second_halfedge(edge_index);
      int second_face_index = he_to_corner[second_halfedge_index].first;
      int second_face_vertex_index = he_to_corner[second_halfedge_index].second;
      update_independent_variable_vector(
        variable_values,
        edge_gradients[second_face_index][second_face_vertex_index],
        start_index);
    }
  }
}

void
build_variable_vertex_indices_map(int num_vertices,
                                  const std::vector<int>& variable_vertices,
                                  std::vector<int>& global_vertex_indices)
{
  // Get variable vertex indices
  global_vertex_indices = std::vector<int>(num_vertices, -1);
  for (size_t i = 0; i < variable_vertices.size(); ++i) {
    global_vertex_indices[variable_vertices[i]] = i;
  }
}

void
build_variable_edge_indices_map(
  int num_faces,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  std::vector<std::array<int, 3>>& global_edge_indices)
{
  global_edge_indices.resize(num_faces);
  for (size_t i = 0; i < variable_edges.size(); ++i) {
    int edge_index = variable_edges[i];
    int h0 = halfedge.edge_to_first_halfedge(edge_index);
    int f0 = he_to_corner[h0].first;
    int f0_vertex_index = he_to_corner[h0].second;
    global_edge_indices[f0][f0_vertex_index] = i;

    if (!halfedge.is_boundary_edge(edge_index)) {
      int h1 = halfedge.edge_to_second_halfedge(edge_index);
      int f1 = he_to_corner[h1].first;
      int f1_vertex_index = he_to_corner[h1].second;
      global_edge_indices[f1][f1_vertex_index] = i;
    }
  }
}
