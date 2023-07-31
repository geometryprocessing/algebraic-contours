// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "affine_manifold.h"
#include "common.h"
#include "differentiable_variable.h"
#include "halfedge.h"
#include "polynomial_function.h"
#include "position_data.h"
#include <Eigen/CholmodSupport>
#include <Eigen/Sparse>

/// \file powell_sabin_local_to_global_indexing.h
///
/// Methods to generate local to global indexing maps for six and twelve split
/// Powell-Sabin spline surfaces.

/// Compute the index of a vertex position variable in a local DOF vector.
///
/// @param[in] face_vertex_index: index of the vertex in the face
/// @param[in] coord: coordinate of the variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the local DOF vector
size_t
generate_local_vertex_position_variable_index(size_t face_vertex_index,
                                              size_t coord,
                                              size_t dimension = 3);

/// Compute the index of a vertex gradient variable in a local DOF vector.
///
/// @param[in] face_vertex_index: index of the vertex in the face
/// @param[in] row: row of the gradient matrix variable
/// @param[in] col: column of the gradient matrix variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the local DOF vector
size_t
generate_local_vertex_gradient_variable_index(size_t face_vertex_index,
                                              size_t row,
                                              size_t col,
                                              size_t dimension = 3);

/// Compute the index of a edge gradient variable in a local DOF vector.
///
/// @param[in] face_vertex_index: index of the edge in the face
/// @param[in] coord: coordinate of the variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the local DOF vector
size_t
generate_local_edge_gradient_variable_index(size_t face_edge_index,
                                            size_t coord,
                                            size_t dimension = 3);

/// Compute the index of a vertex position variable in a global DOF vector.
///
/// @param[in] vertex_index: index of the vertex in the mesh
/// @param[in] coord: coordinate of the variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the global DOF vector
int
generate_global_vertex_position_variable_index(int vertex_index,
                                               int coord,
                                               int dimension = 3);

/// Compute the index of a vertex gradient variable in a global DOF vector.
///
/// @param[in] num_variable_vertices: number of variable vertices for the optimization
/// @param[in] vertex_index: index of the vertex in the mesh
/// @param[in] row: row of the gradient matrix variable
/// @param[in] col: column of the gradient matrix variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the global DOF vector
int
generate_global_vertex_gradient_variable_index(int num_variable_vertices,
                                               int vertex_index,
                                               int row,
                                               int col,
                                               int dimension = 3);

/// Compute the index of an edge gradient variable in a global DOF vector.
///
/// @param[in] num_variable_vertices: number of variable vertices for the optimization
/// @param[in] edge_index: index of the edge in the mesh
/// @param[in] coord: coordinate of the variable
/// @param[in] dimension: number of coordinate dimensions
/// @return index of the variable in the global DOF vector
int
generate_global_edge_gradient_variable_index(int num_variable_vertices,
                                             int edge_index,
                                             int coord,
                                             int dimension = 3);

/// Given vertex positions and gradients and a list of variable vertices, assemble
/// the vector of global vertex degrees of freedom
///
/// This is the complete list of degrees of freedom for the six split, and it is
/// a subset of the degrees of freedom for the twelve split.
///
/// @param[in] vertex_positions: list of vertex position values
/// @param[in] vertex_gradients: list of vertex gradient matrices 
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[out] variable_values: vertex DOF vector
void
generate_six_split_variable_value_vector(
  const std::vector<SpatialVector>& vertex_positions,
  const std::vector<Matrix2x3r>& vertex_gradients,
  const std::vector<int>& variable_vertices,
  VectorXr& variable_values);

/// Given vertex positions and gradients, edge gradients, and lists
/// of variable vertices and edges assemble the vector of global
/// degrees of freedom for the twelve split
///
/// @param[in] vertex_positions: list of vertex position values
/// @param[in] vertex_gradients: list of vertex gradient matrices 
/// @param[in] edge_gradients: list of edge gradient normal vectors 
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[in] variable_edges: list of variable edge indices
/// @param[in] halfedge: halfedge data structure
/// @param[in] he_to_corner: map from halfedges to opposite triangle corners
/// @param[out] variable_values: twelve-split DOF vector
void
generate_twelve_split_variable_value_vector(
  const std::vector<SpatialVector>& vertex_positions,
  const std::vector<Matrix2x3r>& vertex_gradients,
  const std::vector<std::array<SpatialVector, 3>>& edge_gradients,
  const std::vector<int>& variable_vertices,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  VectorXr& variable_values);

/// Given the global vertex indices of a triangle, compute the map from the
/// local DOF vector indices for this triangle to their indices in the global
/// DOF vector for the six-split
///
/// This is used as a subroutine for the twelve-split local to global map.
///
/// @param[in] global_vertex_indices: global indices of the triangle vertices
/// @param[in] num_variable_vertices: number of variable vertices
/// @param[out] local_to_global_map: map from local to global DOF indices
void
generate_six_split_local_to_global_map(
  const std::array<int, 3>& global_vertex_indices,
  int num_variable_vertices,
  std::vector<int>& local_to_global_map);

/// Given the global vertex and edge indices of a triangle, compute the map
/// from the local DOF vector indices for this triangle to their indices in
/// the global DOF vector for the twelve-split.
///
/// @param[in] global_vertex_indices: global indices of the triangle vertices
/// @param[in] global_vertex_indices: global indices of the triangle edges 
/// @param[in] num_variable_vertices: number of variable vertices
/// @param[out] local_to_global_map: map from local to global DOF indices
void
generate_twelve_split_local_to_global_map(
  const std::array<int, 3>& global_vertex_indices,
  const std::array<int, 3>& global_edge_indices,
  int num_variable_vertices,
  std::vector<int>& local_to_global_map);

/// Extract vertex positions from the global DOF vector.
///
/// @param[in] variable_values: twelve-split DOF vector
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[out] vertex_positions: list of vertex position values
void
update_position_variables(const VectorXr& variable_values,
                          const std::vector<int>& variable_vertices,
                          std::vector<SpatialVector>& vertex_positions);

/// Extract vertex gradients from the global DOF vector.
///
/// @param[in] variable_values: twelve-split DOF vector
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[out] vertex_gradients: list of vertex gradient values
void
update_vertex_gradient_variables(const VectorXr& variable_values,
                                 const std::vector<int>& variable_vertices,
                                 std::vector<Matrix2x3r>& vertex_gradients);

/// Extract edge gradients from the global DOF vector.
///
/// @param[in] variable_values: twelve-split DOF vector
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[in] variable_edges: list of variable edge indices
/// @param[in] halfedge: halfedge data structure
/// @param[in] he_to_corner: map from halfedges to opposite triangle corners
/// @param[out] edge_gradients: list of edge gradient values
void
update_edge_gradient_variables(
  const VectorXr& variable_values,
  const std::vector<int>& variable_vertices,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  std::vector<std::array<SpatialVector, 3>>& edge_gradients);

/// Generate a map from all vertices to a list of variable vertices or -1 for
/// vertices that are not variable.
///
/// @param[in] num_vertices: total number of vertices
/// @param[in] variable_vertices: list of variable vertex indices
/// @param[out] global_vertex_indices: map from vertex indices to variable vertices
void
build_variable_vertex_indices_map(int num_vertices,
                                  const std::vector<int>& variable_vertices,
                                  std::vector<int>& global_vertex_indices);


/// Generate a map from all edges to a list of variable edges or -1 for
/// edges that are not variable.
///
/// @param[in] num_faces: total number of faces
/// @param[in] variable_edges: list of variable edge indices
/// @param[in] halfedge: halfedge data structure
/// @param[in] he_to_corner: map from halfedges to opposite triangle corners
/// @param[out] global_edge_indices: map from edge indices to variable edges
void
build_variable_edge_indices_map(
  int num_faces,
  const std::vector<int>& variable_edges,
  const Halfedge& halfedge,
  const std::vector<std::pair<Eigen::Index, Eigen::Index>>& he_to_corner,
  std::vector<std::array<int, 3>>& global_edge_indices);

/// Update global energy, derivatives, and hessian with local per face values
///
/// @param[in] local_energy: local energy value
/// @param[in] local_derivatives: local energy gradient
/// @param[in] local_hessian: local energy Hessian
/// @param[in] local_to_global_map: map from local to global DOF indices
/// @param[out] energy: global energy value
/// @param[out] derivatives: global energy gradient
/// @param[out] hessian: global energy Hessian
template<typename Gradient, typename Hessian>
void
update_energy_quadratic(const double& local_energy,
                        const Gradient& local_derivatives,
                        const Hessian& local_hessian,
                        const std::vector<int>& local_to_global_map,
                        double& energy,
                        VectorXr& derivatives,
                        std::vector<Eigen::Triplet<double>>& hessian_entries)
{
  spdlog::trace("Adding local face energy {}", local_energy);
  SPDLOG_TRACE("Local to global map: {}",
               formatted_vector(local_to_global_map, ", "));
  // Update energy
  energy += local_energy;

  // Update derivatives
  int num_local_indices = local_to_global_map.size();
  for (int local_index = 0; local_index < num_local_indices; ++local_index) {
    int global_index = local_to_global_map[local_index];
    if (global_index < 0)
      continue; // Skip fixed variables with no global index
    derivatives[global_index] += local_derivatives[local_index];
  }

  // Update hessian entries
  for (int local_index_i = 0; local_index_i < num_local_indices;
       ++local_index_i) {
    // Get global row index, skipping fixed variables with no global index
    int global_index_i = local_to_global_map[local_index_i];
    if (global_index_i < 0)
      continue;

    for (int local_index_j = 0; local_index_j < num_local_indices;
         ++local_index_j) {
      // Get global column index, skipping fixed variables with no global index
      int global_index_j = local_to_global_map[local_index_j];
      if (global_index_j < 0)
        continue; 

      // Get Hessian entry value
      double hessian_value = local_hessian(local_index_i, local_index_j);

      // Assemble global Hessian entry
      hessian_entries.push_back(
        Eigen::Triplet<double>(global_index_i, global_index_j, hessian_value));
    }
  }
}

/// Build a triplet of face vertex values from a global array of vertex variables
///
/// @param[in] variables: global variables
/// @param[in] i: first variable index
/// @param[in] j: second variable index
/// @param[in] k: third variable index
/// @param[out] face_variable_vector: variables for face Tijk
template<typename T>
void
build_face_variable_vector(const std::vector<T>& variables,
                           int i,
                           int j,
                           int k,
                           std::array<T, 3>& face_variable_vector)
{
  face_variable_vector[0] = variables[i];
  face_variable_vector[1] = variables[j];
  face_variable_vector[2] = variables[k];
}