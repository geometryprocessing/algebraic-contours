// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "position_data.h"

// #define VIEW_POWELL_SABIN

// View position data at one triangle
void
view_triangle_corner_data(
  std::vector<TriangleCornerFunctionData>& corner_data)
{
  // Organize geometric data for polyscope
  Eigen::MatrixXd V(3, 3);
  Eigen::MatrixXi F(1, 3);
  MatrixXr first_gradient(3, 3);
  MatrixXr second_gradient(3, 3);
  for (size_t j = 0; j < 3; ++j) {
    F(0, j) = j;
    V.row(j) = corner_data[j].function_value;
    first_gradient.row(j) = corner_data[j].first_edge_derivative;
    second_gradient.row(j) = corner_data[j].second_edge_derivative;
  }
  spdlog::info("Position data vertices: {}", V);

  // Build viewer with gradients
  polyscope::registerSurfaceMesh("surface", V, F);
  polyscope::getSurfaceMesh("surface")->addVertexVectorQuantity(
    "first_gradient", first_gradient);
  polyscope::getSurfaceMesh("surface")->addVertexVectorQuantity(
    "second_gradient", second_gradient);

  // Show corner data
  polyscope::show();
}

// Helper function to add corner data for a given chart
void
generate_affine_manifold_chart_corner_data(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const VertexManifoldChart& chart,
  const Matrix2x3r& gradient,
  std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data)
{
  // Build position data for the corners adjacent to the given vertex
  for (size_t face_index = 0; face_index < chart.face_one_ring.size();
       ++face_index) {
    // Get the face and index of the vertex in it
    int f = chart.face_one_ring[face_index];
    int face_vertex_index =
      find_face_vertex_index(F.row(f), chart.vertex_index);

    // Compute position and edge derivatives from layout positions
    SpatialVector pi = V.row(chart.vertex_index);
    OneFormXr uj = chart.one_ring_uv_positions.row(face_index);
    OneFormXr uk = chart.one_ring_uv_positions.row(face_index + 1);
    SpatialVector dij = uj * gradient;
    SpatialVector dik = uk * gradient;
    spdlog::trace("Vertex position: {}", pi);
    spdlog::trace("Vertex gradient:\n{}", gradient);
    spdlog::trace("First edge layout: {}", uj);
    spdlog::trace("Second edge layout: {}", uk);
    spdlog::trace("First edge derivative: {}", dij);
    spdlog::trace("Second edge derivative: {}", dik);

    // Build position data
    corner_data[f][face_vertex_index] =
      (TriangleCornerFunctionData(pi, dij, dik));
  }
}

void
generate_affine_manifold_corner_data(
  const Eigen::MatrixXd& V,
  const AffineManifold& affine_manifold,
  const std::vector<Matrix2x3r>& gradients,
  std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data)
{
  // Resize the gradient and position data
  corner_data.resize(affine_manifold.num_faces());

  // Compute the gradient and position data per vertex
  for (Eigen::Index i = 0; i < V.rows(); ++i) {
    generate_affine_manifold_chart_corner_data(
      V,
      affine_manifold.get_faces(),
      affine_manifold.get_vertex_chart(i),
      gradients[i],
      corner_data);
  }
}

void
generate_affine_manifold_midpoint_data(
  const AffineManifold& affine_manifold,
  const std::vector<std::array<Matrix2x3r, 3>>& edge_gradients,
  std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data)
{
  Eigen::MatrixXi const& F = affine_manifold.get_faces();

  // Set midpoint data per face corner
  midpoint_data.resize(affine_manifold.num_faces());
  for (AffineManifold::Index i = 0; i < affine_manifold.num_faces(); ++i) {
    for (AffineManifold::Index j = 0; j < 3; ++j) {
      // Get local edge chart
      EdgeManifoldChart const& chart = affine_manifold.get_edge_chart(i, j);

      // Get the corner index for the top face
      int f_top = chart.top_face_index;
      int v_top = chart.top_vertex_index;
      int j_top = find_face_vertex_index(F.row(f_top), v_top);

      // Only process top faces of edge charts to prevent redundancy
      if (f_top != i) continue; 
      assert(j_top == j);

      // Compute midpoint to opposite corner derivative for the top face
      PlanarPoint uv_top = chart.top_vertex_uv_position;
      midpoint_data[f_top][j_top].normal_derivative =
        uv_top * edge_gradients[i][j];
      SPDLOG_TRACE("Midpoint data for corner ({}, {}) is {} = {}\n{}",
                   f_top,
                   j_top,
                   midpoint_data[f_top][j_top].normal_derivative.transpose(),
                   uv_top,
                   edge_gradients[i][j]);

      // Only set the bottom vertex if the edge is not on the boundary
      if (!chart.is_boundary) {
        // Get the corner index for the bottom face
        int f_bottom = chart.bottom_face_index;
        int v_bottom = chart.bottom_vertex_index;
        int j_bottom = find_face_vertex_index(F.row(f_bottom), v_bottom);

        // Compute midpoint to opposite corner derivative for the bottom face
        PlanarPoint uv_bottom = chart.bottom_vertex_uv_position;
        midpoint_data[f_bottom][j_bottom].normal_derivative =
          uv_bottom * edge_gradients[i][j];
        SPDLOG_TRACE("Midpoint data for corner ({}, {}) is {}",
                     f_bottom,
                     j_bottom,
                     midpoint_data[f_bottom][j_bottom].normal_derivative);
      }
    }
  }
}

void
compute_edge_midpoint_with_gradient(
  const TriangleCornerFunctionData& edge_origin_corner_data,
  const TriangleCornerFunctionData& edge_dest_corner_data,
  SpatialVector& midpoint,
  SpatialVector& midpoint_edge_gradient)
{
  SpatialVector fi = edge_origin_corner_data.function_value;
  SpatialVector fj = edge_dest_corner_data.function_value;
  SpatialVector dij = edge_origin_corner_data.first_edge_derivative;
  SpatialVector dji = edge_dest_corner_data.second_edge_derivative;
  midpoint = 0.5 * (fi + fj) + 0.125 * (dij + dji);
  midpoint_edge_gradient = 2.0 * (fj - fi) + 0.5 * (dji - dij);
}

void
generate_corner_data_matrices(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  MatrixXr& position_matrix,
  MatrixXr& first_derivative_matrix,
  MatrixXr& second_derivative_matrix)
{
  int num_faces = corner_data.size();

  // Resize the corner data matrices
  position_matrix.resize(3 * num_faces, 3);
  first_derivative_matrix.resize(3 * num_faces, 3);
  second_derivative_matrix.resize(3 * num_faces, 3);

  // Organize position data into matrices
  for (int i = 0; i < num_faces; ++i) {
    for (int j = 0; j < 3; ++j) {
      position_matrix.row(3 * i + j) = corner_data[i][j].function_value;
      first_derivative_matrix.row(3 * i + j) =
        corner_data[i][j].first_edge_derivative;
      second_derivative_matrix.row(3 * i + j) =
        corner_data[i][j].second_edge_derivative;
    }
  }
}


void
generate_midpoint_data_matrices(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  const std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data,
  MatrixXr& position_matrix,
  MatrixXr& tangent_derivative_matrix,
  MatrixXr& normal_derivative_matrix)
{
  int num_faces = corner_data.size();

  // Resize the corner data matrices
  position_matrix.resize(3 * num_faces, 3);
  tangent_derivative_matrix.resize(3 * num_faces, 3);
  normal_derivative_matrix.resize(3 * num_faces, 3);

  // Organize position data into matrices
  for (int i = 0; i < num_faces; ++i) {
    for (int j = 0; j < 3; ++j) {
      SpatialVector midpoint;
      SpatialVector midpoint_edge_gradient;
      compute_edge_midpoint_with_gradient(corner_data[i][(j + 1) % 3],
                                          corner_data[i][(j + 2) % 3],
                                          midpoint,
                                          midpoint_edge_gradient);
      position_matrix.row(3 * i + j) = midpoint;
      tangent_derivative_matrix.row(3 * i + j) = midpoint_edge_gradient;
      normal_derivative_matrix.row(3 * i + j) =
        midpoint_data[i][j].normal_derivative;
    }
  }
}
