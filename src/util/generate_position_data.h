// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.


#pragma once
#include "common.h"
#include "generate_shapes.h"

// ***************************
// Parametric surface functors
// ***************************

class QuadraticPositionFunction{
public:
  QuadraticPositionFunction(
    double _uv_coeff,
    double _uu_coeff,
    double _vv_coeff
  ) {
    uv_coeff = _uv_coeff;
    uu_coeff = _uu_coeff;
    vv_coeff = _vv_coeff;
  }

  SpatialVector operator()(
    double u,
    double v
  ) {
    SpatialVector point;
    point(0) = u;
    point(1) = v;
    point(2) = uv_coeff * u * v + uu_coeff * u * u + vv_coeff * v * v;

    return point;
  }

  double uv_coeff;
  double uu_coeff;
  double vv_coeff;
};


class QuadraticGradientFunction{
public:
  QuadraticGradientFunction(
    double _uv_coeff,
    double _uu_coeff,
    double _vv_coeff
  ) {
    uv_coeff = _uv_coeff;
    uu_coeff = _uu_coeff;
    vv_coeff = _vv_coeff;
  }

  Matrix2x3r operator()(
    double u,
    double v
  ) {
    MatrixXr gradient(2, 3);

    // Generate derivative with respect to u
    gradient(0, 0) = 1.0;
    gradient(0, 1) = 0.0;
    gradient(0, 2) = uv_coeff * v + 2 * uu_coeff * u;

    // Generate derivative with respect to v
    gradient(1, 0) = 0.0;
    gradient(1, 1) = 1.0;
    gradient(1, 2) = uv_coeff * u + 2 * vv_coeff * v;

    return gradient;
  }

  double uv_coeff;
  double uu_coeff;
  double vv_coeff;
};

class TorusPositionFunction {
public:
  TorusPositionFunction(
    double _major_radius,
    double _minor_radius,
    int _resolution,
    double _angle_offset
  ) {
    R = _major_radius;
    r = _minor_radius;
    resolution = _resolution;
    angle_offset = _angle_offset;
  }

  SpatialVector operator()(
    double u,
    double v
  ) {
    SpatialVector point;
    double theta = generate_angle(u, resolution, angle_offset);
    double phi = generate_angle(v, resolution, angle_offset);
    point(0) = (R + r * std::cos(theta)) * std::cos(phi);
    point(1) = (R + r * std::cos(theta)) * std::sin(phi);
    point(2) = r * std::sin(theta);

    return point;
  }

  double R;
  double r;
  int resolution;
  double angle_offset;
};


class TorusGradientFunction {
public:
  TorusGradientFunction(
    double _major_radius,
    double _minor_radius,
    int _resolution,
    double _angle_offset
  ) {
    R = _major_radius;
    r = _minor_radius;
    resolution = _resolution;
    angle_offset = _angle_offset;
  }

  Matrix2x3r operator()(
    double u,
    double v
  ) {
    Matrix2x3r gradient(2, 3);

    // Compute angle derivatives
    double theta = generate_angle(u, resolution, angle_offset);
    double phi = generate_angle(v, resolution, angle_offset);
    double dtheta_du = generate_angle_derivative(resolution);
    double dphi_dv = generate_angle_derivative(resolution);

    // Generate derivative with respect to u
    gradient(0, 0) = -r * std::sin(theta) * std::cos(phi) * dtheta_du;
    gradient(0, 1) = -r * std::sin(theta) * std::sin(phi) * dtheta_du;
    gradient(0, 2) = r * std::cos(theta) * dtheta_du;

    // Generate derivative with respect to v
    gradient(1, 0) = -(R + r * std::cos(theta)) * std::sin(phi) * dphi_dv;
    gradient(1, 1) = (R + r * std::cos(theta)) * std::cos(phi) * dphi_dv;
    gradient(1, 2) = 0.0;

    return gradient;
  }

  double R;
  double r;
  int resolution;
  double angle_offset;
};

template <typename PositionFunction>
void generate_point_grid(
  int resolution,
  PositionFunction position_func,
  std::vector< std::vector<VectorXr> > &point_grid
) {
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      double u = i * 1.0;
      double v = j * 1.0;
      point_grid[i][j] = position_func(u, v);
    }
  }
}


template <typename PositionFunction>
void generate_parametric_affine_manifold_vertex_positions(
  PositionFunction &position_func,
  ParametricAffineManifold &parametric_affine_manifold,
  MatrixXr &vertex_positions
) {
  size_t num_vertices = parametric_affine_manifold.num_vertices();
  vertex_positions.resize(num_vertices, 3);
  for (size_t vi = 0; vi < num_vertices; ++vi)
  {
    // Get global vertex uv coordinates
    PlanarPoint uv;
    parametric_affine_manifold.get_vertex_global_uv(vi, uv);
    double u = uv[0];
    double v = uv[1];

    // Compute gradient at these coordinates
    vertex_positions.row(vi) = position_func(u, v);
  }
}


template <typename GradientFunction>
void generate_parametric_affine_manifold_vertex_gradients(
  GradientFunction &gradient_func,
  ParametricAffineManifold &parametric_affine_manifold,
  std::vector<Matrix2x3r> &vertex_gradients
) {
  size_t num_vertices = parametric_affine_manifold.num_vertices();
  vertex_gradients.resize(num_vertices);
  for (size_t vi = 0; vi < num_vertices; ++vi)
  {
    // Get global vertex uv coordinates
    PlanarPoint uv;
    parametric_affine_manifold.get_vertex_global_uv(vi, uv);
    double u = uv[0];
    double v = uv[1];

    // Compute gradient at these coordinates
    vertex_gradients[vi] = gradient_func(u, v);
  }
}


template <typename GradientFunction>
void generate_parametric_affine_manifold_edge_gradients(
  GradientFunction &gradient_func,
  ParametricAffineManifold &parametric_affine_manifold,
  Halfedge &halfedge,
  std::vector<std::array<SpatialVector, 3>> &edge_gradients
) {
  Eigen::MatrixXi const & F = parametric_affine_manifold.get_faces();
  size_t num_faces = parametric_affine_manifold.num_faces();

  edge_gradients.resize(num_faces);
  for (size_t face_index = 0; face_index < num_faces; ++face_index)
  {
    for (size_t face_vertex_index = 0; face_vertex_index < 3; ++face_vertex_index)
    {
      EdgeManifoldChart const & chart = parametric_affine_manifold.get_edge_chart(
        face_index,
        face_vertex_index
      );
      size_t f_top = chart.top_face_index;
      size_t v_top = chart.top_vertex_index;
      int j_top = find_face_vertex_index(F.row(f_top), v_top);
      if ( f_top != face_index ) continue; // Only process top faces of edge charts to prevent redundancy
      assert( j_top == face_vertex_index );

      // Get face indices
      size_t i = F(face_index, (face_vertex_index + 1) % 3);
      size_t j = F(face_index, (face_vertex_index + 2) % 3);
      size_t k = F(face_index, face_vertex_index);
      assert (k == v_top);

      // Get vertex uv positions
      PlanarPoint uvi, uvj;
      parametric_affine_manifold.get_vertex_global_uv(i, uvi);
      parametric_affine_manifold.get_vertex_global_uv(j, uvj);

      // Get midpoint uv position and gradient
      PlanarPoint uvij = 0.5 * (uvi + uvj);
      Matrix2x3r Gij = gradient_func(uvij[0], uvij[1]);
      spdlog::debug(
        "Gradient for corner {}, {} is: \n{}",
        face_index,
        face_vertex_index,
        Gij
      );

      // Get normal uv directions perpendicular to the edge direction
      Matrix2x2r R;
      R <<
        0, 1,
        -1, 0;
      PlanarPoint dji = (uvj - uvi); 
      PlanarPoint uv_perp_ij = dji * R;

      // Build midpoint normals (indexed opposite the edge)
      edge_gradients[face_index][face_vertex_index] = uv_perp_ij * Gij;
      spdlog::debug(
        "Reduced edge gradient for corner {}, {} is: \n{}",
        face_index,
        face_vertex_index,
        edge_gradients[face_index][face_vertex_index]
      );

      // If the edge isn't on the boundary, set the other face corner corresponding to it
      if (!chart.is_boundary)
      {
        size_t f_bottom = chart.bottom_face_index;
        size_t v_bottom = chart.bottom_vertex_index;
        int j_bottom = find_face_vertex_index(F.row(f_bottom), v_bottom);
        edge_gradients[f_bottom][j_bottom] = edge_gradients[face_index][face_vertex_index];
      }
    }
  }
}


template <typename PositionFunction, typename GradientFunction>
void generate_parametric_affine_manifold_corner_data(
  PositionFunction position_func,
  GradientFunction gradient_func,
  ParametricAffineManifold &parametric_affine_manifold,
  std::vector<std::array<TriangleCornerFunctionData, 3>> &corner_data
) {
  size_t num_faces = parametric_affine_manifold.num_faces();
  corner_data.resize(num_faces);
  for (size_t face_index = 0; face_index < num_faces; ++face_index)
  {
    // Get face vertex indices
    Eigen::MatrixXi const & F = parametric_affine_manifold.get_faces();
    size_t i = F(face_index, 0);
    size_t j = F(face_index, 1);
    size_t k = F(face_index, 2);

    // Get vertex uv positions
    PlanarPoint uvi, uvj, uvk;
    parametric_affine_manifold.get_vertex_global_uv(i, uvi);
    parametric_affine_manifold.get_vertex_global_uv(j, uvj);
    parametric_affine_manifold.get_vertex_global_uv(k, uvk);

    // Get vertex positions
    SpatialVector vi = position_func(uvi[0], uvi[1]);
    SpatialVector vj = position_func(uvj[0], uvj[1]);
    SpatialVector vk = position_func(uvk[0], uvk[1]);

    // Get vertex gradients
    Matrix2x3r Gi = gradient_func(uvi[0], uvi[1]);
    Matrix2x3r Gj = gradient_func(uvj[0], uvj[1]);
    Matrix2x3r Gk = gradient_func(uvk[0], uvk[1]);

    // Get uv directions
    PlanarPoint dij = uvj - uvi;
    PlanarPoint dik = uvk - uvi;
    PlanarPoint djk = uvk - uvj;
    PlanarPoint dji = uvi - uvj;
    PlanarPoint dki = uvi - uvk;
    PlanarPoint dkj = uvj - uvk;

    // Build first corner data 
    corner_data[face_index][0].function_value = vi;
    corner_data[face_index][0].first_edge_derivative = dij * Gi;
    corner_data[face_index][0].second_edge_derivative = dik * Gi;

    // Build second corner data 
    corner_data[face_index][1].function_value = vj;
    corner_data[face_index][1].first_edge_derivative = djk * Gj;
    corner_data[face_index][1].second_edge_derivative = dji * Gj;

    // Build third corner data 
    corner_data[face_index][2].function_value = vk;
    corner_data[face_index][2].first_edge_derivative = dki * Gk;
    corner_data[face_index][2].second_edge_derivative = dkj * Gk;
  }
}


template <typename GradientFunction>
void generate_parametric_affine_manifold_midpoint_data(
  GradientFunction gradient_func,
  ParametricAffineManifold &parametric_affine_manifold,
  std::vector<std::array<TriangleMidpointFunctionData, 3>> &midpoint_data
) {
  size_t num_faces = parametric_affine_manifold.num_faces();
  midpoint_data.resize(num_faces);
  for (size_t face_index = 0; face_index < num_faces; ++face_index)
  {
    // Get face vertex indices
    Eigen::MatrixXi const & F = parametric_affine_manifold.get_faces();
    size_t i = F(face_index, 0);
    size_t j = F(face_index, 1);
    size_t k = F(face_index, 2);

    // Get vertex uv positions
    PlanarPoint uvi, uvj, uvk;
    parametric_affine_manifold.get_vertex_global_uv(i, uvi);
    parametric_affine_manifold.get_vertex_global_uv(j, uvj);
    parametric_affine_manifold.get_vertex_global_uv(k, uvk);

    // Get midpoint uv positions
    PlanarPoint uvij = 0.5 * (uvi + uvj);
    PlanarPoint uvjk = 0.5 * (uvj + uvk);
    PlanarPoint uvki = 0.5 * (uvk + uvi);

    // Get midpoint gradients
    Matrix2x3r Gij = gradient_func(uvij[0], uvij[1]);
    Matrix2x3r Gjk = gradient_func(uvjk[0], uvjk[1]);
    Matrix2x3r Gki = gradient_func(uvki[0], uvki[1]);

    // Get uv directions
    PlanarPoint nij = uvk - uvij;
    PlanarPoint njk = uvi - uvjk;
    PlanarPoint nki = uvj - uvki;

    // Build midpoint normals (indexed opposite the edge)
    midpoint_data[face_index][0].normal_derivative = njk * Gjk;
    midpoint_data[face_index][1].normal_derivative = nki * Gki;
    midpoint_data[face_index][2].normal_derivative = nij * Gij;
  }
}