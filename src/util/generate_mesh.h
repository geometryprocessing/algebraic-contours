// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <Eigen/Core>

#include "common.h"
#include "conic.h"
#include "generate_shapes.h"
#include <vector>
#include <Eigen/Core>


/// Generate a surface mesh for a quadratic spline surface with given dimension.
///
/// @param[in] dimension: dimension of the mesh
/// @return Surface mesh for Zwart Powell patches corresponding to the given dimension
void generate_zwart_powell_mesh(int dimension);

/// Generate a VF representation for a quadratic spline surface.
///
/// @param[in] control_point_grid: Control point for the quadratic surface
/// @param[out] V: Vertices for the mesh
/// @param[out] F: Faces for the mesh
void compute_mesh_from_control_point_grid(
  const std::vector< std::vector<VectorXr> > &control_point_grid,
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F
);