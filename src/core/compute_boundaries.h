// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "halfedge.h"
#include "quadratic_spline_surface.h"

/// \file compute_boundaries.h
///
/// Methods to compute the boundaries of a mesh

/// Given a mesh, compute the edges on the boundary.
///
/// @param[in] F: mesh faces
/// @param[out] face_boundary_edges: edges of the triangles that are boundaries,
/// indexed by face and opposite corner local face vertex index
void
compute_face_boundary_edges(
  const Eigen::MatrixXi& F,
  std::vector<std::pair<int, int>>& face_boundary_edges);

/// Given a mesh, compute the vertices on the boundary.
///
/// @param[in] F: mesh faces
/// @param[out] boundary_vertices: vertices of the mesh on the boundary
void
compute_boundary_vertices(const Eigen::MatrixXi& F,
                          std::vector<int>& boundary_vertices);

