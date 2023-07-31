// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "line_segment.h"

/// Representation of a convex polygon in R^2 that supports containment queries,
/// sampling, boundary segments and vertices computation, triangulation, and
/// boundary parametrization.
class ConvexPolygon
{
public:
  ConvexPolygon() {}

  ConvexPolygon(
    const std::array<Eigen::Matrix<double, 3, 1>, 3>& boundary_segments_coeffs);

  ConvexPolygon(const Eigen::Matrix<double, 3, 2>& vertices);

  // TODO Implement constructor from collection of points

  /// Return true iff point is in the convex polygon
  bool contains(const PlanarPoint& point) const;

  void intersect_patch_boundaries(
    const Eigen::Matrix<double, 3, 1>& first_boundary_segment_coeffs,
    const Eigen::Matrix<double, 3, 1>& second_boundary_segment_coeffs,
    PlanarPoint& intersection) const;

  std::array<Eigen::Matrix<double, 3, 1>, 3> const& get_boundary_segments()
    const;

  Eigen::Matrix<double, 3, 2> const& get_vertices() const;

  void parametrize_patch_boundaries(
    std::array<LineSegment, 3>& patch_boundaries) const;

  // Triangulate domain with
  void triangulate(size_t num_refinements,
                   Eigen::MatrixXd& V,
                   Eigen::MatrixXi& F) const;

  void sample(size_t num_samples,
              std::vector<PlanarPoint>& domain_points) const;

private:
  std::array<Eigen::Matrix<double, 3, 1>, 3> m_boundary_segments_coeffs;
  Eigen::Matrix<double, 3, 2> m_vertices;
};
