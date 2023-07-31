// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"

/// Discretize the given curve segments as a polyline curve network
///
/// @param[in] contour_segments: curve segments
/// @param[in] curve_disc_params: parameters for the contour discretization
/// @param[out] points: points of the curve network
/// @param[out] polylines: polyline indices of the curve network
template<size_t degree, size_t dimension>
void
discretize_curve_segments(
  const std::vector<RationalFunction<degree, dimension>>& curve_segments,
  const CurveDiscretizationParameters& curve_disc_params,
  std::vector<Eigen::Matrix<double, 1, dimension>>& points,
  std::vector<std::vector<int>>& polylines)
{
  points.clear();
  polylines.clear();
  for (size_t k = 0; k < curve_segments.size(); ++k) {
    RationalFunction<degree, dimension> const& curve_segment =
      curve_segments[k];

    // Write curves
    int num_samples = curve_disc_params.num_samples;
    std::vector<Eigen::Matrix<double, 1, dimension>> points_k;
    curve_segment.sample_points(num_samples, points_k);

    // Build polyline for the given curve
    std::vector<int> polyline;
    polyline.resize(points_k.size());
    for (size_t l = 0; l < points_k.size(); ++l) {
      polyline[l] = points.size() + l;
    }

    append(points, points_k);
    polylines.push_back(polyline);
  }
}

/// Discretize the boundaries of a spline surface patches.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[out] points: points of the boundary curve network
/// @param[out] polylines: polyline indices of the boundary curve network
void
discretize_patch_boundaries(const QuadraticSplineSurface& spline_surface,
                            std::vector<SpatialVector>& points,
                            std::vector<std::vector<size_t>>& polylines);

/// Sample the Darboux frames of the contours on the surface.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] contour_domain_curve_segments: local parametric domain contour segments
/// @param[in] contour_segments: surface contour segments
/// @param[in] contour_patch_indices: spline surface patch indices for the contour segments
/// @param[in] curve_disc_params: parameters for the contour discretization
/// @param[out] base_points: sampled points on the contours
/// @param[out] tangents: tangents to the contours at the sampled points
/// @param[out] normals: surface normals at the sampled points
/// @param[out] tangent_normals: cross products of normals and tangents at the sampled points
void
sample_contour_frames(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<Conic>& contour_domain_curve_segments,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  const CurveDiscretizationParameters& curve_disc_params,
  std::vector<SpatialVector>& base_points,
  std::vector<SpatialVector>& tangents,
  std::vector<SpatialVector>& normals,
  std::vector<SpatialVector>& tangent_normals);
