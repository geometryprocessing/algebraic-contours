// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"

/// \file compute_cusps.h
///
/// Methods to compute intersections of a ray and a quadratic surface.

/// Given a quadratic spline surface and a ray, compute the intersections of the
/// ray with each patch.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] ray_mapping_coeffs: coefficients for the linear ray
/// @param[out] patch_indices: indices of the patches with intersections
/// @param[out] surface_intersections: parameters of the intersections on the
/// surface
/// @param[out] ray_intersections: parameters of the intersections on the ray
/// TODO What are the ray_intersection_call and ray_bounding_box_call?
void
compute_spline_surface_ray_intersections(
  const QuadraticSplineSurface& spline_surface,
  const Matrix2x3r& ray_mapping_coeffs,
  std::vector<QuadraticSplineSurface::PatchIndex>& patch_indices,
  std::vector<PlanarPoint>& surface_intersections,
  std::vector<double>& ray_intersections,
  long long& ray_intersections_call,
  long long& ray_bounding_box_call);

/// Given a ray with intersection points and a point on the ray, partition the points into
/// those above and below the given point
///
/// @param[in] ray_mapping_coeffs: coefficients for the linear ray
/// @param[in] comparison_point: point on the ray to compare other intersections
/// against
/// @param[in] ray_intersections: parameters for intersection points on the ray
/// @param[out] ray_intersections_below: intersections below the comparison
/// point
/// @param[out] ray_intersections_above: intersections above the comparison
/// point
void
partition_ray_intersections(const Matrix2x3r& ray_mapping_coeffs,
                            const SpatialVector& comparison_point,
                            const std::vector<double>& ray_intersections,
                            std::vector<double>& ray_intersections_below,
                            std::vector<double>& ray_intersections_above);
