// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"

/// \file ray_intersection_heuristics.h
///
/// Methods to accelerate computation of intersections of a ray and a quadratic
/// surface.

/// Compute the 3D bounding boxes for the spline surface patches.
///
/// @param[in] spline_surface: surface to compute patch bounding boxes for
/// @param[out] min_points: minimum coefficients points of the bounding boxes
/// @param[out] max_points: maximum coefficients points of the bounding boxes
void
compute_spline_surface_bounding_boxes(
  const QuadraticSplineSurface& spline_surface,
  std::vector<SpatialVector>& min_points,
  std::vector<SpatialVector>& max_points);