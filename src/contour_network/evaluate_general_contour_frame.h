// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"

/// \file evaluate_general_contour_frame.h
///
/// Methods to evaluate the frame for general contours anywhere on a surface

/// Evaluate the tangent frame to the curve tau * N = c for some implicit
/// constant c on a quadratic spline surface at some global coordinates.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] frame: camera projection frame
/// @param[in] domain_point: point in the patch domain to evaluate the tangent
/// at
/// @param[in] patch_index: index of the patch of the quadratic spline surface
/// @param[out] tangent: tangent to a contour function level set at the given
/// input point
/// @param[out] normal: surface normal at the given input point
/// @param[out] tangent_normal: tangent plane vector perpendicular to the
/// tangent vector
void
evaluate_spline_surface_general_contour_frame(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const PlanarPoint& domain_point,
  const QuadraticSplineSurface::PatchIndex& patch_index,
  SpatialVector& tangent,
  SpatialVector& normal,
  SpatialVector& tangent_normal);