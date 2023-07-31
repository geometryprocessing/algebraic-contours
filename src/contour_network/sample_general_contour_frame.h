// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"

/// \file sample_general_contour_frame.h
///
/// Methods to sample the frame associated with a general contour function on a
/// given surface.

/// Sample the tangent frame to the curve tau * N = c for some implicit constant
/// c on a quadratic spline surface at point determined by the surface
/// discretization parameters.
///
/// @param[in] control_point_grid: control points defining the spline surface
/// @param[in] frame: camera view frame
/// @param[in] surface_disc_params: parameters defining the discretization of
/// the surface
/// @param[out] tangents: sampled tangents to a contour function level set
/// @param[out] normals: sampled surface normals
/// @param[out] tangent_normals: sampled tangent plane vectors perpendicular to
/// the tangent vectors
// void
// sample_spline_surface_general_contour_frame(
//   const QuadraticSplineSurface& spline_surface,
//   const MatrixXr& frame,
//   const SurfaceDiscretizationParameters& surface_disc_params,
//   MatrixXr& tangents,
//   MatrixXr& normals,
//   MatrixXr& tangent_normals);
