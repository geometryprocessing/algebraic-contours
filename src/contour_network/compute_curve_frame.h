// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"

/// \file compute_curve_frame.h
///
/// Methods to compute a curve aligned frame for quadratic and spline surfaces.

/// @brief Compute the frame aligned to a curve on a quadratic surface.
///
/// @param[in] surface_mapping_coeffs: coefficients for the quadratic surface
/// @param[in] normal_mapping_coeffs: coefficients for the quadratic surface
/// normal
/// @param[in] domain_curve_segment: curve segment in the parametric domain of
/// the surface
/// @param[out] surface_curve_tangent: function for the tangent on the curve
/// @param[out] surface_curve_normal: function for the surface normal on the
/// curve
/// @param[out] surface_curve_tangent_normal: function for the tangent normal on
/// the curve
void
compute_quadratic_surface_curve_frame(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_mapping_coeffs,
  const Conic& domain_curve_segment,
  RationalFunction<8, 3>& surface_curve_tangent,
  RationalFunction<4, 3>& surface_curve_normal,
  RationalFunction<12, 3>& surface_curve_tangent_normal);

/// Compute the frame aligned to a given curve on a quadratic spline surface patch.
///
/// @param[in] spline_surface_patch: quadratic spline surface patch
/// @param[in] domain_curve_segment: curve segment in the parametric domain of
/// the surface
/// @param[in] patch_index: index for the patch containing the curve
/// @param[out] surface_curve_tangent: function for the tangent on the curve
/// @param[out] surface_curve_normal: function for the surface normal on the
/// curve
/// @param[out] surface_curve_tangent_normal: function for the tangent normal on
/// the curve
void
compute_spline_surface_patch_curve_frame(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Conic& domain_curve_segment,
  RationalFunction<8, 3>& surface_curve_tangent,
  RationalFunction<4, 3>& surface_curve_normal,
  RationalFunction<12, 3>& surface_curve_tangent_normal);
