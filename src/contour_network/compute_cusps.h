// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"

/// \file compute_cusps.h
///
/// Methods to compute cusps for a quadratic spline surface.

/// @brief Compute an implicit function for a quadratic surface with cusps at the roots.
///
/// Note that this function is high degree and may have spurious cusps.
/// compute_spline_surface_cusps is preferred for finding cusps.
///
/// @param[in] surface_mapping_coeffs: coefficients for the quadratic surface
/// @param[in] normal_mapping_coeffs: coefficients for the quadratic surface
/// normal
/// @param[in] frame: projection frame
/// @param[in] contour_domain_curve_segment: local parametric domain contour segment
/// @param[out] cusp_functions: implicit cusp function
void
compute_quadratic_surface_cusp_function(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_mapping_coeffs,
  const Matrix3x3r& frame,
  const Conic& contour_domain_curve_segment,
  RationalFunction<12, 1>& cusp_function);

/// @brief Compute implicit functions per patch for a spline surface with cusps at the roots.
///
/// Note that these functions are high degree and may have spurious cusps.
/// compute_spline_surface_cusps is preferred for finding cusps.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] frame: projection frame
/// @param[in] contour_domain_curve_segments: local parametric domain contour segments
/// @param[in] patch_indices: spline surface patch indices for the contour segments
/// @param[out] cusp_functions: implicit cusp functions per contour segment
void
compute_spline_surface_cusp_functions(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const std::vector<Conic>& contour_domain_curve_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& patch_indices,
  std::vector<RationalFunction<12, 1>>& cusp_functions);

/// @brief Compute interior and boundary cusps per patch for a spline surface.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] contour_domain_curve_segments: local parametric domain contour segments
/// @param[in] contour_segments: surface contour segments
/// @param[in] patch_indices: spline surface patch indices for the contour segments
/// @param[in] closed_contours: list of indices of segments for complete surface
/// contours
/// @param[out] interior_cusps: paramater points of interior cusps per contour segment
/// @param[out] boundary_cusps: paramater points of boundary cusps per contour segment
/// @param[out] has_cusp_at_base: boolean per contour segment indicating if a cusp is at the base
/// @param[out] has_cusp_at_tip: boolean per contour segment indicating if a cusp is at the tip
void
compute_spline_surface_cusps(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<Conic>& contour_domain_curve_segments,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& patch_indices,
  const std::vector<std::vector<int>>& closed_contours,
  std::vector<std::vector<double>>& interior_cusps,
  std::vector<std::vector<double>>& boundary_cusps,
  std::vector<bool>& has_cusp_at_base,
  std::vector<bool>& has_cusp_at_tip);
