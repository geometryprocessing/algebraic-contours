// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"
#include "intersection_data.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"

/// \file compute_contours.h
///
/// Methods to compute a contours for quadratic surfaces.

/// Parameters for contour generation
struct ContourParameters
{
  bool closed_surface =
    true; // True if the surface and thus each contour is closed
};

/// @brief Given the coefficients for the normal vector function and a projection
/// frame, compute the coefficients for the implicit contour function.
///
/// @param[in] normal_mapping_coeffs: coefficients for the surface normal
/// function
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[out] contour_equation_coeffs; implicit contour function coefficients
void
compute_contour_equation(const Matrix6x3r& normal_mapping_coeffs,
                         const Matrix3x3r& frame,
                         Vector6r& contour_equation_coeffs);

/// @brief Given a quadratic surface (with normals) defined by the surface and normal
/// mappings, compute the rational contour segments in the parametric domain and
/// on the surface with respect to the given frame
///
/// @param[in] surface_mapping_coeffs: coefficients for the quadratic surface
/// @param[in] normal_mapping_coeffs: coefficients for the quadratic surface
/// normal
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[in] contour_params: parameters for the contour extraction
/// @param[out] contour_domain_curve_segments: parametric domain contour
/// segments
/// @param[out] contour_segments: surface contour segments
void
compute_quadratic_surface_contours(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_mapping_coeffs,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments);

/// @brief Given a quadratic spline surface and a projection frame, compute the
/// rational functions defining the contour segments.
///
/// Both (quadratic) parametric domain contour functions and (higher order)
/// surface contour functions are computed. The indices of the corresponding
/// patches for each segment are also extracted.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[in] contour_params: parameters for the contour extraction
/// @param[out] contour_domain_curve_segments: local parametric domain contour
/// segments
/// @param[out] contour_segments: surface contour segments
/// @param[out] contour_patch_indices: spline surface patch indices for the
/// contour segments
void
compute_spline_surface_contours(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices);
  
/// Given a quadratic spline surface, compute the rational functions defining
/// the surface boundary.
///
/// Both (quadratic) parametric domain boundary functions and (higher order)
/// surface boundary functions are computed. The indices of the corresponding
/// patches for each segment are also extracted.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] patch_boundary_edges: edges of the patch triangle domains that
/// are boundaries
/// @param[out] boundary_domain_curve_segments: local parametric domain boundary
/// segments
/// @param[out] boundary_segments: surface boundary segments
/// @param[out] boundary_patch_indices: spline surface patch indices for the
/// boundary segments
void
compute_spline_surface_boundaries(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<std::pair<int, int>>& patch_boundary_edges,
  std::vector<Conic>& boundary_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& boundary_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& boundary_patch_indices);

/// @brief Given a quadratic spline surface and a projection frame, compute the
/// rational functions defining both the contour segments and the surface
/// boundaries;
///
/// Both (quadratic) parametric domain contour functions and (higher order)
/// surface contour functions are computed. The indices of the corresponding
/// patches for each segment are also extracted and a flag for whether they are
/// boundaries.
///
/// An initial collection of intersections of the boundary and the interior
/// contours are also computed, which exist for any planar projection.
///
/// @param[in] spline_surface: quadratic spline surface
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[out] contour_domain_curve_segments: local parametric domain contour
/// segments
/// @param[out] contour_segments: surface contour segments
/// @param[out] contour_patch_indices: spline surface patch indices for the
/// contour segments
/// @param[out] contour_is_boundary: true iff the patch is a boundary contour
/// @param[out] intersection_data: data for the intersections of the contours
/// and boundaries
/// @param[out] num_intersections: number of intersections
void
compute_spline_surface_contours_and_boundaries(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const std::vector<std::pair<int, int>>& patch_boundary_edges,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  std::vector<bool>& contour_is_boundary,
  std::vector<std::vector<IntersectionData>>& contour_intersections,
  int& num_intersections);

/// @brief Pad contour domains by some given amount in the parametric domain
///
/// @param[out] contour_domain_curve_segments: local parametric domain contour
/// segments
/// @param[out] contour_segments: surface contour segments
/// @param[out] planar_contour_segments: projected planar contour segments
/// @param[in] pad_amount: amount to pad the domain ends by
void
pad_contours(std::vector<Conic>& contour_domain_curve_segments,
             std::vector<RationalFunction<4, 3>>& contour_segments,
             std::vector<RationalFunction<4, 2>>& planar_contour_segments,
             double pad_amount = 0.0);