// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "compute_rational_bezier_curve_intersection.h"

#include "compute_contours.h"
#include "intersection_data.h"
#include "intersection_heuristics.h"
#include "rational_function.h"

/// \file compute_intersections.h
///
/// Methods to compute intersections for quadratic surfaces.

// Parameters for intersection computations
struct IntersectionParameters
{
  bool use_heuristics =
    true; // If true, use heuristics to check if there are no intersections
  double trim_amount = 1e-5; // Amount to trim ends of contour segments by; intersections
    // that are trimmed are clamped to the endpoint.
};

/// Compute the intersections between the two planar curve segments using existing
/// bounding box and control point data.
///
/// @param[in] first_planar_curve: first planar curve segment
/// @param[in] second_planar_curve: second planar curve segment
/// @param[in] intersect_params: parameters for the intersection computation
/// @param[out] first_curve_intersections: list of intersection knots in the
/// first curve
/// @param[out] second_curve_intersections: list of intersection knots in the
/// second curve
/// @param[out] intersection_stats: diagnostic information about the
/// intersection
/// @param[in] first_bounding_box: bounding box for the first curve
/// @param[in] second_bounding_box: bounding box for the second curve
/// @param[in] first_bezier_control_points: control points for the first curve
/// @param[in] second_bezier_control_points: control points for the second curve
void
compute_planar_curve_intersections(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  const IntersectionParameters& intersect_params,
  std::vector<double>& first_curve_intersections,
  std::vector<double>& second_curve_intersections,
  IntersectionStats& intersection_stats,
  const std::pair<PlanarPoint, PlanarPoint>& first_bounding_box,
  const std::pair<PlanarPoint, PlanarPoint>& second_bounding_box,
  const Eigen::Matrix<double, 5, 3>& first_bezier_control_points,
  const Eigen::Matrix<double, 5, 3>& second_bezier_control_points);

/// Compute the intersections between the two planar curve segments.
///
/// @param[in] first_planar_curve: first planar curve segment
/// @param[in] second_planar_curve: second planar curve segment
/// @param[in] intersect_params: parameters for the intersection computation
/// @param[out] first_curve_intersections: list of intersection knots in the
/// first curve
/// @param[out] second_curve_intersections: list of intersection knots in the
/// second curve
/// @param[out] intersection_stats: diagnostic information about the
/// intersection
void
compute_planar_curve_intersections(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  const IntersectionParameters& intersect_params,
  std::vector<double>& first_curve_intersections,
  std::vector<double>& second_curve_intersections,
  IntersectionStats& intersection_stats);

/// Compute all intersections between the planar curve image segments.
///
/// For each intersection, the knot in the curve and the other curve it
/// corresponds to are recorded.
///
/// @param[in] image_segments: list of planar curve segments
/// @param[in] intersect_params: parameters for the intersection computation
/// @param[out] intersections: list of lists of intersection knots
/// @param[out] intersection_indices: list of lists of intersection indices
/// @param[out] contour_intersection: list of lists of full intersection data
/// @param[in, out] num_intersections: total number of intersections to increment
/// TODO What is intersection_call?
void
compute_intersections(
  const std::vector<RationalFunction<4, 2>>& image_segments,
  const IntersectionParameters& intersect_params,
  std::vector<std::vector<double>>& intersections,
  std::vector<std::vector<size_t>>& intersection_indices,
  std::vector<std::vector<IntersectionData>>& contour_intersections,
  int& num_intersections,
  long long& intersection_call);

/// Get planar curve split parameters until the refined curves cannot have self
/// intersections
///
/// @param[in] planar_curves: list of planar curve segments
/// @param[out] split_points: points to split curves at to avoid self
/// intersections
void
split_planar_curves_no_self_intersection(
  const std::vector<RationalFunction<4, 2>>& planar_curves,
  std::vector<std::vector<double>>& split_points);
