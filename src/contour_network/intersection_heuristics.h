// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "rational_function.h"

#include <chrono>

/// \file intersection_heuristics.h
///
/// Methods to quickly determine if two planar curves do not intersect.

struct IntersectionStats
{
  int num_intersection_tests = 0;
  int num_bezier_nonoverlaps = 0;
  long long bounding_box_call = 0;
  long long intersection_call = 0;
};

/// @brief Compute the bezier control points for a curve over a domain
///
/// The domain may differ from the domain of the planar curve itself.
///
/// @param[in] planar_curve: planar curve segment
/// @param[in] t_min: lower bound of the domain
/// @param[in] t_max: upper bound of the domain
/// @param[out] bezier_points: homogeneous control points of the curve
void
compute_homogeneous_bezier_points_over_interval(
  const RationalFunction<4, 2>& planar_curve,
  double t_min,
  double t_max,
  Eigen::Matrix<double, 5, 3>& bezier_points);

/// @brief Compute a bounding box for a planar curve using Bezier control points for the
/// two subcurves split at the middle
///
/// @param[in] planar_curve: planar curve segment
/// @param[out] lower_left_point: lower left point of the bounding box
/// @param[out] upper_right_point: upper right point of the bounding box
void
compute_bezier_bounding_box(const RationalFunction<4, 2>& planar_curve,
                            PlanarPoint& lower_left_point,
                            PlanarPoint& upper_right_point);

/// @brief Compute bounding boxes for a list of planar curve using Bezier control points.
///
/// @param[in] planar_curves: planar curve segments
/// @param[out] lower_left_points: lower left points of the bounding boxes
/// @param[out] upper_right_points: upper right points of the bounding boxes
void
compute_bezier_bounding_boxes(
  const std::vector<RationalFunction<4, 2>>& planar_curves,
  std::vector<PlanarPoint>& lower_left_points,
  std::vector<PlanarPoint>& upper_right_points);

/// @brief Get the bounding box containing two bounding boxes
///
/// @param[in] first_box_lower_left_point: lower left point of the first bounding box
/// @param[in] first_box_upper_right_point: upper right point of the first bounding box
/// @param[in] second_box_lower_left_point: lower left point of the second bounding box
/// @param[in] second_box_upper_right_point: upper right point of the second bounding box
/// @param[out] lower_left_point: lower left point of the combined bounding box
/// @param[out] upper_right_point: upper right point of the combined bounding box
void
combine_bounding_boxes(const PlanarPoint& first_box_lower_left_point,
                       const PlanarPoint& first_box_upper_right_point,
                       const PlanarPoint& second_box_lower_left_point,
                       const PlanarPoint& second_box_upper_right_point,
                       PlanarPoint& lower_left_point,
                       PlanarPoint& upper_right_point);

/// @brief Check if a test point is in a bounding box with extreme corners lower left
/// point and upper right point
///
/// @param[in] test_point: point to test for containment
/// @param[in] lower_left_point: lower left point of the bounding box
/// @param[in] upper_right_point: upper right point of the bounding box
/// @return true iff the test point is in the bounding box
bool
is_in_bounding_box(const PlanarPoint& test_point,
                   const PlanarPoint& lower_left_point,
                   const PlanarPoint& upper_right_point);

/// @brief Determine if two planar curves cannot intersect by heuristics
///
/// This method can only determine if two curves do not intersect and may have
/// false negatives.
///
/// @param[in] first_planar_curve: first planar curve segment
/// @param[in] second_planar_curve: second planar curve segment
/// @param[in, out] intersection_stats: statistics for intersections computation
/// @return true if the two curves do not intersect
bool
are_nonintersecting_by_heuristic(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  IntersectionStats& intersection_stats);

/// @brief Determine if two planar curves cannot intersect by heuristics using just
/// the bounding boxes of the curve
///
/// This method can only determine if two curves do not intersect and may have
/// false negatives.
///
/// @param[in] first_bounding_box: bounding box of the first planar curve segment
/// @param[in] second_bounding_box: bounding box of the second planar curve segment
/// @param[in, out] intersection_stats: statistics for intersections computation
/// @return true if the two curves do not intersect
bool
are_nonintersecting_by_heuristic(
  const std::pair<PlanarPoint, PlanarPoint>& first_bounding_box,
  const std::pair<PlanarPoint, PlanarPoint>& second_bounding_box,
  IntersectionStats& intersection_stats);

/// @brief Hash bounding boxes by x and y coordinates.
///
/// @param[in] bounding_boxes: bounding boxes to hash
/// @param[out] hash_table: lists of bounding boxes per hash region
/// @param[out] reverse_hash_table: map from bounding boxes to ids of hash
/// regions containing them
void
compute_bounding_box_hash_table(
  const std::vector<std::pair<PlanarPoint, PlanarPoint>>& bounding_boxes,
  std::vector<int> hash_table[50][50],
  std::vector<std::vector<int>>& reverse_hash_table);

/// @brief Determine if the line between the endpoints of two curves intersect.
///
/// Note that linear intersection does not imply curve intersection, and the
/// the absence of linear intersection does not imply the curves do not intersect.
/// It just roughly indicates potential intersections for simple curves.
///
/// @param[in] first_planar_curve: first planar curve segment
/// @param[in] second_planar_curve: second planar curve segment
/// @return true iff the lines between the respective endpoints of the curves intersect
bool
has_linear_intersection(const RationalFunction<4, 2>& first_planar_curve,
                        const RationalFunction<4, 2>& second_planar_curve);
