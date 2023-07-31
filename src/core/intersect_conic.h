// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"
#include "convex_polygon.h"

/// @brief Intersect a parametrized conic segment with an implicit line.
///
/// @param[in] C_param: parametrized conic segment
/// @param[in] L_coeffs: implicit line function
/// @param[out] intersections: list of intersection points in the domain of the
/// conic
/// @return true iff TODO
bool
intersect_conic_with_line(const Conic& C_param,
                          const Eigen::Matrix<double, 3, 1>& L_coeffs,
                          std::vector<double>& intersections);

/// @brief Intersect a parametrized conic segment with a convex polygon and
/// split the conic into segments between intersection points.
///
/// The indices of the line segments of the convex polygon bounding the conic
/// segments (or -1 for open segments) are also recorded
///
/// @param[in] conic: parametrized conic segment
/// @param[in] convex_polygon: convex polygon to intersect
/// @param[out] conic_segments: parametrized conic split at intersection
/// @param[out] line_intersection_indices: indices of the bounding polygon edges
void
intersect_conic_with_convex_polygon(
  const Conic& conic,
  const ConvexPolygon& convex_polygon,
  std::vector<Conic>& conic_segments,
  std::vector<std::pair<int, int>>& line_intersection_indices);

/// @brief Intersect a parametrized conic segment in a cone patch with a convex
/// polygon and split the conic into segments between the intersection points.
///
/// The indices of the line segments of the convex polygon bounding the conic
/// segments (or -1 for open segments) are also recorded
///
/// @param[in] conic: parametrized conic segment
/// @param[in] convex_polygon: convex polygon to intersect
/// @param[in] cone_corner_index: index of the corner where the cone is located
/// @param[out] conic_segment: parametrized conic split at intersection
/// @param[out] line_intersection_indices: indices of the bounding polygon edges
/// @return true iff an intersection is found
bool
intersect_conic_in_cone_patch(const Conic& conic,
                              const ConvexPolygon& convex_polygon,
                              size_t cone_corner_index,
                              Conic& conic_segment,
                              std::pair<int, int>& line_intersection_indices);

/// Check if a parameterized line conic segment that passes through a cone
/// corner in a cone patch intersects the interior of the triangle
///
/// @param[in] conic: parametrized conic segment
/// @param[in] convex_polygon: convex polygon to intersect
/// @param[in] cone_corner_index: index of the corner where the cone is located
/// @return true iff the line through the conic stably intersects the interior
/// of the domain
bool
check_if_conic_intersects_cone_patch_domain(const Conic& conic,
                                            const ConvexPolygon& convex_polygon,
                                            size_t cone_corner_index);
