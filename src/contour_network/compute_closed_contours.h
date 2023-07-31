// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "rational_function.h"

/// \file compute_closed_contours.h
///
/// Methods to chain contour segments into closed contours.

/// @brief Given contour segments on a surface, chain them together to generate the
/// full closed contours.
///
/// @param[in] contour_segments: surface contour segments
/// @param[out] contours: list of indices of segments for complete surface
/// contours
/// @param[out] contour_labels: index of the contour corresponding to each
/// segment
void
compute_closed_contours(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<std::vector<int>>& contours,
  std::vector<int>& contour_labels);

/// @brief Given closed contours, write the distance between the endpoints of
/// the segments to file.
///
/// @param[in] contour_segments: surface contour segments
/// @param[in] contours: list indices of segments for complete surface contours
/// @param[in] filename: file to write contour distances to
void
write_contour_closure_distances(
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<int>>& contours,
  const std::string& filename);
