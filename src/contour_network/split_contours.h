// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"

void
split_contours(
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<RationalFunction<4, 2>>& planar_contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  std::vector<bool>& contour_is_boundary,
  std::vector<std::map<std::string, int>> contour_segment_labels);