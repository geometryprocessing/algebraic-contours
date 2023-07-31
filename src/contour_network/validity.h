// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "rational_function.h"

bool
is_valid_spatial_mapping(const Matrix6x3r& spatial_mapping_coeffs);

bool
is_valid_frame(const Matrix3x3r& frame);

bool
are_valid_contour_segments(
  const std::vector<RationalFunction<4, 3>>& contour_segments);
