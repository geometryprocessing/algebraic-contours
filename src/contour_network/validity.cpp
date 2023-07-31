// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "validity.h"

bool
is_valid_spatial_mapping(const Matrix6x3r& spatial_mapping_coeffs)
{
  // Must map to R3
  if (spatial_mapping_coeffs.cols() != 3)
    return false;

  return true;
}

bool
is_valid_frame(const Matrix3x3r& frame)
{
  // Must be a 3x3 matrix
  // if (frame.rows() != 3) return false;
  // if (frame.cols() != 3) return false;

  // Must have determinant 1
  if (!float_equal(frame.determinant(), 1.0))
    return false;

  return true;
}

bool
are_valid_contour_segments(
  const std::vector<RationalFunction<4, 3>>& contour_segments)
{
  for (size_t i = 0; i < contour_segments.size(); ++i) {
    if (!contour_segments[i].domain().is_finite())
      return false;
    if (!contour_segments[i].domain().is_compact())
      return false;
  }

  return true;
}
