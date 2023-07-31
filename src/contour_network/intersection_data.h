// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// Data for intersections
struct IntersectionData
{
  double knot; // Parameter value for the the intersection in the given curve
  size_t intersection_index; // ID of the intersecting curve
  double intersection_knot; // Parameter of the intersection in the intersecting
                            // curve's domain
  bool is_base = false; // True iff the intersection is at the base of the curve
  bool is_tip = false;  // True iff the intersection is at the tip of the curve
  int id;               // Unique identifier for the intersection
  bool is_redundant = false; // Flag for redundant intersections

  /// Check if the knot is the tip of an oriented curve
  ///
  /// @param[in] domain: domain for the curve
  /// @param[in] eps: epsilon tolerance for the check
  void check_if_tip(const Interval& domain, double eps)
  {
    is_tip = (float_equal(domain.get_upper_bound(), knot, eps));
  }

  /// Check if the knot is the base of an oriented curve
  ///
  /// @param[in] domain: domain for the curve
  /// @param[in] eps: epsilon tolerance for the check
  void check_if_base(const Interval& domain, double eps)
  {
    is_base = (float_equal(domain.get_lower_bound(), knot, eps));
  }
};

// Comparator for Intersection data with respect to knot values
struct knot_less_than
{
  inline bool operator()(const IntersectionData& data_1,
                         const IntersectionData& data_2)
  {
    return (data_1.knot < data_2.knot);
  }
};
