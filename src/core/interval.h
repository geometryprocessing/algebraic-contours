// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// @brief Representation of an interval in R that can be bounded or unbounded
/// and closed or open on both ends. Also supports sampling and predicates such
/// as compactness and containment.
class Interval
{
public:
  Interval() { reset_bounds(); }

  Interval(double lower_bound, double upper_bound)
  {
    set_lower_bound(lower_bound);
    set_upper_bound(upper_bound);
  }

  void set_lower_bound(double lower_bound, bool is_open = false);

  void set_upper_bound(double upper_bound, bool is_open = false);

  void trim_lower_bound(double trim_amount);
  void trim_upper_bound(double trim_amount);
  void pad_lower_bound(double pad_amount);
  void pad_upper_bound(double pad_amount);

  double get_lower_bound() const;
  double get_upper_bound() const;
  double get_center() const;
  double get_length() const;

  bool is_bounded_below() const;
  bool is_bounded_above() const;

  bool is_open_below() const;
  bool is_open_above() const;

  bool is_finite() const;
  bool is_compact() const;

  void reset_bounds();

  bool contains(double t) const;
  bool is_in_interior(double t) const;

  std::vector<double> sample_points(int num_points) const;

  std::string formatted_interval() const;

private:
  double m_t0;
  double m_t1;
  bool m_bounded_below;
  bool m_bounded_above;
  bool m_open_below;
  bool m_open_above;
};
