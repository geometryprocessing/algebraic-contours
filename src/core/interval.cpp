// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "interval.h"

void
Interval::set_lower_bound(double lower_bound, bool is_open)
{
  m_t0 = lower_bound;
  m_bounded_below = true;
  m_open_below = is_open;
}

void
Interval::set_upper_bound(double upper_bound, bool is_open)
{
  m_t1 = upper_bound;
  m_bounded_above = true;
  m_open_above = is_open;
}

void
Interval::trim_lower_bound(double trim_amount)
{
  // Don't trim if the interval is not bounded or has length of order trim
  // amount
  if (!is_bounded_below() || (2.0 * std::abs(trim_amount) > get_length()))
    return;
  m_t0 += trim_amount;
}

void
Interval::trim_upper_bound(double trim_amount)
{
  // Don't trim if the interval is not bounded or has length of order trim
  // amount
  if (!is_bounded_above() || (2.0 * std::abs(trim_amount) > get_length()))
    return;
  m_t1 -= trim_amount;
}

void
Interval::pad_lower_bound(double pad_amount)
{
  // Don't trim if the interval is not bounded or has length of order trim
  // amount
  if (!is_bounded_below() || (pad_amount < 0))
    return;
  m_t0 -= pad_amount;
}

void
Interval::pad_upper_bound(double pad_amount)
{
  // Don't trim if the interval is not bounded or has length of order trim
  // amount
  if (!is_bounded_above() || (pad_amount < 0))
    return;
  m_t1 += pad_amount;
}

void
Interval::reset_bounds()
{
  m_t0 = -std::numeric_limits<double>::infinity();
  m_t1 = std::numeric_limits<double>::infinity();
  m_bounded_below = false;
  m_bounded_above = false;
  m_open_below = true;
  m_open_above = true;
}

bool
Interval::is_bounded_above() const
{
  return m_bounded_above;
}

bool
Interval::is_bounded_below() const
{
  return m_bounded_below;
}

bool
Interval::is_open_above() const
{
  return m_open_above;
}

bool
Interval::is_open_below() const
{
  return m_open_below;
}

bool
Interval::is_finite() const
{
  return (is_bounded_above() && is_bounded_below());
}

bool
Interval::is_compact() const
{
  return (is_finite() && !is_open_above() && !is_open_below());
}

double
Interval::get_lower_bound() const
{
  return m_t0;
}

double
Interval::get_upper_bound() const
{
  return m_t1;
}

double
Interval::get_center() const
{
  return 0.5 * (m_t0 + m_t1);
}

double
Interval::get_length() const
{
  if (is_finite()) {
    return get_upper_bound() - get_lower_bound();
  } else {
    return std::numeric_limits<double>::infinity();
  }
}

bool
Interval::contains(double t) const
{
  if ((!m_bounded_below) && (!m_bounded_above)) {
    return true;
  }

  if ((!m_bounded_below) && (t <= m_t1)) {
    return true;
  }

  if ((!m_bounded_above) && (t >= m_t0)) {
    return true;
  }

  if ((t < m_t0) && (!m_open_below)) {
    return false;
  }

  if ((t <= m_t0) && (m_open_below)) {
    return false;
  }

  if ((t > m_t1) && (!m_open_above)) {
    return false;
  }

  if ((t >= m_t1) && (m_open_above)) {
    return false;
  }

  return true;
}

bool
Interval::is_in_interior(double t) const
{
  if (t <= m_t0)
    return false;
  if (t >= m_t1)
    return false;

  return true;
}

std::vector<double>
Interval::sample_points(int num_points) const
{
  std::vector<double> points;
  points.reserve(num_points);
  // FIXME Implement trivial case and infinite case
  if (num_points <= 1) {
    return points;
  }

  // Unbounded
  if ((!is_bounded_above()) && (!is_bounded_below())) {
    for (int i = 0; i < num_points; ++i) {
      points.push_back(-num_points / 20.0 + 0.1 * i);
    }
    return points;
  }
  // Unbounded below, closed above
  else if ((!is_bounded_below()) && (!is_open_above())) {
    for (int i = 0; i < num_points; ++i) {
      points.push_back(get_upper_bound() - 0.1 * i);
    }
    return points;
  }
  // Unbounded below, open above
  else if ((!is_bounded_below()) && (is_open_above())) {
    for (int i = 1; i <= num_points; ++i) {
      points.push_back(get_upper_bound() - 0.1 * i);
    }
    return points;
  }
  // Unbounded above, closed below
  else if ((!is_bounded_above()) && (!is_open_below())) {
    for (int i = 0; i < num_points; ++i) {
      points.push_back(get_lower_bound() + 0.1 * i);
    }
    return points;
  }
  // Unbounded above, open below
  else if ((!is_bounded_above()) && (is_open_below())) {
    SPDLOG_TRACE("Unbounded above ({}, infty)", get_upper_bound());
    for (int i = 1; i <= num_points; ++i) {
      points.push_back(get_lower_bound() + 0.1 * i);
    }
    return points;
  }

  // Bounded cases
  double t0 = get_lower_bound();
  double t1 = get_upper_bound();

  // Closed
  if ((!is_open_below()) && (!is_open_above())) {
    for (int i = 0; i < num_points; ++i) {
      double s = static_cast<float>(i) / static_cast<float>(num_points - 1);
      double t = (1.0 - s) * t0 + s * t1;
      points.push_back(t);
    }
  }
  // Open below, closed above
  else if ((is_open_below()) && (!is_open_above())) {
    for (int i = 1; i <= num_points; ++i) {
      double s = static_cast<float>(i) / static_cast<float>(num_points);
      double t = (1.0 - s) * t0 + s * t1;
      assert(t != t0);
      points.push_back(t);
    }
  }
  // Closed below, open above
  else if ((!is_open_below()) && (is_open_above())) {
    for (int i = 0; i < num_points; ++i) {
      double s = static_cast<float>(i) / static_cast<float>(num_points);
      double t = (1.0 - s) * t0 + s * t1;
      assert(t != t1);
      points.push_back(t);
    }
  }
  // Open
  else if ((is_open_below()) && (is_open_above())) {
    for (int i = 1; i <= num_points; ++i) {
      double s = static_cast<float>(i) / static_cast<float>(num_points + 1);
      double t = (1.0 - s) * t0 + s * t1;
      points.push_back(t);
    }
  }

  return points;
}

std::string
Interval::formatted_interval() const
{
  std::stringstream interval_string;
  // Unbounded case
  if ((!is_bounded_above()) && (!is_bounded_below())) {
    interval_string << "(-infty, infty)";
  } else if ((!is_bounded_below()) && (!is_open_above())) {
    interval_string << "(-infty, " << get_upper_bound() << "]";
  } else if ((!is_bounded_below()) && (is_open_above())) {
    interval_string << "(-infty, " << get_upper_bound() << ")";
  } else if ((!is_bounded_above()) && (!is_open_below())) {
    interval_string << "[" << get_lower_bound() << ", infty)";
  } else if ((!is_bounded_above()) && (is_open_below())) {
    interval_string << "(" << get_lower_bound() << ", infty)";
  }

  // Bounded case
  double t0 = get_lower_bound();
  double t1 = get_upper_bound();
  if ((!is_open_below()) && (!is_open_above())) {
    interval_string << std::fixed << std::setprecision(17) << "[" << t0 << ", "
                    << t1 << "]";
  } else if ((is_open_below()) && (!is_open_above())) {
    interval_string << std::fixed << std::setprecision(17) << "(" << t0 << ", "
                    << t1 << "]";
  } else if ((!is_open_below()) && (is_open_above())) {
    interval_string << std::fixed << std::setprecision(17) << "[" << t0 << ", "
                    << t1 << ")";
  } else if ((is_open_below()) && (is_open_above())) {
    interval_string << std::fixed << std::setprecision(17) << "(" << t0 << ", "
                    << t1 << ")";
  }

  return interval_string.str();
}
