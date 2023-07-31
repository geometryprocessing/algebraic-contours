// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

typedef Eigen::Vector3d Point;

// TODO neeed to rewrite these classes, fake inheritance
class Interval_ink
{
public:
  double begin, end;

  Interval_ink() {}
  Interval_ink(double a, double b)
    : begin(a)
    , end(b)
  {
  }
  // Interval_ink(const Interval_ink& I) = default;

  double extent() const { return end - begin; }

  double max() const { return end; }
  double min() const { return begin; }
  double middle() const { return (begin + end) / 2; }
  void setEnds(double a, double b)
  {
    if (a <= b) {
      begin = a;
      end = b;
    } else {
      begin = b;
      end = a;
    }
  }

  double valueAt(double t) const { return (1 - t) * begin + t * end; }
};

class OptInterval : public Interval_ink
{
public:
  double begin, end;
  OptInterval(){};
  OptInterval(double a, double b)
    : begin(a)
    , end(b)
  {
  }
  OptInterval(const Interval_ink& I)
    : begin(I.begin)
    , end(I.end)
  {
  }
  OptInterval(const OptInterval& I) = default;

  double max() const { return end; }

  double min() const { return begin; }

  double middle() const { return (begin + end) / 2; }

  void setEnds(double a, double b)
  {
    if (a <= b) {
      begin = a;
      end = b;
    } else {
      begin = b;
      end = a;
    }
  }

  double valueAt(double t) const { return (1 - t) * begin + t * end; }

  double extent() const { return end - begin; }
};

// util functions
void
map_to(Interval_ink& J, Interval_ink const& I);
void
map_to(Interval_ink& J, OptInterval const& I);
bool
float_equal_ink(const double& a, const double& b, double eps);
bool
Point_equal(const Point& a, const Point& b, double eps);
Point
lerp(double t, Point const& a, Point const& b);
Point
middle_point(Point const& p1, Point const& p2);
size_t
get_precision(Interval_ink const& I);
bool
is_constant(std::vector<Point> const& A, double precision);
bool
is_constant(double P[5][3]);
void
left_portion(double t, std::vector<Point>& B);
void
right_portion(double t, std::vector<Point>& B);
void
portion(std::vector<Point>& B, Interval_ink const& I);
void
portion(std::vector<Point>& B, OptInterval const& I);

// clip fatline
double
dot3(double v[3], double w[3]);
void
cross3(double v[3], double w[3], double res[3]);
void
fatline(double Q[5][3], double Lmin[3], double Lmax[3]);
void
clipline(double P[5][3], double L[3], double clip_range[2]);
void
clipfatline(double P[5][3], double Q[5][3], double clip_range[2]);

// constants for inkscope code
const double MAX_PRECISION = 1e-8;
const double MIN_CLIPPED_SIZE_THRESHOLD = 0.8;
const Interval_ink UNIT_INTERVAL(0, 1);
const OptInterval EMPTY_INTERVAL;
const Interval_ink H1_INTERVAL(0, 0.5);
const Interval_ink H2_INTERVAL(nextafter(0.5, 1.0), 1.0);

// inkscope code
void
iterate(std::vector<Interval_ink>& domsA,
        std::vector<Interval_ink>& domsB,
        std::vector<Point> const& A,
        std::vector<Point> const& B,
        Interval_ink const& domA,
        Interval_ink const& domB,
        double precision);

void
get_solutions(std::vector<std::pair<double, double>>& xs,
              std::vector<Point> const& A,
              std::vector<Point> const& B,
              double precision);

void
find_intersections_bezier_clipping(std::vector<std::pair<double, double>>& xs,
                                   std::vector<Point> const& A,
                                   std::vector<Point> const& B,
                                   double precision);

// split rational bezier curves into 2

bool
check_split_criteria(Point q[5]);
void
split_bezier_curve(Point p[5], Point q[5], Point r[5], double t);
void
split_bezier_curve_no_self_intersection(
  Point p[5],
  std::vector<std::array<Point, 5>>& result);

void
split_bezier_curve_no_self_intersection(
  Point p[5],
  double start,
  double end,
  std::vector<double>& split_point_params);