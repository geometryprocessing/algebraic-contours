// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_rational_bezier_curve_intersection.h"

void
map_to(Interval_ink& J, Interval_ink const& I)
{
  J.setEnds(J.valueAt(I.min()), J.valueAt(I.max()));
}

void
map_to(Interval_ink& J, OptInterval const& I)
{
  J.setEnds(J.valueAt(I.min()), J.valueAt(I.max()));
}

bool
float_equal_ink(const double& a, const double& b, double eps = 1e-9)
{
  if (abs(a - b) < eps) {
    return true;
  }
  return false;
}

bool
Point_equal(const Point& a, const Point& b, double eps = 1e-6)
{
  // homogeneous points x, y, w
  Eigen::Vector2d aa(a[0] / a[2], a[1] / a[2]), bb(b[0] / b[2], b[1] / b[2]);
  if ((aa - bb).norm() < eps)
    return true;
  return false;
}

Point
lerp(double t, Point const& a, Point const& b)
{
  return (1 - t) * a + t * b;
}

// This is actually position comparison, should get to same weight
Point
middle_point(Point const& p1, Point const& p2)
{
  Point p1s(p1[0] / p1[2], p1[1] / p1[2], 1);
  Point p2s(p2[0] / p2[2], p2[1] / p2[2], 1);
  return lerp(0.5, p1s, p2s);
}

size_t
get_precision(Interval_ink const& I)
{
  double d = I.extent();
  double e = 0.1, p = 10;
  int n = 0;
  while (n < 16 && d < e) {
    p *= 10;
    e = 1 / p;
    ++n;
  }
  return n;
}

bool
is_constant(std::vector<Point> const& A, double precision = 1e-6)
{
  for (unsigned int i = 1; i < A.size(); ++i) {
    if (!Point_equal(A[i], A[0], precision))
      return false;
  }
  return true;
}

bool
is_constant(double P[5][3], double precision = 1e-6)
{
  for (unsigned int i = 1; i < 5; ++i) {
    if (!Point_equal(Point(P[0][0], P[0][1], P[0][2]),
                     Point(P[i][0], P[i][1], P[i][2]),
                     precision))
      return false;
  }
  return true;
}

void
left_portion(double t, std::vector<Point>& B)
{
  size_t n = B.size();
  for (size_t i = 1; i < n; ++i) {
    for (size_t j = n - 1; j > i - 1; --j) {
      B[j] = lerp(t, B[j - 1], B[j]);
    }
  }
}

void
right_portion(double t, std::vector<Point>& B)
{
  size_t n = B.size();
  for (size_t i = 1; i < n; ++i) {
    for (size_t j = 0; j < n - i; ++j) {
      B[j] = lerp(t, B[j], B[j + 1]);
    }
  }
}

void
portion(std::vector<Point>& B, Interval_ink const& I)
{
  if (float_equal_ink(I.min(), 0.0)) {
    if (float_equal_ink(I.max(), 1.0))
      return;
    left_portion(I.max(), B);
    return;
  }
  right_portion(I.min(), B);
  if (float_equal_ink(I.max(), 1.0))
    return;
  double t = I.extent() / (1 - I.min());
  left_portion(t, B);
}

void
portion(std::vector<Point>& B, OptInterval const& I)
{
  if (float_equal_ink(I.min(), 0.0)) {
    if (float_equal_ink(I.max(), 1.0))
      return;
    left_portion(I.max(), B);
    return;
  }
  right_portion(I.min(), B);
  if (float_equal_ink(I.max(), 1.0))
    return;
  double t = I.extent() / (1 - I.min());
  left_portion(t, B);
}

// clip fatline
double
dot3(double v[3], double w[3])
{
  return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]);
}

void
cross3(double v[3], double w[3], double res[3])
{
  res[0] = v[1] * w[2] - v[2] * w[1];
  res[1] = -v[0] * w[2] + v[2] * w[0];
  res[2] = v[0] * w[1] - v[1] * w[0];
  return;
}

void
fatline(double Q[5][3], double Lmin[3], double Lmax[3])
{
  int i;
  double c;
  int nQ;
  double cmin;
  double cmax;
  double L[3];
  nQ = 5;
  L[0] = 0;
  L[1] = 0;
  L[2] = 0;
  cross3(Q[0], Q[nQ - 1], L);
  cmin = L[2];
  cmax = L[2];
  for (i = 2; i <= nQ; i++) {
    c = -Q[i - 1][0] * L[0] - Q[i - 1][1] * L[1];
    cmin = fmin(cmin, c);
    cmax = fmax(cmax, c);
  }
  Lmin[0] = -L[0];
  Lmin[1] = -L[1];
  Lmin[2] = -cmin;
  Lmax[0] = L[0];
  Lmax[1] = L[1];
  Lmax[2] = cmax;
  return;
}

void
clipline(double P[5][3], double L[3], double clip_range[2])
{
  double E0[3];
  double En[3];
  double Ei[3];
  double L0i[3];
  double Lni[3];
  double V0i[3];
  double Vni[3];
  double v0x[5];
  double vnx[5];
  double ey[3];
  double t1;
  double t2;
  int nP;
  double nPf;
  int j;
  int i;
  L0i[0] = 0;
  L0i[1] = 0;
  L0i[2] = 0;
  Lni[0] = 0;
  Lni[1] = 0;
  Lni[2] = 0;
  V0i[0] = 0;
  V0i[1] = 0;
  V0i[2] = 0;
  Vni[0] = 0;
  Vni[1] = 0;
  Vni[2] = 0;
  v0x[0] = 0;
  v0x[1] = 0;
  v0x[2] = 0;
  v0x[3] = 0;
  v0x[4] = 0;
  vnx[0] = 0;
  vnx[1] = 0;
  vnx[2] = 0;
  vnx[3] = 0;
  vnx[4] = 0;
  nP = 5;
  v0x[0] = 0.0e0;
  vnx[nP - 1] = 0.0e0;
  E0[0] = 0;
  E0[1] = dot3(L, P[0]);
  E0[2] = 1;
  En[0] = 1;
  En[1] = dot3(L, P[nP - 1]);
  En[2] = 1;
  ey[0] = 0;
  ey[1] = 1;
  ey[2] = 0;
  nPf = (double)(nP - 1);
  for (i = 1; i <= nP; i++) {
    Ei[0] = (double)(i - 1) / nPf;
    Ei[1] = dot3(L, P[i - 1]);
    Ei[2] = 1;
    cross3(E0, Ei, L0i);
    cross3(En, Ei, Lni);
    if (1 < i) {
      cross3(L0i, ey, V0i);
      v0x[i - 1] = V0i[0] / V0i[2];
    }
    if (i < nP) {
      cross3(Lni, ey, Vni);
      vnx[i - 1] = Vni[0] / Vni[2];
    }
  }
  if (0.0e0 <= E0[1])
    t1 = 0.0e0;
  else {
    t1 = 0.1e1;
    for (j = 2; j <= nP; j++)
      if (0.0e0 < v0x[j - 1])
        t1 = fmin(t1, v0x[j - 1]);
  }
  if (0.0e0 <= En[1])
    t2 = 0.1e1;
  else {
    t2 = 0.0e0;
    for (j = 1; j <= nP - 1; j++)
      if (vnx[j - 1] < 0.1e1)
        t2 = fmax(t2, vnx[j - 1]);
  }
  if (t2 <= t1) {
    clip_range[0] = -0.1e1;
    clip_range[1] = -0.1e1;
  } else {
    clip_range[0] = t1;
    clip_range[1] = t2;
  }
  return;
}

void
clipfatline(double P[5][3],
            double Q[5][3],
            double clip_range[2],
            double precision = 1e-8)
{
  for (int i = 0; i < 5; i++) {
    //        P[i][0]/=P[i][2];
    //        P[i][1]/=P[i][2];
    //        P[i][2] = 1;
    Q[i][0] /= Q[i][2];
    Q[i][1] /= Q[i][2];
    Q[i][2] = 1;
  }
  double Lmin[3];
  double Lmax[3];
  double tmin[2];
  double tmax[2];
  Lmin[0] = 0;
  Lmin[1] = 0;
  Lmin[2] = 0;
  Lmax[0] = 0;
  Lmax[1] = 0;
  Lmax[2] = 0;
  tmin[0] = 0;
  tmin[1] = 0;
  tmax[0] = 0;
  tmax[1] = 0;
  if (is_constant(Q, precision)) {
    if (is_constant(P, precision)) {
      clip_range[0] = 0;
      clip_range[1] = 1;
      return;
    } else {
      Eigen::Vector2d direction(-(P[4][1] / P[4][2] - P[0][1] / P[0][2]),
                                P[4][0] / P[4][2] - P[0][0] / P[0][2]);
      direction = direction.normalized();
      Lmin[0] = direction[1];
      Lmin[1] = -direction[0];
      Lmin[2] = -direction[1] * (Q[0][0] + Q[4][0]) / 2 +
                direction[0] * (Q[0][1] + Q[4][1]) / 2;
      Lmax[0] = -direction[1];
      Lmax[1] = direction[0];
      Lmax[2] = direction[1] * (Q[0][0] + Q[4][0]) / 2 -
                direction[0] * (Q[0][1] + Q[4][1]) / 2;
    }
  } else {
    fatline(Q, Lmin, Lmax);
  }

  clipline(P, Lmin, tmin);
  clipline(P, Lmax, tmax);
  if (tmin[0] == -0.1e1 || tmax[0] == -0.1e1) {
    clip_range[0] = -1;
    clip_range[1] = -1;
    ;
  } else {
    clip_range[0] = fmax(tmax[0], tmin[0]);
    clip_range[1] = fmin(tmin[1], tmax[1]);
    ;
  }
}

// inkscope code
void
iterate(std::vector<Interval_ink>& domsA,
        std::vector<Interval_ink>& domsB,
        std::vector<Point> const& A,
        std::vector<Point> const& B,
        Interval_ink const& domA,
        Interval_ink const& domB,
        double precision)
{
  // in order to limit recursion
  static size_t counter = 0;
  if (float_equal_ink(domA.extent(), 1.0) &&
      float_equal_ink(domB.extent(), 1.0))
    counter = 0;
  if (++counter > 100)
    return;

  if (precision < 1e-8)
    precision = 1e-8;

  std::vector<Point> pA = A;
  std::vector<Point> pB = B;
  std::vector<Point>* C1 = &pA;
  std::vector<Point>* C2 = &pB;

  Interval_ink dompA = domA;
  Interval_ink dompB = domB;
  Interval_ink* dom1 = &dompA;
  Interval_ink* dom2 = &dompB;

  std::shared_ptr<OptInterval> dom;

  if (is_constant(A, precision) && is_constant(B, precision)) {
    Point M1 = middle_point(C1->front(), C1->back());
    Point M2 = middle_point(C2->front(), C2->back());
    if (Point_equal(M1, M2)) {
      domsA.push_back(domA);
      domsB.push_back(domB);
    }
    return;
  }

  size_t iter = 0;
  while (++iter < 100 &&
         (dompA.extent() >= precision || dompB.extent() >= precision)) {

    double P[5][3], Q[5][3], clip_range[2];
    assert((*C1).size() == 5);
    assert((*C2).size() == 5);
    for (int i = 0; i < 5; i++) {
      P[i][0] = (*C1)[i][0];
      P[i][1] = (*C1)[i][1];
      P[i][2] = (*C1)[i][2];
      Q[i][0] = (*C2)[i][0];
      Q[i][1] = (*C2)[i][1];
      Q[i][2] = (*C2)[i][2];
    }

    clipfatline(Q, P, clip_range, precision);

    dom = std::make_shared<OptInterval>(clip_range[0], clip_range[1]);

    if (float_equal_ink(dom->max(), -1.0) &&
        float_equal_ink(dom->min(), -1.0)) {
      return;
    }

    // all other cases where dom[0] > dom[1] are invalid
    assert(dom->min() <= dom->max());

    map_to(*dom2, *dom);
    portion(*C2, *dom);

    if (is_constant(*C2, precision) && is_constant(*C1, precision)) {
      Point M1 = middle_point(C1->front(), C1->back());
      Point M2 = middle_point(C2->front(), C2->back());

      if (Point_equal(M1, M2))
        break; // append the new interval
      else
        return; // exit without appending any new interval
    }

    // if we have clipped less than 20% than we need to subdive the curve
    // with the largest domain into two sub-curves
    if (dom->extent() > MIN_CLIPPED_SIZE_THRESHOLD) {
      std::vector<Point> pC1, pC2;
      Interval_ink dompC1, dompC2;
      if (dompA.extent() > dompB.extent()) {
        pC1 = pC2 = pA;
        portion(pC1, H1_INTERVAL);
        portion(pC2, H2_INTERVAL);
        dompC1 = dompC2 = dompA;
        map_to(dompC1, H1_INTERVAL);
        map_to(dompC2, H2_INTERVAL);
        iterate(domsA, domsB, pC1, pB, dompC1, dompB, precision);
        iterate(domsA, domsB, pC2, pB, dompC2, dompB, precision);
      } else {
        pC1 = pC2 = pB;
        portion(pC1, H1_INTERVAL);
        portion(pC2, H2_INTERVAL);
        dompC1 = dompC2 = dompB;
        map_to(dompC1, H1_INTERVAL);
        map_to(dompC2, H2_INTERVAL);
        iterate(domsB, domsA, pC1, pA, dompC1, dompA, precision);
        iterate(domsB, domsA, pC2, pA, dompC2, dompA, precision);
      }
      return;
    }

    std::swap(C1, C2);
    std::swap(dom1, dom2);
  }
  domsA.push_back(dompA);
  domsB.push_back(dompB);
}

void
get_solutions(std::vector<std::pair<double, double>>& xs,
              std::vector<Point> const& A,
              std::vector<Point> const& B,
              double precision)
{
  std::pair<double, double> ci;
  std::vector<Interval_ink> domsA, domsB;
  iterate(domsA, domsB, A, B, UNIT_INTERVAL, UNIT_INTERVAL, precision);
  if (domsA.size() != domsB.size()) {
    assert(domsA.size() == domsB.size());
  }
  xs.clear();
  xs.reserve(domsA.size());
  for (size_t i = 0; i < domsA.size(); ++i) {
    ci.first = domsA[i].middle();
    ci.second = domsB[i].middle();
    xs.push_back(ci);
  }
}

void
find_intersections_bezier_clipping(std::vector<std::pair<double, double>>& xs,
                                   std::vector<Point> const& A,
                                   std::vector<Point> const& B,
                                   double precision)
{
  get_solutions(xs, A, B, precision);
}

// split rational curve
bool
check_split_criteria(Point q[5])
{
  double degree_sum = 0;
  Point p[5];

  for (int i = 0; i < 5; i++) {
    // standardlize to w=1
    p[i] = q[i] / q[i][2];
  }
  for (int i = 0; i < 3; i++) {
    Eigen::Vector3d v1 = p[i + 1] - p[i];
    Eigen::Vector3d v2 = p[i + 2] - p[i + 1];
    v1 = v1.normalized();
    v2 = v2.normalized();
    degree_sum += acos(v1.dot(v2));
  }
  if (degree_sum > M_PI)
    return false;
  return true;
}

void
split_bezier_curve(Point p[5], Point q[5], Point r[5], double t)
{
  // De Casteljau's Algo
  // 5th level
  q[0] = p[0];
  r[4] = p[4];
  // 4th level
  q[1] = (1 - t) * p[0] + t * p[1];
  Point m1 = (1 - t) * p[1] + t * p[2];
  Point m2 = (1 - t) * p[2] + t * p[3];
  r[3] = (1 - t) * p[3] + t * p[4];
  // 3rd level
  q[2] = (1 - t) * q[1] + t * m1;
  Point m3 = (1 - t) * m1 + t * m2;
  r[2] = (1 - t) * m2 + t * r[3];
  // 2nd level
  q[3] = (1 - t) * q[2] + t * m3;
  r[1] = (1 - t) * m3 + t * r[2];
  // 1st level
  q[4] = (1 - t) * q[3] + t * r[1];
  r[0] = q[4];
}

void
split_bezier_curve_no_self_intersection(
  Point p[5],
  std::vector<std::array<Point, 5>>& result)
{
  if (check_split_criteria(p)) {
    std::array<Point, 5> q{ { p[0], p[1], p[2], p[3], p[4] } };
    result.push_back(q);
    return;
  }
  Point q[5], r[5];
  split_bezier_curve(p, q, r, 0.5);
  split_bezier_curve_no_self_intersection(q, result);
  split_bezier_curve_no_self_intersection(r, result);
}

// void split_bezier_curve_no_self_intersection(
//     Point p[5], double start, double end, std::vector<double>
//     &split_point_params) {
//   if (check_split_criteria(p)) {
//     return;
//   }
//   Point q[5], r[5];
//   split_bezier_curve(p, q, r, 0.5);
//   split_point_params.push_back((start+end)/2.0);
//   split_bezier_curve_no_self_intersection(q, start, (start+end)/2.0, result);
//   split_bezier_curve_no_self_intersection(r, (start+end)/2.0, end, result);
// }

void
split_bezier_curve_no_self_intersection(Point p[5],
                                        double start,
                                        double end,
                                        std::vector<double>& split_point_params)
{
  if (check_split_criteria(p)) {
    return;
  }
  if (float_equal(start, end)) {
    spdlog::error("Split until interval length zero");
    return;
  }
  Point q[5], r[5];
  split_bezier_curve(p, q, r, 0.5);
  split_point_params.push_back((start + end) / 2.0);
  split_bezier_curve_no_self_intersection(
    q, start, (start + end) / 2.0, split_point_params);
  split_bezier_curve_no_self_intersection(
    r, (start + end) / 2.0, end, split_point_params);
}
