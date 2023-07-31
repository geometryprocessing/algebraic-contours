// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "intersect_conic.h"

static std::vector<double> intersections_poly;

typedef std::pair<double, int> indexed_intersection;
static std::vector<indexed_intersection> indexed_intersections;

void
swap_line_coordinates(Eigen::Matrix<double, 3, 1>& L_coeffs)
{
  std::swap(L_coeffs(1), L_coeffs(2));
}

void
swap_conic_coordinates(Vector6r& C_coeffs)
{
  std::swap(C_coeffs(1), C_coeffs(2));
  std::swap(C_coeffs(4), C_coeffs(5));
}

void
eliminate_x(Vector6r& C_coeffs, Eigen::Matrix<double, 3, 1>& L_coeffs)
{
  assert(!float_equal(L_coeffs(1), 0.0));

  // Get initial coefficients
  double L0 = L_coeffs(0);
  double Lx = L_coeffs(1);
  double Ly = L_coeffs(2);
  double C0 = C_coeffs(0);
  double Cx = C_coeffs(1);
  double Cy = C_coeffs(2);
  double Cxy = C_coeffs(3);
  double Cxx = C_coeffs(4);
  double Cyy = C_coeffs(5);

  // Perform explicit substitutions
  C_coeffs(0) = C0 - Cx * (L0 / Lx) + Cxx * std::pow((L0 / Lx), 2);
  C_coeffs(1) = 0.0;
  C_coeffs(2) =
    Cy - Cx * (Ly / Lx) - Cxy * (L0 / Lx) - 2 * Cxx * ((L0 * Ly) / (Lx * Lx));
  C_coeffs(3) = 0.0;
  C_coeffs(4) = 0.0;
  C_coeffs(5) = Cyy - Cxy * (Ly / Lx) + Cxx * pow((Ly / Lx), 2);
}

void
eliminate_y(Vector6r& C_coeffs, Eigen::Matrix<double, 3, 1>& L_coeffs)
{
  assert(!float_equal(L_coeffs(2), 0.0));

  // Swap x and y coordinates
  swap_conic_coordinates(C_coeffs);
  swap_line_coordinates(L_coeffs);

  // Eliminate x coordinate in swapped equations
  eliminate_x(C_coeffs, L_coeffs);

  // Swap x and y again
  swap_conic_coordinates(C_coeffs);
  swap_line_coordinates(L_coeffs);
}

double
solve_for_x(Eigen::Matrix<double, 3, 1>& L_coeffs, double y)
{
  assert(!float_equal(L_coeffs(1), 0.0));

  return (-L_coeffs(0) - L_coeffs(2) * y) / L_coeffs(1);
}

// FIXME make intersections fixed size
void
intersect_conic_with_line(Vector6r& C_coeffs,
                          Eigen::Matrix<double, 3, 1>& L_coeffs,
                          MatrixXr& intersections)
{
  // Swap x and y coordinates if necessary so x can be eliminated
  bool swapped = false;
  // if (!float_equal(L_coeffs(1), 0.0)) {
  if (std::abs(L_coeffs(1)) < std::abs(L_coeffs(2))) {
    swap_conic_coordinates(C_coeffs);
    swap_line_coordinates(L_coeffs);
    swapped = true;
  }

  // Eliminate x and solve for TODO Remove non stable
  eliminate_x(C_coeffs, L_coeffs);
  bool use_stable_quadratic = true;
  if (use_stable_quadratic) {
    Eigen::Matrix<double, 3, 1> intersection_quadratic = { C_coeffs(0),
                                                           C_coeffs(2),
                                                           C_coeffs(5) };
    std::array<double, 2> solutions;
    int num_solutions;
    quadratic_real_roots(intersection_quadratic, solutions, num_solutions);
    if (num_solutions == 1) {
      double y0 = solutions[0];
      double x0 = solve_for_x(L_coeffs, y0);
      intersections.resize(1, 2);
      intersections << x0, y0;
    } else if (num_solutions == 2) {
      double y0 = solutions[0];
      double y1 = solutions[1];
      double x0 = solve_for_x(L_coeffs, y0);
      double x1 = solve_for_x(L_coeffs, y1);
      intersections.resize(2, 2);
      intersections << x0, y0, x1, y1;
    }
  } else {
    double discriminant =
      compute_discriminant(C_coeffs(0), C_coeffs(2), C_coeffs(5));
    if (float_equal(discriminant, 0.0)) {
      double y0 = -C_coeffs(2) / (2.0 * C_coeffs(5));
      double x0 = solve_for_x(L_coeffs, y0);
      intersections.resize(1, 2);
      intersections << x0, y0;
    } else if (discriminant > 0.0) {
      double y0 =
        (-C_coeffs(2) - std::sqrt(discriminant)) / (2.0 * C_coeffs(5));
      double y1 =
        (-C_coeffs(2) + std::sqrt(discriminant)) / (2.0 * C_coeffs(5));
      double x0 = solve_for_x(L_coeffs, y0);
      double x1 = solve_for_x(L_coeffs, y1);
      intersections.resize(2, 2);
      intersections << x0, y0, x1, y1;
    }
  }

  // Swap back x and y coordinates if necessary
  if (swapped) {
    swap_conic_coordinates(C_coeffs);
    swap_line_coordinates(L_coeffs);
  }
}

bool
indexed_intersection_comparator(indexed_intersection& intersection_1,
                                indexed_intersection& intersection_2)
{
  return intersection_1.first < intersection_2.first;
}

bool
intersect_conic_with_line(const Conic& C_param,
                          const Eigen::Matrix<double, 3, 1>& L_coeffs,
                          std::vector<double>& intersections)
{
  Matrix3x2r P_coeffs = C_param.get_numerators();
  Eigen::Matrix<double, 3, 1> Q_coeffs = C_param.get_denominator();

  // Get equation for intersection
  Eigen::Matrix<double, 3, 1> X_coeffs = P_coeffs.col(0);
  Eigen::Matrix<double, 3, 1> Y_coeffs = P_coeffs.col(1);
  double a = L_coeffs(0);
  double b = L_coeffs(1);
  double c = L_coeffs(2);
  Eigen::Matrix<double, 3, 1> I_coeffs =
    a * Q_coeffs + b * X_coeffs + c * Y_coeffs;

  bool use_stable_quadratic = true;
  if (use_stable_quadratic) {
    std::array<double, 2> solutions;
    int num_solutions;
    quadratic_real_roots(I_coeffs, solutions, num_solutions);
    if (num_solutions == 0) {
      return false;
    } else if (num_solutions == 1) {
      double t0 = solutions[0];
      if (C_param.is_in_domain(t0)) {
        intersections.push_back(t0);
      }
      return true;
    } else if (num_solutions == 2) {
      double t0 = std::min(solutions[0], solutions[1]);
      double t1 = std::max(solutions[0], solutions[1]);
      assert(t0 <= t1);
      if (C_param.is_in_domain(t0)) {
        intersections.push_back(t0);
      }
      if (C_param.is_in_domain(t1)) {
        intersections.push_back(t1);
      }
      return true;
    }

    return false;
  }

  if (float_equal(I_coeffs(2), 0.0)) {
    if (float_equal(I_coeffs(1), 0.0)) {
      SPDLOG_TRACE("Degenerate intersection polynomial {}",
                   formatted_polynomial(I_coeffs));
      return false;
    }

    double t = -I_coeffs(0) / I_coeffs(1);
    if (C_param.is_in_domain(t)) {
      SPDLOG_TRACE("Linear intersection polynomial {}",
                   formatted_polynomial(I_coeffs));
      SPDLOG_TRACE("Intersection point {}", t);
      intersections.push_back(t);
      return true;
    } else {
      return false;
    }
  }

  double discriminant =
    compute_discriminant(I_coeffs(0), I_coeffs(1), I_coeffs(2));
  if (float_equal(discriminant, 0.0)) {
    double t = -I_coeffs(1) / (2.0 * I_coeffs(2));

    if (C_param.is_in_domain(t)) {
      intersections.push_back(t);
      SPDLOG_TRACE("Singular quadratic intersection polynomial {}",
                   formatted_polynomial(I_coeffs));
      SPDLOG_TRACE("Intersection point {}", t);
      return true;
    } else {
      return false;
    }
  } else if (discriminant > 0.0) {
    double t0 = (-I_coeffs(1) - std::sqrt(discriminant)) / (2.0 * I_coeffs(2));
    double t1 = (-I_coeffs(1) + std::sqrt(discriminant)) / (2.0 * I_coeffs(2));
    if (C_param.is_in_domain(t0)) {
      intersections.push_back(t0);
    }
    if (C_param.is_in_domain(t1))
      intersections.push_back(t1);
    return true;
  }

  return false;
}

bool
are_in_polygon(const std::vector<Conic> conic_segments,
               const ConvexPolygon& convex_polygon)
{
  for (size_t i = 0; i < conic_segments.size(); ++i) {
    if (!convex_polygon.contains(conic_segments[i].mid_point()))
      return false;
  }

  return true;
}

void
intersect_conic_with_convex_polygon(
  const Conic& conic,
  const ConvexPolygon& convex_polygon,
  std::vector<Conic>& conic_segments,
  std::vector<std::pair<int, int>>& line_intersection_indices)
{
  intersections_poly.reserve(100000);
  indexed_intersections.reserve(100000);
  const std::array<Eigen::Matrix<double, 3, 1>, 3>& polygon_boundary_coeffs =
    convex_polygon.get_boundary_segments();

  conic_segments.clear();
  line_intersection_indices.clear();
  double t0;
  double t1;
  double t_sample;
  PlanarPoint p_sample;

  // Get all intersections of the conic with the lines bounding the polygon
  // std::vector<double> intersections;
  intersections_poly.clear();
  indexed_intersections.clear();
  // std::vector<indexed_intersection> indexed_intersections;
  // intersections.reserve(10000);
  for (size_t i = 0; i < polygon_boundary_coeffs.size(); ++i) {
    Eigen::Matrix<double, 3, 1> L_coeffs = polygon_boundary_coeffs[i];
    intersect_conic_with_line(conic, L_coeffs, intersections_poly);

    // Record intersections with the line
    for (size_t j = indexed_intersections.size(); j < intersections_poly.size();
         ++j) {
      indexed_intersections.push_back(std::make_pair(intersections_poly[j], i));
    }
    assert(indexed_intersections.size() == intersections_poly.size());
  }
  SPDLOG_TRACE("Intersections: {}", formatted_vector(intersections_poly, ", "));
  std::sort(indexed_intersections.begin(), indexed_intersections.end());

  if (indexed_intersections.empty()) {
    t0 = conic.domain().get_lower_bound();
    t1 = conic.domain().get_upper_bound();
    if ((conic.domain().is_bounded_above()) &&
        (conic.domain().is_bounded_below())) {
      t_sample = 0.5 * (t0 + t1);
    } else if (conic.domain().is_bounded_below()) {
      t_sample = t0 + 1.0;
    } else if (conic.domain().is_bounded_above()) {
      t_sample = t1 - 1.0;
    } else {
      t_sample = 0.0;
    }

    p_sample = conic(t_sample);
    if (convex_polygon.contains(p_sample)) {
      conic_segments.push_back(conic);
      line_intersection_indices.push_back(std::make_pair(-1, -1));
    }
    return;
  }

  // First segment
  t0 = conic.domain().get_lower_bound();
  t1 = indexed_intersections.front().first;
  t_sample = std::max<double>(0.5 * (t0 + t1), t1 - 1);
  p_sample = conic(t_sample);
  SPDLOG_TRACE("Sampling at {}", t_sample);
  if (convex_polygon.contains(p_sample)) {
    Conic conic_segment(conic);
    conic_segment.domain().set_upper_bound(t1, false);
    conic_segments.push_back(conic_segment);
    line_intersection_indices.push_back(
      std::make_pair(-1, indexed_intersections.front().second));
  }

  // Determine whether conic segments are contained in the polygon
  for (size_t i = 0; i < indexed_intersections.size() - 1; ++i) {
    t0 = indexed_intersections[i].first;
    t1 = indexed_intersections[i + 1].first;
    t_sample = 0.5 * (t0 + t1);
    p_sample = conic(t_sample);
    SPDLOG_TRACE("Sampling at {}", t_sample);
    if (convex_polygon.contains(p_sample)) {
      Conic conic_segment(conic);
      conic_segment.domain().set_lower_bound(t0, false);
      conic_segment.domain().set_upper_bound(t1, false);
      conic_segments.push_back(conic_segment);
      line_intersection_indices.push_back(std::make_pair(
        indexed_intersections[i].second, indexed_intersections[i + 1].second));
    }
  }

  // Last segment
  t0 = indexed_intersections.back().first;
  t1 = conic.domain().get_upper_bound();
  t_sample = std::min<double>(0.5 * (t0 + t1), t0 + 1);
  SPDLOG_TRACE("Sampling at {}", t_sample);
  p_sample = conic(t_sample);
  if (convex_polygon.contains(p_sample)) {
    Conic conic_segment(conic);
    conic_segment.domain().set_lower_bound(t0, false);
    conic_segments.push_back(conic_segment);
    line_intersection_indices.push_back(
      std::make_pair(indexed_intersections.back().second, -1));
  }

  assert(are_in_polygon(conic_segments, convex_polygon));
}

bool
intersect_conic_in_cone_patch(const Conic& conic,
                              const ConvexPolygon& convex_polygon,
                              size_t cone_corner_index,
                              Conic& conic_segment,
                              std::pair<int, int>& line_intersection_indices)
{
  // Get boundary for the edge opposite the cone corner
  const std::array<Eigen::Matrix<double, 3, 1>, 3>& polygon_boundaries =
    convex_polygon.get_boundary_segments();
  size_t polygon_boundary_index = cone_corner_index;
  Eigen::Matrix<double, 3, 1> L_coeffs =
    polygon_boundaries[polygon_boundary_index];

  // Get the intersection with the line
  std::vector<double> intersections;
  intersections.reserve(10000);
  intersect_conic_with_line(conic, L_coeffs, intersections);

  // Check that there are precisely two intersections
  if (intersections.size() > 1) {
    SPDLOG_ERROR("More than two intersections found in cone patch");
    return false;
  }

  if (intersections.size() == 0) {
    SPDLOG_TRACE("No intersection of the ray with the opposing line");
    return false;
  }

  // Split the conic at the intersection (depending on what kind of ray it is)
  conic_segment = Conic(conic);
  double t = intersections[0];
  if (conic_segment.domain().is_bounded_above()) {
    conic_segment.domain().set_lower_bound(t, false);
    line_intersection_indices.first = polygon_boundary_index;
    line_intersection_indices.second = -1;
  } else if (conic_segment.domain().is_bounded_below()) {
    conic_segment.domain().set_upper_bound(t, false);
    line_intersection_indices.first = -1;
    line_intersection_indices.second = polygon_boundary_index;
  } else {
    SPDLOG_ERROR("Cone conic is a line, not the expected ray");
    return false;
  }

  return true;
}

bool
check_if_conic_intersects_cone_patch_domain(const Conic& conic,
                                            const ConvexPolygon& convex_polygon,
                                            size_t cone_corner_index)
{
  // Get implicit line equation L0 + Lx x + Ly y for the conic with 0 constant
  // term
  Eigen::Matrix<double, 3, 2> P_coeffs = conic.get_numerators();
  double Lx = P_coeffs(1, 1);
  double Ly = -P_coeffs(1, 0);
  double L0 = -(Lx * P_coeffs(0, 0) + Ly * P_coeffs(0, 1));

  // Get the two domain points that the conic does not pass through
  Eigen::Matrix<double, 3, 2> const& uv = convex_polygon.get_vertices();
  PlanarPoint domain_point_0 = uv.row((cone_corner_index + 1) % 3);
  PlanarPoint domain_point_1 = uv.row((cone_corner_index + 2) % 3);

  // Get sign of the two points relative to the line
  double domain_point_0_sign =
    L0 + Lx * domain_point_0[0] + Ly * domain_point_0[1];
  double domain_point_1_sign =
    L0 + Lx * domain_point_1[0] + Ly * domain_point_1[1];

  // Check for sign difference
  return ((domain_point_0_sign * domain_point_1_sign) < 0.0);
}