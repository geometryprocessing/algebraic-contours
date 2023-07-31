// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_ray_intersections_pencil_method.h"

#include "intersection_heuristics.h"
#include "polynomial_function.h"
#include "ray_intersection_heuristics.h"

int
solve_quadr(double q[3], double thres, double solution[2]);
int
solve_lin_quadr(double a_in[6],
                double b_in[6],
                double thres,
                double solution[2][2]);

int
solve_quadr(double q[3], double thres, double solution[2])
{
  double discr;
  int num_sol;
  if (thres <= fabs(q[0])) {
    discr = -0.4e1 * q[0] * q[2] + q[1] * q[1];
    if (thres * thres <= discr) {
      if (0.0e0 < q[1]) {
        solution[0] = 0.2e1 * q[2] / (-q[1] - sqrt(discr));
        solution[1] = (-q[1] - sqrt(discr)) / q[0] / 0.2e1;
      } else {
        solution[0] = (-q[1] + sqrt(discr)) / q[0] / 0.2e1;
        solution[1] = 0.2e1 * q[2] / (-q[1] + sqrt(discr));
      }
      num_sol = 2;
    } else if (0.0e0 <= discr) {
      solution[0] = -q[1] / q[0] / 0.2e1;
      num_sol = 1;
    } else
      num_sol = 0;
  } else if (thres <= fabs(q[1])) {
    solution[0] = -q[2] / q[1];
    num_sol = 1;
  } else
    num_sol = 0;
  return (num_sol);
}

// uu vv uv 1 u v
int
solve_lin_quadr(double a_in[6],
                double b_in[6],
                double thres,
                double solution[2][2])
{
  double quadr_coeff[3];
  double v_sol[2];
  double u_sol[2];
  int num_sol;
  int i;
  double a[6];
  double b[6];
  double ma;
  double mb;
  ma = a_in[0];
  mb = b_in[0];
  for (i = 2; i <= 6; i++) {
    ma = fmax(fabs(a_in[i - 1]), ma);
    mb = fmax(fabs(b_in[i - 1]), mb);
  }
  if (ma < thres || mb < thres)
    return (0);
  else
    for (i = 1; i <= 6; i++) {
      a[i - 1] = a_in[i - 1] / ma;
      b[i - 1] = b_in[i - 1] / mb;
    }
  quadr_coeff[0] = 0;
  quadr_coeff[1] = 0;
  quadr_coeff[2] = 0;
  if (fabs(b[5]) <= fabs(b[4]) && thres <= fabs(b[4])) {
    quadr_coeff[0] =
      a[0] * b[5] * b[5] - a[2] * b[5] * b[4] + a[1] * b[4] * b[4];
    quadr_coeff[1] = 0.2e1 * a[0] * b[3] * b[5] - a[2] * b[3] * b[4] -
                     a[4] * b[5] * b[4] + a[5] * b[4] * b[4];
    quadr_coeff[2] =
      b[3] * b[3] * a[0] - a[4] * b[3] * b[4] + a[3] * b[4] * b[4];
    v_sol[0] = 0;
    v_sol[1] = 0;
    num_sol = solve_quadr(quadr_coeff, thres, v_sol);
    if (num_sol == 2) {
      solution[0][0] = -(v_sol[0] * b[5] + b[3]) / b[4];
      solution[0][1] = v_sol[0];
      solution[1][0] = -(v_sol[1] * b[5] + b[3]) / b[4];
      solution[1][1] = v_sol[1];
    } else if (num_sol == 1) {
      solution[0][0] = -(v_sol[0] * b[5] + b[3]) / b[4];
      solution[0][1] = v_sol[0];
    }
  } else if (fabs(b[4]) < fabs(b[5]) && thres <= fabs(b[5])) {
    quadr_coeff[0] =
      a[0] * b[5] * b[5] - a[2] * b[5] * b[4] + a[1] * b[4] * b[4];
    quadr_coeff[1] = -a[2] * b[3] * b[5] + 0.2e1 * a[1] * b[3] * b[4] +
                     a[4] * b[5] * b[5] - a[5] * b[4] * b[5];
    quadr_coeff[2] =
      a[1] * b[3] * b[3] - a[5] * b[3] * b[5] + a[3] * b[5] * b[5];
    u_sol[0] = 0;
    u_sol[1] = 0;
    num_sol = solve_quadr(quadr_coeff, thres, u_sol);
    if (num_sol == 2) {
      solution[0][0] = u_sol[0];
      solution[0][1] = -(u_sol[0] * b[4] + b[3]) / b[5];
      solution[1][0] = u_sol[1];
      solution[1][1] = -(u_sol[1] * b[4] + b[3]) / b[5];
    } else if (num_sol == 1) {
      solution[0][0] = u_sol[0];
      solution[0][1] = -(u_sol[0] * b[4] + b[3]) / b[5];
    }
  } else
    num_sol = 0;
  return (num_sol);
}

bool
pencil_first_part(
  double coeff_F[6],
  double coeff_G[6],
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& intersection_points)
{
  num_intersections = 0;
  double coeff_threshold = 1e-10;

  bool F_linear = false;
  bool G_linear = false;
  if (abs(coeff_F[0]) < coeff_threshold && abs(coeff_F[1]) < coeff_threshold &&
      abs(coeff_F[3]) < coeff_threshold) {
    F_linear = true;
  }
  if (abs(coeff_G[0]) < coeff_threshold && abs(coeff_G[1]) < coeff_threshold &&
      abs(coeff_G[3]) < coeff_threshold) {
    G_linear = true;
  }

  // different cases

  if (F_linear && G_linear) {
    // case both linear

    // std::cout << "in 2 linear" << std::endl;
    double d = coeff_F[4] * coeff_G[5] - coeff_F[5] * coeff_G[4];
    if (abs(d) > coeff_threshold) {
      double u = -(coeff_F[2] * coeff_G[5] - coeff_F[5] * coeff_G[2]) / d;
      double v = (coeff_F[2] * coeff_G[4] - coeff_F[4] * coeff_G[2]) / d;

      if (0.0 <= u && 0.0 <= v && u + v <= 1.0) {
        intersection_points[num_intersections] = PlanarPoint(u, v);
        num_intersections++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else if (F_linear && !G_linear) {
    // F linear   G quadr

    double solution[2][2];
    double coeff_F_sol[6] = { coeff_F[0], coeff_F[1], coeff_F[3],
                              coeff_F[2], coeff_F[4], coeff_F[5] };
    double coeff_G_sol[6] = { coeff_G[0], coeff_G[1], coeff_G[3],
                              coeff_G[2], coeff_G[4], coeff_G[5] };
    int num_solution =
      solve_lin_quadr(coeff_G_sol, coeff_F_sol, coeff_threshold, solution);

    for (int i = 0; i < num_solution; i++) {
      if (0.0 <= solution[i][0] && 0.0 <= solution[i][1] &&
          solution[i][0] + solution[i][1] <= 1.0) {
        intersection_points[num_intersections] =
          PlanarPoint(solution[i][0], solution[i][1]);
        num_intersections++;
      }
    }
    if (num_intersections > 0) {
      return true;
    } else {
      return false;
    }
  } else if (!F_linear && G_linear) {
    // F quadr   G linear

    // std::cout << "in F quadr G linear" << std::endl;
    double solution[2][2];
    double coeff_F_sol[6] = { coeff_F[0], coeff_F[1], coeff_F[3],
                              coeff_F[2], coeff_F[4], coeff_F[5] };
    double coeff_G_sol[6] = { coeff_G[0], coeff_G[1], coeff_G[3],
                              coeff_G[2], coeff_G[4], coeff_G[5] };
    int num_solution =
      solve_lin_quadr(coeff_F_sol, coeff_G_sol, coeff_threshold, solution);
    for (int i = 0; i < num_solution; i++) {
      if (0.0 <= solution[i][0] && 0.0 <= solution[i][1] &&
          solution[i][0] + solution[i][1] <= 1.0) {
        intersection_points[num_intersections] =
          PlanarPoint(solution[i][0], solution[i][1]);
        num_intersections++;
      }
    }
    if (num_intersections > 0) {
      return true;
    } else {
      return false;
    }
  }

  // std::cout << "in F G quadr" << std::endl;

  // F and G both quadr

  // input F = a uu + b vv + c + d uv + e u + f v
  // input G = l uu + m vv + n + o uv + p u + q v

  // code convention
  // F = a uu + b vv + c + 2d uv + 2e u + 2f v
  // G = l uu + m vv + n + 2o uv + 2p u + 2q v

  // so the each last 3 inputs should be divided by 2

  int intersection_flag = false;

  double a, b, c, d, e, f, l, m, n, o, p, q;
  a = coeff_F[0];
  b = coeff_F[1];
  c = coeff_F[2];
  d = coeff_F[3] / 2.0;
  e = coeff_F[4] / 2.0;
  f = coeff_F[5] / 2.0;
  l = coeff_G[0];
  m = coeff_G[1];
  n = coeff_G[2];
  o = coeff_G[3] / 2.0;
  p = coeff_G[4] / 2.0;
  q = coeff_G[5] / 2.0;

  // cubic equation
  double a0, a1, a2, a3;
  a0 = (l * m * n + 2.0 * o * p * q) - (l * q * q + m * p * p + n * o * o);
  a1 = (a * m * n + l * b * n + l * m * c +
        2.0 * (d * p * q + o * e * q + o * p * f)) -
       (a * q * q + b * p * p + c * o * o +
        2.0 * (l * f * q + m * e * p + n * d * o));
  a2 = (a * b * n + a * m * c + l * b * c +
        2.0 * (o * e * f + d * e * q + d * p * f)) -
       (l * f * f + m * e * e + n * d * d +
        2.0 * (a * f * q + b * e * p + c * d * o));
  a3 = (a * b * c + 2.0 * d * e * f) - (a * f * f + b * e * e + c * d * d);

  //   std::cout << "cubic coeffs: a3 a2 a1 a0:" << std::endl;
  //   std::cout << a3 << " " << a2 << " " << a1 << " " << a0 << std::endl;

  int num_cubic_real_roots = 0;
  std::array<double, 3> cubic_real_roots;
  // check if a3 = 0;
  if (abs(a3) < coeff_threshold) {
    // quadr equation
    double q[3] = { a2, a1, a0 };
    double solution[2];
    int num_solution = solve_quadr(q, coeff_threshold, solution);
    for (int i = 0; i < num_solution; i++) {
      cubic_real_roots[num_cubic_real_roots] = solution[i];
      num_cubic_real_roots++;
    }
  } else {
    // solve cubic
    Eigen::Vector4d cubic_coeffs(a0, a1, a2, a3);
    Eigen::PolynomialSolver<double, 3> cubic_solver;
    cubic_solver.compute(cubic_coeffs);
    Eigen::PolynomialSolver<double, 3>::RootsType cubic_roots =
      cubic_solver.roots();

    // check real roots
    double imag_threshold = 1e-12;
    // std::cout << "cubic roots: " << std::endl;
    for (int i = 0; i < 3; i++) {
      //   std::cout << cubic_roots[i] << std::endl;
      if (abs(cubic_roots[i].imag()) < imag_threshold) {
        cubic_real_roots[num_cubic_real_roots] = cubic_roots[i].real();
        num_cubic_real_roots++;
      }
    }
  }

  if (num_cubic_real_roots == 0) {
    // no real root
    // std::cout << "no real root" << std::endl;
    return false;
  }

  double x;
  double determinant = std::numeric_limits<double>::infinity();
  for (int i = 0; i < num_cubic_real_roots; i++) {
    double A = a * cubic_real_roots[i] + l;
    double B = b * cubic_real_roots[i] + m;
    double D = d * cubic_real_roots[i] + o;
    if (determinant > (D * D - A * B)) {
      determinant = D * D - A * B;
      x = cubic_real_roots[i];
    }
  }

  //   std::cout << "x: " << x << std::endl;
  //   std::cout << "determinant: " << determinant << std::endl;

  double A = a * x + l;
  double B = b * x + m;
  double D = d * x + o;
  double C, E, F;

  if (determinant > 1e-10) {
    C = c * x + n;
    E = e * x + p;
    F = f * x + q;

    if (abs(A) < abs(B)) {
      A /= B;
      C /= B;
      D /= B;
      E /= B;
      F /= B;
      B = 1.0;
      double sqrtA = sqrt(D * D - A);
      double sqrtC = sqrt(F * F - C);
      double la1 = D + sqrtA;
      double la2 = D - sqrtA;
      double lc1 = F + sqrtC;
      double lc2 = F - sqrtC;

      if (abs((2.0 * E) - (la1 * lc1 + la2 * lc2)) <
          abs((2.0 * E) - (la1 * lc2 + la2 * lc1))) {
        double tmp = lc1;
        lc1 = lc2;
        lc2 = tmp;
      }

      // quadratic equation: c0 uu + c1 u + c2 = 0  c0 c1 c2 g h
      for (int i = 0; i < 2; i++) {
        double g = (i == 0) ? -la1 : -la2;
        double h = (i == 0) ? -lc1 : -lc2;
        double c0 = a + (2.0 * d + b * g) * g;
        double c1 = 2.0 * ((d + b * g) * h + e + f * g);
        double c2 = (b * h + 2.0 * f) * h + c;

        if ((0.0 < c0) && (((0.0 < c1) && (0.0 < c2)) ||
                           ((0.0 > c2) && (0.0 > c0 + c1 + c2)))) {
          continue;
        }
        if ((0.0 > c0) && (((0.0 > c1) && (0.0 > c2)) ||
                           ((0.0 < c2) && (0.0 < c0 + c1 + c2)))) {
          continue;
        }

        if (0.00000000001 > abs(c0)) {
          if (0.00000000001 < abs(c1)) {
            double u = -c2 / c1;
            double v = g * u + h;
            double w = 1.0 - (u + v);

            if (0.0 <= u && 0.0 <= v && 0.0 <= w) {
              intersection_points[num_intersections] = PlanarPoint(u, v);
              num_intersections++;
              intersection_flag = true;
            }
          }
        } else {

          double discriminant = (c1 * c1) - 4.0 * (c0 * c2);
          if (0.0 <= discriminant) {
            Eigen::PolynomialSolver<double, 2> quadratic_solver;
            Eigen::Vector3d quadratic_coeffs(c2, c1, c0);
            quadratic_solver.compute(quadratic_coeffs);
            Eigen::PolynomialSolver<double, 2>::RootsType quadratic_roots =
              quadratic_solver.roots();

            // with line i
            for (int i = 0; i < 2; i++) {
              double u = quadratic_roots[i].real();
              double v = g * u + h;
              double w = 1.0 - (u + v);

              if (0.0 <= u && 0.0 <= v && 0.0 <= w) {
                intersection_points[num_intersections] = PlanarPoint(u, v);
                num_intersections++;
                intersection_flag = true;
              }
            }
          }
        }
      }

    } else {
      B /= A;
      C /= A;
      D /= A;
      E /= A;
      F /= A;
      A = 1.0;
      double sqrtB = sqrt(D * D - B);
      double sqrtC = sqrt(E * E - C);
      double lb1 = D + sqrtB;
      double lb2 = D - sqrtB;
      double lc1 = E + sqrtC;
      double lc2 = E - sqrtC;

      if (abs((2.0 * F) - (lb1 * lc1 + lb2 * lc2)) <
          abs((2.0 * F) - (lb1 * lc2 + lb2 * lc1))) {
        double tmp = lc1;
        lc1 = lc2;
        lc2 = tmp;
      }

      // quadratic equation: c0 vv + c1 v + c2 = 0  c0 c1 c2 g h
      for (int i = 0; i < 2; i++) {
        double g = (i == 0) ? -lb1 : -lb2;
        double h = (i == 0) ? -lc1 : -lc2;
        double c0 = b + (2.0 * d + a * g) * g;
        double c1 = 2.0 * ((d + a * g) * h + f + e * g);
        double c2 = (a * h + 2.0 * e) * h + c;

        if ((0.0 < c0) && (((0.0 < c1) && (0.0 < c2)) ||
                           ((0.0 > c2) && (0.0 > c0 + c1 + c2)))) {
          continue;
        }
        if ((0.0 > c0) && (((0.0 > c1) && (0.0 > c2)) ||
                           ((0.0 < c2) && (0.0 < c0 + c1 + c2)))) {
          continue;
        }

        if (0.00000000001 > abs(c0)) {
          if (0.00000000001 < abs(c1)) {
            double v = -c2 / c1;
            double u = g * v + h;
            double w = 1.0 - (u + v);
            if (0.0 <= u && 0.0 <= v && 0.0 <= w) {
              intersection_points[num_intersections] = PlanarPoint(u, v);
              num_intersections++;
              intersection_flag = true;
            }
          }
        } else {

          double discriminant = (c1 * c1) - 4.0 * (c0 * c2);
          if (0.0 <= discriminant) {
            Eigen::PolynomialSolver<double, 2> quadratic_solver;
            Eigen::Vector3d quadratic_coeffs(c2, c1, c0);
            quadratic_solver.compute(quadratic_coeffs);
            Eigen::PolynomialSolver<double, 2>::RootsType quadratic_roots =
              quadratic_solver.roots();

            // with line i
            for (int i = 0; i < 2; i++) {
              double v = quadratic_roots[i].real();
              double u = g * v + h;
              double w = 1.0 - (u + v);

              if (0.0 <= u && 0.0 <= v && 0.0 <= w) {
                intersection_points[num_intersections] = PlanarPoint(u, v);
                num_intersections++;
                intersection_flag = true;
              }
            }
          }
        }
      }
    }

  } else {
    // ellipsoid with zero area
    return false;
  }

  return intersection_flag;
}

// 1 u v uv uu vv call
void
solve_quadratic_quadratic_equation_pencil_method(
  double a[6],
  double b[6],
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& intersection_points)
{

  // divide by max coeff to get better precision

  double max = abs(a[0]);
  for (int i = 0; i < 6; i++) {
    if (max < abs(a[i]))
      max = abs(a[i]);
    if (max < abs(b[i]))
      max = abs(b[i]);
  }

  double F[6] = { a[4] / max, a[5] / max, a[0] / max,
                  a[3] / max, a[1] / max, a[2] / max };
  double G[6] = { b[4] / max, b[5] / max, b[0] / max,
                  b[3] / max, b[1] / max, b[2] / max };
  pencil_first_part(F, G, num_intersections, intersection_points);
}

void
compute_spline_surface_patch_ray_intersections_pencil_method(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix2x3r& ray_mapping_coeffs,
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& surface_intersections,
  std::array<double, MAX_PATCH_RAY_INTERSECTIONS>& ray_intersections,
  long long& ray_intersections_call,
  long long& ray_bounding_box_call)
{
  num_intersections = 0;
  const ConvexPolygon& domain = spline_surface_patch.get_domain();
  spdlog::trace("Domain: {}", domain.get_vertices());

  SpatialVector ray_origin = ray_mapping_coeffs.row(0);
  PlanarPoint ray_plane_point = ray_origin.head(2);
  spdlog::trace(
    "Computing intersections for ray origin {} with planar projection {}",
    ray_origin,
    ray_plane_point);

  // Check if the bounding box of the projected patch boundaries contains the
  // ray
  const SpatialVector& min_point =
    spline_surface_patch.get_bounding_box_min_point();
  const SpatialVector& max_point =
    spline_surface_patch.get_bounding_box_max_point();
  PlanarPoint lower_left_point(min_point[0], min_point[1]);
  PlanarPoint upper_right_point(max_point[0], max_point[1]);
  ray_bounding_box_call++;
  if (!is_in_bounding_box(
        ray_plane_point, lower_left_point, upper_right_point)) {
    SPDLOG_TRACE("Skipping intersection test for patch {} and ray {} with "
                 "bounding box ({}, {})",
                 spline_surface_patch,
                 formatted_polynomial(ray_mapping_coeffs),
                 lower_left_point.transpose(),
                 upper_right_point.transpose());
    return;
  }

  ray_intersections_call++;

  SPDLOG_TRACE("Computing intersections for patch {} and ray {}",
               spline_surface_patch,
               formatted_polynomial(ray_mapping_coeffs));

  // Normalize the spline surface patch to have domain triangle u + v <= 1 in
  // [0, 1]^2
  const Matrix6x3r& normalized_surface_mapping_coeffs =
    spline_surface_patch.get_normalized_surface_mapping();

  SPDLOG_TRACE(
    "Coefficients for the surface mapping with normalized domain: {}",
    formatted_bivariate_quadratic_mapping(normalized_surface_mapping_coeffs));

  // 1 u v uv u2 v2
  // to
  //  a uu + b vv + c + d uv + e u + f v

  double coeff_F[6] = { normalized_surface_mapping_coeffs(4, 0),
                        normalized_surface_mapping_coeffs(5, 0),
                        normalized_surface_mapping_coeffs(0, 0) - ray_origin[0],
                        normalized_surface_mapping_coeffs(3, 0),
                        normalized_surface_mapping_coeffs(1, 0),
                        normalized_surface_mapping_coeffs(2, 0) };

  double coeff_G[6] = { normalized_surface_mapping_coeffs(4, 1),
                        normalized_surface_mapping_coeffs(5, 1),
                        normalized_surface_mapping_coeffs(0, 1) - ray_origin[1],
                        normalized_surface_mapping_coeffs(3, 1),
                        normalized_surface_mapping_coeffs(1, 1),
                        normalized_surface_mapping_coeffs(2, 1) };

  // renormalize coeffs
  double max_c = abs(coeff_F[0]);
  for (int i = 0; i < 6; i++) {
    if (max_c < abs(coeff_F[i]))
      max_c = abs(coeff_F[i]);
    if (max_c < abs(coeff_G[i]))
      max_c = abs(coeff_G[i]);
  }

  for (int i = 0; i < 6; i++) {
    coeff_F[i] /= max_c;
    coeff_G[i] /= max_c;
  }

  // Get all intersections of the quadratic normalized surface and the ray
  int num_intersections_all;
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>
    normalized_surface_intersections_all;
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>
    normalized_surface_intersections;
  std::array<double, MAX_PATCH_RAY_INTERSECTIONS> normalized_ray_intersections;
  pencil_first_part(coeff_F,
                    coeff_G,
                    num_intersections_all,
                    normalized_surface_intersections_all);
  if (num_intersections_all > MAX_PATCH_RAY_INTERSECTIONS) {
    spdlog::error(
      "More than the maximum possible number of patch ray intersections found");
  }
  assert(num_intersections_all <= MAX_PATCH_RAY_INTERSECTIONS);
  SPDLOG_TRACE("{} intersections found before pruning", num_intersections_all);

  // Get intersections that are in the domain
  int num_intersections_normalized = 0;
  for (int i = 0; i < num_intersections_all; i++) {
    double t = (normalized_surface_mapping_coeffs(0, 2) +
                normalized_surface_mapping_coeffs(1, 2) *
                  normalized_surface_intersections_all[i][0] +
                normalized_surface_mapping_coeffs(2, 2) *
                  normalized_surface_intersections_all[i][1] +
                normalized_surface_mapping_coeffs(3, 2) *
                  normalized_surface_intersections_all[i][0] *
                  normalized_surface_intersections_all[i][1] +
                normalized_surface_mapping_coeffs(4, 2) *
                  normalized_surface_intersections_all[i][0] *
                  normalized_surface_intersections_all[i][0] +
                normalized_surface_mapping_coeffs(5, 2) *
                  normalized_surface_intersections_all[i][1] *
                  normalized_surface_intersections_all[i][1] -
                ray_origin[2]) /
               ray_mapping_coeffs(1, 2);
    if (t > 0 and t <= 1) {
      normalized_ray_intersections[num_intersections_normalized] = t;
      normalized_surface_intersections[num_intersections_normalized] =
        normalized_surface_intersections_all[i];
      num_intersections_normalized++;
    }
  }

  // Invert the normalization and prune intersections outside of the triangle
  // domain
  for (int i = 0; i < num_intersections_normalized; ++i) {
    PlanarPoint normalized_domain_point = normalized_surface_intersections[i];
    PlanarPoint surface_intersection =
      spline_surface_patch.denormalize_domain_point(normalized_domain_point);
    if (domain.contains(surface_intersection)) {
      surface_intersections[num_intersections] = surface_intersection;
      ray_intersections[num_intersections] = normalized_ray_intersections[i];
      num_intersections++;
    }
    spdlog::trace("Normalized domain point: {}", normalized_domain_point);
    spdlog::trace("Domain point: {}", surface_intersection);
  }

  // Log intersections
  if (num_intersections >= 0) {
    SPDLOG_TRACE("Intersections for patch {} and ray {}",
                 spline_surface_patch,
                 formatted_polynomial(ray_mapping_coeffs));
    SPDLOG_TRACE(
      "Coefficients for the surface mapping with normalized domain: {}",
      formatted_bivariate_quadratic_mapping(normalized_surface_mapping_coeffs));
  }
}