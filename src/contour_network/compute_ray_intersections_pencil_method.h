// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

#include "quadratic_spline_surface.h"

bool
pencil_first_part(
  double coeff_F[6],
  double coeff_G[6],
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& intersection_points);

void
compute_spline_surface_patch_ray_intersections_pencil_method(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix2x3r& ray_mapping_coeffs,
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& surface_intersections,
  std::array<double, MAX_PATCH_RAY_INTERSECTIONS>& ray_intersections,
  long long& ray_intersections_call,
  long long& ray_bounding_box_call);

void
solve_quadratic_quadratic_equation_pencil_method(
  double a[6],
  double b[6],
  int& num_intersections,
  std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>& intersection_points);