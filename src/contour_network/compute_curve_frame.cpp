// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_curve_frame.h"

// Compute tangent of a curve on a quadratic surface
void
compute_quadratic_surface_curve_tangent(
  const Matrix6x3r& surface_mapping_coeffs,
  const Conic& domain_curve,
  RationalFunction<8, 3>& surface_curve_tangent)
{
  // Lift the domain curve to a surface curve
  spdlog::trace("Computing quadratic surface curve tangent");
  spdlog::trace("Domain curve: {}", domain_curve);
  RationalFunction<4, 3> surface_curve;
  domain_curve.pullback_quadratic_function<3>(surface_mapping_coeffs,
                                              surface_curve);
  spdlog::trace("Surface curve: {}", surface_curve);

  // Compute the tangent of the surface curve directly
  surface_curve.compute_derivative(surface_curve_tangent);
}

// Compute surface normal along a curve on a quadratic surface
void
compute_quadratic_surface_curve_normal(
  const Matrix6x3r& normal_mapping_coeffs,
  const Conic& domain_curve,
  RationalFunction<4, 3>& surface_curve_normal)
{
  domain_curve.pullback_quadratic_function<3>(normal_mapping_coeffs,
                                              surface_curve_normal);
}

// General method to compute the surface curve tangent normal from the tangent
// and normal
void
compute_surface_curve_tangent_normal(
  const RationalFunction<8, 3>& surface_curve_tangent,
  const RationalFunction<4, 3>& surface_curve_normal,
  RationalFunction<12, 3>& surface_curve_tangent_normal)
{
  // Compute the cross product of the (vector valued) numerators
  Eigen::Matrix<double, 9, 3> tangent_numerators =
    surface_curve_tangent.get_numerators();
  Eigen::Matrix<double, 5, 3> normal_numerators =
    surface_curve_normal.get_numerators();
  Eigen::Matrix<double, 13, 3> tangent_normal_numerators;
  compute_polynomial_mapping_cross_product<8, 4>(
    tangent_numerators, normal_numerators, tangent_normal_numerators);

  // Compute the scalar product of the denominators
  Eigen::Matrix<double, 9, 1> tangent_denominator =
    surface_curve_tangent.get_denominator();
  Eigen::Matrix<double, 5, 1> normal_denominator =
    surface_curve_normal.get_denominator();
  Eigen::Matrix<double, 13, 1> tangent_normal_denominator;
  compute_polynomial_mapping_product<8, 4, 1>(
    tangent_denominator, normal_denominator, tangent_normal_denominator);

  // Build rational function from the numerators, denominators, and the common
  // domain
  surface_curve_tangent_normal =
    RationalFunction<12, 3>(tangent_normal_numerators,
                            tangent_normal_denominator,
                            surface_curve_tangent.domain());
}

void
compute_quadratic_surface_curve_frame(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_mapping_coeffs,
  const Conic& domain_curve_segment,
  RationalFunction<8, 3>& surface_curve_tangent,
  RationalFunction<4, 3>& surface_curve_normal,
  RationalFunction<12, 3>& surface_curve_tangent_normal)
{
  compute_quadratic_surface_curve_tangent(
    surface_mapping_coeffs, domain_curve_segment, surface_curve_tangent);
  compute_quadratic_surface_curve_normal(
    normal_mapping_coeffs, domain_curve_segment, surface_curve_normal);
  compute_surface_curve_tangent_normal(
    surface_curve_tangent, surface_curve_normal, surface_curve_tangent_normal);
}

void
compute_spline_surface_patch_curve_frame(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Conic& domain_curve_segment,
  RationalFunction<8, 3>& surface_curve_tangent,
  RationalFunction<4, 3>& surface_curve_normal,
  RationalFunction<12, 3>& surface_curve_tangent_normal)
{
  spdlog::trace("Computing spline surface patch curve frame");
  // Get surface and normal mappings for the given patch
  Matrix6x3r surface_mapping_coeffs =
    spline_surface_patch.get_surface_mapping();
  Matrix6x3r normal_mapping_coeffs = spline_surface_patch.get_normal_mapping();

  // Compute frame for quadratic surface patch
  compute_quadratic_surface_curve_frame(surface_mapping_coeffs,
                                        normal_mapping_coeffs,
                                        domain_curve_segment,
                                        surface_curve_tangent,
                                        surface_curve_normal,
                                        surface_curve_tangent_normal);

  spdlog::trace("Surface curve tangent: {}", surface_curve_tangent);
  spdlog::trace("Surface curve normal: {}", surface_curve_normal);
  spdlog::trace("Surface curve tangent normal: {}",
                surface_curve_tangent_normal);
}
