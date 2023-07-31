// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "evaluate_general_contour_frame.h"

#include "compute_contours.h"
#include "evaluate_surface_normal.h"
#include "polynomial_function.h"

// Compute the coefficients for the polynomial gradient function of the general
// contour function tau * N.
void
compute_general_contour_function_gradient(const Matrix6x3r& normal_coeffs,
                                          const Matrix3x3r& frame,
                                          Matrix3x2r& gradient_coeffs)
{
  // Compute the coefficients of the contour function
  Vector6r contour_coeffs;
  compute_contour_equation(normal_coeffs, frame, contour_coeffs);

  // Compute the u, v derivatives
  Eigen::Matrix<double, 3, 6> D_u = u_derivative_matrix();
  Eigen::Matrix<double, 3, 6> D_v = v_derivative_matrix();
  gradient_coeffs.col(0) = D_u * contour_coeffs;
  gradient_coeffs.col(1) = D_v * contour_coeffs;
}

// Compute the coefficients of the polynomial mapping defining the tangents to
// the curves tau * N = c for any constant c on a quadratic surface.
void
compute_quadratic_surface_general_contour_tangent(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_coeffs,
  const Matrix3x3r& frame,
  Matrix6x3r& tangent_coeffs)
{
  // Get surface derivative functions
  Eigen::Matrix<double, 3, 6> D_u = u_derivative_matrix();
  Eigen::Matrix<double, 3, 6> D_v = v_derivative_matrix();
  Matrix3x3r surface_u_derivative = D_u * surface_mapping_coeffs;
  Matrix3x3r surface_v_derivative = D_v * surface_mapping_coeffs;

  // Get contour gradient
  Matrix3x2r gradient_coeffs;
  compute_general_contour_function_gradient(
    normal_coeffs, frame, gradient_coeffs);

  // Compute tangent coefficients
  tangent_coeffs.col(0) =
    -compute_linear_product(surface_u_derivative.col(0),
                            gradient_coeffs.col(1)) +
    compute_linear_product(surface_v_derivative.col(0), gradient_coeffs.col(0));
  tangent_coeffs.col(1) =
    -compute_linear_product(surface_u_derivative.col(1),
                            gradient_coeffs.col(1)) +
    compute_linear_product(surface_v_derivative.col(1), gradient_coeffs.col(0));
  tangent_coeffs.col(2) =
    -compute_linear_product(surface_u_derivative.col(2),
                            gradient_coeffs.col(1)) +
    compute_linear_product(surface_v_derivative.col(2), gradient_coeffs.col(0));
}

// Compute the coefficients of the polynomial mapping defining the tangents to
// the curves tau * N = c for any constant c on quadratic spline surface patch.
void
compute_spline_surface_patch_general_contour_tangent(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix3x3r& frame,
  Matrix6x3r& tangent_coeffs)
{
  // Get patch surface and normal coefficients
  Matrix6x3r surface_mapping_coeffs =
    spline_surface_patch.get_surface_mapping();
  Matrix6x3r normal_coeffs = spline_surface_patch.get_normal_mapping();

  // Compute quadratic patch tangents
  compute_quadratic_surface_general_contour_tangent(
    surface_mapping_coeffs, normal_coeffs, frame, tangent_coeffs);
}

//  Evaluate the tangent to the curve tau * N = c for some implicit constant c
//  on
// a quadratic spline surface patch at some local coordinates.
SpatialVector
evaluate_spline_surface_patch_general_contour_tangent(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix3x3r& frame,
  const PlanarPoint& domain_point)
{
  // Compute the tangent polynomial mapping
  Matrix6x3r tangent_coeffs;
  compute_spline_surface_patch_general_contour_tangent(
    spline_surface_patch, frame, tangent_coeffs);

  // Evaluate the tangent polynomial mapping at the given point
  Eigen::Matrix<double, 1, 6> w = generate_quadratic_monomials(domain_point);
  return w * tangent_coeffs;
}

/// Evaluate the tangent to the curve tau * N = c for some implicit
/// constant c on a quadratic spline surface at some global coordinates.
SpatialVector
evaluate_spline_surface_general_contour_tangent(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const PlanarPoint& domain_point,
  const QuadraticSplineSurface::PatchIndex& patch_index)
{
  // Evaluate the tangent for the surface patch
  return evaluate_spline_surface_patch_general_contour_tangent(
    spline_surface.get_patch(patch_index), frame, domain_point);
}

void
evaluate_spline_surface_general_contour_frame(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const PlanarPoint& domain_point,
  const QuadraticSplineSurface::PatchIndex& patch_index,
  SpatialVector& tangent,
  SpatialVector& normal,
  SpatialVector& tangent_normal)
{
  // Compute the tangent
  tangent = evaluate_spline_surface_general_contour_tangent(
    spline_surface, frame, domain_point, patch_index);

  // Compute the normal
  spline_surface.evaluate_patch_normal(patch_index, domain_point, normal);

  // Compute the tangent normal
  tangent_normal = cross_product<double>(normal, tangent);
}