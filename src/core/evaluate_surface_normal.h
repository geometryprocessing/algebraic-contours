// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// \file evaluate_surface.h
///
/// Methods to evaluate normals to a quadratic surface with Zwart-Powell basis
/// coefficients

/// FIXME Rename this file

/// Compute the quadratic coefficients of the normal vector to a quadratic
/// surface.
///
/// @param[in] surface_mapping_coeffs: coefficients for the quadratic surface
///
/// @return Coefficients for the quadratic polynomial defining the normal vector
/// on the surface
void
generate_quadratic_surface_normal_coeffs(
  const Matrix6x3r& surface_mapping_coeffs,
  Matrix6x3r& normal_mapping_coeffs);
