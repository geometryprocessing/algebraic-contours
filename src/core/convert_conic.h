// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

// Given a symmetric 2x2 matrix A, compute the eigenvalues and a rotation matrix
// U of eigenvectors so that A = U * diag(eigenvalues) * U^T
//
// param[in] A: symmetric matrix to decompose
// param[out] eigenvalues: length 2 vector of eigenvalues with the largest first
// param[out] rotation: rotation matrix such the the columns are the eigenvalues
void
compute_symmetric_matrix_eigen_decomposition(const Matrix2x2r& A,
                                             std::array<double, 2>& eigenvalues,
                                             Matrix2x2r& rotation);

/// Given coefficients conic_coeffs for a conic C, rotate and translate the
/// conic so that it is centered at the origin and axis aligned. A point r in
/// the original conic is mapped to U(r - r_0) in the standard form conic.
///
/// @param[in] conic_coeffs: original coefficients for the conic C
/// @param[out] conic_standard_form: coefficients for the standard form conic C
/// @param[out] rotation: rotation U to convert to standard form
/// @param[out] translation: translation r_0 to convert to standard form.
void
convert_conic_to_standard_form(const Vector6r& conic_coeffs,
                               Vector6r& conic_standard_form,
                               Matrix2x2r& rotation,
                               PlanarPoint& translation);

/// Given coefficients conic_coeffs for a conic C, rotate and translate the
/// conic so that it is centered at the origin and axis aligned.
///
/// @param[in] conic_coeffs: original coefficients for the conic C
/// @return; coefficients for the standard form conic C
Vector6r
convert_conic_to_standard_form(const Vector6r& conic_coeffs);
