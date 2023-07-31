// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "convert_conic.h"

// Given a symmetric 2x2 matrix A, compute the eigenvalues and a rotation matrix
// U of eigenvectors so that A = U * diag(eigenvalues) * U^T
//
// param[in] A: symmetric matrix to decompose
// param[out] eigenvalues: length 2 vector of eigenvalues with the largest first
// param[out] rotation: rotation matrix such the the columns are the eigenvalues
void
compute_symmetric_matrix_eigen_decomposition(const Matrix2x2r& A,
                                             std::array<double, 2>& eigenvalues,
                                             Matrix2x2r& rotation)
{
  assert(float_equal(A(0, 1), A(1, 0)));

  // Compute eigenvalues sigma_1, sigma_2
  // FIXME Remove pow everywhere trace^2 - 4 det
  double discriminant = power(A(0, 0) - A(1, 1), 2) + 4 * A(0, 1) * A(1, 0);
  double sigma_1 = 0.5 * (A(0, 0) + A(1, 1) + std::sqrt(discriminant));
  double sigma_2 = 0.5 * (A(0, 0) + A(1, 1) - std::sqrt(discriminant));
  assert(sigma_1 >= sigma_2);
  eigenvalues = { sigma_1, sigma_2 };

  // Compute rotation matrix U such that A = U diag(sigma_1, sigma_2) U^T
  // First row of rotation matrix is the first eigenvector
  // FIXME Check if discriminant zero instead
  // Otherwise, if a00 < a11, use current formula, if not use other (with
  // A(1,1))
  PlanarPoint eigenvector_1;
  if (!float_equal(A(0, 1), 0.0)) {
    eigenvector_1 << A(0, 1), sigma_1 - A(0, 0);
    assert(!float_equal(eigenvector_1.norm(), 0.0));
    eigenvector_1 /= eigenvector_1.norm();
  }
  // This can be removed
  else if (A(0, 0) < A(1, 1)) {
    eigenvector_1 << 0, 1;
  } else {
    eigenvector_1 << 1, 0;
  }
  rotation.row(0) = eigenvector_1;

  // Second column of rotation matrix is the first eigenvector rotated 90
  // degrees
  PlanarPoint eigenvector_2;
  eigenvector_2 << -eigenvector_1(1), eigenvector_1(0);
  rotation.row(1) = eigenvector_2;
  assert(float_equal(rotation.determinant(), 1.0));
  assert(matrix_equal(
    A,
    rotation.transpose() *
      Eigen::DiagonalMatrix<double, 2>(eigenvalues[0], eigenvalues[1]) *
      rotation));
}

// Given a conic C represented by coefficients conic_coeffs corresponding to
// 1, u, v, uv, u^2, v^2, express the quadratic equation in the form
// 0.5 r^T A r + b^T r + c
//
// param[out] conic_coeffs: coefficients for the conic C in terms of r
// param[out] A: quadratic terms symmetric matrix
// param[out] b: linear terms vector
// param[out] c: constant term
void
convert_conic_to_matrix_form(const Vector6r conic_coeffs,
                             Matrix2x2r& A,
                             Eigen::Matrix<double, 2, 1>& b,
                             double& c)
{
  // Compute A
  A(0, 0) = 2.0 * conic_coeffs(4);
  A(0, 1) = conic_coeffs(3);
  A(1, 0) = conic_coeffs(3);
  A(1, 1) = 2.0 * conic_coeffs(5);

  // Compute b
  b(0) = conic_coeffs(1);
  b(1) = conic_coeffs(2);

  // Compute c
  c = conic_coeffs(0);
}

// Given coefficients conic_coeffs for a conic C, rotate and translate the conic
// so that it is centered at the origin and axis aligned. A point r in the
// original conic is mapped to U(r - r_0) in the standard form conic.
//
// param[out] conic_coeffs: original coefficients for the conic C
// param[out] conic_standard_form: coefficients for the standard form conic C
// param[out] rotation: rotation U to convert to standard form
// param[out] translation: translation r_0 to convert to standard form.
void
convert_conic_to_standard_form(const Vector6r& conic_coeffs,
                               Vector6r& conic_standard_form,
                               Matrix2x2r& rotation,
                               PlanarPoint& translation)
{
  translation.setZero();
  conic_standard_form.setZero();

  // Put contour coefficients in quadratic form (1/2)r^T A r + b^T r + c
  Matrix2x2r A;
  Eigen::Matrix<double, 2, 1> b;
  double c;
  convert_conic_to_matrix_form(conic_coeffs, A, b, c);
  // FIXME Make sure XT A X is identity

  // Get rotation and singular values
  std::array<double, 2> singular_values;
  compute_symmetric_matrix_eigen_decomposition(A, singular_values, rotation);
  double det = singular_values[0] * singular_values[1];

  // Singular conic (parabola, line, etc.)
  if ((float_equal(singular_values[0], 0.0)) ||
      (float_equal(singular_values[1], 0.0))) {
    // FIXME Add warning so we know this case is occurring
    spdlog::warn("Singular conic with det {} and singular values {}, {}",
                 det,
                 singular_values[0],
                 singular_values[1]);

    // Ensure the y coordinate is singular
    if (float_equal(singular_values[0], 0.0)) {
      std::swap(singular_values[0], singular_values[1]);
      std::swap(rotation(0, 0), rotation(1, 0));
      std::swap(rotation(0, 1), rotation(1, 1));

      // Ensure the rotation is still orientation preserving
      rotation.row(1) *= -1;
      assert(float_equal(rotation.determinant(), 1.0));
    }

    // The conic is a line or plane if A = 0
    if (float_equal(singular_values[0], 0.0)) {
      // Normalize the equation
      double normalization_factor = b.norm();
      b /= normalization_factor;
      c /= normalization_factor;

      // Translate b to c
      translation = -c * b;

      // Rotate e1 to the linear term b
      rotation.row(0) = b;
      rotation(1, 0) = -b(1);
      rotation(1, 1) = b(0);
      assert(float_equal(rotation.determinant(), 1.0));
      // assert(matrix_equal(rotation * rotation.transpose(), Matrix2x2r(1.0,
      // 0.0, 0.0, 1.0)));

      // The conic equation is just x = 0
      conic_standard_form(1) = 1.0;
    } else {
      conic_standard_form(4) = 0.5 * singular_values[0];
      Eigen::Matrix<double, 2, 1> Ub = rotation * b;

      // Quadratic in a single variable
      if (float_equal(Ub(1), 0.0)) {
        conic_standard_form(0) = c;
        conic_standard_form(1) = Ub(0);
      }
      // Translate a parabola so that its vertex is at the origin
      else {
        // translation(0) = Ub(0) / std::sqrt(0.5 * abs(singular_values(0)));
        // translation(1) = c / Ub(1);
        // translation = rotation * translation;
        // conic_standard_form(2) = 1.0;
        spdlog::debug("Parabola is {} u^2 + {} u + {} v + {}",
                      0.5 * singular_values[0],
                      Ub(0),
                      Ub(1),
                      c);
        translation(0) = -Ub(0) / singular_values[0];
        translation(1) =
          -(c - 0.5 * singular_values[0] * Ub(0) * Ub(0)) / (Ub(1));
        spdlog::debug("Using translation {} for parabola", translation);
        translation = translation * rotation;
        spdlog::debug("Rotating to translation {}", translation);
        conic_standard_form(2) = Ub(1);
      }

      // Ensure nonzero quadratic term is positive
      if (conic_standard_form(4) < 0.0) {
        conic_standard_form *= -1.0;
      }
    }

  }
  // Nonsingular conics (ellipse, hyperbola, etc.)
  else {
    // Get translation r_0 = -A^{-1} b
    translation = -A.inverse() * b;

    // Compute standard form coefficients for the conic
    conic_standard_form[0] =
      c - 0.5 * (translation * A) * translation.transpose();
    conic_standard_form[4] = 0.5 * singular_values[0];
    conic_standard_form[5] = 0.5 * singular_values[1];
  }

  SPDLOG_TRACE("Standard form: {}",
               formatted_bivariate_quadratic_mapping(conic_standard_form));
  spdlog::trace("Rotation:\n{}", rotation);
  spdlog::trace("Translation:\n{}", translation);
}

// Given coefficients conic_coeffs for a conic C, rotate and translate the conic
// so that it is centered at the origin and axis aligned.
//
// param[out] conic_coeffs: original coefficients for the conic C
// return; coefficients for the standard form conic C
Vector6r
convert_conic_to_standard_form(const Vector6r& conic_coeffs)
{
  Vector6r conic_standard_form;

  Matrix2x2r rotation;
  PlanarPoint translation;
  convert_conic_to_standard_form(
    conic_coeffs, conic_standard_form, rotation, translation);

  return conic_standard_form;
}
