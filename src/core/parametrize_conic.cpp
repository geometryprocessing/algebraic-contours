// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "parametrize_conic.h"

#include "bivariate_quadratic_function.h"
#include "convert_conic.h"

// Identify the type of conic C represented by C_standard form, which is assumed
// to be in standard form.
ConicType
identify_standard_form_conic(const Vector6r& conic_standard_form)
{
  // Extract standard form coefficients
  double c = conic_standard_form(0);
  double b_1 = conic_standard_form(1);
  double b_2 = conic_standard_form(2);
  double sigma_1 = conic_standard_form(4);
  double sigma_2 = conic_standard_form(5);
  double det = sigma_1 * sigma_2;

  // No quadratic terms: Conic is an empty set, line, or plane equation
  if (float_equal(sigma_1, 0.0) && float_equal(sigma_2, 0.0)) {
    if (!float_equal(b_1, 0.0) || !float_equal(b_2, 0.0)) {
      return ConicType::line;
    }
    if (!float_equal(c, 0.0)) {
      return ConicType::empty;
    } else {
      return ConicType::plane;
    }
  }

  // Singular conic: sigma_1 z_1^2 + b_1 z_1 + b_2 z_2 + c = 0
  if (float_equal(sigma_1, 0.0) || float_equal(sigma_2, 0.0)) {
    if (float_equal(b_2, 0.0)) {
      return ConicType::parallel_lines;
    } else {
      return ConicType::parabola;
    }
  }

  // Nonsingular conic with no constant term: sigma_1 z_1^2 + sigma_2 z_2^2 = 0
  if (float_equal(c, 0.0)) {
    if (det > 0) {
      return ConicType::point;
    } else {
      return ConicType::intersecting_lines;
    }
  }

  // Nonsingular conic with constant term: sigma_1 z_1^2 + sigma_2 z_2^2 + c = 0
  if ((sigma_1 > 0) && (sigma_2 > 0) && (c > 0)) {
    return ConicType::empty;
  } else if ((sigma_1 < 0) && (sigma_2 < 0) && (c < 0)) {
    return ConicType::empty;
  } else if ((sigma_1 > 0) && (sigma_2 > 0) && (c < 0)) {
    return ConicType::ellipse;
  } else if ((sigma_1 < 0) && (sigma_2 < 0) && (c > 0)) {
    return ConicType::ellipse;
  } else {
    return ConicType::hyperbola;
  }
}

// Identify the type of conic C represented by conic_coeffs.
//
// param[in] conic_coeffs: quadratic coefficients for the conic
// return: type identifier of the conic
ConicType
identify_conic(const Vector6r& conic_coeffs)
{
  Vector6r conic_standard_form = convert_conic_to_standard_form(conic_coeffs);
  return identify_standard_form_conic(conic_standard_form);
}

void
parametrize_ellipse(const Vector6r& conic_standard_form,
                    std::vector<Conic>& conics)
{
  assert(identify_conic(conic_standard_form) == ConicType::ellipse);

  // Get axes lengths k1, k2
  double c = conic_standard_form(0);
  double sigma_1 = conic_standard_form(4);
  double sigma_2 = conic_standard_form(5);
  double k1 = std::sqrt(abs(c / sigma_1));
  double k2 = std::sqrt(abs(c / sigma_2));

  // Determine sign of c used for orientation
  double sign = sgn(c);

  // Parametrize ellipse as (1/(1 + t^2)) * [2 k1 t, k2 (1 - t^2)]
  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  Eigen::Matrix<double, 3, 1> elliptic_basis_1(3);
  Eigen::Matrix<double, 3, 1> elliptic_basis_2(3);
  elliptic_basis_1 << 1, 0, -1;
  elliptic_basis_2 << 0, 2, 0;
  Q_coeffs << 1.0, 0.0, 1.0;
  if (sign < 0) {
    P_coeffs.col(0) = k1 * elliptic_basis_1;
    P_coeffs.col(1) = k2 * elliptic_basis_2;
  } else {
    P_coeffs.col(0) = -k1 * elliptic_basis_1;
    P_coeffs.col(1) = k2 * elliptic_basis_2;
  }

  // Parametrize upper half over finite interval [-1,1]
  Conic semi_ellipse(P_coeffs, Q_coeffs, ConicType::ellipse);
  semi_ellipse.domain().set_lower_bound(-1.0, false);
  semi_ellipse.domain().set_upper_bound(1.0, false);
  conics.push_back(semi_ellipse);

  // Parametrize lower half (also with closed endpoints)
  semi_ellipse.set_numerators(-P_coeffs);
  semi_ellipse.domain().set_lower_bound(-1.0, false);
  semi_ellipse.domain().set_upper_bound(1.0, false);
  conics.push_back(semi_ellipse);
}

void
parametrize_hyperbola(const Vector6r& conic_standard_form,
                      std::vector<Conic>& conics)
{
  // assert( identify_conic(conic_coeffs) == ConicType::hyperbola );

  // Get axes lengths
  double c = conic_standard_form(0);
  double sigma_1 = conic_standard_form(4);
  double sigma_2 = conic_standard_form(5);
  double k1 = std::sqrt(abs(c / sigma_1));
  double k2 = std::sqrt(abs(c / sigma_2));

  // Parametrize ellipse as (1/(1 - t^2)) * [k1 (1 + t^2), 2 k2 t]
  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  double sign = sgn(c);
  // FIXME Clean up
  Eigen::Matrix<double, 3, 1> hyperbolic_basis_1(3);
  Eigen::Matrix<double, 3, 1> hyperbolic_basis_2(3);
  hyperbolic_basis_1 << 1, 0, 1;
  hyperbolic_basis_2 << 0, 2, 0;
  Q_coeffs << -1.0, 0.0, 1.0;
  if ((sigma_1 < 0) && (sign < 0)) {
    P_coeffs.col(0) = -k1 * hyperbolic_basis_2;
    P_coeffs.col(1) = k2 * hyperbolic_basis_1;
  } else if ((sigma_1 > 0) && (sign < 0)) {
    P_coeffs.col(0) = k1 * hyperbolic_basis_1;
    P_coeffs.col(1) = k2 * hyperbolic_basis_2;
  } else if ((sigma_1 < 0) && (sign > 0)) {
    P_coeffs.col(0) = k1 * hyperbolic_basis_1;
    P_coeffs.col(1) = -k2 * hyperbolic_basis_2;
  } else if ((sigma_1 > 0) && (sign > 0)) {
    P_coeffs.col(0) = k1 * hyperbolic_basis_2;
    P_coeffs.col(1) = k2 * hyperbolic_basis_1;
  }

  // Parametrize one branch over finite interval (-1,1)
  Conic hyperbola_branch(P_coeffs, Q_coeffs, ConicType::hyperbola);
  hyperbola_branch.domain().set_lower_bound(-1.0, true);
  hyperbola_branch.domain().set_upper_bound(1.0, true);
  conics.push_back(hyperbola_branch);

  // Parametrize other branch
  hyperbola_branch.set_numerators(-P_coeffs);
  conics.push_back(hyperbola_branch);
}

void
parametrize_parabola(const Vector6r& conic_standard_form,
                     std::vector<Conic>& conics)
{
  assert(identify_conic(conic_standard_form) == ConicType::parabola);

  // Get parabolic coefficient
  double sigma = conic_standard_form(4);

  // Parametrize ellipse as (1/(1 - t^2)) * [k1 (1 + t^2), 2 k2 t]
  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  P_coeffs << 0.0, 0.0, -1.0, 0.0, 0.0, -sigma;
  Q_coeffs << 1.0, 0.0, 0.0;

  conics.push_back(Conic(P_coeffs, Q_coeffs, ConicType::parabola));
}

// Add the intersecting lines as four rays with consistent orientation
void
parametrize_intersecting_lines(const Vector6r& conic_standard_form,
                               std::vector<Conic>& conics)
{
  assert(identify_conic(conic_standard_form) == ConicType::intersecting_lines);
  conics.clear();

  // Get absolute value of the slope of the lines
  double sigma_1 = conic_standard_form(4);
  double sigma_2 = conic_standard_form(5);

  // Get intervals for the positive and negative rays
  Interval nonpositives, nonnegatives;
  nonpositives.set_upper_bound(0.0, false);
  nonnegatives.set_lower_bound(0.0, false);

  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  P_coeffs << 0.0, 0.0, std::sqrt(abs(sigma_2)), std::sqrt(abs(sigma_1)), 0.0,
    0.0;
  Q_coeffs << 1.0, 0.0, 0.0;
  conics.push_back(
    Conic(P_coeffs, Q_coeffs, nonnegatives, ConicType::intersecting_lines));
  conics.push_back(Conic(
    -1.0 * P_coeffs, Q_coeffs, nonnegatives, ConicType::intersecting_lines));

  P_coeffs(1, 1) *= -1.0;
  conics.push_back(
    Conic(P_coeffs, Q_coeffs, nonpositives, ConicType::intersecting_lines));
  conics.push_back(Conic(
    -1.0 * P_coeffs, Q_coeffs, nonpositives, ConicType::intersecting_lines));
}

void
parametrize_parallel_lines(const Vector6r& conic_standard_form,
                           std::vector<Conic>& conics)
{
  assert(identify_conic(conic_standard_form) == ConicType::parallel_lines);
  conics.clear();

  // Get absolute value of the slope of the lines
  double c = conic_standard_form(0);
  double b1 = conic_standard_form(1);
  double sigma = conic_standard_form(4);
  double discriminant = compute_discriminant(sigma, b1, c);
  double x0 = (0.5 / sigma) * (-b1 - std::sqrt(discriminant));
  double x1 = (0.5 / sigma) * (-b1 + std::sqrt(discriminant));

  // Parametrize ellipse as (1/(1 - t^2)) * [k1 (1 + t^2), 2 k2 t]
  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  P_coeffs << x0, 0.0, 0.0, 1.0, 0.0, 0.0;
  Q_coeffs << 1.0, 0.0, 0.0;
  conics.push_back(Conic(P_coeffs, Q_coeffs, ConicType::parallel_lines));

  P_coeffs(0, 0) = x1;
  P_coeffs(1, 1) = -1.0;
  conics.push_back(Conic(P_coeffs, Q_coeffs, ConicType::parallel_lines));
}

void
parametrize_line(std::vector<Conic>& conics)
{
  conics.clear();

  // All lines are x=0 in standard form
  Matrix3x2r P_coeffs(3, 2);
  Eigen::Matrix<double, 3, 1> Q_coeffs(3);
  P_coeffs << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  Q_coeffs << 1.0, 0.0, 0.0;
  conics.push_back(Conic(P_coeffs, Q_coeffs, ConicType::line));
}

void
parametrize_standard_form_conic(const Vector6r& conic_standard_form,
                                std::vector<Conic>& conics)
{
  assert(is_conic_standard_form(conic_standard_form));
  conics.clear();

  // Parametrize conic based on type
  ConicType conic_type = identify_standard_form_conic(conic_standard_form);
  assert(conic_type == identify_standard_form_conic(-1 * conic_standard_form));

  if (conic_type == ConicType::ellipse) {
    spdlog::trace("Parametrizing ellipse");
    parametrize_ellipse(conic_standard_form, conics);
  } else if (conic_type == ConicType::hyperbola) {
    spdlog::trace("Parametrizing hyperbola");
    parametrize_hyperbola(conic_standard_form, conics);
  } else if (conic_type == ConicType::parabola) {
    spdlog::trace("Parametrizing parabola");
    parametrize_parabola(conic_standard_form, conics);
  } else if (conic_type == ConicType::parallel_lines) {
    spdlog::trace("Parametrizing parallel lines");
    parametrize_parallel_lines(conic_standard_form, conics);
  } else if (conic_type == ConicType::intersecting_lines) {
    spdlog::trace("Parametrizing intersecting lines");
    parametrize_intersecting_lines(conic_standard_form, conics);
  } else if (conic_type == ConicType::line) {
    spdlog::trace("Parametrizing line");
    parametrize_line(conics);
  } else if (conic_type == ConicType::point) {
    spdlog::trace("Skipping degenerate point");
  } else if (conic_type == ConicType::empty) {
    spdlog::trace("Skipping degenerate empty set");
  } else if (conic_type == ConicType::plane) {
    spdlog::error("Entire plane is the solution");
  } else if (conic_type == ConicType::error) {
    spdlog::error("Error in conic");
  } else {
    spdlog::error("Unknown conic");
  }
}

bool
check_standard_form(const Vector6r& conic_coeffs,
                    const Vector6r& conic_standard_form,
                    const Matrix2x2r& rotation,
                    const PlanarPoint& translation)
{
  SPDLOG_TRACE("Checking standard form {} for implicit function {}",
               formatted_bivariate_quadratic_mapping(conic_standard_form, 17),
               formatted_bivariate_quadratic_mapping(conic_coeffs, 17));
  Eigen::Matrix<double, 6, 6> change_of_basis_matrix;
  generate_quadratic_coordinate_affine_transformation_matrix<double>(
    rotation, translation, change_of_basis_matrix);
  Vector6r test_conic_coeffs = change_of_basis_matrix * conic_coeffs;
  SPDLOG_TRACE("Implicit function after rotation is {}",
               formatted_bivariate_quadratic_mapping(test_conic_coeffs, 17));
  return (column_vector_equal(test_conic_coeffs, conic_standard_form, 1e-4));
}

bool
check_parametrized_conic(const Conic& conic, const Vector6r& conic_coeffs)
{
  spdlog::trace(
    "Checking parametrization for conic {} with implicit function {}",
    conic,
    conic_coeffs);
  RationalFunction<4, 1> pullback;
  conic.pullback_quadratic_function<1>(conic_coeffs, pullback);
  spdlog::trace("Pullback by implicit function is {}", pullback);
  return float_equal(pullback.get_numerators().norm(), 0.0, 1e-4);
}

// FIXME Rename
bool
check_orientation(const Conic& conic, const Vector6r& conic_coeffs)
{
  double t = 0.1;
  if (!conic.is_in_domain(t)) {
    t = -0.1;
  }
  spdlog::trace(
    "Checking consistency for conic {} with implicit function {} at {}",
    conic,
    conic_coeffs,
    t);
  RationalFunction<4, 2> tangent;
  conic.compute_derivative(tangent);
  Eigen::Matrix<double, 3, 1> u_derivative =
    u_derivative_matrix() * conic_coeffs;
  Eigen::Matrix<double, 3, 1> v_derivative =
    v_derivative_matrix() * conic_coeffs;
  PlanarPoint point = conic(t);
  double perp_u = evaluate_line(u_derivative, point);
  double perp_v = evaluate_line(v_derivative, point);
  spdlog::trace("Contour gradient: [{}, {}]", perp_u, perp_v);
  PlanarPoint point_tangent = tangent(t);
  double tu = point_tangent[0];
  double tv = point_tangent[1];
  spdlog::trace("Contour parametric tangent: [{}, {}]", tu, tv);

  return ((-tv * perp_u + tu * perp_v) < 0);
}

void
parametrize_conic(const Vector6r& conic_coeffs, std::vector<Conic>& conics)
{
  conics.clear();
  SPDLOG_TRACE("Parametrizing conic with equation: {}",
               formatted_bivariate_quadratic_mapping(conic_coeffs));

  // Get standard form with rotation and translation
  Matrix2x2r rotation;
  PlanarPoint translation;
  Vector6r conic_standard_form;
  convert_conic_to_standard_form(
    conic_coeffs, conic_standard_form, rotation, translation);
  assert(check_standard_form(
    conic_coeffs, conic_standard_form, rotation, translation));

  // Get parametrization of conic
  std::vector<Conic> standard_form_conics;
  parametrize_standard_form_conic(conic_standard_form, standard_form_conics);
  for (size_t i = 0; i < standard_form_conics.size(); ++i) {
    Conic conic = standard_form_conics[i];

    // spdlog::warn(conic.pullback_quadratic_function(conic_standard_form).get_numerators().norm());
    assert(check_parametrized_conic(conic, conic_standard_form));
    assert(check_orientation(conic, conic_standard_form));

    conic.transform(rotation, translation);
    spdlog::trace("Parametrized conic:\n{}", conic);

    assert(check_parametrized_conic(conic, conic_coeffs));
    assert(check_orientation(conic, conic_coeffs));

    conics.push_back(conic);
  }
}

void
parametrize_cone_patch_conic(const Vector6r& conic_coeffs,
                             std::vector<Conic>& conics)
{
  SPDLOG_TRACE("Parametrizing conic from cone patch with equation: {}",
               formatted_bivariate_quadratic_mapping(conic_coeffs));

  // Get standard form with rotation and translation
  Matrix2x2r rotation;
  PlanarPoint translation;
  Vector6r conic_standard_form;
  convert_conic_to_standard_form(
    conic_coeffs, conic_standard_form, rotation, translation);

  // Get singular values
  double sigma_1 = conic_standard_form(4);
  double sigma_2 = conic_standard_form(5);
  double det = sigma_1 * sigma_2;

  // If the determinant is positive, the contour is a point and can be skipped
  if ((det > 0.0) || (float_equal_zero(det))) {
    spdlog::trace(
      "Skipping degenerate point contour in cone patch with determinant {}",
      det);
    return;
  }

  // Get parametrization of conic
  spdlog::trace("Parametrizing intersecting lines in cone patch");
  std::vector<Conic> standard_form_conics;
  parametrize_intersecting_lines(conic_standard_form, standard_form_conics);
  for (size_t i = 0; i < standard_form_conics.size(); ++i) {
    Conic conic = standard_form_conics[i];
    // spdlog::warn(conic.pullback_quadratic_function(conic_standard_form).get_numerators().norm());
    assert(check_parametrized_conic(conic, conic_standard_form));
    assert(check_orientation(conic, conic_standard_form));
    conic.transform(rotation, translation);
    if (!check_parametrized_conic(conic, conic_coeffs)) {
      spdlog::error("Did not parametrize implicit cone patch conic {} with {}",
                    conic_coeffs,
                    conic);
    }
    assert(check_orientation(conic, conic_coeffs));
    spdlog::trace("Parametrized conic:\n{}", conic);
    conics.push_back(conic);
  }
}
