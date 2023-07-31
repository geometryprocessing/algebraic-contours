// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "polynomial_function.h"

// Remove any zeros at the end of the polynomial coefficient vector
// FIXME make constant size if possible
void
remove_polynomial_trailing_coefficients(const VectorXr& A_coeffs,
                                        VectorXr& reduced_coeffs)
{
  // Find last nonzero entry and remove all zero entries after it
  int last_zero = A_coeffs.size();
  while (last_zero > 0) {
    if (!float_equal(A_coeffs[last_zero - 1], 0.0)) {
      reduced_coeffs = A_coeffs.head(last_zero);
      return;
    }

    --last_zero;
  }

  // If zero vector, set reduced coefficients to be empty
  // reduced_coeffs = VectorXr();
  reduced_coeffs = A_coeffs.head(1);
}

void
quadratic_real_roots(const Eigen::Matrix<double, 3, 1>& quadratic_coeffs,
                     std::array<double, 2>& solutions,
                     int& num_solutions,
                     double eps)
{
  double discriminant;

  if (eps <= std::abs(quadratic_coeffs[2])) {
    discriminant = -4 * quadratic_coeffs[0] * quadratic_coeffs[2] +
                   quadratic_coeffs[1] * quadratic_coeffs[1];
    if (eps * eps <= discriminant) {
      if (0.0 < quadratic_coeffs[1]) {
        solutions[0] = 2.0 * quadratic_coeffs[0] /
                       (-quadratic_coeffs[1] - sqrt(discriminant));
        solutions[1] = (-quadratic_coeffs[1] - sqrt(discriminant)) /
                       (2.0 * quadratic_coeffs[2]);
      } else {
        solutions[0] = (-quadratic_coeffs[1] + sqrt(discriminant)) /
                       (2.0 * quadratic_coeffs[2]);
        solutions[1] = 2.0 * quadratic_coeffs[0] /
                       (-quadratic_coeffs[1] + sqrt(discriminant));
      }
      num_solutions = 2;
    } else if (0.0 <= discriminant) {
      solutions[0] = -quadratic_coeffs[1] / (2.0 * quadratic_coeffs[2]);
      num_solutions = 1;
    } else
      num_solutions = 0;
  } else if (eps <= std::abs(quadratic_coeffs[1])) {
    solutions[0] = -quadratic_coeffs[0] / quadratic_coeffs[1];
    num_solutions = 1;
  } else
    num_solutions = 0;
}

std::vector<double>
polynomial_real_roots(const VectorXr& A_coeffs)
{
  SPDLOG_TRACE("Finding roots for {}", formatted_polynomial(A_coeffs));

  // Ensure the polynomial coefficients have no trailing zeros
  VectorXr reduced_coeffs;
  spdlog::trace("Full coefficient vector: {}", A_coeffs);
  remove_polynomial_trailing_coefficients(A_coeffs, reduced_coeffs);
  spdlog::trace("Reduced coefficient vector: {}", reduced_coeffs);

  std::vector<double> roots;
  // check if reduced coeff is 0
  if (reduced_coeffs.size() == 1 and float_equal_zero(reduced_coeffs[0]))
    return roots;

  // Compute the complex roots
  Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
  solver.compute(reduced_coeffs);
  auto solver_roots = solver.roots();
  spdlog::trace("Complex roots: {}", solver_roots);

  // Find the real roots
  // WARNING: Involves a floating point threshold test
  roots.reserve(reduced_coeffs.size());
  for (int i = 0; i < solver_roots.size(); ++i) {
    if (float_equal(solver_roots[i].imag(), 0.0)) {
      roots.push_back(solver_roots[i].real());
    }
  }
  SPDLOG_TRACE("Real roots: {}", formatted_vector(roots));

  return roots;
}

std::string
formatted_monomial(const std::string& variable, int degree)
{
  // Handle degree 0 case
  if (degree < 1)
    return "";

  // Format as "<variable>^<degree>"
  std::stringstream monomial_string;
  monomial_string << variable << "^" << degree;
  return monomial_string.str();
}

std::string
formatted_term(double coefficient, std::string variable, int precision)
{
  std::stringstream term_string;

  // Zero case
  if (float_equal(coefficient, 0.0)) {
    return "";
  }
  // Negative case
  else if (coefficient < 0) {
    term_string << " - " << std::fixed << std::setprecision(precision)
                << abs(coefficient) << " " << variable;
    return term_string.str();
  }
  // Positive case
  else {
    term_string << " + " << std::fixed << std::setprecision(precision)
                << coefficient << " " << variable;
    return term_string.str();
  }
}
