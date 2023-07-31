// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "conic.h"

ConicType
Conic::get_type() const
{
  return m_type;
}

// Assumes row vector points
void
Conic::transform(const Matrix2x2r& rotation, const PlanarPoint& translation)
{
  Matrix3x2r P_rot_coeffs =
    get_numerators() * rotation + get_denominator() * translation;
  // FIXME Remove
  // Matrix3x2r P_coeffs = get_numerators();
  // Eigen::Matrix<double, 3, 1> Px_coeffs = get_numerators().col(0);
  // Eigen::Matrix<double, 3, 1> Py_coeffs = get_numerators().col(1);
  // Eigen::Matrix<double, 3, 1> Q_coeffs = get_denominator();

  // P_rot_coeffs.col(0) = rotation(0, 0) * Px_coeffs + rotation(0, 1) *
  // Py_coeffs; P_rot_coeffs.col(1) = rotation(1, 0) * Px_coeffs + rotation(1,
  // 1) * Py_coeffs; P_rot_coeffs.col(0) += translation(0) * Q_coeffs;
  // P_rot_coeffs.col(1) += translation(1) * Q_coeffs;

  set_numerators(P_rot_coeffs);
}

bool
Conic::is_valid() const
{
  if (m_numerator_coeffs.cols() == 0)
    return false;
  if (m_denominator_coeffs.size() == 0)
    return false;

  return true;
}

// Generated formatted conic string
std::string
Conic::formatted_conic() const
{
  std::stringstream conic_string;
  conic_string << "1/(";
  conic_string << formatted_polynomial<2, 1>(m_denominator_coeffs);
  conic_string << ") [\n  ";
  for (int i = 0; i < m_numerator_coeffs.cols(); ++i) {
    conic_string << formatted_polynomial<2, 1>(m_numerator_coeffs.col(i));
    conic_string << ",\n  ";
  }
  conic_string << "], t in " << m_domain.formatted_interval();

  return conic_string.str();
}

std::ostream&
operator<<(std::ostream& out, const Conic& F)
{
  out << F.formatted_conic();
  return out;
}