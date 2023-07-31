// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "project_curves.h"

#include "validity.h"

void
project_curve(const RationalFunction<4, 3>& spatial_curve,
              const Matrix3x3r& frame,
              RationalFunction<4, 2>& planar_curve)
{
  assert(is_valid_frame(frame));

  // Get coordinates for spatial curve
  spdlog::trace("Getting spatial curve coefficients");
  Eigen::Matrix<double, 5, 3> spatial_P_coeffs = spatial_curve.get_numerators();
  Eigen::Matrix<double, 5, 1> Q_coeffs = spatial_curve.get_denominator();

  // Build planar curve coefficients
  spdlog::trace("Getting planar curve coefficients");
  Eigen::Matrix<double, 5, 2> planar_P_coeffs;
  planar_P_coeffs.col(0) = spatial_P_coeffs * frame.row(0).transpose();
  planar_P_coeffs.col(1) = spatial_P_coeffs * frame.row(1).transpose();

  // Build planar curve
  spdlog::trace("Getting planar curve");
  planar_curve =
    RationalFunction<4, 2>(planar_P_coeffs, Q_coeffs, spatial_curve.domain());
}

void
project_curves(const std::vector<RationalFunction<4, 3>>& spatial_curves,
               const Matrix3x3r& frame,
               std::vector<RationalFunction<4, 2>>& planar_curves)
{
  assert(is_valid_frame(frame));

  size_t num_curves = spatial_curves.size();
  planar_curves.resize(num_curves);
  for (size_t i = 0; i < num_curves; ++i) {
    project_curve(spatial_curves[i], frame, planar_curves[i]);
  }
}