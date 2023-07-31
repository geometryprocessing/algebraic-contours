// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "rational_function.h"

/// \file project_curves.h
///
/// Methods to compute spatial curves to planar curves in image space.

/// @brief Given a rational spatial curve and view direction, project the curve to a
/// planar curve
///
/// @param[in] spatial_curve: curve in R^3 to project
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[out] planar_curve: projected planar curve in R^2
void
project_curve(const RationalFunction<4, 3>& spatial_curve,
              const Matrix3x3r& frame,
              RationalFunction<4, 2>& planar_curve);

/// @brief Given a list of rational spatial curves and view direction, project the
/// curves to a list of planar curves.
///
/// @param[in] spatial_curves: curves in R^3 to project
/// @param[in] frame: 3x3 matrix defining the projection
/// @param[out] planar_curves: projected planar curves in R^2
void
project_curves(const std::vector<RationalFunction<4, 3>>& spatial_curves,
               const Matrix3x3r& frame,
               std::vector<RationalFunction<4, 2>>& planar_curves);
