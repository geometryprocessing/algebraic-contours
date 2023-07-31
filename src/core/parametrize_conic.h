// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "conic.h"

/// @brief Identify the type of conic C represented by conic_coeffs.
///
/// @param[in] conic_coeffs: quadratic coefficients for the conic
/// @return type identifier of the conic
ConicType
identify_conic(const Vector6r& conic_coeffs);

/// @brief Given an implicit quadratic function for a conic, parameterize the
/// conic with explicit rational mappings.
///
/// @param[in] conic_coeffs: implicit conic equation
/// @param[out] conics: segments of the parametrized conic
void
parametrize_conic(const Vector6r& conic_coeffs, std::vector<Conic>& conics);

/// @brief Given an implicit quadratic function for a conic from a singular cone
/// patch, parametrize the conic with explicit rational mappings.
///
/// In this case, the conic can only be a point or intersecting lines.
///
/// @param[in] conic_coeffs: implicit conic equation
/// @param[out] conics: segments of the parametrized conic
void
parametrize_cone_patch_conic(const Vector6r& conic_coeffs,
                             std::vector<Conic>& conics);
