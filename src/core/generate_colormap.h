// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// Given a diverging scalar function with a given center, linearly map the
/// range [center - range, center + range] to the colormap domain [0, 1] with
/// center 0.5.
///
/// @param[in] diverging_function: input function to map to the colormap domain
/// @param[in] range: range of domain values to send to the colormap domain
/// @param[in] center: center of the diverging function to map to 0.5
/// @param[out] colormap_function: output function with the proper colormap
/// domain
void
convert_diverging_function_to_colormap_function(
  const VectorXr& diverging_function,
  double range,
  double center,
  VectorXr& colormap_function);

/// Given a scalar function, linearly map the range [range_min, range_max] to
/// the
// colormap domain [0, 1].
///
/// @param[in] scalar_function: input function to map to the colormap domain
/// @param[in] range_min: minimum value to use for the range
/// @param[in] range_max: maximum value to use for the range
/// @param[out] colormap_function: output function with the proper colormap
/// domain
void
convert_function_to_colormap_function(const VectorXr& scalar_function,
                                      double range_min,
                                      double range_max,
                                      VectorXr& colormap_function);

/// Given a diverging colormap function, generate a colormap that is constantly
/// one color for values below 0.5 and constantly another color for values
/// above.
///
/// @param[in] colormap_function: diverging colormap function with range [0, 1]
/// @param[in] below_color: color to use for values below the center
/// @param[in] above_color: color to use for values above the center
/// @param[out] colormap: RGB matrix of color values
void
generate_binary_colormap(const VectorXr& colormap_function,
                         const Eigen::Matrix<double, 3, 1>& below_color,
                         const Eigen::Matrix<double, 3, 1>& above_color,
                         MatrixXr& colormap);

/// Given a diverging colormap function, generate a colormap that continuously
/// blends from one color for values below 0.5 to another color for values
/// above.
///
/// The C^n smoothness of the colormap is set with the smoothness parameter.
///
/// @param[in] colormap_function: diverging colormap function with range [0, 1]
/// @param[in] below_color: color to use for values below the center
/// @param[in] above_color: color to use for values above the center
/// @param[in] smoothness: order of derivatives to ensure continuity for
/// @param[out] colormap: RGB matrix of color values
void
generate_smooth_diverging_colormap(
  const VectorXr& colormap_function,
  const Eigen::Matrix<double, 3, 1>& below_color,
  const Eigen::Matrix<double, 3, 1>& above_color,
  int smoothness,
  MatrixXr& colormap);

/// Given a vector of category labels, generate a colormap that assigns a random
/// color to each label.
///
/// @param[in] category_labels: labels for discrete categories
/// @param[out] colormap: RGB matrix of color values
void
generate_random_category_colormap(const std::vector<int>& category_labels,
                                  MatrixXr& colormap);

/// Given a colormap function with range [0,1], generate a colormap that
/// continuously blends from one color to another.
///
/// @param[in] colormap_function: colormap function with range [0, 1]
/// @param[in] below_color: color to use for value 0
/// @param[in] above_color: color to use for value 1
/// @param[out] colormap: RGB matrix of color values
void
generate_interpolating_colormap(const VectorXr& colormap_function,
                                const Eigen::Matrix<double, 3, 1>& below_color,
                                const Eigen::Matrix<double, 3, 1>& above_color,
                                MatrixXr& colormap);

Eigen::Matrix<double, 3, 1>
generate_random_color();