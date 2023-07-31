// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "generate_colormap.h"

void
convert_diverging_function_to_colormap_function(
  const VectorXr& diverging_function,
  double range,
  double center,
  VectorXr& colormap_function)
{
  colormap_function.resize(diverging_function.size());

  for (int i = 0; i < diverging_function.size(); ++i) {
    // Linearly map the specified domain to [0, 1]
    colormap_function[i] = ((diverging_function(i) - center) / range) + 0.5;

    // Clip the colormap function to [0, 1]
    colormap_function[i] = std::max(std::min(colormap_function(i), 1.0), 0.0);
  }
}

void
convert_function_to_colormap_function(const VectorXr& scalar_function,
                                      double range_min,
                                      double range_max,
                                      VectorXr& colormap_function)
{
  colormap_function.resize(scalar_function.size());

  for (int i = 0; i < scalar_function.size(); ++i) {
    // Linearly map the specified domain to [0, 1]
    double range = range_max - range_min;
    colormap_function[i] = ((scalar_function(i) - range_min) / range);

    // Clip the colormap function to [0, 1]
    colormap_function[i] = std::max(std::min(colormap_function(i), 1.0), 0.0);
  }
}

void
convert_bounded_function_to_colormap_function(const VectorXr& bounded_function,
                                              double min,
                                              double max,
                                              VectorXr& colormap_function)
{
  colormap_function.resize(bounded_function.size());

  for (int i = 0; i < bounded_function.size(); ++i) {
    // Linearly map the specified domain to [0, 1]
    double range = max - min;
    colormap_function[i] = (bounded_function(i) - min) / range;

    // Clip the colormap function to [0, 1]
    colormap_function[i] = std::max(std::min(colormap_function(i), 1.0), 0.0);
  }
}

void
generate_binary_colormap(const VectorXr& colormap_function,
                         const Eigen::Matrix<double, 3, 1>& below_color,
                         const Eigen::Matrix<double, 3, 1>& above_color,
                         MatrixXr& colormap)
{
  colormap.resize(colormap_function.size(), 3);

  for (int i = 0; i < colormap_function.size(); ++i) {
    if (colormap_function(i) < 0.5) {
      colormap.row(i) = below_color;
    } else {
      colormap.row(i) = above_color;
    }
  }
}

void
generate_smooth_diverging_colormap(
  const VectorXr& colormap_function,
  const Eigen::Matrix<double, 3, 1>& below_color,
  const Eigen::Matrix<double, 3, 1>& above_color,
  int smoothness,
  MatrixXr& colormap)
{
  colormap.resize(colormap_function.size(), 3);
  Eigen::Matrix<double, 3, 1> white_color(1, 1, 1);

  for (int i = 0; i < colormap_function.size(); ++i) {
    // Determine scaling based on smoothness and distance from the center
    double color_scaling =
      std::pow(std::abs(colormap_function[i] - 0.5), smoothness + 1);

    // Set the color based on scaling and sign relative to the center
    if (colormap_function(i) < 0.5) {
      // colormap.row(i) = color_scaling * below_color;
      colormap.row(i) =
        color_scaling * below_color + (0.5 - color_scaling) * white_color;
    } else {
      colormap.row(i) =
        color_scaling * above_color + (0.5 - color_scaling) * white_color;
    }
  }
}

Eigen::Matrix<double, 3, 1>
generate_random_color()
{
  Eigen::Matrix<double, 3, 1> color;

  // Create random number generator
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(0, 1);

  // Generate random color
  color << dist(e2), dist(e2), dist(e2);

  return color;
}

void
generate_random_category_colormap(const std::vector<int>& category_labels,
                                  MatrixXr& colormap)
{
  colormap.setZero(category_labels.size(), 3);

  // Return if empty label set
  if (category_labels.empty())
    return;

  // Get array of random colors
  int max_index =
    *std::max_element(category_labels.begin(), category_labels.end());
  std::vector<Eigen::Matrix<double, 3, 1>> category_colors(max_index + 2);
  for (int i = 0; i <= max_index + 1; ++i) {
    category_colors[i] = generate_random_color();
  }

  // Assign colors corresponding to indices
  for (size_t i = 0; i < category_labels.size(); ++i) {
    colormap.row(i) = category_colors[category_labels[i] + 1];
  }
}

void
generate_interpolating_colormap(const VectorXr& colormap_function,
                                const Eigen::Matrix<double, 3, 1>& below_color,
                                const Eigen::Matrix<double, 3, 1>& above_color,
                                MatrixXr& colormap)
{
  colormap.resize(colormap_function.size(), 3);

  for (int i = 0; i < colormap_function.size(); ++i) {
    double s = colormap_function[i];
    colormap.row(i) = (1 - s) * below_color + s * above_color;
  }
}