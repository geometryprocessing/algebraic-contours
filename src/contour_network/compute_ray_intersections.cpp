// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_ray_intersections.h"

#include "compute_ray_intersections_pencil_method.h"
#include "intersection_heuristics.h"
#include "polynomial_function.h"
#include "ray_intersection_heuristics.h"

void
compute_spline_surface_ray_intersections(
  const QuadraticSplineSurface& spline_surface,
  const Matrix2x3r& ray_mapping_coeffs,
  std::vector<QuadraticSplineSurface::PatchIndex>& patch_indices,
  std::vector<PlanarPoint>& surface_intersections,
  std::vector<double>& ray_intersections,
  long long& ray_intersections_call,
  long long& ray_bounding_box_call)
{
  surface_intersections.clear();
  ray_intersections.clear();
  surface_intersections.reserve(100);
  ray_intersections.reserve(100);
  SPDLOG_TRACE(
    "Computing intersections for spline surface with {} patches and ray {}",
    num_patches,
    formatted_polynomial(ray_mapping_coeffs));

  PlanarPoint ray_plane_point = ray_mapping_coeffs.block(0, 0, 1, 2);
  std::pair<int, int> hash_indices =
    spline_surface.compute_hash_indices(ray_plane_point);

  for (size_t i :
       spline_surface.hash_table[hash_indices.first][hash_indices.second]) {
    int num_intersections;
    std::array<PlanarPoint, MAX_PATCH_RAY_INTERSECTIONS>
      patch_surface_intersections;
    std::array<double, MAX_PATCH_RAY_INTERSECTIONS> patch_ray_intersections;

    compute_spline_surface_patch_ray_intersections_pencil_method(
      spline_surface.get_patch(i),
      ray_mapping_coeffs,
      num_intersections,
      patch_surface_intersections,
      patch_ray_intersections,
      ray_intersections_call,
      ray_bounding_box_call);

    // Add patch intersections to surface intersections arrays
    if (num_intersections > MAX_PATCH_RAY_INTERSECTIONS) {
      spdlog::error("More than four intersections found of a ray with a patch");
    }
    for (int j = 0; j < num_intersections; ++j) {
      patch_indices.push_back(i);
      surface_intersections.push_back(patch_surface_intersections[j]);
      ray_intersections.push_back(patch_ray_intersections[j]);
      spdlog::trace("Patch ray intersection at t={} found, out of {}",
                    patch_ray_intersections[j],
                    num_intersections);
    }
  }
  SPDLOG_TRACE("{} surface ray intersections found",
               surface_intersections.size());
  SPDLOG_TRACE("Spline surface intersection points: {}",
               formatted_vector(surface_intersections));
  SPDLOG_TRACE("Ray intersection points: {}",
               formatted_vector(ray_intersections));
}

void
partition_ray_intersections(const Matrix2x3r& ray_mapping_coeffs,
                            const SpatialVector& comparison_point,
                            const std::vector<double>& ray_intersections,
                            std::vector<double>& ray_intersections_below,
                            std::vector<double>& ray_intersections_above)
{
  size_t num_intersections = ray_intersections.size();
  ray_intersections_below.clear();
  ray_intersections_above.clear();
  ray_intersections_below.reserve(num_intersections);
  ray_intersections_above.reserve(num_intersections);
  SPDLOG_TRACE("Partitioning intersections {} on ray {} around point {}",
               formatted_vector(ray_intersections),
               formatted_polynomial(ray_mapping_coeffs),
               comparison_point);

  // Compute parameter for point
  SpatialVector point_0 = ray_mapping_coeffs.row(0);
  SpatialVector direction = ray_mapping_coeffs.row(1);
  SpatialVector difference = comparison_point - point_0;
  double point_parameter = difference.norm() / direction.norm();
  spdlog::trace("Point parameter: {}", point_parameter);

  // Partition intersections into points above and below the comparison point
  for (size_t i = 0; i < num_intersections; ++i) {
    double ray_intersection = ray_intersections[i];
    if (float_equal(ray_intersection, point_parameter, 1e-7)) {
      continue;
    } else if (ray_intersection < point_parameter) {
      ray_intersections_below.push_back(ray_intersection);
    } else {
      ray_intersections_above.push_back(ray_intersection);
    }
  }
  SPDLOG_TRACE("Intersections below: {}",
               formatted_vector(ray_intersections_below, ", "));
  SPDLOG_TRACE("Intersections above: {}",
               formatted_vector(ray_intersections_above, ", "));
}
