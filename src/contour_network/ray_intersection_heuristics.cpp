// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "ray_intersection_heuristics.h"

#include "intersection_heuristics.h"
#include "polynomial_function.h"
#include "project_curves.h"

// Compute the projection of the patch boundaries to the plane
// void compute_spline_surface_patch_projected_patch_boundaries(
//   const QuadraticSplineSurfacePatch &spline_surface_patch,
//   const MatrixXr &frame,
//   std::array<RationalFunction, 3> &projected_patch_boundaries
// ) {
//   // Compute the spatial patch boundaries
//   std::array<RationalFunction, 3> patch_boundaries;
//   spline_surface_patch.get_patch_boundaries(patch_boundaries);

//   // Project the patch boundaries to the plane
//   for (size_t i = 0; i < patch_boundaries.size(); ++i)
//   {
//     project_curve(
//       patch_boundaries[i],
//       frame,
//       projected_patch_boundaries[i]
//     );
//   }
// }

// void compute_spline_surface_patch_boundaries_bounding_box(
//   const QuadraticSplineSurfacePatch &spline_surface_patch,
//   const MatrixXr &frame,
//   VectorXr &lower_left_point,
//   VectorXr &upper_right_point
// ) {
//   // Compute the projected patch boundaries
//   std::array<RationalFunction, 3> projected_patch_boundaries;
//   compute_spline_surface_patch_projected_patch_boundaries(
//     spline_surface_patch,
//     frame,
//     projected_patch_boundaries
//   );

//   // Compute the bounding boxes of the boundaries and combine them
//   VectorXr current_lower_left_point, prev_lower_left_point;
//   VectorXr current_upper_right_point, prev_upper_right_point;
//   compute_bezier_bounding_box(
//     projected_patch_boundaries[0],
//     lower_left_point,
//     upper_right_point
//   );
//   for (size_t i = 1; i < projected_patch_boundaries.size(); ++i)
//   {
//     prev_lower_left_point = lower_left_point;
//     prev_upper_right_point = upper_right_point;
//     compute_bezier_bounding_box(
//       projected_patch_boundaries[i],
//       current_lower_left_point,
//       current_upper_right_point
//     );
//     combine_bounding_boxes(
//       prev_lower_left_point,
//       prev_upper_right_point,
//       current_lower_left_point,
//       current_upper_right_point,
//       lower_left_point,
//       upper_right_point
//     );
//   }
// }

// void compute_spline_surface_boundaries_bounding_boxes(
//   const QuadraticSplineSurface &spline_surface,
//   const MatrixXr &frame,
//   std::vector<VectorXr> &lower_left_points,
//   std::vector<VectorXr> &upper_right_points
// ) {
//   size_t num_patches = spline_surface.num_patches();
//   lower_left_points.resize(num_patches);
//   upper_right_points.resize(num_patches);
//   for (size_t i = 0; i < num_patches; ++i)
//   {
//     compute_spline_surface_patch_boundaries_bounding_box(
//       spline_surface.get_patch(i),
//       frame,
//       lower_left_points[i],
//       upper_right_points[i]
//     );
//   }
// }

void
compute_spline_surface_bounding_boxes(
  const QuadraticSplineSurface& spline_surface,
  std::vector<SpatialVector>& min_points,
  std::vector<SpatialVector>& max_points)
{
  size_t num_patches = spline_surface.num_patches();
  min_points.resize(num_patches);
  max_points.resize(num_patches);
  for (size_t i = 0; i < num_patches; ++i) {
    spline_surface.get_patch(i).get_bounding_box(min_points[i], max_points[i]);
  }
}
