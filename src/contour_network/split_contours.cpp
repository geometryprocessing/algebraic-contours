// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "split_contours.h"

#include "compute_intersections.h"
#include "compute_rational_bezier_curve_intersection.h"

void
split_contours(
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<RationalFunction<4, 2>>& planar_contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  std::vector<bool>& contour_is_boundary,
  std::vector<std::map<std::string, int>> contour_segment_labels)
{
  std::vector<std::vector<double>> split_knots(0);
  split_planar_curves_no_self_intersection(planar_contour_segments,
                                           split_knots);
  for (size_t i = 0; i < split_knots.size(); ++i) {
    if (split_knots[i].empty())
      continue;

    spdlog::warn("Splitting curve due to potential self intersections");
    Conic lower_domain_curve;
    RationalFunction<4, 3> lower_spatial_curve;
    RationalFunction<4, 2> lower_planar_curve;
    Conic upper_domain_curve;
    RationalFunction<4, 3> upper_spatial_curve;
    RationalFunction<4, 2> upper_planar_curve;
    for (size_t j = 0; j < split_knots[i].size(); ++i) {
      double knot = split_knots[i][j];
      contour_domain_curve_segments[i].split_at_knot(
        knot, lower_domain_curve, upper_domain_curve);
      planar_contour_segments[i].split_at_knot(
        knot, lower_planar_curve, upper_planar_curve);
      contour_segments[i].split_at_knot(
        knot, lower_spatial_curve, upper_spatial_curve);

      // Replace the original index with the upper segment for further
      // splitting
      contour_domain_curve_segments[i] = upper_domain_curve;
      planar_contour_segments[i] = upper_planar_curve;
      contour_segments[i] = upper_spatial_curve;

      // Push back the lower segment and the other necessary data in other
      // arrays
      contour_domain_curve_segments.push_back(lower_domain_curve);
      planar_contour_segments.push_back(lower_planar_curve);
      contour_segments.push_back(lower_spatial_curve);
      contour_patch_indices.push_back(contour_patch_indices[i]);
      contour_is_boundary.push_back(contour_is_boundary[i]);
      contour_segment_labels.push_back(contour_segment_labels[i]);
    }
  }
}
