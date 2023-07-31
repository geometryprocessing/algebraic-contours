// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "discretize.h"

#include "compute_contours.h"
#include "compute_curve_frame.h"
#include "compute_cusps.h"
#include "evaluate_general_contour_frame.h"


// Check that a parametrized line segment satisfies a given implicit equation
bool
satisfies_implicit_line_equation(
  const LineSegment& line_segment,
  const Eigen::Matrix<double, 3, 1>& implicit_line_coeffs)
{
  RationalFunction<1, 1> composition;
  line_segment.pullback_linear_function<1>(implicit_line_coeffs, composition);
  Eigen::Matrix<double, 2, 1> const& numerators = composition.get_numerators();
  if (!float_equal(numerators(0), 0.0))
    return false;
  if (!float_equal(numerators(1), 0.0))
    return false;

  return true;
}

void
discretize_patch_boundaries(const QuadraticSplineSurface& spline_surface,
                            std::vector<SpatialVector>& points,
                            std::vector<std::vector<size_t>>& polylines)
{
  points.clear();
  polylines.clear();

  for (QuadraticSplineSurface::PatchIndex patch_index = 0;
       patch_index < spline_surface.num_patches();
       ++patch_index) {

    std::array<LineSegment, 3> patch_boundaries;
    auto& spline_surface_patch = spline_surface.get_patch(patch_index);
    spline_surface_patch.get_domain().parametrize_patch_boundaries(
      patch_boundaries);
    for (size_t k = 0; k < patch_boundaries.size(); ++k) {
      //  Get points on the boundary curve
      std::vector<PlanarPoint> parameter_points_k;
      patch_boundaries[k].sample_points(5, parameter_points_k);
      std::vector<SpatialVector> points_k(parameter_points_k.size());
      for (size_t i = 0; i < parameter_points_k.size(); ++i) {
        spline_surface_patch.evaluate(parameter_points_k[i], points_k[i]);
      }

      // Build polyline for the given curve
      std::vector<size_t> polyline;
      polyline.resize(points_k.size());
      for (size_t l = 0; l < points_k.size(); ++l) {
        polyline[l] = points.size() + l;
      }

      append(points, points_k);
      polylines.push_back(polyline);
    }
  }
}

void
sample_contour_frames(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<Conic>& contour_domain_curve_segments,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  const CurveDiscretizationParameters& curve_disc_params,
  std::vector<SpatialVector>& base_points,
  std::vector<SpatialVector>& tangents,
  std::vector<SpatialVector>& normals,
  std::vector<SpatialVector>& tangent_normals)
{
  base_points.clear();
  tangents.clear();
  normals.clear();
  tangent_normals.clear();

  // Sample a given number of frames from each contour segment
  int num_points = curve_disc_params.num_tangents_per_segment;
  for (size_t k = 0; k < contour_segments.size(); ++k) {
    Conic parameter_contour_segment = contour_domain_curve_segments[k];
    RationalFunction<4, 3> contour_segment = contour_segments[k];
    QuadraticSplineSurface::PatchIndex patch_indices = contour_patch_indices[k];

    // Compute the frame functions for the given segment
    RationalFunction<8, 3> contour_segment_tangent;
    RationalFunction<4, 3> contour_segment_normal;
    RationalFunction<12, 3> contour_segment_tangent_normal;
    compute_spline_surface_patch_curve_frame(
      spline_surface.get_patch(patch_indices),
      parameter_contour_segment,
      contour_segment_tangent,
      contour_segment_normal,
      contour_segment_tangent_normal);

    // Add the frame samples for the contour segment to the lists
    std::vector<SpatialVector> segment_base_points, segment_tangents,
      segment_normals, segment_tangent_normals;
    contour_segment.sample_points(num_points, segment_base_points);
    contour_segment_tangent.sample_points(num_points, segment_tangents);
    contour_segment_normal.sample_points(num_points, segment_normals);
    contour_segment_tangent_normal.sample_points(num_points,
                                                 segment_tangent_normals);

    append(base_points, segment_base_points);
    append(tangents, segment_tangents);
    append(normals, segment_normals);
    append(tangent_normals, segment_tangent_normals);
  }
}