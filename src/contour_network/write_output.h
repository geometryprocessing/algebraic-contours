// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"
#include "svg.h"

void
sample_surface(Eigen::MatrixXd& V,
               const Matrix3x3r& frame,
               svg::SVG& svgWriter,
               int scale = 800,
               int offset = 400);

void
add_curve_to_svg(const std::vector<PlanarPoint>& points,
                 const std::vector<int>& polyline,
                 svg::SVG& svgWriter,
                 int scale,
                 int offset,
                 Color color = Color(0, 0, 0, 1));

void
add_curve_network_to_svg(const std::vector<SpatialVector>& points,
                         const std::vector<std::vector<int>>& polylines,
                         const Matrix3x3r& frame,
                         svg::SVG& svgWriter,
                         int scale = 800,
                         int offset = 400);

void
add_parametric_curve_to_svg(const std::vector<PlanarPoint> points,
                            svg::SVG& svgWriter,
                            int scale = 800,
                            int offset = 400);

void
add_vectors_to_svg(const std::vector<SpatialVector> base_points,
                   const std::vector<SpatialVector> vectors,
                   const Matrix3x3r& frame,
                   svg::SVG& svgWriter,
                   bool normalize = false,
                   int scale = 800,
                   int offset = 400);

void
add_parametric_vectors_to_svg(const std::vector<PlanarPoint> base_points,
                              const std::vector<PlanarPoint> vectors,
                              svg::SVG& svgWriter,
                              bool normalize,
                              int scale = 800,
                              int offset = 400);

void
write_mesh(std::string filename,
           const Eigen::MatrixXd& V,
           const Eigen::MatrixXi& F);

void
write_mesh(std::string filename,
           const Eigen::MatrixXd& V,
           const Eigen::MatrixXi& F,
           const Eigen::MatrixXd& N);

void
add_point_to_svg(const SpatialVector point,
                 const Matrix3x3r& frame,
                 svg::SVG& svgWriter,
                 int scale = 800,
                 int offset = 400,
                 Color color = Color(1, 0, 0, 1));

void
write_planar_point(const PlanarPoint point,
                   svg::SVG& svgWriter,
                   int scale = 800,
                   int offset = 400,
                   Color color = Color(1, 0, 0, 1));

void
write_contours(const Matrix3x3r& frame,
               const std::vector<RationalFunction<4, 3>>& contour_segments,
               const CurveDiscretizationParameters curve_disc_params,
               svg::SVG& svgWriter);

void
write_planar_curve_segment(
  const RationalFunction<4, 2>& planar_curve_segment,
  const CurveDiscretizationParameters curve_disc_params,
  svg::SVG& svgWriter,
  int scale = 800,
  int offset = 400,
  Color color = Color(0, 0, 0, 1));

void
write_contour_tangents(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const std::vector<Conic>& parameter_contour_segments,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  const CurveDiscretizationParameters& curve_disc_params,
  svg::SVG& svgWriter);

void
write_contour_points(const Matrix3x3r& frame,
                     const RationalFunction<4, 3>& contour_segment,
                     const std::vector<double>& parameter_points,
                     svg::SVG& svgWriter,
                     int scale,
                     int offset,
                     Color color);

void
write_contour_points(
  const Matrix3x3r& frame,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<double>>& parameter_points,
  svg::SVG& svgWriter,
  int scale = 800,
  int offset = 400,
  Color color = Color(1, 0, 0, 1));

void
write_contours_with_annotations(
  const Matrix3x3r& frame,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<double>> interior_cusps,
  const std::vector<std::vector<double>> boundary_cusps,
  const std::vector<std::vector<double>> intersections,
  const CurveDiscretizationParameters curve_disc_params,
  svg::SVG& svgWriter,
  int scale = 800,
  int offset = 400);
