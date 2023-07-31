// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "write_output.h"

#include "compute_curve_frame.h"
#include "discretize.h"

void
add_random_polyline_to_svg(
  std::vector<Eigen::Vector2f>& polyline,
  std::vector<svg::SVG::PathVertexProperty>& vertex_properties,
  svg::SVG& svgWriter)
{
  // Get color
  int r_color = (1 * 739) % 256;
  int g_color = (1 * 181) % 256;
  int b_color = (1 * 607) % 256;

  svgWriter.writePolyline(polyline,
                          0.25,
                          Color(r_color, g_color, b_color, 255),
                          false,
                          vertex_properties);
}

Eigen::Vector2f
project_point(SpatialVector p, const Matrix3x3r& frame, int scale, int offset)
{
  SpatialVector p_3d = p * frame;

  p_3d(0) = -scale * p_3d(0) + offset; // Inverted to account for SVG orientation
  p_3d(1) = scale * p_3d(1) + offset;
  p_3d(2) = scale * p_3d(2) + offset;

  Eigen::Vector3f p_3d_float = p_3d.cast<float>();
  Eigen::Vector2f p_2d;
  p_2d(0) = p_3d_float(0);
  p_2d(1) = p_3d_float(1);

  return p_2d;
}

Eigen::Vector2f
transform_point(PlanarPoint p, int scale, int offset)
{
  Eigen::Vector2f p_float = p.cast<float>();
  Eigen::Vector2f p_transform;
  p_transform(0) = -scale * p_float(0) + offset; // Inverted to account for SVG orientation
  p_transform(1) = scale * p_float(1) + offset;

  return p_transform;
}

void
add_curve_network_to_svg(const std::vector<SpatialVector>& points,
                         const std::vector<std::vector<int>>& polylines,
                         const Matrix3x3r& frame,
                         svg::SVG& svgWriter,
                         int scale,
                         int offset)
{
  for (size_t i = 0; i < polylines.size(); ++i) {
    std::vector<int> polyline = polylines[i];
    std::vector<Eigen::Vector2f> polyline_points(polyline.size());
    std::vector<svg::SVG::PathVertexProperty> vertex_properties;
    vertex_properties.reserve(polyline.size());

    for (size_t j = 0; j < polyline.size(); ++j) {
      SpatialVector point = points[polyline[j]];
      polyline_points[j] = project_point(point, frame, scale, offset);
      vertex_properties.push_back(svg::SVG::PathVertexProperty());
    }
    svgWriter.writePolyline(
      polyline_points, 1.0, Color(0, 0, 0, 1), false, vertex_properties);
  }
}

void
add_curve_to_svg(const std::vector<PlanarPoint>& points,
                 const std::vector<int>& polyline,
                 svg::SVG& svgWriter,
                 int scale,
                 int offset,
                 Color color)
{
  std::vector<Eigen::Vector2f> polyline_points(polyline.size());
  std::vector<svg::SVG::PathVertexProperty> vertex_properties;
  vertex_properties.reserve(polyline.size());

  for (size_t j = 0; j < polyline.size(); ++j) {
    PlanarPoint point = points[polyline[j]];
    polyline_points[j] = transform_point(point, scale, offset);
    vertex_properties.push_back(svg::SVG::PathVertexProperty());
  }
  svgWriter.writePolyline(
    polyline_points, 1.0, color, false, vertex_properties);
}

void
add_parametric_curve_to_svg(const std::vector<PlanarPoint> points,
                            svg::SVG& svgWriter,
                            int scale,
                            int offset)
{
  std::vector<Eigen::Vector2f> polyline;
  std::vector<svg::SVG::PathVertexProperty> vertex_properties;
  for (size_t i = 0; i < points.size(); ++i) {
    Eigen::Vector2f p_2d = points[i].cast<float>();
    p_2d(0) = scale * p_2d(0) + offset;
    p_2d(1) = scale * p_2d(1) + offset;
    polyline.push_back(p_2d);
    vertex_properties.push_back(svg::SVG::PathVertexProperty());
  }

  svgWriter.writePolyline(
    polyline, 0.1, Color(1, 0, 0, 1), false, vertex_properties);
}

void
write_planar_point(const PlanarPoint point,
                   svg::SVG& svgWriter,
                   int scale,
                   int offset,
                   Color color)
{
  Eigen::Vector2f point_2d = transform_point(point, scale, offset);
  svgWriter.writeDot(point_2d, 2.0, color);
}

void
add_point_to_svg(const SpatialVector point,
                 const Matrix3x3r& frame,
                 svg::SVG& svgWriter,
                 int scale,
                 int offset,
                 Color color)
{
  SpatialVector point_3d = point * frame;

  Eigen::Vector3f point_3d_float = point_3d.cast<float>();
  Eigen::Vector2f point_2d;
  point_2d(0) = point_3d_float(0);
  point_2d(1) = point_3d_float(1);

  point_2d(0) = scale * point_2d(0) + offset;
  point_2d(1) = scale * point_2d(1) + offset;

  svgWriter.writeDot(point_2d, 0.25, color);
}

void
add_vectors_to_svg(const std::vector<SpatialVector> base_points,
                   const std::vector<SpatialVector> vectors,
                   const Matrix3x3r& frame,
                   svg::SVG& svgWriter,
                   bool normalize,
                   int scale,
                   int offset)
{
  assert(base_points.size() == vectors.size());
  for (size_t i = 0; i < base_points.size(); ++i) {
    std::vector<Eigen::Vector2f> vector_line;
    std::vector<svg::SVG::PathVertexProperty> vertex_properties;

    SpatialVector base_point_3d = base_points[i] * frame;
    SpatialVector vector_tip_3d = vectors[i] * frame;

    Eigen::Vector3f base_point_3d_float = base_point_3d.cast<float>();
    Eigen::Vector2f base_point_2d;
    base_point_2d(0) = base_point_3d_float(0);
    base_point_2d(1) = base_point_3d_float(1);

    Eigen::Vector3f vector_tip_3d_float = vector_tip_3d.cast<float>();
    Eigen::Vector2f vector_tip_2d;
    vector_tip_2d(0) = vector_tip_3d_float(0);
    vector_tip_2d(1) = vector_tip_3d_float(1);

    base_point_2d(0) = scale * base_point_2d(0) + offset;
    base_point_2d(1) = scale * base_point_2d(1) + offset;

    if (normalize) {
      vector_tip_2d /= vector_tip_2d.norm();
    }
    vector_tip_2d = 0.5 * scale * vector_tip_2d + base_point_2d;

    vector_line.push_back(base_point_2d);
    vector_line.push_back(vector_tip_2d);
    vertex_properties.push_back(svg::SVG::PathVertexProperty());
    vertex_properties.push_back(svg::SVG::PathVertexProperty());

    svgWriter.writePolyline(
      vector_line, 0.1, Color(0, 0, 1, 1), false, vertex_properties);
  }
}

void
add_parametric_vectors_to_svg(const std::vector<PlanarPoint> base_points,
                              const std::vector<PlanarPoint> vectors,
                              svg::SVG& svgWriter,
                              bool normalize,
                              int scale,
                              int offset)
{
  assert(base_points.size() == vectors.size());
  for (size_t i = 0; i < base_points.size(); ++i) {
    std::vector<Eigen::Vector2f> vector_line;
    std::vector<svg::SVG::PathVertexProperty> vertex_properties;

    Eigen::Vector2f base_point_2d = base_points[i].cast<float>();

    Eigen::Vector2f vector_tip_2d = vectors[i].cast<float>();

    base_point_2d(0) = scale * base_point_2d(0) + offset;
    base_point_2d(1) = scale * base_point_2d(1) + offset;

    if (normalize) {
      vector_tip_2d /= vector_tip_2d.norm();
    }
    vector_tip_2d = 0.1 * scale * vector_tip_2d + base_point_2d;

    vector_line.push_back(base_point_2d);
    vector_line.push_back(vector_tip_2d);
    vertex_properties.push_back(svg::SVG::PathVertexProperty());
    vertex_properties.push_back(svg::SVG::PathVertexProperty());

    svgWriter.writePolyline(
      vector_line, 0.1, Color(0, 0, 1, 1), false, vertex_properties);
  }
}

void
sample_surface(Eigen::MatrixXd& V,
               const Matrix3x3r& frame,
               svg::SVG& svgWriter,
               int scale,
               int offset)
{
  for (int i = 0; i < V.rows(); ++i) {
    SpatialVector p = V.row(i);
    p = p * frame;
    Eigen::Vector2f dot(2);
    dot(0) = scale * p(0) + offset;
    dot(1) = scale * p(1) + offset;
    svgWriter.writeDot(dot, 0.1, Color(0, 1, 0, 0.1));
  }
}

void
write_svg(std::string filename)
{
  Eigen::Vector2i viewport(600, 400);
  svg::SVG svgWriter(filename, viewport);
}

void
write_mesh(std::string filename,
           const Eigen::MatrixXd& V,
           const Eigen::MatrixXi& F)
{
  // Open file
  std::ofstream file(filename);

  // Write vertices to the file
  for (int i = 0; i < V.rows(); ++i) {
    file << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
  }

  // Write faces to the file
  for (int i = 0; i < F.rows(); ++i) {
    file << "f " << F(i, 0) + 1 << " " << F(i, 1) + 1 << " " << F(i, 2) + 1
         << std::endl;
  }
  // Close the file
  file.close();
}

void
write_mesh(std::string filename,
           const Eigen::MatrixXd& V,
           const Eigen::MatrixXi& F,
           const Eigen::MatrixXd& N)
{
  // Open file
  std::ofstream file(filename);

  // Write vertices to the file
  for (int i = 0; i < V.rows(); ++i) {
    file << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
  }

  // Write vertex normals to the file
  for (int i = 0; i < N.rows(); ++i) {
    file << "vn " << N(i, 0) << " " << N(i, 1) << " " << N(i, 2) << std::endl;
  }

  // Write faces to the file
  for (int i = 0; i < F.rows(); ++i) {
    file << "f " << F(i, 0) + 1 << "//" << F(i, 0) + 1 << " " << F(i, 1) + 1
         << "//" << F(i, 1) + 1 << " " << F(i, 2) + 1 << "//" << F(i, 2) + 1
         << std::endl;
  }
  // Close the file
  file.close();
}

void
write_contour_tangents(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const std::vector<Conic>& parameter_contour_segments,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  const CurveDiscretizationParameters& curve_disc_params,
  svg::SVG& svgWriter)
{
  std::vector<SpatialVector> tangent_points;
  std::vector<SpatialVector> tangents;
  std::vector<SpatialVector> normals;
  std::vector<SpatialVector> tangent_normals;
  sample_contour_frames(spline_surface,
                        parameter_contour_segments,
                        contour_segments,
                        contour_patch_indices,
                        curve_disc_params,
                        tangent_points,
                        tangents,
                        normals,
                        tangent_normals);

  add_vectors_to_svg(tangent_points, tangents, frame, svgWriter, false);

  if (!tangent_points.empty()) {
    add_point_to_svg(tangent_points.front(), frame, svgWriter);
    add_point_to_svg(tangent_points.back(), frame, svgWriter);
  }
}

void
write_contours(const Matrix3x3r& frame,
               const std::vector<RationalFunction<4, 3>>& contour_segments,
               const CurveDiscretizationParameters curve_disc_params,
               svg::SVG& svgWriter)
{
  std::vector<SpatialVector> points;
  std::vector<std::vector<int>> polylines;

  discretize_curve_segments<4, 3>(
    contour_segments, curve_disc_params, points, polylines);

  add_curve_network_to_svg(points, polylines, frame, svgWriter);
}

void
write_planar_curve_segment(
  const RationalFunction<4, 2>& planar_curve_segment,
  const CurveDiscretizationParameters curve_disc_params,
  svg::SVG& svgWriter,
  int scale,
  int offset,
  Color color)
{
  std::vector<PlanarPoint> points;
  std::vector<int> polyline;
  planar_curve_segment.discretize(curve_disc_params, points, polyline);
  add_curve_to_svg(points, polyline, svgWriter, scale, offset, color);
}

void
write_contour_points(const Matrix3x3r& frame,
                     const RationalFunction<4, 3>& contour_segment,
                     const std::vector<double>& parameter_points,
                     svg::SVG& svgWriter,
                     int scale,
                     int offset,
                     Color color)
{
  // Write points
  for (size_t i = 0; i < parameter_points.size(); ++i) {
    SpatialVector point = contour_segment(parameter_points[i]);
    add_point_to_svg(point, frame, svgWriter, scale, offset, color);
  }
}

void
write_contour_points(
  const Matrix3x3r& frame,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<double>>& parameter_points,
  svg::SVG& svgWriter,
  int scale,
  int offset,
  Color color)
{
  // Write points
  for (size_t i = 0; i < parameter_points.size(); ++i) {
    write_contour_points(frame,
                         contour_segments[i],
                         parameter_points[i],
                         svgWriter,
                         scale,
                         offset,
                         color);
  }
}

void
write_contours_with_annotations(
  const Matrix3x3r& frame,
  const std::vector<RationalFunction<4, 3>>& contour_segments,
  const std::vector<std::vector<double>> interior_cusps,
  const std::vector<std::vector<double>> boundary_cusps,
  const std::vector<std::vector<double>> intersections,
  const CurveDiscretizationParameters curve_disc_params,
  svg::SVG& svgWriter,
  int scale,
  int offset)
{
  // Write contours
  write_contours(frame, contour_segments, curve_disc_params, svgWriter);

  // Write points
  write_contour_points(frame,
                       contour_segments,
                       interior_cusps,
                       svgWriter,
                       scale,
                       offset,
                       Color(1, 0, 0, 1) // Red
  );
  write_contour_points(frame,
                       contour_segments,
                       boundary_cusps,
                       svgWriter,
                       scale,
                       offset,
                       Color(1, 0.75, 0.75, 1) // Pink
  );
  write_contour_points(frame,
                       contour_segments,
                       intersections,
                       svgWriter,
                       scale,
                       offset,
                       Color(0, 1, 0, 1) // Green
  );
}