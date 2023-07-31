// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "convex_polygon.h"

// Compute the implicit form of a line between two points
void
compute_line_between_points(const PlanarPoint& point_0,
                            const PlanarPoint& point_1,
                            Eigen::Matrix<double, 3, 1>& line_coeffs)
{
  double x0 = point_0[0];
  double y0 = point_0[1];
  double x1 = point_1[0];
  double y1 = point_1[1];

  line_coeffs[0] = x0 * y1 - x1 * y0;
  line_coeffs[1] = y0 - y1;
  line_coeffs[2] = x1 - x0;
}

// Compute the parametric form of a line between two points
void
compute_parametric_line_between_points(const PlanarPoint& point_0,
                                       const PlanarPoint& point_1,
                                       LineSegment& line_segment)
{
  // Set numerator
  Matrix2x2r numerators;
  numerators(0, 0) = point_0(0);
  numerators(0, 1) = point_0(1);
  numerators(1, 0) = point_1(0) - point_0(0);
  numerators(1, 1) = point_1(1) - point_0(1);

  // Set domain interval [0, 1]
  Interval domain;
  domain.set_lower_bound(0, false);
  domain.set_upper_bound(1, false);

  line_segment = LineSegment(numerators, domain);
}

// Refine a mesh with midpoint subdivision
void
refine_triangles(const Eigen::MatrixXd& V,
                 const Eigen::MatrixXi& F,
                 Eigen::MatrixXd& V_refined,
                 Eigen::MatrixXi& F_refined)
{
  Eigen::Index num_faces = F.rows();
  V_refined.resize(num_faces * 6, 2);
  F_refined.resize(num_faces * 4, 3);
  for (Eigen::Index i = 0; i < num_faces; ++i) {
    Eigen::VectorXd v0 = V.row(F(i, 0));
    Eigen::VectorXd v1 = V.row(F(i, 1));
    Eigen::VectorXd v2 = V.row(F(i, 2));

    // Add vertices for refined face
    V_refined.row(6 * i + 0) = v0;
    V_refined.row(6 * i + 1) = v1;
    V_refined.row(6 * i + 2) = v2;
    V_refined.row(6 * i + 3) = (v0 + v1) / 2.0;
    V_refined.row(6 * i + 4) = (v1 + v2) / 2.0;
    V_refined.row(6 * i + 5) = (v2 + v0) / 2.0;

    // Add refined faces
    F_refined.row(4 * i + 0) = Eigen::Vector3i(6 * i + 0, 6 * i + 3, 6 * i + 5);
    F_refined.row(4 * i + 1) = Eigen::Vector3i(6 * i + 1, 6 * i + 4, 6 * i + 3);
    F_refined.row(4 * i + 2) = Eigen::Vector3i(6 * i + 2, 6 * i + 5, 6 * i + 4);
    F_refined.row(4 * i + 3) = Eigen::Vector3i(6 * i + 3, 6 * i + 4, 6 * i + 5);
  }
}

ConvexPolygon::ConvexPolygon(
  const std::array<Eigen::Matrix<double, 3, 1>, 3>& boundary_segments_coeffs)
  : m_boundary_segments_coeffs(boundary_segments_coeffs)
{
  PlanarPoint v0, v1, v2;
  intersect_patch_boundaries(
    boundary_segments_coeffs[1], boundary_segments_coeffs[2], v0);
  intersect_patch_boundaries(
    boundary_segments_coeffs[2], boundary_segments_coeffs[0], v1);
  intersect_patch_boundaries(
    boundary_segments_coeffs[0], boundary_segments_coeffs[1], v2);
  m_vertices.row(0) = v0;
  m_vertices.row(1) = v1;
  m_vertices.row(2) = v2;
}

ConvexPolygon::ConvexPolygon(const Eigen::Matrix<double, 3, 2>& vertices)
  : m_vertices(vertices)
{
  size_t num_vertices = vertices.rows();
  for (size_t i = 0; i < num_vertices; ++i) {
    compute_line_between_points(vertices.row(i),
                                vertices.row((i + 1) % num_vertices),
                                m_boundary_segments_coeffs[i]);
  }
}

bool
ConvexPolygon::contains(const PlanarPoint& point) const
{
  for (size_t i = 0; i < m_boundary_segments_coeffs.size(); i++) {
    Eigen::Matrix<double, 3, 1> L_coeffs = m_boundary_segments_coeffs[i];
    if ((L_coeffs(0) + L_coeffs(1) * point(0) + L_coeffs(2) * point(1)) < 0.0) {
      return false;
    }
  }

  return true;
}

void
ConvexPolygon::triangulate(size_t num_refinements,
                           Eigen::MatrixXd& V,
                           Eigen::MatrixXi& F) const
{
  // TODO Can generalize to arbitrary domain if needed
  V = m_vertices;
  F.resize(1, 3);
  F << 0, 1, 2;
  for (size_t i = 0; i < num_refinements; ++i) {
    Eigen::MatrixXd V_refined;
    Eigen::MatrixXi F_refined;
    refine_triangles(V, F, V_refined, F_refined);
    V = V_refined;
    F = F_refined;
  }
}

void
ConvexPolygon::sample(size_t num_samples,
                      std::vector<PlanarPoint>& domain_points) const
{
  domain_points.clear();

  // TODO Make actual bounding box
  PlanarPoint lower_left_corner(-1, -1);
  PlanarPoint upper_right_corner(1, 1);
  double x0 = lower_left_corner[0];
  double y0 = lower_left_corner[1];
  double x1 = upper_right_corner[0];
  double y1 = upper_right_corner[1];

  // Compute points
  VectorXr x_axis;
  VectorXr y_axis;
  generate_linspace(x0, x1, num_samples, x_axis);
  generate_linspace(y0, y1, num_samples, y_axis);
  for (size_t i = 0; i < num_samples; ++i) {
    for (size_t j = 0; j < num_samples; ++j) {
      PlanarPoint point(x_axis[i], y_axis[j]);
      if (contains(point)) {
        domain_points.push_back(point);
      }
    }
  }
}

std::array<Eigen::Matrix<double, 3, 1>, 3> const&
ConvexPolygon::get_boundary_segments() const
{
  return m_boundary_segments_coeffs;
}

Eigen::Matrix<double, 3, 2> const&
ConvexPolygon::get_vertices() const
{
  return m_vertices;
}

void
ConvexPolygon::parametrize_patch_boundaries(
  std::array<LineSegment, 3>& patch_boundaries) const
{
  size_t num_vertices = m_vertices.rows();
  for (size_t i = 0; i < num_vertices; ++i) {
    compute_parametric_line_between_points(
      m_vertices.row(i),
      m_vertices.row((i + 1) % num_vertices),
      patch_boundaries[i]);
  }
}

void
ConvexPolygon::intersect_patch_boundaries(
  const Eigen::Matrix<double, 3, 1>& first_boundary_segment_coeffs,
  const Eigen::Matrix<double, 3, 1>& second_boundary_segment_coeffs,
  PlanarPoint& intersection) const
{
  double a00 = first_boundary_segment_coeffs[0];
  double a10 = first_boundary_segment_coeffs[1];
  double a01 = first_boundary_segment_coeffs[2];
  double b00 = second_boundary_segment_coeffs[0];
  double b10 = second_boundary_segment_coeffs[1];
  double b01 = second_boundary_segment_coeffs[2];

  double x, y;
  // Solve for y in terms of x first
  if (!float_equal(a01, 0.0)) {
    double my = -a10 / a01;
    double by = -a00 / a01;
    assert(!float_equal(b10 + b01 * my, 0.0));
    x = -(b00 + b01 * by) / (b10 + b01 * my);
    y = my * x + by;
  }
  // Solve for x in terms of y first
  else if (!float_equal(a10, 0.0)) {
    double mx = -a01 / a10;
    double bx = -a00 / a10;
    assert(!float_equal(b01 + b10 * mx, 0.0));
    y = -(b00 + b10 * bx) / (b01 + b10 * mx);
    x = mx * y + bx;
  } else {
    spdlog::error("Degenerate line");
    return;
  }

  // Build intersection
  intersection << x, y;
}
