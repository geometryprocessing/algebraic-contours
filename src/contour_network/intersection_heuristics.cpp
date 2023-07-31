// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "intersection_heuristics.h"

#include "polynomial_function.h"

// Return true iff the two boxes defined by the extreme points are disjoint
bool
are_disjoint_bounding_boxes(const PlanarPoint& first_box_lower_left_point,
                            const PlanarPoint& first_box_upper_right_point,
                            const PlanarPoint& second_box_lower_left_point,
                            const PlanarPoint& second_box_upper_right_point)
{
  // Extract points
  double min_x_0 = first_box_lower_left_point[0];
  double min_y_0 = first_box_lower_left_point[1];
  double max_x_0 = first_box_upper_right_point[0];
  double max_y_0 = first_box_upper_right_point[1];
  double min_x_1 = second_box_lower_left_point[0];
  double min_y_1 = second_box_lower_left_point[1];
  double max_x_1 = second_box_upper_right_point[0];
  double max_y_1 = second_box_upper_right_point[1];

  // Check if boxes are disjoint
  if (min_x_0 > max_x_1)
    return true;
  if (min_x_1 > max_x_0)
    return true;
  if (min_y_0 > max_y_1)
    return true;
  if (min_y_1 > max_y_0)
    return true;

  // Otherwise, there is overlap
  return false;
}

// Check if the bounding box contains the curve over interval
// WARNING: May have false positives
bool
is_valid_bounding_box(const RationalFunction<4, 2>& planar_curve,
                      double t_min,
                      double t_max,
                      PlanarPoint& lower_left_point,
                      PlanarPoint& upper_right_point)
{
  double t_avg = (t_min + t_max) / 2.0;

  // Build test points
  PlanarPoint test_point_1 = planar_curve(t_min + 1e-6);
  PlanarPoint test_point_2 = planar_curve(t_avg);
  PlanarPoint test_point_3 = planar_curve(t_max - 1e-6);
  spdlog::trace("Testing points on curve at {}, {}, {}: {}, {}, {}",
                t_min,
                t_avg,
                t_max,
                test_point_1,
                test_point_2,
                test_point_3);

  // Check all points
  if (!is_in_bounding_box(test_point_1, lower_left_point, upper_right_point)) {
    spdlog::warn("Point {} at {} not in bounding box {}, {}",
                 test_point_1,
                 t_min,
                 lower_left_point,
                 upper_right_point);
    // return false;
  }
  if (!is_in_bounding_box(test_point_2, lower_left_point, upper_right_point)) {
    spdlog::warn("Point {} at {} not in bounding box {}, {}",
                 test_point_2,
                 t_avg,
                 lower_left_point,
                 upper_right_point);
    // return false;
  }
  if (!is_in_bounding_box(test_point_3, lower_left_point, upper_right_point)) {
    spdlog::warn("Point {} at {} not in bounding box {}, {}",
                 test_point_3,
                 t_max,
                 lower_left_point,
                 upper_right_point);
    // return false;
  }

  // Valid otherwise
  return true;
}

// Pad the bounding box by some epsilon
void
pad_bounding_box(PlanarPoint& lower_left_point,
                 PlanarPoint& upper_right_point,
                 double padding = 0)
{
  lower_left_point[0] -= padding;
  lower_left_point[1] -= padding;
  upper_right_point[0] += padding;
  upper_right_point[1] += padding;
}

// Compute bezier control points for a rational curve from a given change
// of coordinates matrix
void
compute_homogeneous_bezier_points_from_matrix(
  const RationalFunction<4, 2>& planar_curve,
  const Eigen::Matrix<double, 5, 5>& monomial_to_bezier_matrix,
  Eigen::Matrix<double, 5, 3>& bezier_points)
{
  // Get planar curve numerator polynomial coefficients in the monomial basis
  Eigen::Matrix<double, 5, 2> P_coeffs = planar_curve.get_numerators();
  Eigen::Matrix<double, 5, 1> Q_coeffs = planar_curve.get_denominator();
  Eigen::Matrix<double, 5, 1> x_coeffs, y_coeffs, w_coeffs;
  x_coeffs.setZero();
  x_coeffs.head(P_coeffs.rows()) = P_coeffs.col(0);
  y_coeffs.setZero();
  y_coeffs.head(P_coeffs.rows()) = P_coeffs.col(1);
  w_coeffs.setZero();
  w_coeffs.head(Q_coeffs.size()) = Q_coeffs;

  // Compute bezier homogeneous points
  // bezier_points.resize(5, 3);
  bezier_points.col(0) = monomial_to_bezier_matrix * x_coeffs;
  bezier_points.col(1) = monomial_to_bezier_matrix * y_coeffs;
  bezier_points.col(2) = monomial_to_bezier_matrix * w_coeffs;
}

// Compute the bezier points for a planar curve over the full domain [-1, 1]
void
compute_homogeneous_bezier_points(const RationalFunction<4, 2>& planar_curve,
                                  Eigen::Matrix<double, 5, 3>& bezier_points)
{
  spdlog::trace("Computing Bezier coefficients");

  // Compute matrix to go from monomial coefficients to Bezier coefficients
  Eigen::Matrix<double, 5, 5> monomial_to_bezier_matrix;
  monomial_to_bezier_matrix << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.0, -0.5,
    -1.0, 1.0, 0.0, -(1.0 / 3.0), 0.0, 1.0, 1.0, -0.5, 0.0, 0.5, -1.0, 1.0,
    -1.0, 1.0, -1.0, 1.0;

  // Compute bezier points from the matrix
  compute_homogeneous_bezier_points_from_matrix(
    planar_curve, monomial_to_bezier_matrix, bezier_points);
}

void
compute_homogeneous_bezier_points_over_interval(
  const RationalFunction<4, 2>& planar_curve,
  double t_min,
  double t_max,
  Eigen::Matrix<double, 5, 3>& bezier_points)
{
  spdlog::trace(
    "Computing Bezier coefficients for interval [{}, {}]", t_min, t_max);
  double r = t_max - t_min;

  // Compute matrix to go from monomial coefficients to Bezier coefficients
  Eigen::Matrix<double, 5, 5> monomial_to_bezier_matrix;

  // First row
  monomial_to_bezier_matrix(0, 0) = power(r + t_min, 0);
  monomial_to_bezier_matrix(0, 1) = power(r + t_min, 1);
  monomial_to_bezier_matrix(0, 2) = power(r + t_min, 2);
  monomial_to_bezier_matrix(0, 3) = power(r + t_min, 3);
  monomial_to_bezier_matrix(0, 4) = power(r + t_min, 4);

  // Second row
  monomial_to_bezier_matrix(1, 0) = 1.0;
  monomial_to_bezier_matrix(1, 1) = 0.75 * r + t_min;
  monomial_to_bezier_matrix(1, 2) = 0.5 * ((r + 2 * t_min) * (r + t_min));
  monomial_to_bezier_matrix(1, 3) =
    0.25 * ((r + 4 * t_min) * power(r + t_min, 2));
  monomial_to_bezier_matrix(1, 4) = t_min * power(r + t_min, 3);

  // Third row
  monomial_to_bezier_matrix(2, 0) = 1.0;
  monomial_to_bezier_matrix(2, 1) = 0.5 * r + t_min;
  monomial_to_bezier_matrix(2, 2) =
    power(r, 2) / 6.0 + r * t_min + power(t_min, 2);
  monomial_to_bezier_matrix(2, 3) =
    0.5 * (t_min * (r + 2 * t_min) * (r + t_min));
  monomial_to_bezier_matrix(2, 4) = power(t_min, 2) * power(r + t_min, 2);

  // Fourth row
  monomial_to_bezier_matrix(3, 0) = 1.0;
  monomial_to_bezier_matrix(3, 1) = 0.25 * r + t_min;
  monomial_to_bezier_matrix(3, 2) = 0.5 * (t_min * (r + 2 * t_min));
  monomial_to_bezier_matrix(3, 3) =
    0.25 * (power(t_min, 2) * (3 * r + 4 * t_min));
  monomial_to_bezier_matrix(3, 4) = power(t_min, 3) * (r + t_min);

  // Fifth row
  monomial_to_bezier_matrix(4, 0) = power(t_min, 0);
  monomial_to_bezier_matrix(4, 1) = power(t_min, 1);
  monomial_to_bezier_matrix(4, 2) = power(t_min, 2);
  monomial_to_bezier_matrix(4, 3) = power(t_min, 3);
  monomial_to_bezier_matrix(4, 4) = power(t_min, 4);

  // Compute bezier points from the matrix
  compute_homogeneous_bezier_points_from_matrix(
    planar_curve, monomial_to_bezier_matrix, bezier_points);
}

// Compute a bounding box for a planar curve over a domain using Bezier control
// points
void
compute_bezier_bounding_box_over_domain(
  const RationalFunction<4, 2>& planar_curve,
  double t_min,
  double t_max,
  PlanarPoint& lower_left_point,
  PlanarPoint& upper_right_point)
{
  spdlog::trace("Computing bezier bounding box for {} over [{}, {}]",
                planar_curve,
                t_min,
                t_max);

  // Convert to bezier coordinates
  Eigen::Matrix<double, 5, 3> bezier_points;
  compute_homogeneous_bezier_points_over_interval(
    planar_curve, t_min, t_max, bezier_points);

  // Normalize homogeneous coordinates
  Eigen::Matrix<double, 5, 1> bezier_x_coords;
  Eigen::Matrix<double, 5, 1> bezier_y_coords;
  for (Eigen::Index i = 0; i < bezier_points.rows(); ++i) {
    bezier_x_coords[i] = bezier_points(i, 0) / bezier_points(i, 2);
    bezier_y_coords[i] = bezier_points(i, 1) / bezier_points(i, 2);
  }

  // Bezier points should interpolate the endpoints
  assert(float_equal(planar_curve(t_min)[0], bezier_x_coords[4]));
  assert(float_equal(planar_curve(t_min)[1], bezier_y_coords[4]));
  assert(float_equal(planar_curve(t_max)[0], bezier_x_coords[0]));
  assert(float_equal(planar_curve(t_max)[1], bezier_y_coords[0]));

  // Get the max and min x values from the points
  double x_min = column_vector_min(bezier_x_coords);
  double x_max = column_vector_max(bezier_x_coords);

  // Get the max and min y values from the points
  double y_min = column_vector_min(bezier_y_coords);
  double y_max = column_vector_max(bezier_y_coords);

  // Build lower left point
  lower_left_point.resize(2);
  lower_left_point[0] = x_min;
  lower_left_point[1] = y_min;
  spdlog::trace("Lower left point: {}", lower_left_point);

  // Build upper right point
  upper_right_point.resize(2);
  upper_right_point[0] = x_max;
  upper_right_point[1] = y_max;
  spdlog::trace("Upper right point: {}", upper_right_point);

  // Check validity
  assert(is_valid_bounding_box(planar_curve, t_min, t_max, lower_left_point,
                               upper_right_point));
}

void
compute_bezier_bounding_box(const RationalFunction<4, 2>& planar_curve,
                            PlanarPoint& lower_left_point,
                            PlanarPoint& upper_right_point)
{
  // Get domain and domain midpoint
  double t_min = planar_curve.domain().get_lower_bound();
  double t_max = planar_curve.domain().get_upper_bound();
  double t_avg = (t_max + t_min) / 2.0;

  // Get first bezier bounding box
  PlanarPoint first_box_lower_left_point;
  PlanarPoint first_box_upper_right_point;
  compute_bezier_bounding_box_over_domain(planar_curve,
                                          t_min,
                                          t_avg,
                                          first_box_lower_left_point,
                                          first_box_upper_right_point);

  // Get second bezier bounding box
  PlanarPoint second_box_lower_left_point;
  PlanarPoint second_box_upper_right_point;
  compute_bezier_bounding_box_over_domain(planar_curve,
                                          t_avg,
                                          t_max,
                                          second_box_lower_left_point,
                                          second_box_upper_right_point);

  // Combine two bounding boxes
  combine_bounding_boxes(first_box_lower_left_point,
                         first_box_upper_right_point,
                         second_box_lower_left_point,
                         second_box_upper_right_point,
                         lower_left_point,
                         upper_right_point);

  // Inversely pad (that is, trim) the bounding box slightly to prevent
  // common endpoint intersections
  pad_bounding_box(lower_left_point, upper_right_point, -1e-10);

  // Check validity
  assert(is_valid_bounding_box(planar_curve, t_min, t_max, lower_left_point,
    upper_right_point));
}

void
compute_bezier_bounding_boxes(
  const std::vector<RationalFunction<4, 2>>& planar_curves,
  std::vector<PlanarPoint>& lower_left_points,
  std::vector<PlanarPoint>& upper_right_points)
{
  lower_left_points.resize(planar_curves.size());
  upper_right_points.resize(planar_curves.size());
  for (size_t i = 0; i < planar_curves.size(); ++i) {
    compute_bezier_bounding_box(
      planar_curves[i], lower_left_points[i], upper_right_points[i]);
  }
}

void
combine_bounding_boxes(const PlanarPoint& first_box_lower_left_point,
                       const PlanarPoint& first_box_upper_right_point,
                       const PlanarPoint& second_box_lower_left_point,
                       const PlanarPoint& second_box_upper_right_point,
                       PlanarPoint& lower_left_point,
                       PlanarPoint& upper_right_point)
{
  // Build point cloud from the input bounding box points
  Eigen::Matrix<double, 4, 2> points;
  points.row(0) = first_box_lower_left_point;
  points.row(1) = first_box_upper_right_point;
  points.row(2) = second_box_lower_left_point;
  points.row(3) = second_box_upper_right_point;

  // Get the bounding box of the point cloud
  compute_point_cloud_bounding_box(points, lower_left_point, upper_right_point);
}

bool
is_in_bounding_box(const PlanarPoint& test_point,
                   const PlanarPoint& lower_left_point,
                   const PlanarPoint& upper_right_point)
{
  // Equivalent to checking if the trivial bounding box at the test point
  // overlaps the bounding box
  return !are_disjoint_bounding_boxes(
    test_point, test_point, lower_left_point, upper_right_point);
}

// Check if two curves have disjoint Bezier bounding boxes
bool
are_disjoint_bezier_bounding_boxes(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve)
{
  // Compute first bezier bounding box
  PlanarPoint first_box_lower_left_point;
  PlanarPoint first_box_upper_right_point;
  compute_bezier_bounding_box(first_planar_curve,
                              first_box_lower_left_point,
                              first_box_upper_right_point);

  // Compute second bezier bounding box
  PlanarPoint second_box_lower_left_point;
  PlanarPoint second_box_upper_right_point;
  compute_bezier_bounding_box(second_planar_curve,
                              second_box_lower_left_point,
                              second_box_upper_right_point);

  // Compare bezier bounding boxes
  return are_disjoint_bounding_boxes(first_box_lower_left_point,
                                     first_box_upper_right_point,
                                     second_box_lower_left_point,
                                     second_box_upper_right_point);
}

// Check if two bounding boxes are disjoint
bool
are_disjoint_bezier_bounding_boxes(
  const std::pair<PlanarPoint, PlanarPoint>& first_bounding_box,
  const std::pair<PlanarPoint, PlanarPoint>& second_bounding_box)
{

  PlanarPoint first_box_lower_left_point = first_bounding_box.first;
  PlanarPoint first_box_upper_right_point = first_bounding_box.second;
  PlanarPoint second_box_lower_left_point = second_bounding_box.first;
  PlanarPoint second_box_upper_right_point = second_bounding_box.second;

  // Compare bezier bounding boxes
  return are_disjoint_bounding_boxes(first_box_lower_left_point,
                                     first_box_upper_right_point,
                                     second_box_lower_left_point,
                                     second_box_upper_right_point);
}

bool
are_nonintersecting_by_heuristic(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  IntersectionStats& intersection_stats)
{
  // Check bezier bounding boxes
  if (are_disjoint_bezier_bounding_boxes(first_planar_curve,
                                         second_planar_curve)) {
    spdlog::trace("Bezier bounding boxes do not overlap");
    intersection_stats.num_bezier_nonoverlaps++;
    return true;
  }

  // Otherwise, we cannot determine if there are intersections
  return false;
}

bool
are_nonintersecting_by_heuristic(
  const std::pair<PlanarPoint, PlanarPoint>& first_bounding_box,
  const std::pair<PlanarPoint, PlanarPoint>& second_bounding_box,
  IntersectionStats& intersection_stats)
{
  // Check bezier bounding boxes
  if (are_disjoint_bezier_bounding_boxes(first_bounding_box,
                                         second_bounding_box)) {
    spdlog::trace("Bezier bounding boxes do not overlap");
    intersection_stats.num_bezier_nonoverlaps++;
    return true;
  }

  // Otherwise, we cannot determine if there are intersections
  return false;
}

void
compute_bounding_box_hash_table(
  const std::vector<std::pair<PlanarPoint, PlanarPoint>>& bounding_boxes,
  std::vector<int> hash_table[50][50],
  std::vector<std::vector<int>>& reverse_hash_table)
{
  reverse_hash_table.clear();
  int num_interval = 50;
  int num_segments = bounding_boxes.size();
  if (num_segments == 0)
    return;
  for (int i = 0; i < num_segments; i++) {
    std::vector<int> a;
    reverse_hash_table.push_back(a);
  }

  double segments_bbox_x_min = bounding_boxes[0].first[0];
  double segments_bbox_y_min = bounding_boxes[0].first[1];
  double segments_bbox_x_max = bounding_boxes[0].second[0];
  double segments_bbox_y_max = bounding_boxes[0].second[1];

  for (int i = 1; i < num_segments; i++) {
    if (segments_bbox_x_min > bounding_boxes[i].first[0])
      segments_bbox_x_min = bounding_boxes[i].first[0];
    if (segments_bbox_y_min > bounding_boxes[i].first[1])
      segments_bbox_y_min = bounding_boxes[i].first[1];
    if (segments_bbox_x_max < bounding_boxes[i].second[0])
      segments_bbox_x_max = bounding_boxes[i].second[0];
    if (segments_bbox_y_max < bounding_boxes[i].second[1])
      segments_bbox_y_max = bounding_boxes[i].second[1];
  }

  double x_interval =
    (segments_bbox_x_max - segments_bbox_x_min) / num_interval;
  double y_interval =
    (segments_bbox_y_max - segments_bbox_y_min) / num_interval;

  double eps = PLANAR_BOUNDING_BOX_PRECISION;

  for (int i = 0; i < num_segments; i++) {
    int left_x =
      (bounding_boxes[i].first[0] - eps - segments_bbox_x_min) / x_interval;
    int right_x =
      num_interval -
      int((segments_bbox_x_max - eps - bounding_boxes[i].second[0]) /
          x_interval) -
      1;
    int left_y =
      (bounding_boxes[i].first[1] - eps - segments_bbox_y_min) / y_interval;
    int right_y =
      num_interval -
      int((segments_bbox_y_max - eps - bounding_boxes[i].second[1]) /
          y_interval) -
      1;

    for (int j = left_x; j <= right_x; j++) {
      for (int k = left_y; k <= right_y; k++) {
        hash_table[j][k].push_back(i);
        reverse_hash_table[i].push_back(j * num_interval + k);
      }
    }
  }
}

// helper for linear intersections
bool
ccw(PlanarPoint p0, PlanarPoint p1, PlanarPoint p2)
{
  return ((p2[1] - p0[1]) * (p1[0] - p0[0]) >
          (p1[1] - p0[1]) * (p2[0] - p0[0]));
}

bool
has_linear_intersection(const RationalFunction<4, 2>& first_planar_curve,
                        const RationalFunction<4, 2>& second_planar_curve)
{
  PlanarPoint l1_0 = first_planar_curve.start_point();
  PlanarPoint l1_1 = first_planar_curve.end_point();
  PlanarPoint l2_0 = second_planar_curve.start_point();
  PlanarPoint l2_1 = second_planar_curve.end_point();

  return ((ccw(l1_0, l2_0, l2_1) != ccw(l1_1, l2_0, l2_1)) &&
          (ccw(l1_0, l1_1, l2_0) != ccw(l1_0, l1_1, l2_1)));
}
