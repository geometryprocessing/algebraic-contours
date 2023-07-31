// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_intersections.h"

#include "polynomial_function.h"
#include <chrono>

// Map from the uniform domain [0, 1] to the planar curve domain
double
convert_spline_to_planar_curve_parameter(
  const RationalFunction<4, 2>& planar_curve,
  double t_spline,
  double epsilon = 0)
{
  double t_min = planar_curve.domain().get_lower_bound() + epsilon;
  double t_max = planar_curve.domain().get_upper_bound() - epsilon;
  return interval_lerp(0, 1, t_max, t_min, t_spline);
}

// Compute planar curve intersections with Bezier clipping
void
compute_bezier_clipping_planar_curve_intersections(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  std::vector<std::pair<double, double>>& intersection_points,
  const Eigen::Matrix<double, 5, 3>& first_bezier_control_points,
  const Eigen::Matrix<double, 5, 3>& second_bezier_control_points,
  double epsilon = 0)
{
  intersection_points.clear();

  // FIXME This is inefficient
  Point curve1[5], curve2[5];
  for (int i = 0; i < 5; i++) {
    curve1[i] = first_bezier_control_points.row(i);
    curve2[i] = second_bezier_control_points.row(i);
  }

  // FIXME
  if (!check_split_criteria(curve1)) {
    // std::cout << "potential self intersection p1!" << std::endl;
  }

  if (!check_split_criteria(curve2)) {
    // std::cout << "potential self intersection p2!" << std::endl;
  }

  std::vector<std::pair<double, double>> intersection_param_inkscope;
  std::vector<Point> P1(5);
  std::vector<Point> P2(5);
  for (int i = 0; i < 5; i++) {
    P1[i] = first_bezier_control_points.row(i);
    P2[i] = second_bezier_control_points.row(i);
  }

  // // std::cout << "in" << std::endl;
  find_intersections_bezier_clipping(
    intersection_param_inkscope,
    P1,
    P2,
    FIND_INTERSECTIONS_BEZIER_CLIPPING_PRECISION); // 1e-7
  // // std::cout << intersection_param_inkscope.size() << std::endl;

  for (size_t i = 0; i < intersection_param_inkscope.size(); i++) {
    double t_spline = intersection_param_inkscope[i].first;
    double s_spline = intersection_param_inkscope[i].second;
    double t = convert_spline_to_planar_curve_parameter(
      first_planar_curve, t_spline, epsilon);
    double s = convert_spline_to_planar_curve_parameter(
      second_planar_curve, s_spline, epsilon);
    intersection_points.push_back(std::make_pair(t, s));
  }
}

// Prune curve intersection points to the proper domains
void
prune_intersection_points(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  const std::vector<std::pair<double, double>>& intersection_points,
  std::vector<double>& first_curve_intersections,
  std::vector<double>& second_curve_intersections)
{
  for (size_t i = 0; i < intersection_points.size(); ++i) {
    double t = intersection_points[i].first;
    double s = intersection_points[i].second;

    // Trim points entirely out of domain of one of the two curves
    if (!first_planar_curve.is_in_domain_interior(t))
      continue;
    if (!second_planar_curve.is_in_domain_interior(s))
      continue;

    // Add points otherwise
    first_curve_intersections.push_back(t);
    second_curve_intersections.push_back(s);
  }
}

void
compute_planar_curve_intersections(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  const IntersectionParameters& intersect_params,
  std::vector<double>& first_curve_intersections,
  std::vector<double>& second_curve_intersections,
  IntersectionStats& intersection_stats,
  const std::pair<PlanarPoint, PlanarPoint>& first_bounding_box,
  const std::pair<PlanarPoint, PlanarPoint>& second_bounding_box,
  const Eigen::Matrix<double, 5, 3>& first_bezier_control_points,
  const Eigen::Matrix<double, 5, 3>& second_bezier_control_points)
{
  intersection_stats.num_intersection_tests++;
  auto t1 = std::chrono::high_resolution_clock::now();
  spdlog::trace("Finding intersections for {} and {}",
                first_planar_curve,
                second_planar_curve);

  intersection_stats.bounding_box_call++;

  if (intersect_params.use_heuristics &&
      are_nonintersecting_by_heuristic(
        first_bounding_box, second_bounding_box, intersection_stats)) {
    return;
  }

  intersection_stats.intersection_call++;

  // Compute intersection points by Bezier clipping
  std::vector<std::pair<double, double>> intersection_points;
  try {
    compute_bezier_clipping_planar_curve_intersections(
      first_planar_curve,
      second_planar_curve,
      intersection_points,
      first_bezier_control_points,
      second_bezier_control_points);
  } catch (std::runtime_error) {
    intersection_points.clear();
    spdlog::error("Failed to find intersection points");
  }

  // Prune the computed intersections to ensure they are in the correct domain
  prune_intersection_points(first_planar_curve,
                            second_planar_curve,
                            intersection_points,
                            first_curve_intersections,
                            second_curve_intersections);

  // Record the time spent finding intersections
  auto t2 = std::chrono::high_resolution_clock::now();
  double total_time =
    std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  spdlog::trace("Finding intersections took {} ms", total_time);
}

void
compute_planar_curve_intersections(
  const RationalFunction<4, 2>& first_planar_curve,
  const RationalFunction<4, 2>& second_planar_curve,
  const IntersectionParameters& intersect_params,
  std::vector<double>& first_curve_intersections,
  std::vector<double>& second_curve_intersections,
  IntersectionStats& intersection_stats)
{
  // Compute bounding boxes
  PlanarPoint lower_left_point, upper_right_point;
  std::pair<PlanarPoint, PlanarPoint> first_bounding_box;
  std::pair<PlanarPoint, PlanarPoint> second_bounding_box;
  compute_bezier_bounding_box(
    first_planar_curve, lower_left_point, upper_right_point);
  first_bounding_box = std::make_pair(lower_left_point, upper_right_point);
  compute_bezier_bounding_box(
    second_planar_curve, lower_left_point, upper_right_point);
  second_bounding_box = std::make_pair(lower_left_point, upper_right_point);

  // Compute Bezier points
  Eigen::Matrix<double, 5, 3> first_bezier_control_points;
  Eigen::Matrix<double, 5, 3> second_bezier_control_points;
  compute_homogeneous_bezier_points_over_interval(
    first_planar_curve,
    first_planar_curve.domain().get_lower_bound(),
    first_planar_curve.domain().get_upper_bound(),
    first_bezier_control_points);
  compute_homogeneous_bezier_points_over_interval(
    second_planar_curve,
    second_planar_curve.domain().get_lower_bound(),
    second_planar_curve.domain().get_upper_bound(),
    second_bezier_control_points);

  // Compute intersections with computed bounding boxes and Bezier points
  compute_planar_curve_intersections(first_planar_curve,
                                     second_planar_curve,
                                     intersect_params,
                                     first_curve_intersections,
                                     second_curve_intersections,
                                     intersection_stats,
                                     first_bounding_box,
                                     second_bounding_box,
                                     first_bezier_control_points,
                                     second_bezier_control_points);
}

void
compute_intersections(
  const std::vector<RationalFunction<4, 2>>& image_segments,
  const IntersectionParameters& intersect_params,
  std::vector<std::vector<double>>& intersections,
  std::vector<std::vector<size_t>>& intersection_indices,
  std::vector<std::vector<IntersectionData>>& contour_intersections,
  int& num_intersections,
  long long& intersection_call)
{
  intersections.clear();
  intersections.resize(image_segments.size());
  intersection_indices.clear();
  intersection_indices.resize(image_segments.size());

  // Setup intersection diagnostic tools
  IntersectionStats intersection_stats;

  // Compute all rational bezier control points
  std::vector<Eigen::Matrix<double, 5, 3>> image_segments_bezier_control_points;
  image_segments_bezier_control_points.reserve(image_segments.size());
  for (size_t i = 0; i < image_segments.size(); i++) {
    Eigen::Matrix<double, 5, 3> bezier_control_points;
    compute_homogeneous_bezier_points_over_interval(
      image_segments[i],
      image_segments[i].domain().get_lower_bound(),
      image_segments[i].domain().get_upper_bound(),
      bezier_control_points);
    image_segments_bezier_control_points.push_back(bezier_control_points);
  }

  // Compute all bounding boxes
  std::vector<std::pair<PlanarPoint, PlanarPoint>> image_segments_bounding_box;
  image_segments_bounding_box.reserve(image_segments.size());
  for (size_t i = 0; i < image_segments.size(); i++) {
    PlanarPoint lower_left_point, upper_right_point;
    compute_bezier_bounding_box(
      image_segments[i], lower_left_point, upper_right_point);
    image_segments_bounding_box.push_back(
      std::make_pair(lower_left_point, upper_right_point));
  }

  // Hash by uv
  int num_interval = 50;
  std::vector<int>
    hash_table[50][50]; // FIXME Make global: change both here and num_interval
  std::vector<std::vector<int>> reverse_hash_table;
  compute_bounding_box_hash_table(
    image_segments_bounding_box, hash_table, reverse_hash_table);

  // Compute intersections
  int num_segments = image_segments.size();
  for (int image_segment_index = 0; image_segment_index < num_segments;
       ++image_segment_index) {
    const std::vector<int>& cells = reverse_hash_table[image_segment_index];
    std::vector<bool> visited(image_segment_index, false);

    for (int cell : cells) {
      int j = cell / num_interval;
      int k = cell % num_interval;
      for (int i : hash_table[j][k]) {
        if (i >= image_segment_index || visited[i])
          continue;
        visited[i] = true;

        // Iterate over image segments with lower indices
        spdlog::trace("Computing segments {}, {} out of {}",
                      image_segment_index,
                      i,
                      image_segments.size());

        // Compute intersections between the two image segments
        std::vector<double> current_segment_intersections;
        std::vector<double> other_segment_intersections;
        compute_planar_curve_intersections(
          image_segments[image_segment_index],
          image_segments[i],
          intersect_params,
          current_segment_intersections,
          other_segment_intersections,
          intersection_stats,
          image_segments_bounding_box[image_segment_index],
          image_segments_bounding_box[i],
          image_segments_bezier_control_points[image_segment_index],
          image_segments_bezier_control_points[i]);
        append(intersections[image_segment_index],
               current_segment_intersections);
        append(intersections[i], other_segment_intersections);

        // Record the respective indices corresponding to the intersections
        for (size_t k = 0; k < other_segment_intersections.size(); ++k) {
          intersection_indices[image_segment_index].push_back(i);
          intersection_indices[i].push_back(image_segment_index);
        }
        // }

        // Build full intersection data
        for (size_t k = 0; k < other_segment_intersections.size(); ++k) {
          IntersectionData current_intersection_data;
          current_intersection_data.knot = current_segment_intersections[k];
          current_intersection_data.intersection_index = i;
          current_intersection_data.intersection_knot =
            other_segment_intersections[k];
          current_intersection_data.id = num_intersections;
          current_intersection_data.check_if_tip(
            image_segments[image_segment_index].domain(),
            intersect_params.trim_amount);
          current_intersection_data.check_if_base(
            image_segments[image_segment_index].domain(),
            intersect_params.trim_amount);
          contour_intersections[image_segment_index].push_back(
            current_intersection_data);

          // Build complementary boundary intersection data
          IntersectionData other_intersection_data;
          other_intersection_data.knot = other_segment_intersections[k];
          other_intersection_data.intersection_index = image_segment_index;
          other_intersection_data.intersection_knot =
            current_segment_intersections[k];
          other_intersection_data.id = num_intersections;
          other_intersection_data.check_if_tip(image_segments[i].domain(),
                                               intersect_params.trim_amount);
          other_intersection_data.check_if_base(image_segments[i].domain(),
                                                intersect_params.trim_amount);
          contour_intersections[i].push_back(other_intersection_data);
          num_intersections++;
        }
      }
    }
  }

  // Record intersection information
  intersection_call = intersection_stats.intersection_call;
  SPDLOG_INFO("Number of intersection tests: {}",
              intersection_stats.num_intersection_tests);
  SPDLOG_INFO("Number of nonoverlapping Bezier boxes: {}",
              intersection_stats.num_bezier_nonoverlaps);
}

void
split_planar_curves_no_self_intersection(
  const std::vector<RationalFunction<4, 2>>& planar_curves,
  std::vector<std::vector<double>>& split_points)
{
  split_points.resize(planar_curves.size());
  for (size_t i = 0; i < planar_curves.size(); i++) {
    if (float_equal_zero(planar_curves[i].domain().get_length(), 1e-6)) {
      spdlog::error("Splitting curve of length {}",
                    planar_curves[i].domain().get_length());
    }
    // Get Bezier points
    Eigen::Matrix<double, 5, 3> bezier_control_points;
    compute_homogeneous_bezier_points_over_interval(
      planar_curves[i],
      planar_curves[i].domain().get_lower_bound(),
      planar_curves[i].domain().get_upper_bound(),
      bezier_control_points);
    Point curve[5];
    for (int i = 0; i < 5; i++) {
      curve[i] = bezier_control_points.row(i);
    }

    // Get splits in the Bezier domain
    std::vector<double> split_points_bezier(0);
    split_bezier_curve_no_self_intersection(curve, 0, 1, split_points_bezier);

    // Get splits in the planar curve domain
    split_points[i].resize(split_points_bezier.size());
    for (size_t j = 0; j < split_points_bezier.size(); j++) {
      split_points[i][j] = convert_spline_to_planar_curve_parameter(
        planar_curves[i], split_points_bezier[j]);
    }
  }
}
