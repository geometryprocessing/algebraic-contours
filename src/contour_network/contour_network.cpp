// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "contour_network.h"

#include "compute_closed_contours.h"
#include "compute_cusps.h"
#include "igl/Timer.h"
#include "project_curves.h"
#include "split_contours.h"
#include "write_output.h"

// Build contour labels mapping label names to integer data per contour segment
void
build_contour_labels(
  const std::vector<size_t>& contour_patch_indices,
  const std::vector<bool>& contour_is_boundary,
  std::vector<std::map<std::string, int>>& contour_segment_labels)
{
  size_t num_segments = contour_patch_indices.size();
  contour_segment_labels.resize(num_segments);
  for (size_t i = 0; i < num_segments; ++i) {
    contour_segment_labels[i]["surface_patch"] = contour_patch_indices[i];
    contour_segment_labels[i]["is_boundary"] = contour_is_boundary[i] ? 0 : 1;
  }
}

ContourNetwork::ContourNetwork()
{
  clear();
}

ContourNetwork::ContourNetwork(
  const QuadraticSplineSurface& spline_surface,
  const IntersectionParameters& intersect_params,
  const InvisibilityParameters& invisibility_params,
  std::vector<std::pair<int, int>> patch_boundary_edges)
{
  init_contour_network(spline_surface,
                       intersect_params,
                       invisibility_params,
                       patch_boundary_edges);
}

void
ContourNetwork::clear()
{
  ProjectedCurveNetwork::clear();
}

void
ContourNetwork::add_contour_network_to_viewer() const
{
  polyscope::init();

  // Add curve network to viewer
  add_spatial_network_to_viewer();
}

void
ContourNetwork::view(const QuadraticSplineSurface& spline_surface) const
{
  polyscope::init();
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  spline_surface.add_surface_to_viewer(OFF_WHITE, DISCRETIZATION_LEVEL);
  add_contour_network_to_viewer();
  polyscope::show();
  polyscope::removeAllStructures();
}

void
ContourNetwork::view_contours() const
{
  polyscope::init();
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  add_spatial_network_to_viewer();
  polyscope::show();
  polyscope::removeAllStructures();
}

void
ContourNetwork::screenshot(const std::string& filename,
                           const QuadraticSplineSurface& spline_surface,
                           SpatialVector camera_position,
                           SpatialVector camera_target,
                           bool use_orthographic) const
{
  // Add the contour network to the surface
  polyscope::init();
  spline_surface.add_surface_to_viewer(OFF_WHITE, DISCRETIZATION_LEVEL);
  add_contour_network_to_viewer();

  // Build the cameras for the viewer
  glm::vec3 glm_camera_position = { camera_position[0],
                                    camera_position[1],
                                    camera_position[2] };
  glm::vec3 glm_camera_target = { camera_target[0],
                                  camera_target[1],
                                  camera_target[2] };

  // Set up the cameras
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  if (use_orthographic) {
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  }

  // Take the screenshot
  polyscope::screenshot(filename);
  SPDLOG_INFO("Screenshot saved to {}", filename);
  polyscope::removeAllStructures();
}

void
ContourNetwork::write_rasterized_contours(const std::string& filename,
                                          SpatialVector camera_position,
                                          SpatialVector camera_target) const
{
  polyscope::init();
  add_spatial_network_to_viewer();

  // Build the cameras for the viewer
  glm::vec3 glm_camera_position = { camera_position[0],
                                    camera_position[1],
                                    camera_position[2] };
  glm::vec3 glm_camera_target = { camera_target[0],
                                  camera_target[1],
                                  camera_target[2] };

  // Set up the cameras
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;

  // Take the screenshot
  polyscope::screenshot(filename);
  SPDLOG_INFO("Screenshot saved to {}", filename);
  polyscope::removeAllStructures();
}

// Initialize the contour network
void
ContourNetwork::init_contour_network(
  const QuadraticSplineSurface& spline_surface,
  const IntersectionParameters& intersect_params,
  const InvisibilityParameters& invisibility_params,
  std::vector<std::pair<int, int>> patch_boundary_edges)
{
  Matrix3x3r frame = Matrix3x3r::Identity();

  // Compute contours
  std::vector<RationalFunction<4, 3>> contour_segments;
  std::vector<Conic> contour_domain_curve_segments;
  std::vector<QuadraticSplineSurface::PatchIndex> contour_patch_indices;
  std::vector<bool> contour_is_boundary;

  igl::Timer timer;
  timer.start();
  std::vector<std::vector<IntersectionData>> contour_intersections;
  int num_intersections;
  compute_spline_surface_contours_and_boundaries(spline_surface,
                                                 frame,
                                                 patch_boundary_edges,
                                                 contour_domain_curve_segments,
                                                 contour_segments,
                                                 contour_patch_indices,
                                                 contour_is_boundary,
                                                 contour_intersections,
                                                 num_intersections);

  // Build contour labels for boundary contours and patch locations
  std::vector<std::map<std::string, int>> contour_segment_labels;
  build_contour_labels(
    contour_patch_indices, contour_is_boundary, contour_segment_labels);

  // Project contours to the plane
  std::vector<RationalFunction<4, 2>> planar_contour_segments;
  project_curves(contour_segments, frame, planar_contour_segments);

  // FIXME Split contours at self intersections
  // split_contours(
  //  contour_domain_curve_segments, contour_segments, planar_contour_segments,
  //  contour_patch_indices, contour_is_boundary, contour_segment_labels);

  // Chain the contour segments into closed contours
  std::vector<std::vector<int>> contours;
  std::vector<int> contour_labels;
  compute_closed_contours(contour_segments, contours, contour_labels);

  // Pad contour domains by an epsilon
  pad_contours(contour_domain_curve_segments,
               contour_segments,
               planar_contour_segments,
               invisibility_params.pad_amount);

  compute_contour_time = timer.getElapsedTime();
  timer.stop();

  // Get cusp points and intersections if needed
  size_t num_segments = contour_segments.size();
  std::vector<std::vector<double>> interior_cusps(num_segments);
  std::vector<std::vector<double>> boundary_cusps(num_segments);
  std::vector<bool> has_cusp_at_base(num_segments);
  std::vector<bool> has_cusp_at_tip(num_segments);
  timer.start();
  compute_spline_surface_cusps(
    spline_surface,
    contour_domain_curve_segments,
    contour_segments,
    contour_patch_indices,
    contours,
    interior_cusps,
    boundary_cusps,
    has_cusp_at_base,
    has_cusp_at_tip);

  compute_cusp_time = timer.getElapsedTime();
  timer.stop();
  SPDLOG_TRACE("Found {} interior cusps", nested_vector_size(interior_cusps));
  SPDLOG_TRACE("Found {} boundary cusps", nested_vector_size(boundary_cusps));

  // Compute planar contour intersections
  timer.start();
  std::vector<std::vector<double>> intersection_knots(num_segments);
  std::vector<std::vector<size_t>> intersection_indices(num_segments);
  compute_intersections(planar_contour_segments,
                        intersect_params,
                        intersection_knots,
                        intersection_indices,
                        contour_intersections,
                        num_intersections,
                        intersection_call);

  compute_intersection_time = timer.getElapsedTime();
  timer.stop();
  SPDLOG_TRACE("Found {} intersections",
               nested_vector_size(intersection_knots));

  // Optionally write contours before any graph construction
  if (invisibility_params.write_contour_soup) {
    CurveDiscretizationParameters curve_disc_params;
    Eigen::Vector2i viewport = Eigen::Vector2i(800, 800);
    svg::SVG svgWriter("contour_soup.svg", viewport);
    write_contours_with_annotations(frame,
                                    contour_segments,
                                    interior_cusps,
                                    boundary_cusps,
                                    intersection_knots,
                                    curve_disc_params,
                                    svgWriter);
  }

  // Build the curve network from the contours
  timer.start();
  init_projected_curve_network(contour_domain_curve_segments,
                               contour_segments,
                               planar_contour_segments,
                               contour_segment_labels,
                               contours,
                               contour_labels,
                               interior_cusps,
                               has_cusp_at_base,
                               intersection_knots,
                               intersection_indices,
                               contour_intersections,
                               num_intersections);
  compute_projected_time = timer.getElapsedTime();
  timer.stop();

  // Check the validity of the topological graph structure
#if CHECK_VALIDITY
  if (!is_valid_abstract_curve_network()) {
    spdlog::error("Invalid abstract curve network made");
    return;
  }

  // Check the validity of the geometric graph structure
  if (!is_valid_projected_curve_network()) {
    spdlog::error("Invalid projected curve network made");
    return;
  }
#endif

  // Compute the quantitative invisibility
  timer.start();
  compute_quantitative_invisibility(spline_surface, invisibility_params);
  compute_visibility_time = timer.getElapsedTime();
  timer.stop();
}

// *********************
// Direct QI Computation
// *********************

// Compute the ray mapping coeffs for a given sample point
void
ContourNetwork::generate_ray_mapping_coeffs(
  const SpatialVector& sample_point,
  Matrix2x3r& ray_mapping_coeffs) const
{
  SpatialVector view_direction(0, 0, 1);
  ray_mapping_coeffs.row(0) = sample_point - 20 * view_direction;
  ray_mapping_coeffs.row(1) = 40 * view_direction;
}

// Compute the QI from a given point from the list of intersection parameters
// with the surface
int
ContourNetwork::compute_quantitative_invisibility_from_ray_intersections(
  const Matrix2x3r& ray_mapping_coeffs,
  const SpatialVector& point,
  const std::vector<double>& ray_intersections) const
{
  // Partition intersections into points above and below the sample point
  if (ray_intersections.empty()) {
    return 0;
  } else {
    std::vector<double> ray_intersections_below;
    std::vector<double> ray_intersections_above;
    partition_ray_intersections(ray_mapping_coeffs,
                                point,
                                ray_intersections,
                                ray_intersections_below,
                                ray_intersections_above);

    // Set QI as the number of intersection points occluding the sample point
    return ray_intersections_below.size();
  }
}

// Compute the QI for a single segment with robust polling of three sample
// points
int
ContourNetwork::compute_segment_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex segment_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check segment validity
  if (!(is_valid_segment_index(segment_index))) {
    spdlog::debug("Attempting to propagate QI at invalid segment {}",
                  segment_index);
    return -1;
  }

  segment_number++;

  long long ray_int_call = 0;
  long long ray_bbox_call = 0;
  std::array<int, 3> qi_poll;
  RationalFunction<4, 3> const& spatial_curve =
    segment_spatial_curve(segment_index);

  std::array<double, 3> sample_parameters = { 0.5, 0.25, 0.75 };
  for (size_t i = 0; i < 3; ++i) {
    // Sample a point on the contour
    SpatialVector sample_point;
    spatial_curve.evaluate_normalized_coordinate(sample_parameters[i],
                                                 sample_point);
    spdlog::trace("Sample point: {}", sample_point);

    // Build ray mapping coefficients
    Matrix2x3r ray_mapping_coeffs;
    generate_ray_mapping_coeffs(sample_point, ray_mapping_coeffs);
    spdlog::trace("Ray mapping coefficients: {}", ray_mapping_coeffs);

    // Compute intersections of the ray with the surface
    std::vector<QuadraticSplineSurface::PatchIndex> patch_indices;
    std::vector<PlanarPoint> surface_intersections;
    std::vector<double> ray_intersections;

    patch_indices.reserve(10000);
    surface_intersections.reserve(10000);
    ray_intersections.reserve(10000);

    ray_number++;

    // Compute the intersections of the ray with the spline surface
    compute_spline_surface_ray_intersections(spline_surface,
                                             ray_mapping_coeffs,
                                             patch_indices,
                                             surface_intersections,
                                             ray_intersections,
                                             ray_int_call,
                                             ray_bbox_call);

    // Compute the QI from the ray intersections
    qi_poll[i] = compute_quantitative_invisibility_from_ray_intersections(
      ray_mapping_coeffs, sample_point, ray_intersections);

    // If no polling, return after first computation
    if (!invisibility_params.poll_segment_points) {
      return qi_poll[i];
    }
  }
  SPDLOG_TRACE(
    "QI poll values: {}, {}, {}", qi_poll[0], qi_poll[1], qi_poll[2]);

  ray_intersection_call += ray_int_call;
  ray_bounding_box_call += ray_bbox_call;

  // Poll for a majority
  if (qi_poll[0] == qi_poll[1])
    return qi_poll[0];
  if (qi_poll[1] == qi_poll[2])
    return qi_poll[1];
  if (qi_poll[2] == qi_poll[0])
    return qi_poll[2];

  // Arbitrarily choose the midpoint if no majority
  spdlog::warn("Could not compute consistent segment qi amongst {}, {}, {}",
               qi_poll[0],
               qi_poll[1],
               qi_poll[2]);
  return qi_poll[1];
}

// Compute the QI for a chain of segments with robust polling of three segments
int
ContourNetwork::compute_chain_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex start_segment_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check segment validity
  if (!(is_valid_segment_index(start_segment_index))) {
    spdlog::debug("Attempting to compute QI at invalid segment {}",
                  start_segment_index);
    return -1;
  }
  // Count number of segments
  SegmentChainIterator counter_iter =
    get_segment_chain_iterator(start_segment_index);
  size_t num_segments = 0;
  for (; !counter_iter.at_end_of_chain(); ++counter_iter) {
    ++num_segments;
  }

  chain_number++;

  // Determine three segments to sample
  SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
  int jump_size = 0;
  if (num_segments >= 5) {
    jump_size = num_segments / 5; // Warning: possibly slow

    // Skip forward one jump
    for (int j = 0; j < jump_size; ++j) {
      iter++;
    }
  }
  // Small number of segment case
  else if (num_segments >= 3) {
    jump_size = 1;
  }
  // Less than 3 case (just sample one segment)
  else {
    jump_size = 0;
  }

  // Sample at three segments
  std::array<int, 3> qi_poll;
  for (size_t i = 0; i < 3; ++i) {
    SegmentIndex segment_index = *iter;
    qi_poll[i] = compute_segment_quantitative_invisibility(
      spline_surface, segment_index, invisibility_params);

    // If no polling or trivial jump size, return immediately
    if ((!invisibility_params.poll_chain_segments) || (jump_size == 0)) {
      return qi_poll[i];
    }

    // Jump to next query point otherwise
    for (int j = 0; j < jump_size; ++j) {
      iter++;
    }
  }

  // Poll for majority
  if (qi_poll[0] == qi_poll[1])
    return qi_poll[0];
  if (qi_poll[1] == qi_poll[2])
    return qi_poll[1];
  if (qi_poll[2] == qi_poll[0])
    return qi_poll[2];

  // Choose the center segment if no majority found
  spdlog::warn("Could not compute consistent chain qi amongst {}, {}, {}",
               qi_poll[0],
               qi_poll[1],
               qi_poll[2]);
  return qi_poll[1];
}

// ****************
// Chain QI Methods
// ****************

// Propagate QI from the start segment forward to the next special node and
// return the final node index
ContourNetwork::NodeIndex
ContourNetwork::chain_quantitative_invisibility_forward(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex start_segment_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check segment validity
  if (!(is_valid_segment_index(start_segment_index))) {
    spdlog::debug("Attempting to propagate QI at invalid segment {}",
                  start_segment_index);
    return -1;
  }

  // Get QI of the start segment
  int quantitative_invisibility =
    get_segment_quantitative_invisibility(start_segment_index);
  if (quantitative_invisibility < 0) {
    spdlog::warn("Attempted to copy negative QI forward");
    return -1;
  }

  // Propagate the QI along the chain
  SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
  NodeIndex from_node_index = from(start_segment_index);
  NodeIndex to_node_index = to(start_segment_index);
  ++iter;
  for (; !iter.at_end_of_chain(); ++iter) {
    SegmentIndex segment_index = *iter;
    from_node_index = from(segment_index);
    to_node_index = to(segment_index);
    set_segment_quantitative_invisibility(segment_index,
                                          quantitative_invisibility);
    set_node_quantitative_invisibility(from_node_index,
                                       quantitative_invisibility);

    // Check against direct computation
    if (invisibility_params.check_chaining) {
      SPDLOG_DEBUG("Checking QI chaining at knot node");
      check_quantitative_invisibility_propagation(
        spline_surface, segment_index, invisibility_params);
    }
  }

  return to_node_index;
}

// Propagate QI from the start segment backward to the next special node and
// return the final node index
ContourNetwork::NodeIndex
ContourNetwork::chain_quantitative_invisibility_backward(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex start_segment_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check segment validity
  if (!(is_valid_segment_index(start_segment_index))) {
    spdlog::debug("Attempting to propagate QI at invalid segment {}",
                  start_segment_index);
    return -1;
  }

  // Get QI of the start segment
  int quantitative_invisibility =
    get_segment_quantitative_invisibility(start_segment_index);
  if (quantitative_invisibility < 0) {
    spdlog::debug("Attempted to propagate negative QI backward");
    return -1;
  }

  // Propagate the QI along the chain
  SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
  NodeIndex from_node_index = from(start_segment_index);
  NodeIndex to_node_index = to(start_segment_index);
  --iter;
  for (; !iter.at_reverse_end_of_chain(); --iter) {
    SegmentIndex segment_index = *iter;
    from_node_index = from(segment_index);
    to_node_index = to(segment_index);
    set_segment_quantitative_invisibility(segment_index,
                                          quantitative_invisibility);
    set_node_quantitative_invisibility(to_node_index,
                                       quantitative_invisibility);

    // Check against direct computation
    if (invisibility_params.check_chaining) {
      SPDLOG_DEBUG("Checking QI chaining at knot node");
      check_quantitative_invisibility_propagation(
        spline_surface, segment_index, invisibility_params);
    }
  }

  return from_node_index;
}

// ****************************
// Local QI Propagation Methods
// ****************************

// Chain QI from the start segment forward to the next special node and continue
// propagation
void
ContourNetwork::propagate_quantitative_invisibility_forward(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex start_segment_index,
  const InvisibilityParameters& invisibility_params)
{
  if (!(is_valid_segment_index(start_segment_index))) {
    spdlog::debug("Attempting to propagate QI at invalid segment {}",
                  start_segment_index);
    return;
  }

  NodeIndex chain_end_node_index = chain_quantitative_invisibility_forward(
    spline_surface, start_segment_index, invisibility_params);
  propagate_quantitative_invisibility_at_node(
    spline_surface, chain_end_node_index, invisibility_params);
}

// Chain QI from the start segment backward to the next special node and
// continue propagation
void
ContourNetwork::propagate_quantitative_invisibility_backward(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex start_segment_index,
  const InvisibilityParameters& invisibility_params)
{
  if (!(is_valid_segment_index(start_segment_index))) {
    spdlog::debug("Attempting to propagate QI at invalid segment {}",
                  start_segment_index);
    return;
  }

  NodeIndex chain_start_node_index = chain_quantitative_invisibility_backward(
    spline_surface, start_segment_index, invisibility_params);
  propagate_quantitative_invisibility_at_node(
    spline_surface, chain_start_node_index, invisibility_params);
}

// Propagate QI at a node with case work for different node types
void
ContourNetwork::propagate_quantitative_invisibility_at_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index,
  const InvisibilityParameters& invisibility_params)
{
  spdlog::debug("Propagating QI at node {}", node_index);

  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return;
  }

  // Stop processing the node if node QI already set, and mark the node as
  // processed otherwise
  if (node_quantitative_invisibility_is_set(node_index))
    return;
  mark_node_quantitative_invisibility_as_set(node_index);

  // Attempt to propagate QI at the node
  bool success = true;
  if (is_intersection_node(node_index)) {
    SPDLOG_INFO("Propagating QI at intersection node");
    success = propagate_quantitative_invisibility_at_intersection_node(
      spline_surface, node_index, invisibility_params);
  } else if (is_interior_cusp_node(node_index)) {
    SPDLOG_INFO("Propagating QI at interior cusp node");
    success = propagate_quantitative_invisibility_at_cusp_node(
      spline_surface, node_index, invisibility_params);
  } else if (is_boundary_cusp_node(node_index)) {
    SPDLOG_INFO("Propagating QI at boundary cusp node");
    propagate_quantitative_invisibility_at_boundary_cusp_node(
      spline_surface, node_index, invisibility_params);
  } else if (is_marked_knot_node(node_index)) {
    SPDLOG_INFO("Propagating QI at marked node");
    propagate_quantitative_invisibility_at_marked_knot_node(
      spline_surface, node_index, invisibility_params);
  }

  // Fallback to direct computation and chaining if not successful
  if (!success) {
    int qi_out = compute_chain_quantitative_invisibility(
      spline_surface, out(node_index), invisibility_params);
    int qi_in = compute_chain_quantitative_invisibility(
      spline_surface, in(node_index), invisibility_params);
    set_segment_quantitative_invisibility(out(node_index), qi_out);
    set_segment_quantitative_invisibility(in(node_index), qi_in);
    set_node_quantitative_invisibility(node_index, qi_out);
    propagate_quantitative_invisibility_forward(
      spline_surface, out(node_index), invisibility_params);
    propagate_quantitative_invisibility_backward(
      spline_surface, in(node_index), invisibility_params);
  }
}

// Propagate QI if it is known for the in or out segment at a node to the other
// node and return the QI of the out node
int
ContourNetwork::propagate_quantitative_invisibility_across_node(
  NodeIndex node_index,
  int change_in_quantitative_invisibility)
{
  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return -1;
  }

  // Check if boundary node (no work can be done)
  SegmentIndex in_segment = in(node_index);
  SegmentIndex out_segment = out(node_index);
  if (!is_valid_segment_index(in_segment) ||
      !is_valid_segment_index(out_segment)) {
    spdlog::debug("Cannot propagate QI across terminal boundary node");
    return -1;
  }
  int qi_in = get_segment_quantitative_invisibility(in_segment);
  int qi_out = get_segment_quantitative_invisibility(out_segment);

  // Do nothing if both segments already have QI
  if ((qi_in >= 0) && (qi_out >= 0)) {
    int new_qi = qi_in + change_in_quantitative_invisibility;
    if (qi_out != new_qi) {
      spdlog::warn(
        "Inconsistent QI propagation of {} = {} + {} to known value {}",
        new_qi,
        qi_in,
        change_in_quantitative_invisibility,
        qi_out);
      return -1;
    }
    return qi_out;
  }
  // Propagate to the out segment if the in segment has known QI
  else if (qi_in >= 0) {
    int new_qi = qi_in + change_in_quantitative_invisibility;
    if (new_qi < 0) {
      spdlog::warn("Attempted to propagate negative QI across node");
      return -1;
    }
    set_segment_quantitative_invisibility(out_segment, new_qi);
    set_node_quantitative_invisibility(node_index, new_qi);
    return new_qi;
  }
  // Propagate to the out segment if the in segment has known QI
  else if (qi_out >= 0) {
    int new_qi = qi_out - change_in_quantitative_invisibility;
    if (new_qi < 0) {
      spdlog::warn("Attempted to propagate negative QI across node");
      return -1;
    }

    set_segment_quantitative_invisibility(in_segment, new_qi);
    set_node_quantitative_invisibility(node_index, qi_out);
    return qi_out;
  }
  // Cannot propagate if both are unassigned
  else {
    spdlog::warn(
      "Cannot propagate QI across node if both segments are unassigned");
    return -1;
  }
}

// Propagate the QI at an intersection node
bool
ContourNetwork::propagate_quantitative_invisibility_at_intersection_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index,
  const InvisibilityParameters& invisibility_params)
{
  spdlog::debug("Propagating QI at intersection node {}", node_index);
  Eigen::Matrix<double, 3, 1> tau(0, 0, 1);

  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return false;
  }

  // Get the intersecting node and check it is valid
  NodeIndex intersecting_node_index = intersection(node_index);
  if (!(is_valid_node_index(intersecting_node_index))) {
    spdlog::error(
      "Attempting to propagate QI with intersection rule at invalid node {}",
      intersecting_node_index);
    return false;
  }

  // Handle boundary case separately
  if (is_tnode(node_index)) {
    return propagate_quantitative_invisibility_at_boundary_intersection_node(
      spline_surface, node_index, invisibility_params);
  }

  // Determine if the node is above the intersection node
  SpatialVector node_point = node_spatial_point(node_index);
  SpatialVector intersecting_node_point =
    node_spatial_point(intersecting_node_index);
  double node_tau_projection = node_point * tau;
  double intersecting_node_tau_projection = intersecting_node_point * tau;
  bool is_above = (node_tau_projection < intersecting_node_tau_projection);
  spdlog::debug("Node has view direction projection {}", node_tau_projection);
  spdlog::debug("Intersecting node has view direction projection {}",
                intersecting_node_tau_projection);

  // Determine if the contour crosses the intersection with positive
  // orientation, meaning the intersecting tangent is ccw from the node tangent.
  // This means the intersecting tangent is outward normal to the projected
  // surface along the node contour.
  PlanarPoint node_tangent = node_planar_tangent(node_index);
  PlanarPoint intersecting_node_tangent =
    node_planar_tangent(intersecting_node_index);
  PlanarPoint node_tangent_normal(2);
  node_tangent_normal[0] = -node_tangent[1];
  node_tangent_normal[1] = node_tangent[0];
  double orientation_ind =
    dot_product<double, 2>(node_tangent_normal, intersecting_node_tangent);
  bool is_positively_orientated = (orientation_ind > 0);
  SPDLOG_DEBUG("Node tangent {} with planar projection {}",
               node_spatial_tangent(node_index),
               node_tangent);
  SPDLOG_DEBUG("Intersecting node tangent {} with planar projection {}",
               node_spatial_tangent(intersecting_node_index),
               intersecting_node_tangent);

  // Propagate QI at the intersecting according to casework
  int new_qi;
  // Above with positive orientation
  if ((is_above) && (is_positively_orientated)) {
    spdlog::debug(
      "Input node is above, and the intersection is positively oriented");
    new_qi = propagate_quantitative_invisibility_across_node(node_index, 0);
    if (new_qi < 0) {
      spdlog::error("Failed to propagate QI across intersection node");
      return false;
    }
    set_segment_quantitative_invisibility(out(intersecting_node_index), new_qi);
    set_segment_quantitative_invisibility(in(intersecting_node_index),
                                          new_qi + 2);
    set_node_quantitative_invisibility(intersecting_node_index, new_qi);
  }
  // Above with negative orientation
  else if ((is_above) && (!is_positively_orientated)) {
    spdlog::debug(
      "Input node is above, and the intersection is negatively oriented");
    new_qi = propagate_quantitative_invisibility_across_node(node_index, 0);
    if (new_qi < 0) {
      spdlog::error("Failed to propagate QI across intersection node");
      return false;
    }
    set_segment_quantitative_invisibility(out(intersecting_node_index),
                                          new_qi + 2);
    set_segment_quantitative_invisibility(in(intersecting_node_index), new_qi);
    set_node_quantitative_invisibility(intersecting_node_index, new_qi + 2);
  }
  // Below with positive orientation
  else if ((!is_above) && (is_positively_orientated)) {
    spdlog::debug(
      "Input node is below, and the intersection is positively oriented");
    new_qi = propagate_quantitative_invisibility_across_node(node_index, 2);
    if (new_qi < 0) {
      spdlog::error("Failed to propagate QI across intersection node");
      return false;
    }
    set_segment_quantitative_invisibility(out(intersecting_node_index),
                                          new_qi - 2);
    set_segment_quantitative_invisibility(in(intersecting_node_index),
                                          new_qi - 2);
    set_node_quantitative_invisibility(intersecting_node_index, new_qi - 2);
  }
  // Below with negative orientation
  else if ((!is_above) && (!is_positively_orientated)) {
    spdlog::debug(
      "Input node is below, and the intersection is negatively oriented");
    new_qi = propagate_quantitative_invisibility_across_node(node_index, -2);
    if (new_qi < 0) {
      spdlog::error("Failed to propagate QI across intersection node");
      return false;
    }
    set_segment_quantitative_invisibility(out(intersecting_node_index), new_qi);
    set_segment_quantitative_invisibility(in(intersecting_node_index), new_qi);
    set_node_quantitative_invisibility(intersecting_node_index, new_qi);
  }

  // Get the QI for adjacent segments
  int qi_in = get_segment_quantitative_invisibility(in(node_index));
  int qi_out = get_segment_quantitative_invisibility(out(node_index));
  int qi_intersecting_in =
    get_segment_quantitative_invisibility(in(intersecting_node_index));
  int qi_intersecting_out =
    get_segment_quantitative_invisibility(out(intersecting_node_index));
  spdlog::debug("Propagated QI across intersection node {} -> {} with "
                "intersecting node {} -> {}",
                qi_in,
                qi_out,
                qi_intersecting_in,
                qi_intersecting_out);

  // Check the QI against direct computation
  if (invisibility_params.check_propagation) {
    SPDLOG_DEBUG("Checking QI propagation at intersection node");
    check_quantitative_invisibility_propagation(
      spline_surface, in(node_index), invisibility_params);
    check_quantitative_invisibility_propagation(
      spline_surface, out(node_index), invisibility_params);
    check_quantitative_invisibility_propagation(
      spline_surface, in(intersecting_node_index), invisibility_params);
    check_quantitative_invisibility_propagation(
      spline_surface, out(intersecting_node_index), invisibility_params);
  }

  // View the intersection locally
  if (invisibility_params.view_intersections) {
    view_intersection_node(spline_surface, node_index);
  }

  // Propagate all nodes
  propagate_quantitative_invisibility_forward(
    spline_surface, out(node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(node_index), invisibility_params);
  propagate_quantitative_invisibility_forward(
    spline_surface, out(intersecting_node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(intersecting_node_index), invisibility_params);

  return true;
}

// Propagate the QI at a boundary intersection node
bool
ContourNetwork::
  propagate_quantitative_invisibility_at_boundary_intersection_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params)
{
  spdlog::debug("Propagating QI at boundary intersection node {}", node_index);

  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return false;
  }

  // Get the intersecting node and check it is valid
  NodeIndex intersecting_node_index = intersection(node_index);
  if (!(is_valid_node_index(intersecting_node_index))) {
    spdlog::error(
      "Attempting to propagate QI with intersection rule at invalid node {}",
      intersecting_node_index);
    return false;
  }

  // Node out segment QI
  int qi_out = compute_chain_quantitative_invisibility(
    spline_surface, out(node_index), invisibility_params);
  if (qi_out >= 0) {
    set_segment_quantitative_invisibility(out(node_index), qi_out);
  }

  // Node in segment QI
  int qi_in = compute_chain_quantitative_invisibility(
    spline_surface, in(node_index), invisibility_params);
  if (qi_in >= 0) {
    set_segment_quantitative_invisibility(in(node_index), qi_in);
  }

  // Node QI
  set_node_quantitative_invisibility(node_index, std::max<int>(qi_in, qi_out));

  // Intersection node out segment
  int qi_intersection_out = compute_chain_quantitative_invisibility(
    spline_surface, out(intersecting_node_index), invisibility_params);
  if (qi_intersection_out >= 0) {
    set_segment_quantitative_invisibility(out(intersecting_node_index),
                                          qi_intersection_out);
  }

  // Intersection node in segment QI
  int qi_intersection_in = compute_chain_quantitative_invisibility(
    spline_surface, in(intersecting_node_index), invisibility_params);
  if (qi_intersection_in >= 0) {
    set_segment_quantitative_invisibility(in(intersecting_node_index),
                                          qi_intersection_in);
  }

  // Intersection node QI
  set_node_quantitative_invisibility(
    intersecting_node_index,
    std::max<int>(qi_intersection_in, qi_intersection_out));

  // Propagate QI directly at the boundary intersection. For the segment that
  // does not exist, nothing is done in these function calls
  spdlog::debug("Propagated QI across boundary intersection node {} -> {} with "
                "intersecting node {} -> {}",
                qi_in,
                qi_out,
                qi_intersection_in,
                qi_intersection_out);
  propagate_quantitative_invisibility_forward(
    spline_surface, out(node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(node_index), invisibility_params);
  propagate_quantitative_invisibility_forward(
    spline_surface, out(intersecting_node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(intersecting_node_index), invisibility_params);

  return true;
}

// Propagate QI at a cusp
bool
ContourNetwork::propagate_quantitative_invisibility_at_cusp_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return false;
  }

  // Determine if the tangent is pointing toward or away from the camera
  Eigen::Matrix<double, 3, 1> tau(0, 0, 1);
  SpatialVector node_tangent = node_spatial_tangent(node_index);
  double node_tangent_tau_projection = node_tangent * tau;
  bool is_reversed = (node_tangent_tau_projection < 0);
  SPDLOG_TRACE("Tangent at cusp node", node_tangent.transpose());
  SPDLOG_TRACE("In segment midpoint is {} and out midpoint segment is {}",
               segment_spatial_curve(in(node_index)).mid_point(),
               segment_spatial_curve(out(node_index)).mid_point());

  // Propagate the QI
  int new_qi;
  if (is_reversed) {
    new_qi = propagate_quantitative_invisibility_across_node(node_index, -1);
  } else {
    new_qi = propagate_quantitative_invisibility_across_node(node_index, 1);
  }
  if (new_qi < 0) {
    spdlog::error("Failed to propagate QI across cusp node");
    return false;
  }

  int qi_in = get_segment_quantitative_invisibility(in(node_index));
  int qi_out = get_segment_quantitative_invisibility(out(node_index));
  SPDLOG_DEBUG("Propagated QI across cusp node {} -> {}", qi_in, qi_out);

  // Check the QI against direct computation
  if (invisibility_params.check_propagation) {
    SPDLOG_DEBUG("Checking QI propagation at cusp node");
    check_quantitative_invisibility_propagation(
      spline_surface, in(node_index), invisibility_params);
    check_quantitative_invisibility_propagation(
      spline_surface, out(node_index), invisibility_params);
  }

  // View the cusp locally
  if (invisibility_params.view_cusps) {
    view_cusp_node(spline_surface, node_index);
  }

  // Propagate the QI
  propagate_quantitative_invisibility_forward(
    spline_surface, out(node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(node_index), invisibility_params);

  return true;
}

// Propagate QI at a boundary cusp node by just computing QI directly and
// propagating
void
ContourNetwork::propagate_quantitative_invisibility_at_boundary_cusp_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index,
  const InvisibilityParameters& invisibility_params)
{
  // Check node validity
  if (!(is_valid_node_index(node_index))) {
    spdlog::debug("Attempting to propagate QI at invalid node {}", node_index);
    return;
  }

  int qi_out = compute_chain_quantitative_invisibility(
    spline_surface, out(node_index), invisibility_params);
  int qi_in = compute_chain_quantitative_invisibility(
    spline_surface, in(node_index), invisibility_params);
  set_segment_quantitative_invisibility(out(node_index), qi_out);
  set_segment_quantitative_invisibility(in(node_index), qi_in);
  set_node_quantitative_invisibility(node_index, qi_out);
  propagate_quantitative_invisibility_forward(
    spline_surface, out(node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(node_index), invisibility_params);
  return;
}

// Propagate QI at a marked knot node by propagating known QI across it without
// change
void
ContourNetwork::propagate_quantitative_invisibility_at_marked_knot_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index,
  const InvisibilityParameters& invisibility_params)
{
  int new_qi = propagate_quantitative_invisibility_across_node(node_index, 0);
  if (new_qi < 0) {
    spdlog::error("Failed to propagate QI across marked node");
    return;
  }
  propagate_quantitative_invisibility_forward(
    spline_surface, out(node_index), invisibility_params);
  propagate_quantitative_invisibility_backward(
    spline_surface, in(node_index), invisibility_params);
}

// Check QI propagation at segment with direct computation
void
ContourNetwork::check_quantitative_invisibility_propagation(
  const QuadraticSplineSurface& spline_surface,
  SegmentIndex segment_index,
  const InvisibilityParameters& invisibility_params)
{
  int qi = get_segment_quantitative_invisibility(segment_index);
  int direct_qi = compute_segment_quantitative_invisibility(
    spline_surface, segment_index, invisibility_params);
  if (qi != direct_qi) {
    spdlog::warn(
      "Propagated QI {} differs from direct computation {}", qi, direct_qi);
    // set_segment_quantitative_invisibility(segment_index, direct_qi);
  } else {
    SPDLOG_DEBUG("Propagated QI {} correctly", qi);
  }
}

// Set QI to 0 for all segments
void
ContourNetwork::compute_default_quantitative_invisibility()
{
  // Set QI to 0 for all segments independently
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    set_segment_quantitative_invisibility(i, 0);
  }
}

// Compute the quantitative invisibility for each segment of the contour
// network.
// Note: This method is extremely slow and should only be used for
// validation on small examples.
void
ContourNetwork::compute_direct_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  const InvisibilityParameters& invisibility_params)
{
  // Compute QI for all segments independently
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    int qi = compute_segment_quantitative_invisibility(
      spline_surface, i, invisibility_params);
    set_segment_quantitative_invisibility(i, qi);
  }
}

// Compute the quantitative invisibility for each chain of the contour network
// and propagate it to all segments
void
ContourNetwork::compute_chained_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  const InvisibilityParameters& invisibility_params)
{
  // Compute QI for all segment chains independently
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    // Compute QI for the first segment
    NodeIndex start_node_index = chain_start_nodes[i];
    SegmentIndex start_segment_index = out(start_node_index);
    int qi_start = compute_chain_quantitative_invisibility(
      spline_surface, start_segment_index, invisibility_params);
    set_segment_quantitative_invisibility(start_segment_index, qi_start);

    // Propagate the QI to the other chain segments
    chain_quantitative_invisibility_forward(
      spline_surface, start_segment_index, invisibility_params);
  }
}

// Compute the quantitative invisibility by computing it per chain and propagating
// across nodes with local rules when possible
void
ContourNetwork::compute_propagated_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  const InvisibilityParameters& invisibility_params)
{
  // Compute QI for all segment chains independently
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();

  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    NodeIndex start_node_index = chain_start_nodes[i];
    SegmentIndex start_segment_index = out(start_node_index);
    SPDLOG_TRACE("Attempting to continue QI propagation start node {}",
                 start_node_index);

    // Compute QI for the first segment and propagate it if it doesn't already
    // exist
    if (!node_quantitative_invisibility_is_set(start_node_index)) {
      int qi_start = compute_chain_quantitative_invisibility(
        spline_surface, start_segment_index, invisibility_params);
      set_node_quantitative_invisibility(start_node_index, qi_start);
      set_segment_quantitative_invisibility(start_segment_index, qi_start);
      SPDLOG_INFO("Propagating QI of {} from start node {}",
                  get_segment_quantitative_invisibility(start_segment_index),
                  start_segment_index);
      propagate_quantitative_invisibility_forward(
        spline_surface, start_segment_index, invisibility_params);
    }
  }
}

// Compute the quantitative invisibility with a method of choice
void
ContourNetwork::compute_quantitative_invisibility(
  const QuadraticSplineSurface& spline_surface,
  const InvisibilityParameters& invisibility_params)
{

  // compute hash table here
  // Perform chosen invisibility method
  // Optionally skip QI computation
  if (invisibility_params.invisibility_method == InvisibilityMethod::none) {
    compute_default_quantitative_invisibility();
  }
  // Compute QI per segment
  else if (invisibility_params.invisibility_method ==
           InvisibilityMethod::direct) {
    igl::Timer timer;
    SPDLOG_INFO("start compute direct QI");
    timer.start();
    compute_direct_quantitative_invisibility(spline_surface,
                                             invisibility_params);
    SPDLOG_INFO("Direct QI compute Time cost: {}s", timer.getElapsedTime());
  }
  // Compute QI per chain
  else if (invisibility_params.invisibility_method ==
           InvisibilityMethod::chaining) {
    igl::Timer timer;
    SPDLOG_INFO("start compute chaining QI");
    timer.start();
    compute_chained_quantitative_invisibility(spline_surface,
                                              invisibility_params);
    SPDLOG_INFO("Chaining QI compute Time cost: {}s", timer.getElapsedTime());
  }
  // Compute full QI propagation
  else if (invisibility_params.invisibility_method ==
           InvisibilityMethod::propagation) {
    igl::Timer timer;
    SPDLOG_INFO("start compute propagation QI");
    timer.start();
    compute_propagated_quantitative_invisibility(spline_surface,
                                                 invisibility_params);
    SPDLOG_INFO("Propagation QI compute Time cost: {}s",
                timer.getElapsedTime());
  }

  // Check for errors in the QI
  std::vector<int> quantitative_invisibility;
  enumerate_quantitative_invisibility(quantitative_invisibility);
  SPDLOG_TRACE(formatted_vector(quantitative_invisibility, ", "));
  if (vector_contains(quantitative_invisibility, -1)) {
    spdlog::error("Negative QI present in final values");
  }
}

// View all local features in the contour network
// Warning: Should only be used for small meshes
void
ContourNetwork::view_local_features(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<size_t>& patch_indices,
  const std::vector<SegmentIndex>& segment_indices,
  const std::vector<NodeIndex>& node_indices) const
{

  // View patches
  for (size_t i = 0; i < patch_indices.size(); ++i) {
    size_t patch_index = patch_indices[i];
    std::string patch_name = "local_patch_" + std::to_string(i);
    spline_surface.get_patch(patch_index).add_patch_to_viewer(patch_name);
  }

  // View segments
  for (size_t i = 0; i < segment_indices.size(); ++i) {
    size_t segment_index = segment_indices[i];
    std::string segment_name = "local_segment_" + std::to_string(i);
    segment_spatial_curve(segment_index).add_curve_to_viewer(segment_name);
  }

  // View Nodes
  MatrixXr node_points(node_indices.size(), 3);
  MatrixXr node_in_tangents(node_indices.size(), 3);
  MatrixXr node_out_tangents(node_indices.size(), 3);
  for (size_t i = 0; i < node_indices.size(); ++i) {
    size_t node_index = node_indices[i];
    node_points.row(i) = node_spatial_point(node_index);
    node_in_tangents.row(i) = node_spatial_in_tangent(node_index);
    node_out_tangents.row(i) = node_spatial_out_tangent(node_index);
  }
  polyscope::registerPointCloud("node_points", node_points);
  polyscope::getPointCloud("node_points")
    ->addVectorQuantity("node_in_tangents", node_in_tangents);
  polyscope::getPointCloud("node_points")
    ->addVectorQuantity("node_out_tangents", node_out_tangents);

  // Show viewer with orthographic projection
  polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;

  auto viewer_callback = [&]() {
    if (ImGui::Button("Serialize")) {
      QuadraticSplineSurface local_surface =
        spline_surface.subsurface(patch_indices);
      local_surface.save_obj("local_features.obj");
      local_surface.write_spline("local_features");
    }
  };
  polyscope::state::userCallback = viewer_callback;

  polyscope::show();
  polyscope::removeAllStructures();
}

// View an intersection node in the contour network
// Warning: Should only be used for small meshes
void
ContourNetwork::view_intersection_node(
  const QuadraticSplineSurface& spline_surface,
  NodeIndex node_index) const
{
  // Get intersection nodes
  NodeIndex intersecting_node_index = intersection(node_index);
  std::vector<NodeIndex> node_indices = { node_index, intersecting_node_index };

  // Get intersection segments
  SegmentIndex in_segment = in(node_index);
  SegmentIndex out_segment = out(node_index);
  SegmentIndex intersecting_in_segment = in(intersecting_node_index);
  SegmentIndex intersecting_out_segment = out(intersecting_node_index);
  std::vector<SegmentIndex> segment_indices = {
    in_segment, out_segment, intersecting_in_segment, intersecting_out_segment
  };

  // Get intersection patches
  size_t in_patch = get_segment_label(in_segment, "surface_patch");
  size_t out_patch = get_segment_label(out_segment, "surface_patch");
  size_t intersecting_in_patch =
    get_segment_label(intersecting_in_segment, "surface_patch");
  size_t intersecting_out_patch =
    get_segment_label(intersecting_out_segment, "surface_patch");
  std::vector<QuadraticSplineSurface::PatchIndex> patch_indices = {
    in_patch, out_patch, intersecting_in_patch, intersecting_out_patch
  };

  spline_surface.add_surface_to_viewer();
  add_spatial_network_to_viewer();

  // View with local features
  view_local_features(
    spline_surface, patch_indices, segment_indices, node_indices);
}

// View a cusp node in the contour network
// Warning: Should only be used for small meshes
void
ContourNetwork::view_cusp_node(const QuadraticSplineSurface& spline_surface,
                               NodeIndex node_index) const
{
  // Get cusp nodes
  std::vector<NodeIndex> node_indices = { node_index };

  // Get cusp segments
  SegmentIndex in_segment = in(node_index);
  SegmentIndex out_segment = out(node_index);
  std::vector<SegmentIndex> segment_indices = { in_segment, out_segment };

  // Get cusp patches
  size_t in_patch = get_segment_label(in_segment, "surface_patch");
  size_t out_patch = get_segment_label(out_segment, "surface_patch");
  std::vector<QuadraticSplineSurface::PatchIndex> patch_indices = { in_patch,
                                                                    out_patch };
  spline_surface.add_surface_to_viewer(); // FIXME
  add_spatial_network_to_viewer();        // FIXME

  view_local_features(
    spline_surface, patch_indices, segment_indices, node_indices);
}
