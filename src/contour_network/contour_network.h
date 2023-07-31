// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

#include "compute_contours.h"
#include "compute_intersections.h"
#include "compute_ray_intersections.h"
#include "projected_curve_network.h"
#include "quadratic_spline_surface.h"
#include "rational_function.h"

/// \file contour_network.h
///
/// Methods to compute the contour curve network for a spline surface with view
/// frame.

enum class InvisibilityMethod
{
  none,       // Set all QI to 0
  direct,     // Ray test per segment
  chaining,   // Ray test per chain of segments between features
  propagation // Ray test for connected components with local propagation
};

/// @brief Parameters for the invisibility computation
struct InvisibilityParameters
{
  double pad_amount = 1e-9; /// Padding for contour domains
  bool write_contour_soup =
    false; /// Option to write contours before graph construction for diagnostics

  /// Method for computing quantitative visibility
  InvisibilityMethod invisibility_method = InvisibilityMethod::chaining;

  /// Options to view each local propagation step during computation for debugging
  bool view_intersections = false;
  bool view_cusps = false;

  /// Options for redundancy checks
  bool poll_chain_segments =
    true; // Sample and poll 3 segments for majority per chain QI
  bool poll_segment_points =
    false; // Sample and poll 3 points for majority per segment QI

  /// Consistency checks
  bool check_chaining = false;
  bool check_propagation = false;
};

/// @brief Class to compute the projected contours of a quadratic spline surface
/// and represent them as a curve network. Also computes the quantitative
/// invisibility for the contours.
class ContourNetwork : public ProjectedCurveNetwork
{
public:
  /// Default constructor
  ContourNetwork();

  /// Constructor that takes a spline surface and computes the full projected
  /// contour curve network with standard viewing frame along the z axis.
  ///
  /// @param[in] spline_surface: quadratic spline surface to build contours for
  /// @param[in] intersect_params: parameters for the intersection methods
  /// @param[in] intersect_params: parameters for the invisibility methods
  /// @param[in] patch_boundary_edges: patch boundary edge indices (default
  /// none)
  ContourNetwork(const QuadraticSplineSurface& spline_surface,
                 const IntersectionParameters& intersect_params,
                 const InvisibilityParameters& invisibility_params,
                 std::vector<std::pair<int, int>> patch_boundary_edges =
                   std::vector<std::pair<int, int>>(0));

  /// Reset contour network to a default state.
  void clear();

  /// Add the spatial network to the viewer
  void add_contour_network_to_viewer() const;

  /// View the contour network surface and spatial network
  ///
  /// @param[in] spline_surface: underlying quadratic spline surface
  void view(const QuadraticSplineSurface& spline_surface) const;

  /// View the contour network without the underlying surface
  void view_contours() const;

  /// Save a screenshot of the contour network to file
  ///
  /// @param[in] filename: file to save the screenshot to
  /// @param[in] spline_surface: underlying quadratic spline surface
  /// @param[in] camera_position: camera position for the screenshot
  /// @param[in] camera_target: camera target for the screenshot
  /// @param[in] use_orthographic: use orthographic perspective if true
  void screenshot(const std::string& filename,
                  const QuadraticSplineSurface& spline_surface,
                  SpatialVector camera_position = SpatialVector(0, 0, 2),
                  SpatialVector camera_target = SpatialVector(0, 0, 0),
                  bool use_orthographic = false) const;

  /// Save a screenshot of just the contours to file
  ///
  /// @param[in] filename: file to save the screenshot to
  /// @param[in] camera_position: camera position for the screenshot
  /// @param[in] camera_target: camera target for the screenshot
  void write_rasterized_contours(
    const std::string& filename,
    SpatialVector camera_position = SpatialVector(0, 0, 2),
    SpatialVector camera_target = SpatialVector(0, 0, 0)) const;

  /// Reset timing information
  void reset_counter()
  {
    ray_intersection_call = 0;
    ray_bounding_box_call = 0;
    ray_number = 0;
    chain_number = 0;
    segment_number = 0;
    interior_cusp_number = 0;
    boundary_cusp_number = 0;
    intersection_call = 0;

    surface_update_position_time = 0;
    compute_contour_time = 0;
    compute_cusp_time = 0;
    compute_intersection_time = 0;
    compute_visibility_time = 0;
    compute_projected_time = 0;
  }

  long long ray_intersection_call = 0;
  long long ray_bounding_box_call = 0;
  long long ray_number = 0;
  long long chain_number = 0;
  long long segment_number = 0;
  long long interior_cusp_number = 0;
  long long boundary_cusp_number = 0;
  long long intersection_call = 0;

  double surface_update_position_time = 0;
  double compute_contour_time = 0;
  double compute_cusp_time = 0;
  double compute_intersection_time = 0;
  double compute_visibility_time = 0;
  double compute_projected_time = 0;

private:
  void init_contour_network(
    const QuadraticSplineSurface& spline_surface,
    const IntersectionParameters& intersect_params,
    const InvisibilityParameters& invisibility_params,
    std::vector<std::pair<int, int>> patch_boundary_edges);

  // *********************
  // Direct QI Computation
  // *********************

  void generate_ray_mapping_coeffs(const SpatialVector& sample_point,
                                   Matrix2x3r& ray_mapping_coeffs) const;

  int compute_quantitative_invisibility_from_ray_intersections(
    const Matrix2x3r& ray_mapping_coeffs,
    const SpatialVector& point,
    const std::vector<double>& ray_intersections) const;

  int compute_segment_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex segment_index,
    const InvisibilityParameters& invisibility_params);

  int compute_chain_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex start_segment_index,
    const InvisibilityParameters& invisibility_params);

  // ****************
  // Chain QI Methods
  // ****************

  NodeIndex chain_quantitative_invisibility_forward(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex start_segment_index,
    const InvisibilityParameters& invisibility_params);

  NodeIndex chain_quantitative_invisibility_backward(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex start_segment_index,
    const InvisibilityParameters& invisibility_params);

  // ****************************
  // Local QI Propagation Methods
  // ****************************

  void propagate_quantitative_invisibility_forward(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex start_segment_index,
    const InvisibilityParameters& invisibility_params);

  void propagate_quantitative_invisibility_backward(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex start_segment_index,
    const InvisibilityParameters& invisibility_params);

  void propagate_quantitative_invisibility_at_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  int propagate_quantitative_invisibility_across_node(
    NodeIndex node_index,
    int change_in_quantitative_invisibility);

  bool propagate_quantitative_invisibility_at_intersection_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  bool propagate_quantitative_invisibility_at_boundary_intersection_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  bool propagate_quantitative_invisibility_at_cusp_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  void propagate_quantitative_invisibility_at_boundary_cusp_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  void propagate_quantitative_invisibility_at_marked_knot_node(
    const QuadraticSplineSurface& spline_surface,
    NodeIndex node_index,
    const InvisibilityParameters& invisibility_params);

  void check_quantitative_invisibility_propagation(
    const QuadraticSplineSurface& spline_surface,
    SegmentIndex segment_index,
    const InvisibilityParameters& invisibility_params);

  // *****************
  // Global QI Methods
  // *****************

  void compute_default_quantitative_invisibility();

  void compute_direct_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    const InvisibilityParameters& invisibility_params);

  void compute_chained_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    const InvisibilityParameters& invisibility_params);

  void compute_propagated_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    const InvisibilityParameters& invisibility_params);

  void compute_quantitative_invisibility(
    const QuadraticSplineSurface& spline_surface,
    const InvisibilityParameters& invisibility_params);

  // *******
  // Viewers
  // *******

  void view_ray_intersection(const Matrix2x3r& ray_mapping_coeffs) const;

  void view_local_features(const QuadraticSplineSurface& spline_surface,
                           const std::vector<size_t>& patch_indices,
                           const std::vector<SegmentIndex>& segment_indices,
                           const std::vector<NodeIndex>& node_indices) const;

  void view_intersection_node(const QuadraticSplineSurface& spline_surface,
                              NodeIndex node_index) const;

  void view_cusp_node(const QuadraticSplineSurface& spline_surface,
                      NodeIndex node_index) const;
};
