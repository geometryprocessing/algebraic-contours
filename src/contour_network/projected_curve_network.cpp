// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "projected_curve_network.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include "discretize.h"
#include "generate_colormap.h"
#include "write_output.h"

// Initialize the curve network without intersections directly from the input
// curve data, marking cusps between boundaries.
// FIXME Rename to indicate no cusps either
void
build_projected_curve_network_without_intersections(
  const std::vector<Conic>& parameter_segments,
  const std::vector<RationalFunction<4, 3>>& spatial_segments,
  const std::vector<RationalFunction<4, 2>>& planar_segments,
  const std::vector<std::map<std::string, int>>& segment_labels,
  const std::vector<std::vector<int>>& chains,
  const std::vector<bool>& has_cusp_at_base,
  std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<SegmentGeometry>& segments,
  std::vector<NodeGeometry>& nodes)
{
  AbstractCurveNetwork::SegmentIndex num_segments = spatial_segments.size();

  // Initialize segments without intersections
  spdlog::info("Initializing segment geometry without intersections");
  segments.resize(num_segments);
  for (AbstractCurveNetwork::SegmentIndex i = 0; i < num_segments; ++i) {
    segments[i] = SegmentGeometry(parameter_segments[i],
                                  spatial_segments[i],
                                  planar_segments[i],
                                  segment_labels[i]);
  }

  // Initialize node bases without intersections
  spdlog::info("Initializing node geometry without intersections");
  nodes.resize(num_segments);
  for (AbstractCurveNetwork::NodeIndex i = 0; i < num_segments; ++i) {
    nodes[i] = NodeGeometry();
    if (has_cusp_at_base[i]) {
      spdlog::info("Marking node {} as boundary cusp", i);
      nodes[i].mark_as_boundary_cusp();
    }
  }

  // Initialize one node for the base of each segment. Since we assume all
  // vertices are degree 2 or 1, there are always as many nodes as segments.
  arange(num_segments, out_array);

  // Connect the tips of segments to base nodes at the same location in space or
  // cap them with new nodes if no such base node exists. WARNING: Assumes the
  // base of the segment has the same index
  to_array.assign(num_segments, -1);
  for (size_t i = 0; i < chains.size(); ++i) {
    for (size_t j = 1; j < chains[i].size(); ++j) {
      to_array[chains[i][j - 1]] = chains[i][j];
    }

    // Close vector if the curve is closed or add a cap node otherwise
    AbstractCurveNetwork::SegmentIndex first_segment = chains[i].front();
    AbstractCurveNetwork::SegmentIndex last_segment = chains[i].back();
    if (vector_equal(spatial_segments[first_segment].start_point(),
                     spatial_segments[last_segment].end_point(),
                     1e-6)) {
      to_array[last_segment] = first_segment;
      out_array[to_array[last_segment]] = first_segment;
    } else {
      // Add a node with no outgoing vector
      to_array[last_segment] = out_array.size();
      out_array.push_back(-1);
      nodes.push_back(NodeGeometry());
      nodes.back().mark_as_path_end_node();
    }
  }
}

// Record the start of open chains and also mark an arbitrary node on each
// closed contour
void
mark_open_chain_endpoints(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  const std::vector<std::vector<int>>& chains,
  std::vector<NodeGeometry>& nodes)
{
  // Build from array from the network topology
  std::vector<AbstractCurveNetwork::NodeIndex> from_array;
  build_from_array(to_array, out_array, from_array);

  for (size_t i = 0; i < chains.size(); ++i) {
    // Get the first and last segments in the chain
    AbstractCurveNetwork::SegmentIndex first_segment = chains[i].front();
    AbstractCurveNetwork::SegmentIndex last_segment = chains[i].back();

    if (to_array[last_segment] != from_array[first_segment]) {
      AbstractCurveNetwork::NodeIndex start_node = from_array[first_segment];
      AbstractCurveNetwork::NodeIndex end_node = to_array[last_segment];
      nodes[start_node].mark_as_path_start_node();
      nodes[end_node].mark_as_path_end_node();
    }
  }
}

// Helper function to split a given segment at a knot, updating the curve
// network data and getting the indices of the new node and segments in them
void
split_segment_at_knot(
  AbstractCurveNetwork::SegmentIndex original_segment_index,
  double knot,
  std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<SegmentGeometry>& segments,
  std::vector<NodeGeometry>& nodes,
  AbstractCurveNetwork::SegmentIndex& lower_segment_index,
  AbstractCurveNetwork::SegmentIndex& upper_segment_index,
  AbstractCurveNetwork::NodeIndex& knot_node_index)
{
  // Check segment validity
  if (to_array.size() != segments.size()) {
    spdlog::error("Segment geometry and topology mismatch");
    return;
  }

  // Check node validity
  if (out_array.size() != nodes.size()) {
    spdlog::error("Node geometry and topology mismatch");
    return;
  }

  // Split the segment at the knot
  SegmentGeometry lower_segment;
  SegmentGeometry upper_segment;
  segments[original_segment_index].split_at_knot(
    knot, lower_segment, upper_segment);
  spdlog::trace("Segment {} split into {} and {} at {}",
                segments[original_segment_index],
                lower_segment,
                upper_segment,
                knot);

  // Add the new split segments to the array
  AbstractCurveNetwork::SegmentIndex new_segment_index =
    static_cast<AbstractCurveNetwork::SegmentIndex>(segments.size());
  segments[original_segment_index] = lower_segment;
  segments.push_back(upper_segment);
  nodes.push_back(NodeGeometry());

  // Update output segment indices
  lower_segment_index = original_segment_index;
  upper_segment_index = new_segment_index;

  // Add node for the knot to the topology
  AbstractCurveNetwork::NodeIndex original_segment_to_index =
    to_array[original_segment_index];
  knot_node_index =
    static_cast<AbstractCurveNetwork::NodeIndex>(out_array.size());
  out_array.push_back(new_segment_index);
  to_array[original_segment_index] = knot_node_index;
  to_array.push_back(original_segment_to_index);
}

// Remove intersections between adjacent segments at the tip/base of the
// contours and one of two intersections if it is at the respective tip and base
// of two adjacent segments
void
remove_redundant_intersections(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  int num_intersections,
  std::vector<std::vector<IntersectionData>>& intersection_data)
{
  // Build map from intersection ids to their index in intersection data
  std::vector<std::vector<std::pair<size_t, size_t>>> intersection_segments(
    num_intersections);
  for (int i = 0; i < num_intersections; ++i) {
    intersection_segments[i].clear();
  }
  for (size_t i = 0; i < intersection_data.size(); ++i) {
    for (size_t j = 0; j < intersection_data[i].size(); ++j) {
      int intersection_id = intersection_data[i][j].id;
      intersection_segments[intersection_id].push_back(std::make_pair(i, j));
    }
  }

  for (size_t i = 0; i < intersection_segments.size(); ++i) {
    if (intersection_segments[i].size() != 2) {
      spdlog::error("Should have two intersection per id");
      continue;
    }

    // Get intersection data per segment
    std::pair<size_t, size_t> first_segment_indices =
      intersection_segments[i][0];
    std::pair<size_t, size_t> second_segment_indices =
      intersection_segments[i][1];
    AbstractCurveNetwork::SegmentIndex first_segment_index =
      first_segment_indices.first;
    AbstractCurveNetwork::SegmentIndex second_segment_index =
      second_segment_indices.first;
    AbstractCurveNetwork::SegmentIndex first_intersection_index =
      first_segment_indices.second;
    AbstractCurveNetwork::SegmentIndex second_intersection_index =
      second_segment_indices.second;
    IntersectionData const& first_segment_data =
      intersection_data[first_segment_index][first_intersection_index];
    IntersectionData const& second_segment_data =
      intersection_data[second_segment_index][second_intersection_index];

    // Don't process intersections already marked as redundant
    if (intersection_data[first_segment_index][first_intersection_index]
          .is_redundant)
      continue;
    if (intersection_data[second_segment_index][second_intersection_index]
          .is_redundant)
      continue;

    // Remove intersections at endpoints between adjacent segments
    if (out_array[to_array[first_segment_index]] == second_segment_index) {
      if ((first_segment_data.is_tip) || (second_segment_data.is_base)) {
        intersection_data[first_segment_index][first_intersection_index]
          .is_redundant = true;
        intersection_data[second_segment_index][second_intersection_index]
          .is_redundant = true;
      }
    }
    if (out_array[to_array[second_segment_index]] == first_segment_index) {
      if ((second_segment_data.is_tip) || (first_segment_data.is_base)) {
        intersection_data[first_segment_index][first_intersection_index]
          .is_redundant = true;
        intersection_data[second_segment_index][second_intersection_index]
          .is_redundant = true;
      }
    }

    // If a contour intersects the tip of a contour and there is an intersection
    // with next contour that is not already marked as redundant, mark this one
    // as redundant
    if (first_segment_data.is_tip) {
      for (size_t j = 0; j < intersection_data[second_segment_index].size();
           ++j) {
        if (intersection_data[second_segment_index][j].is_redundant)
          continue;

        AbstractCurveNetwork::SegmentIndex third_segment_index =
          intersection_data[second_segment_index][j].intersection_index;
        if (out_array[to_array[first_segment_index]] == third_segment_index) {
          intersection_data[first_segment_index][first_intersection_index]
            .is_redundant = true;
          intersection_data[second_segment_index][second_intersection_index]
            .is_redundant = true;
        }
      }
    }
    if (second_segment_data.is_tip) {
      for (size_t j = 0; j < intersection_data[first_segment_index].size();
           ++j) {
        if (intersection_data[first_segment_index][j].is_redundant)
          continue;

        AbstractCurveNetwork::SegmentIndex third_segment_index =
          intersection_data[first_segment_index][j].intersection_index;
        if (out_array[to_array[second_segment_index]] == third_segment_index) {
          intersection_data[first_segment_index][first_intersection_index]
            .is_redundant = true;
          intersection_data[second_segment_index][second_intersection_index]
            .is_redundant = true;
        }
      }
    }
  }
}

// Split segments at all of the intersection points, and create maps from the
// new segments to their original indices and from the original indices to their
// corresponding split segment indices
void
split_segments_at_intersections(
  const std::vector<std::vector<IntersectionData>>& intersection_data,
  int num_intersections,
  std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<SegmentGeometry>& segments,
  std::vector<NodeGeometry>& nodes,
  std::vector<AbstractCurveNetwork::SegmentIndex>& original_segment_indices,
  std::vector<std::vector<AbstractCurveNetwork::SegmentIndex>>&
    split_segment_indices,
  std::vector<std::vector<AbstractCurveNetwork::NodeIndex>>& intersection_nodes)
{
  spdlog::info("Splitting segments at intersections");
  size_t num_segments = segments.size();
  size_t num_nodes = nodes.size();
  intersection_nodes.assign(num_intersections,
                            std::vector<AbstractCurveNetwork::NodeIndex>(0));
  // size_t num_intersections = nested_vector_size(intersections);

  // Check segment validity
  if (to_array.size() != segments.size()) {
    spdlog::error("Segment geometry and topology mismatch");
    return;
  }

  // Check node validity
  if (out_array.size() != nodes.size()) {
    spdlog::error("Node geometry and topology mismatch");
    return;
  }

  // Initially the map is the identity
  arange(num_segments, original_segment_indices);
  split_segment_indices.resize(num_segments);
  for (size_t i = 0; i < num_segments; ++i) {
    split_segment_indices[i].clear();
    split_segment_indices[i].push_back(i);
  }

  // Mark all base nodes before modifying connectivity
  std::vector<AbstractCurveNetwork::NodeIndex> from_array;
  build_from_array(to_array, out_array, from_array);
  for (size_t i = 0; i < intersection_data.size(); ++i) {
    AbstractCurveNetwork::SegmentIndex segment_index = i;
    SegmentGeometry segment = segments[segment_index];
    std::vector<IntersectionData> segment_intersection_data =
      intersection_data[segment_index];
    for (size_t j = 0; j < segment_intersection_data.size(); ++j) {
      // Skip redundant intersections
      if (segment_intersection_data[j].is_redundant) {
        continue;
      }

      // Mark intersections at the tip of the contour
      if (segment_intersection_data[j].is_base) {
        // Add node to global registry
        size_t intersection_id = segment_intersection_data[j].id;
        AbstractCurveNetwork::NodeIndex intersection_node_index =
          from_array[segment_index];
        intersection_nodes[intersection_id].push_back(intersection_node_index);
      }
    }
  }

  for (size_t i = 0; i < intersection_data.size(); ++i) {
    // Get segment i and sorted intersections
    AbstractCurveNetwork::SegmentIndex segment_index = i;
    SegmentGeometry segment = segments[segment_index];
    // std::vector<double> segment_intersections = intersections[segment_index];
    // std::sort(segment_intersections.begin(), segment_intersections.end());
    // split_segment_indices[segment_index].resize(segment_intersections.size()
    // + 1);
    std::vector<IntersectionData> segment_intersection_data =
      intersection_data[segment_index];
    std::sort(segment_intersection_data.begin(),
              segment_intersection_data.end(),
              knot_less_than());

    split_segment_indices[segment_index].reserve(
      segment_intersection_data.size() + 1);
    for (size_t j = 0; j < segment_intersection_data.size(); ++j) {
      // Skip redundant intersections
      if (segment_intersection_data[j].is_redundant) {
        continue;
      }

      // Skip already handled base points
      if (segment_intersection_data[j].is_base) {
        continue;
      }

      // Mark intersections at the tip of the contour
      if (segment_intersection_data[j].is_tip) {
        // Mark base as intersection
        AbstractCurveNetwork::NodeIndex intersection_node_index =
          to_array[segment_index];

        // Add node to global registry
        size_t intersection_id = segment_intersection_data[j].id;
        intersection_nodes[intersection_id].push_back(intersection_node_index);
      }
      // Split segment at intersection and mark if interior
      else {
        // FIXME The base is currently not tracked and can't be removed
        // As a hack, length zero intersections are removed

        // double intersection = segment_intersections[j];
        double intersection = segment_intersection_data[j].knot;
        AbstractCurveNetwork::SegmentIndex lower_segment_index;
        AbstractCurveNetwork::SegmentIndex upper_segment_index;
        AbstractCurveNetwork::NodeIndex intersection_node_index;
        split_segment_at_knot(segment_index,
                              intersection,
                              to_array,
                              out_array,
                              segments,
                              nodes,
                              lower_segment_index,
                              upper_segment_index,
                              intersection_node_index);

        // Add node to global registry
        size_t intersection_id = segment_intersection_data[j].id;
        intersection_nodes[intersection_id].push_back(intersection_node_index);

        // Update segment indices record
        original_segment_indices.push_back(-1);
        split_segment_indices[i].back() = lower_segment_index;
        split_segment_indices[i].push_back(upper_segment_index);
        original_segment_indices[lower_segment_index] = i;
        original_segment_indices[upper_segment_index] = i;

        // Continue splitting upper segment
        segment_index = upper_segment_index;
      }
    }
  }

  spdlog::debug(
    "Split {} segments with {} nodes into {} segments with {} nodes",
    num_segments,
    num_nodes,
    segments.size(),
    nodes.size());
}

// Create a map from nodes to their corresponding intersection point if they are
// an intersection node and -1 otherwise.
void
connect_segment_intersections(
  const std::vector<SegmentGeometry>& segments,
  const std::vector<std::vector<IntersectionData>>& intersection_data,
  const std::vector<std::vector<AbstractCurveNetwork::NodeIndex>>&
    intersection_nodes,
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  const std::vector<std::vector<AbstractCurveNetwork::SegmentIndex>>&
    split_segment_indices,
  std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array,
  std::vector<NodeGeometry>& nodes)
{
  spdlog::info("Connecting segments intersections");

  // Initialize all intersection indices to -1
  AbstractCurveNetwork::NodeIndex num_nodes(out_array.size());
  intersection_array.resize(num_nodes);
  std::fill(intersection_array.begin(), intersection_array.end(), -1);

  bool use_combinatorial_data = true;
  if (use_combinatorial_data) {
    for (size_t i = 0; i < intersection_nodes.size(); ++i) {
      if (intersection_nodes[i].size() == 0)
        continue;

      if (intersection_nodes[i].size() != 2) {
        spdlog::error("Node intersection is not a pair");
        continue;
      }
      AbstractCurveNetwork::NodeIndex first_intersection_node =
        intersection_nodes[i][0];
      AbstractCurveNetwork::NodeIndex second_intersection_node =
        intersection_nodes[i][1];

      // Ignore intersections that already exist
      if (intersection_array[first_intersection_node] >= 0) {
        continue;
      } else if (intersection_array[second_intersection_node] >= 0) {
        continue;
      } else {
        intersection_array[first_intersection_node] = second_intersection_node;
        intersection_array[second_intersection_node] = first_intersection_node;
        nodes[first_intersection_node].mark_as_intersection();
        nodes[second_intersection_node].mark_as_intersection();
      }
    }

    return;
  }

  AbstractCurveNetwork::SegmentIndex num_original_segments(
    split_segment_indices.size());
  for (AbstractCurveNetwork::SegmentIndex i = 0; i < num_original_segments;
       ++i) {
    AbstractCurveNetwork::SegmentIndex num_original_segment_partitions(
      split_segment_indices[i].size());
    for (AbstractCurveNetwork::SegmentIndex j = 0;
         j < num_original_segment_partitions - 1;
         ++j) {
      AbstractCurveNetwork::SegmentIndex original_segment_partition =
        split_segment_indices[i][j];
      AbstractCurveNetwork::NodeIndex new_node_index =
        to_array[original_segment_partition];

      // Skip nodes that have already been processed
      if (intersection_array[new_node_index] >= 0)
        continue;

      // AbstractCurveNetwork::SegmentIndex
      // intersecting_segment_index(intersection_indices[i][j]);
      AbstractCurveNetwork::SegmentIndex intersecting_segment_index =
        intersection_data[i][j].intersection_index;
      AbstractCurveNetwork::SegmentIndex num_intersecting_segment_partitions(
        split_segment_indices[intersecting_segment_index].size());
      for (AbstractCurveNetwork::SegmentIndex k = 0;
           k < num_intersecting_segment_partitions;
           ++k) {
        AbstractCurveNetwork::SegmentIndex intersecting_segment_partition =
          split_segment_indices[intersecting_segment_index][k];
        if (vector_equal(segments[original_segment_partition]
                           .get_planar_curve()
                           .end_point(),
                         segments[intersecting_segment_partition]
                           .get_planar_curve()
                           .end_point())) {
          intersection_array[to_array[original_segment_partition]] =
            to_array[intersecting_segment_partition];
          intersection_array[to_array[intersecting_segment_partition]] =
            to_array[original_segment_partition];
        }
      }
    }
  }

  // Check consistency of out and intersection arrays
  if (out_array.size() != intersection_array.size()) {
    spdlog::error("Inconsistent number of intersections and nodes after "
                  "intersections are split");
    return;
  }
}

// Split segments at all of the cusp points, and create maps from the new
// segments to their original indices and from the original indices to their
// corresponding split segment indices
void
split_segments_at_cusps(
  const std::vector<std::vector<double>>& interior_cusps,
  std::vector<AbstractCurveNetwork::SegmentIndex>& original_segment_indices,
  std::vector<std::vector<AbstractCurveNetwork::SegmentIndex>>&
    split_segment_indices,
  std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array,
  std::vector<SegmentGeometry>& segments,
  std::vector<NodeGeometry>& nodes)
{
  spdlog::info("Splitting segments at cusps");

  size_t num_segments = segments.size();
  size_t num_nodes = nodes.size();
  AbstractCurveNetwork::SegmentIndex num_original_segments =
    interior_cusps.size();

  // Check segment validity
  if (to_array.size() != segments.size()) {
    spdlog::error("Segment geometry and topology mismatch");
    return;
  }

  // Check node validity
  if (out_array.size() != nodes.size()) {
    spdlog::error("Node geometry and topology mismatch");
    return;
  }

  for (AbstractCurveNetwork::SegmentIndex i = 0; i < num_original_segments;
       ++i) {
    for (size_t j = 0; j < interior_cusps[i].size(); ++j) {
      double cusp = interior_cusps[i][j];
      for (size_t k = 0; k < split_segment_indices[i].size(); ++k) {
        AbstractCurveNetwork::SegmentIndex segment_partition_index =
          split_segment_indices[i][k];
        SegmentGeometry& segment_partition = segments[segment_partition_index];
        if (segment_partition.get_parameter_curve().domain().is_in_interior(
              cusp)) {
          AbstractCurveNetwork::SegmentIndex lower_segment_index;
          AbstractCurveNetwork::SegmentIndex upper_segment_index;
          AbstractCurveNetwork::NodeIndex knot_node_index;
          split_segment_at_knot(segment_partition_index,
                                cusp,
                                to_array,
                                out_array,
                                segments,
                                nodes,
                                lower_segment_index,
                                upper_segment_index,
                                knot_node_index);

          // Mark new node as interior cusp
          spdlog::info("Marking node {} as interior cusp", knot_node_index);
          nodes[knot_node_index].mark_as_interior_cusp();

          // Add trivial intersection information
          AbstractCurveNetwork::NodeIndex num_nodes = intersection_array.size();
          if (knot_node_index == num_nodes) {
            intersection_array.push_back(-1);
          } else {
            spdlog::error("Invalid node index used for knot split");
          }

          // Update segment indices record
          split_segment_indices[i][k] = lower_segment_index;
          split_segment_indices[i].push_back(upper_segment_index);
          assert((static_cast<size_t>(lower_segment_index) ==
                  original_segment_indices.size()) ||
                 (static_cast<size_t>(upper_segment_index) ==
                  original_segment_indices.size()));
          original_segment_indices.push_back(i); // FIXME Risky assumption
          original_segment_indices[lower_segment_index] = i;
          original_segment_indices[upper_segment_index] = i;

          break;
        }
      }
    }
  }

  spdlog::debug(
    "Split {} segments with {} nodes into {} segments with {} nodes",
    num_segments,
    num_nodes,
    segments.size(),
    nodes.size());

  // Check consistency of out and intersection arrays
  if (out_array.size() != intersection_array.size()) {
    spdlog::error(
      "Inconsistent number of intersections and nodes after cusps are split");
    return;
  }
}

/// @brief View a polygon with boolean function in polyscope
void
view_polygon(std::vector<PlanarPoint>& points,
             std::vector<int>& polyline,
             std::vector<bool>& is_marked)
{
  polyscope::init();

  // Build points
  MatrixXr points_mat_2d = convert_nested_vector_to_matrix(points);
  MatrixXr points_mat(points_mat_2d.rows(), 3);
  points_mat.setZero();
  points_mat.block(0, 0, points_mat_2d.rows(), 2) = points_mat_2d;

  // Build edges
  std::vector<Edge> edges(polyline.size());
  edges[0] = { polyline.back(), polyline.front() };
  for (size_t j = 1; j < polyline.size(); ++j) {
    edges[j] = { polyline[j - 1], polyline[j] };
  }

  // Build vertex function
  VectorXr vertex_function_vector(is_marked.size());
  for (size_t j = 0; j < is_marked.size(); ++j) {
    vertex_function_vector[j] = is_marked[j] ? 1 : 0;
  }

  polyscope::registerCurveNetwork("polygon", points_mat, edges);
  polyscope::registerPointCloud("polygon_vertices", points_mat);
  polyscope::getCurveNetwork("polygon")->setRadius(0.0025);
  polyscope::getPointCloud("polygon_vertices")
    ->addScalarQuantity("is_marked", vertex_function_vector);
  polyscope::show();
  polyscope::removeAllStructures();
}

// ***********************
// Projected Curve Network
// ***********************

ProjectedCurveNetwork::ProjectedCurveNetwork()
{
  clear_topology();
  clear_geometry();
}

ProjectedCurveNetwork::ProjectedCurveNetwork(
  const std::vector<Conic>& parameter_segments,
  const std::vector<RationalFunction<4, 3>>& spatial_segments,
  const std::vector<RationalFunction<4, 2>>& planar_segments,
  const std::vector<std::map<std::string, int>>& segment_labels,
  const std::vector<std::vector<int>>& chains,
  const std::vector<int>& chain_labels,
  const std::vector<std::vector<double>> interior_cusps,
  const std::vector<bool>& has_cusp_at_base,
  const std::vector<std::vector<double>>& intersections,
  const std::vector<std::vector<size_t>>& intersection_indices,
  std::vector<std::vector<IntersectionData>>& intersection_data,
  int num_intersections)
{
  init_projected_curve_network(parameter_segments,
                               spatial_segments,
                               planar_segments,
                               segment_labels,
                               chains,
                               chain_labels,
                               interior_cusps,
                               has_cusp_at_base,
                               intersections,
                               intersection_indices,
                               intersection_data,
                               num_intersections);

#if CHECK_VALIDITY
  if (!is_valid_abstract_curve_network()) {
    spdlog::error("Invalid abstract curve network made");
    return;
  }

  if (!is_valid_projected_curve_network()) {
    spdlog::error("Invalid projected curve network made");
    return;
  }
#endif
}

// ******
// Counts
// ******

AbstractCurveNetwork::SegmentIndex
ProjectedCurveNetwork::num_segments() const
{
  return SegmentIndex(m_segments.size());
}

AbstractCurveNetwork::NodeIndex
ProjectedCurveNetwork::num_nodes() const
{
  return NodeIndex(m_nodes.size());
}

// ****************
// Segment geometry
// ****************

Conic const&
ProjectedCurveNetwork::segment_parameter_curve(SegmentIndex segment_index) const
{
  assert(is_valid_segment_index(segment_index));
  return m_segments[segment_index].get_parameter_curve();
}

RationalFunction<4, 2> const&
ProjectedCurveNetwork::segment_planar_curve(SegmentIndex segment_index) const
{
  assert(is_valid_segment_index(segment_index));
  return m_segments[segment_index].get_planar_curve();
}

RationalFunction<4, 3> const&
ProjectedCurveNetwork::segment_spatial_curve(SegmentIndex segment_index) const
{
  assert(is_valid_segment_index(segment_index));
  return m_segments[segment_index].get_spatial_curve();
}

int
ProjectedCurveNetwork::get_segment_label(SegmentIndex segment_index,
                                         const std::string& label_name) const
{
  return m_segments[segment_index].get_segment_label(label_name);
}

int
ProjectedCurveNetwork::get_segment_quantitative_invisibility(
  SegmentIndex segment_index) const
{
  return m_segments[segment_index].get_quantitative_invisibility();
}

void
ProjectedCurveNetwork::set_segment_quantitative_invisibility(
  SegmentIndex segment_index,
  int new_quantitative_invisibility)
{
  if (new_quantitative_invisibility < 0) {
    spdlog::error("Cannot set negative segment QI");
    return;
  }
  m_segments[segment_index].set_quantitative_invisibility(
    new_quantitative_invisibility);
}

int
ProjectedCurveNetwork::get_node_quantitative_invisibility(
  NodeIndex node_index) const
{
  return m_nodes[node_index].get_quantitative_invisibility();
}

void
ProjectedCurveNetwork::set_node_quantitative_invisibility(
  NodeIndex node_index,
  int new_quantitative_invisibility)
{
  if (new_quantitative_invisibility < 0) {
    spdlog::error("Cannot set negative node QI");
    return;
  }
  m_nodes[node_index].set_quantitative_invisibility(
    new_quantitative_invisibility);
}

void
ProjectedCurveNetwork::enumerate_parameter_curves(
  std::vector<Conic>& parameter_curves) const
{
  parameter_curves.resize(num_segments());
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    parameter_curves[i] = m_segments[i].get_parameter_curve();
  }
}

void
ProjectedCurveNetwork::enumerate_planar_curves(
  std::vector<RationalFunction<4, 2>>& planar_curves) const
{
  planar_curves.resize(num_segments());
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    planar_curves[i] = m_segments[i].get_planar_curve();
  }
}

void
ProjectedCurveNetwork::enumerate_spatial_curves(
  std::vector<RationalFunction<4, 3>>& spatial_curves) const
{
  spatial_curves.resize(num_segments());
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    spatial_curves[i] = m_segments[i].get_spatial_curve();
  }
}

void
ProjectedCurveNetwork::enumerate_spatial_nodes(
  std::vector<SpatialVector>& spatial_knot_nodes,
  std::vector<SpatialVector>& spatial_marked_knot_nodes,
  std::vector<SpatialVector>& spatial_intersection_nodes,
  std::vector<SpatialVector>& spatial_interior_cusp_nodes,
  std::vector<SpatialVector>& spatial_boundary_cusp_nodes,
  std::vector<SpatialVector>& spatial_path_start_nodes,
  std::vector<SpatialVector>& spatial_path_end_nodes) const
{
  spatial_knot_nodes.reserve(num_nodes());
  spatial_marked_knot_nodes.reserve(num_nodes());
  spatial_intersection_nodes.reserve(num_nodes());
  spatial_interior_cusp_nodes.reserve(num_nodes());
  spatial_boundary_cusp_nodes.reserve(num_nodes());
  spatial_path_start_nodes.reserve(num_nodes());
  spatial_path_end_nodes.reserve(num_nodes());

  // Iterate over nodes, sorting into the appropriate categories
  for (NodeIndex ni = 0; ni < num_nodes(); ++ni) {
    SpatialVector spatial_point = node_spatial_point(ni);
    if (is_knot_node(ni)) {
      spatial_knot_nodes.push_back(spatial_point);
    } else if (is_marked_knot_node(ni)) {
      spatial_marked_knot_nodes.push_back(spatial_point);
    } else if (is_intersection_node(ni)) {
      spatial_intersection_nodes.push_back(spatial_point);
    } else if (is_interior_cusp_node(ni)) {
      spatial_interior_cusp_nodes.push_back(spatial_point);
    } else if (is_boundary_cusp_node(ni)) {
      spatial_boundary_cusp_nodes.push_back(spatial_point);
    } else if (is_path_start_node(ni)) {
      spatial_path_start_nodes.push_back(spatial_point);
    } else if (is_path_end_node(ni)) {
      spatial_path_end_nodes.push_back(spatial_point);
    }
  }

  spdlog::debug("Enumerated {} spatial intersection nodes",
                spatial_intersection_nodes.size());
}

void
ProjectedCurveNetwork::enumerate_cusp_spatial_tangents(
  std::vector<SpatialVector>& spatial_interior_cusp_in_tangents,
  std::vector<SpatialVector>& spatial_interior_cusp_out_tangents,
  std::vector<SpatialVector>& spatial_boundary_cusp_in_tangents,
  std::vector<SpatialVector>& spatial_boundary_cusp_out_tangents) const
{
  spatial_interior_cusp_in_tangents.reserve(num_nodes());
  spatial_interior_cusp_out_tangents.reserve(num_nodes());
  spatial_boundary_cusp_in_tangents.reserve(num_nodes());
  spatial_boundary_cusp_out_tangents.reserve(num_nodes());

  // Iterate over nodes, sorting into the appropriate categories
  for (NodeIndex ni = 0; ni < num_nodes(); ++ni) {
    if (is_interior_cusp_node(ni)) {
      SpatialVector spatial_in_tangent = node_spatial_in_tangent(ni);
      SpatialVector spatial_out_tangent = node_spatial_out_tangent(ni);
      spatial_interior_cusp_in_tangents.push_back(spatial_in_tangent);
      spatial_interior_cusp_out_tangents.push_back(spatial_out_tangent);
    } else if (is_boundary_cusp_node(ni)) {
      SpatialVector spatial_in_tangent = node_spatial_in_tangent(ni);
      SpatialVector spatial_out_tangent = node_spatial_out_tangent(ni);
      spatial_boundary_cusp_in_tangents.push_back(spatial_in_tangent);
      spatial_boundary_cusp_out_tangents.push_back(spatial_out_tangent);
    }
  }
}

void
ProjectedCurveNetwork::enumerate_planar_nodes(
  std::vector<PlanarPoint>& planar_knot_nodes,
  std::vector<PlanarPoint>& planar_marked_knot_nodes,
  std::vector<PlanarPoint>& planar_intersection_nodes,
  std::vector<PlanarPoint>& planar_interior_cusp_nodes,
  std::vector<PlanarPoint>& planar_boundary_cusp_nodes,
  std::vector<PlanarPoint>& planar_path_start_nodes,
  std::vector<PlanarPoint>& planar_path_end_nodes) const
{
  planar_knot_nodes.reserve(num_nodes());
  planar_marked_knot_nodes.reserve(num_nodes());
  planar_intersection_nodes.reserve(num_nodes());
  planar_interior_cusp_nodes.reserve(num_nodes());
  planar_boundary_cusp_nodes.reserve(num_nodes());
  planar_path_start_nodes.reserve(num_nodes());
  planar_path_end_nodes.reserve(num_nodes());

  // Iterate over nodes, sorting into the appropriate categories
  for (NodeIndex ni = 0; ni < num_nodes(); ++ni) {
    PlanarPoint planar_point = node_planar_point(ni);
    if (is_knot_node(ni)) {
      planar_knot_nodes.push_back(planar_point);
    } else if (is_marked_knot_node(ni)) {
      planar_marked_knot_nodes.push_back(planar_point);
    } else if (is_intersection_node(ni)) {
      planar_intersection_nodes.push_back(planar_point);
    } else if (is_interior_cusp_node(ni)) {
      planar_interior_cusp_nodes.push_back(planar_point);
    } else if (is_boundary_cusp_node(ni)) {
      planar_boundary_cusp_nodes.push_back(planar_point);
    } else if (is_path_start_node(ni)) {
      planar_path_start_nodes.push_back(planar_point);
    } else if (is_path_end_node(ni)) {
      planar_path_end_nodes.push_back(planar_point);
    }
  }
}

void
ProjectedCurveNetwork::add_spatial_network_to_viewer() const
{
  add_spatial_segments_to_viewer();
  add_spatial_nodes_to_viewer();
}

void
ProjectedCurveNetwork::spatial_network_viewer() const
{
  polyscope::init();
  add_spatial_network_to_viewer();
  polyscope::show();
  polyscope::removeAllStructures();
}

void
ProjectedCurveNetwork::get_closed_curves(
  std::vector<NodeIndex>& closed_curve_start_nodes,
  std::vector<NodeIndex>& closed_curve_end_nodes) const
{
  closed_curve_start_nodes.clear();
  closed_curve_end_nodes.clear();

  // Maintain record of nodes that have been covered by a closed curve
  NodeIndex num_nodes = get_num_nodes();
  std::vector<bool> is_covered_node(num_nodes, false);

  // Find closed curve start and end nodes
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  closed_curve_start_nodes.reserve(chain_start_nodes.size());
  closed_curve_end_nodes.reserve(chain_start_nodes.size());
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    // Get a chain start node
    NodeIndex start_node_index = chain_start_nodes[i];
    if (is_covered_node[start_node_index])
      continue;
    is_covered_node[start_node_index] = true;

    // Go backward along the chain until the start is found or a closed loop is
    // obtained
    NodeIndex node_index = start_node_index;
    while (true) {
      // Get incoming segment
      SegmentIndex in_segment_index = in(node_index);

      // Add node if it is an initial endpoint
      if (!is_valid_segment_index(in_segment_index)) {
        closed_curve_start_nodes.push_back(node_index);
        break;
      }

      // Get previous node (always valid)
      node_index = from(in_segment_index);
      is_covered_node[node_index] = true;

      // Add node if it is the start node and thus a closed loop
      if (node_index == start_node_index) {
        closed_curve_start_nodes.push_back(node_index);
        break;
      }
    }

    // Go forward along the curve until an endpoint is reached
    node_index = start_node_index;
    while (true) {
      // Get outgoing segment
      SegmentIndex out_segment_index = out(node_index);

      // Add node if it is a terminal endpoint
      if (!is_valid_segment_index(out_segment_index)) {
        closed_curve_end_nodes.push_back(node_index);
        break;
      }

      // Get next node (always valid)
      node_index = to(out_segment_index);
      is_covered_node[node_index] = true;

      // Add node if it is the start node and thus a closed loop
      if (node_index == start_node_index) {
        closed_curve_end_nodes.push_back(node_index);
        break;
      }
    }
  }
}

void
ProjectedCurveNetwork::get_visible_curves(
  std::vector<NodeIndex>& visible_curve_start_nodes,
  std::vector<NodeIndex>& visible_curve_end_nodes) const
{
  visible_curve_start_nodes.clear();
  visible_curve_end_nodes.clear();

  // Maintain record of segments that have been covered by a visible curve
  NodeIndex num_segments = get_num_segments();
  std::vector<bool> is_covered_segment(num_segments, false);

  // Find visible curve start and end nodes
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  visible_curve_start_nodes.reserve(chain_start_nodes.size());
  visible_curve_end_nodes.reserve(chain_start_nodes.size());
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    // Get a chain start node
    NodeIndex start_node_index = chain_start_nodes[i];

    // Skip if the outgoing segment is not visible
    if ((!is_valid_segment_index(out(start_node_index))) ||
        (get_segment_quantitative_invisibility(out(start_node_index)) != 0)) {
      continue;
    }

    // Skip already covered segments
    if (is_covered_segment[out(start_node_index)])
      continue;

    // Go backward along the chain until an invalid/invisible segment is found
    // or a closed loop is obtained
    NodeIndex node_index = start_node_index;
    while (true) {
      // Get incoming segment
      SegmentIndex in_segment_index = in(node_index);

      // Add node if it is an initial endpoint or the previous segment is
      // invisible
      if ((!is_valid_segment_index(in_segment_index)) ||
          (get_segment_quantitative_invisibility(in_segment_index) != 0)) {
        visible_curve_start_nodes.push_back(node_index);
        break;
      }

      // Get previous node (always valid)
      is_covered_segment[in_segment_index] = true;
      node_index = from(in_segment_index);

      // Add node if it is the start node and thus a closed loop
      if (node_index == start_node_index) {
        visible_curve_start_nodes.push_back(node_index);
        break;
      }
    }

    // Go forward along the curve until an endpoint is reached
    node_index = start_node_index;
    while (true) {
      // Get outgoing segment
      SegmentIndex out_segment_index = out(node_index);

      // Add node if it is a terminal endpoint
      if ((!is_valid_segment_index(out_segment_index)) ||
          (get_segment_quantitative_invisibility(out_segment_index) != 0)) {
        visible_curve_end_nodes.push_back(node_index);
        break;
      }

      // Get next node (always valid)
      is_covered_segment[out_segment_index] = true;
      node_index = to(out_segment_index);

      // Add node if it is the start node and thus a closed loop
      if (node_index == start_node_index) {
        visible_curve_end_nodes.push_back(node_index);
        break;
      }
    }
  }
}

// *************
// Node geometry
// *************

bool
ProjectedCurveNetwork::is_knot_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_knot();
}

bool
ProjectedCurveNetwork::is_marked_knot_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_marked_knot();
}

bool
ProjectedCurveNetwork::is_intersection_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_intersection();
}

bool
ProjectedCurveNetwork::is_interior_cusp_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_interior_cusp();
}

bool
ProjectedCurveNetwork::is_boundary_cusp_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_boundary_cusp();
}

bool
ProjectedCurveNetwork::is_path_start_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_path_start_node();
}

bool
ProjectedCurveNetwork::is_path_end_node(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));
  if (!is_valid_node_index(node_index)) {
    spdlog::error("Invalid node query");
    return false;
  }
  return m_nodes[node_index].is_path_end_node();
}

PlanarPoint
ProjectedCurveNetwork::node_planar_point(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));

  if (is_valid_segment_index(in(node_index))) {
    return segment_planar_curve(in(node_index)).end_point();
  } else if (is_valid_segment_index(out(node_index))) {
    return segment_planar_curve(out(node_index)).start_point();
  } else {
    return PlanarPoint();
  }
}

SpatialVector
ProjectedCurveNetwork::node_spatial_point(NodeIndex node_index) const
{
  assert(is_valid_node_index(node_index));

  if (is_valid_segment_index(in(node_index))) {
    return segment_spatial_curve(in(node_index)).end_point();
  } else if (is_valid_segment_index(out(node_index))) {
    return segment_spatial_curve(out(node_index)).start_point();
  } else {
    return SpatialVector();
  }
}

PlanarPoint
ProjectedCurveNetwork::node_planar_in_tangent(NodeIndex node_index) const
{
  RationalFunction<8, 2> tangent_curve;
  segment_planar_curve(in(node_index)).compute_derivative(tangent_curve);
  return tangent_curve.end_point();
}

PlanarPoint
ProjectedCurveNetwork::node_planar_out_tangent(NodeIndex node_index) const
{
  RationalFunction<8, 2> tangent_curve;
  segment_planar_curve(out(node_index)).compute_derivative(tangent_curve);
  return tangent_curve.start_point();
}

PlanarPoint
ProjectedCurveNetwork::node_planar_tangent(NodeIndex node_index) const
{
  if (is_valid_segment_index(in(node_index))) {
    return node_planar_in_tangent(node_index);
  } else if (is_valid_segment_index(out(node_index))) {
    return node_planar_out_tangent(node_index);
  } else {
    spdlog::error("Isolated node");
    return PlanarPoint(0, 0);
  }
}

SpatialVector
ProjectedCurveNetwork::node_spatial_in_tangent(NodeIndex node_index) const
{
  RationalFunction<8, 3> tangent_curve;
  segment_spatial_curve(in(node_index)).compute_derivative(tangent_curve);
  return tangent_curve.end_point();
}

SpatialVector
ProjectedCurveNetwork::node_spatial_out_tangent(NodeIndex node_index) const
{
  RationalFunction<8, 3> tangent_curve;
  segment_spatial_curve(out(node_index)).compute_derivative(tangent_curve);
  return tangent_curve.start_point();
}

SpatialVector
ProjectedCurveNetwork::node_spatial_tangent(NodeIndex node_index) const
{
  if (is_valid_segment_index(in(node_index))) {
    return node_spatial_in_tangent(node_index);
  } else if (is_valid_segment_index(out(node_index))) {
    return node_spatial_out_tangent(node_index);
  } else {
    spdlog::error("Isolated node");
    return SpatialVector(0, 0, 0);
  }
}

// Main constructor implementation
void
ProjectedCurveNetwork::init_projected_curve_network(
  const std::vector<Conic>& parameter_segments,
  const std::vector<RationalFunction<4, 3>>& spatial_segments,
  const std::vector<RationalFunction<4, 2>>& planar_segments,
  const std::vector<std::map<std::string, int>>& segment_labels,
  const std::vector<std::vector<int>>& chains,
  const std::vector<int>& chain_labels,
  const std::vector<std::vector<double>> interior_cusps,
  const std::vector<bool>& has_cusp_at_base,
  const std::vector<std::vector<double>>& intersections,
  const std::vector<std::vector<size_t>>& intersection_indices,
  std::vector<std::vector<IntersectionData>>& intersection_data,
  int num_intersections)
{
  size_t num_segments = planar_segments.size();
  // size_t num_split_segments = num_segments + 2 * num_intersections; // FIXME
  // not correct with snapping and boundary

  if (parameter_segments.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (spatial_segments.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (planar_segments.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (segment_labels.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (chain_labels.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (interior_cusps.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (has_cusp_at_base.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (intersections.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  if (intersection_indices.size() != num_segments) {
    spdlog::error("Inconsistent number of segments");
  }
  spdlog::info("Building projected curve network for {} segments",
               num_segments);

  // Connect segments into chains before splitting at intersections
  std::vector<NodeIndex> to_array;
  std::vector<SegmentIndex> out_array;
  build_projected_curve_network_without_intersections(parameter_segments,
                                                      spatial_segments,
                                                      planar_segments,
                                                      segment_labels,
                                                      chains,
                                                      has_cusp_at_base,
                                                      to_array,
                                                      out_array,
                                                      m_segments,
                                                      m_nodes);
  assert(is_valid_curve_data(to_array, out_array));
  mark_open_chain_endpoints(to_array, out_array, chains, m_nodes);

  // Remove intersections that are redundant
  remove_redundant_intersections(
    to_array, out_array, num_intersections, intersection_data);
  assert(is_valid_curve_data(to_array, out_array));

  // Split segments at intersections while maintaining a record of the original
  // segments
  std::vector<SegmentIndex> original_segment_indices;
  std::vector<std::vector<SegmentIndex>> split_segment_indices;
  std::vector<std::vector<NodeIndex>> intersection_nodes;
  split_segments_at_intersections(intersection_data,
                                  num_intersections,
                                  to_array,
                                  out_array,
                                  m_segments,
                                  m_nodes,
                                  original_segment_indices,
                                  split_segment_indices,
                                  intersection_nodes);
  assert(is_valid_curve_data(to_array, out_array));

  for (size_t i = 0; i < intersection_nodes.size(); ++i) {
    spdlog::trace(
      "Intersection {}: {}", i, formatted_vector(intersection_nodes[i]));
    if ((intersection_nodes[i].size() != 2) &&
        (intersection_nodes[i].size() != 0)) {
      SPDLOG_WARN("Intersection {} does not have two nodes: {}",
                  i,
                  formatted_vector(intersection_nodes[i]));
    }
  }

  // Link intersection nodes
  std::vector<NodeIndex> intersection_array(out_array.size(), -1);
  connect_segment_intersections(m_segments,
                                intersection_data,
                                intersection_nodes,
                                to_array,
                                out_array,
                                split_segment_indices,
                                intersection_array,
                                m_nodes);
  assert(is_valid_minimal_curve_network_data(
    to_array, out_array, intersection_array));

  // Further split segments at cusps
  split_segments_at_cusps(interior_cusps,
                          original_segment_indices,
                          split_segment_indices,
                          to_array,
                          out_array,
                          intersection_array,
                          m_segments,
                          m_nodes);
  assert(is_valid_minimal_curve_network_data(
    to_array, out_array, intersection_array));

  // Rebuild topology with intersection and cusp splits
  update_topology(to_array, out_array, intersection_array);

  // Record chain start points
  init_chain_start_nodes();

  // Check validity
  for (NodeIndex node_index = 0; node_index < num_nodes(); ++node_index) {
    if ((is_intersection_node(node_index)) &&
        (!is_valid_node_index(intersection(node_index)))) {
      spdlog::error("Intersection node {} does not have a valid intersection",
                    node_index);
    }
  }
}

// Add all special nodes except the path end nodes to the list of chain start
// nodes
// WARNING: This method is a little dangerous; it modifies the segments as it
// iterates over them
void
ProjectedCurveNetwork::init_chain_start_nodes()
{
  size_t num_nodes = m_nodes.size();

  // Get all nodes that are special (and not path end nodes)
  m_chain_start_nodes.clear();
  for (size_t ni = 0; ni < num_nodes; ++ni) {
    if ((!m_nodes[ni].is_knot()) && (!m_nodes[ni].is_path_end_node())) {
      if (out(ni) < 0)
        continue; // Hack to skip intersection path end nodes
      m_chain_start_nodes.push_back(ni);
    }
  }

  // Mark any missed chains
  bool all_nodes_covered = false;
  while (!all_nodes_covered) {
    // Get list of all covered nodes
    std::vector<bool> is_covered_node(num_nodes, false);
    for (size_t i = 0; i < m_chain_start_nodes.size(); ++i) {
      AbstractCurveNetwork::NodeIndex ni = m_chain_start_nodes[i];
      is_covered_node[ni] = true;
      AbstractCurveNetwork::SegmentIndex start_si = out(ni);
      if (!is_valid_segment_index(start_si)) {
        spdlog::error("Start node is an end point");
        return;
      }

      // Check chain from start node
      ProjectedCurveNetwork::SegmentChainIterator iter =
        get_segment_chain_iterator(start_si);
      for (; !iter.at_end_of_chain(); ++iter) {
        AbstractCurveNetwork::SegmentIndex si = *iter;
        is_covered_node[to(si)] = true;
      }
    }

    // Determine if all nodes covered
    all_nodes_covered = true;
    for (size_t ni = 0; ni < num_nodes; ++ni) {
      if (!is_covered_node[ni]) {
        spdlog::debug("Marking node {} on closed featureless contour", ni);
        m_nodes[ni].mark_as_marked_knot();
        m_chain_start_nodes.push_back(ni);
        all_nodes_covered = false;
        break;
      }
    }
  }
}

ProjectedCurveNetwork::SegmentChainIterator
ProjectedCurveNetwork::get_segment_chain_iterator(
  SegmentIndex segment_index) const
{
  assert(is_valid_segment_index(segment_index));
  return SegmentChainIterator(*this, segment_index);
}

void
ProjectedCurveNetwork::discretize_segment_chain(
  SegmentChainIterator& iter,
  std::vector<PlanarPoint>& points,
  std::vector<int>& polyline) const
{
  CurveDiscretizationParameters curve_disc_params;
  if (iter.at_end_of_chain())
    return;

  // Start the chain polyline
  points.clear();
  std::vector<int> start_polyline;
  RationalFunction<4, 2> const& planar_curve_start =
    m_segments[*iter].get_planar_curve();
  planar_curve_start.discretize(curve_disc_params, points, start_polyline);

  // Add other chain segment points, skipping shared endpoints
  while (!iter.at_end_of_chain()) {
    RationalFunction<4, 2> const& planar_curve_segment =
      m_segments[*iter].get_planar_curve();
    std::vector<PlanarPoint> segment_points;
    std::vector<int> segment_polyline;
    planar_curve_segment.discretize(
      curve_disc_params, segment_points, segment_polyline);
    for (size_t i = 1; i < segment_points.size(); ++i) {
      points.push_back(segment_points[i]);
    }
    iter++;
  }

  // Build trivial polyline
  polyline.resize(points.size());
  for (size_t l = 0; l < points.size(); ++l) {
    polyline[l] = l;
  }
}

void
ProjectedCurveNetwork::discretize_curve(
  NodeIndex start_node_index,
  NodeIndex end_node_index,
  const CurveDiscretizationParameters& curve_disc_params,
  std::vector<PlanarPoint>& points,
  std::vector<int>& polyline,
  std::vector<bool>& is_cusp) const
{

  if (!is_valid_node_index(start_node_index))
    return;
  if (!is_valid_node_index(end_node_index))
    return;

  // Get the points of the curve
  SegmentIndex node_index = start_node_index;
  points.clear();
  is_cusp.clear();
  points.push_back(node_planar_point(start_node_index));
  is_cusp.push_back(
    (is_interior_cusp_node(node_index) || (is_boundary_cusp_node(node_index))));
  while (true) {
    // Get point of the outgoing segment
    SegmentIndex segment_index = out(node_index);
    if (!is_valid_segment_index(segment_index))
      break;
    RationalFunction<4, 2> const& planar_curve_segment =
      m_segments[segment_index].get_planar_curve();
    std::vector<PlanarPoint> segment_points;
    std::vector<int> segment_polyline;
    planar_curve_segment.discretize(
      curve_disc_params, segment_points, segment_polyline);
    for (size_t i = 1; i < segment_points.size(); ++i) {
      points.push_back(segment_points[i]);
      is_cusp.push_back(false);
    }

    // Get the next node in the curve, check if it's a cusp and break if it's
    // the end
    node_index = to(segment_index);
    is_cusp.back() = (is_interior_cusp_node(node_index) ||
                      (is_boundary_cusp_node(node_index)));
    if (node_index == end_node_index)
      break;
  }

  // Build trivial polyline
  polyline.resize(points.size());
  for (size_t l = 0; l < points.size(); ++l) {
    polyline[l] = l;
  }

  // If the segment is closed, make the loop closed
  if (start_node_index == end_node_index) {
    points.pop_back();
    is_cusp.pop_back();
    polyline.back() = polyline.front();
  }
}

// Check if next/prev are consistent
bool
is_valid_next_prev_pair(
  const std::vector<int>& next,
  const std::vector<int>& prev
) {
  int next_size = next.size();
  int prev_size = prev.size();

  // Check for consistent sizes
  if (next_size != prev_size)
  {
    spdlog::error("Inconsistent prev/next sizes");
    return false;
  }

  for (int i = 0; i < next_size; ++i)
  {
    // Check prev[next] is the identity where it is defined
    if ((next[i] != -1) && (prev[next[i]] != i))
    {
      spdlog::error("prev[next] is not the identity");
      return false;
    }

    // Check next[prev] is the identity where it is defined
    if ((prev[i] != -1) && (next[prev[i]] != i))
    {
      spdlog::error("next[prev] is not the identity");
      return false;
    }
  }

  return true;
}

void
ProjectedCurveNetwork::simplify_curves(
  const std::vector<NodeIndex>& visible_curve_start_nodes,
  const std::vector<NodeIndex>& visible_curve_end_nodes,
  const std::vector<std::vector<PlanarPoint>>& all_points,
  std::vector<std::vector<PlanarPoint>>& simplified_points,
  std::vector<std::vector<int>>& simplified_polylines) const
{
  int num_curves = visible_curve_start_nodes.size();

  // Do an O(n^2) search for joined curves
  std::vector<int> prev(num_curves, -1);
  std::vector<int> next(num_curves, -1);
  for (int i = 0; i < num_curves; ++i) {
    NodeIndex end_node_i = visible_curve_end_nodes[i];
    for (int j = 0; j < num_curves; ++j) {
      // Don't check for self overlap
      if (i == j)
        continue;

      NodeIndex start_node_j = visible_curve_start_nodes[j];

      // Check if the two curves are adjacent up to some threshold
      PlanarPoint pi = node_planar_point(end_node_i);
      PlanarPoint pj = node_planar_point(start_node_j);
      PlanarPoint difference = pj - pi;
      if (difference.dot(difference) > 1e-7)
        continue;

      // Get next visible node after the end node i
      int node_index = end_node_i;
      while (true) {
        // Get outgoing segment
        SegmentIndex out_segment_index = out(node_index);

        // Add node if it is a terminal endpoint
        if ((!is_valid_segment_index(out_segment_index)) ||
            (get_segment_quantitative_invisibility(out_segment_index) == 0)) {
          break;
        }

        // Get next node (always valid)
        node_index = to(out_segment_index);

        // Check if closed loop
        if (node_index == end_node_i) {
          break;
        }
      }

      // Check if the two curves are adjacent in the connectivity
      if (node_index != start_node_j)
        continue;

      // Mark connectivity
      next[i] = j;
      prev[j] = i;
    }
  }

#if CHECK_VALIDITY
  // Check if the connectivity has consistent next/prev pairs
  if (!is_valid_next_prev_pair(next, prev))
  {
    spdlog::error("Invalid next/prev pair found");
    simplified_points.clear();
    simplified_polylines.clear();
    return;
  }
#endif

  // Combine the sorted lines
  std::vector<bool> covered(num_curves, false);
  for (int i = 0; i < num_curves; ++i) {
    // Skip already covered curves
    if (covered[i])
      continue;

    // Go back until the start of the chain is found
    int start_index = i;
    SPDLOG_INFO("Processing start index {}", start_index);
    while (prev[start_index] != -1) {
      start_index = prev[start_index];

      // Avoid infinite loop
      if (start_index == i)
        break;
    }

    // Iterate until end of chain found
    std::vector<PlanarPoint> points(0);
    std::vector<int> polyline(0);
    int current_index = start_index;
    int prev_index;
    do {
      covered[current_index] = true;
      for (size_t j = 1; j < all_points[current_index].size(); ++j) {
        polyline.push_back(points.size());
        points.push_back(all_points[current_index][j - 1]);
      }
      prev_index = current_index;
      current_index = next[current_index];
      SPDLOG_INFO("Processing index {}", current_index);

      // Break if current index is invalid or check for consistency
      if (current_index == -1) break;
      PlanarPoint end_point = all_points[prev_index].back();
      PlanarPoint start_point = all_points[current_index].front();
      PlanarPoint difference = start_point - end_point;
      if (difference.dot(difference) > 2e-4)
      {
        spdlog::error("Points {} and {} are distant", end_point, start_point);
      }

    } while (current_index != start_index);

    // Close loop or add last point
    if (current_index == start_index) {
      polyline.push_back(0);
    } else {
      polyline.push_back(points.size());
      points.push_back(all_points[prev_index].back());
    }

    // Add polyline to the global list
    simplified_points.push_back(points);
    simplified_polylines.push_back(polyline);
  }
}

void
ProjectedCurveNetwork::write_uniform_segments(svg::SVG& svgWriter) const
{
  CurveDiscretizationParameters curve_disc_params;
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    RationalFunction<4, 2> const& planar_curve =
      m_segments[i].get_planar_curve();
    write_planar_curve_segment(
      planar_curve, curve_disc_params, svgWriter, 800, 400);
  }
}

void
ProjectedCurveNetwork::write_uniform_visible_segments(svg::SVG& svgWriter) const
{
  CurveDiscretizationParameters curve_disc_params;
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    if (m_segments[i].get_quantitative_invisibility() > 0)
      continue;
    RationalFunction<4, 2> const& planar_curve =
      m_segments[i].get_planar_curve();
    write_planar_curve_segment(
      planar_curve, curve_disc_params, svgWriter, 800, 400);
  }
}

void
ProjectedCurveNetwork::write_contrast_invisible_segments(
  svg::SVG& svgWriter) const
{
  CurveDiscretizationParameters curve_disc_params;
  Color invisible_color = Color(1, 0, 0, 1);
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    RationalFunction<4, 2> const& planar_curve =
      m_segments[i].get_planar_curve();
    if (m_segments[i].get_quantitative_invisibility() > 0) {
      write_planar_curve_segment(
        planar_curve, curve_disc_params, svgWriter, 800, 400, invisible_color);
    } else {
      write_planar_curve_segment(
        planar_curve, curve_disc_params, svgWriter, 800, 400);
    }
  }
}

void
ProjectedCurveNetwork::write_random_chains(svg::SVG& svgWriter) const
{
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    // Get random chain color
    Eigen::Matrix<double, 3, 1> color_without_alpha = generate_random_color();
    Color color(color_without_alpha[0],
                color_without_alpha[1],
                color_without_alpha[2],
                1.0);

    // Write chain
    NodeIndex start_node_index = chain_start_nodes[i];
    SegmentIndex start_segment_index = out(start_node_index);
    SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    discretize_segment_chain(iter, points, polyline);
    add_curve_to_svg(points, polyline, svgWriter, 800, 400, color);
  }
}

void
ProjectedCurveNetwork::write_uniform_chains(svg::SVG& svgWriter) const
{
  Color color(0, 0, 0, 1);
  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    NodeIndex start_node_index = chain_start_nodes[i];
    SegmentIndex start_segment_index = out(start_node_index);
    SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    discretize_segment_chain(iter, points, polyline);
    add_curve_to_svg(points, polyline, svgWriter, 800, 400, color);
  }
}

void
ProjectedCurveNetwork::write_uniform_visible_chains(svg::SVG& svgWriter) const
{
  Color color(0, 0, 0, 1);

  std::vector<int> const& chain_start_nodes = get_chain_start_nodes();
  for (size_t i = 0; i < chain_start_nodes.size(); ++i) {
    NodeIndex start_node_index = chain_start_nodes[i];
    SegmentIndex start_segment_index = out(start_node_index);
    if (m_segments[start_segment_index].get_quantitative_invisibility() > 0)
      continue;
    SegmentChainIterator iter = get_segment_chain_iterator(start_segment_index);
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    discretize_segment_chain(iter, points, polyline);
    add_curve_to_svg(points, polyline, svgWriter, 800, 400, color);
  }
}

void
ProjectedCurveNetwork::write_uniform_visible_curves(svg::SVG& svgWriter) const
{
  Color color(0, 0, 0, 1);
  CurveDiscretizationParameters curve_disc_params;

  // Get visible curves indices
  std::vector<NodeIndex> visible_curve_start_nodes;
  std::vector<NodeIndex> visible_curve_end_nodes;
  get_visible_curves(visible_curve_start_nodes, visible_curve_end_nodes);
  spdlog::info("{} start points and {} end points found",
               visible_curve_start_nodes.size(),
               visible_curve_end_nodes.size());
  SPDLOG_INFO("Start nodes are {}",
              formatted_vector(visible_curve_start_nodes, ", "));
  SPDLOG_INFO("End nodes are {}",
              formatted_vector(visible_curve_end_nodes, ", "));

  // Write all curves
  for (size_t i = 0; i < visible_curve_start_nodes.size(); ++i) {
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    std::vector<bool> is_cusp;
    discretize_curve(visible_curve_start_nodes[i],
                     visible_curve_end_nodes[i],
                     curve_disc_params,
                     points,
                     polyline,
                     is_cusp);
    add_curve_to_svg(points, polyline, svgWriter, 800, 400, color);
  }
}

void
ProjectedCurveNetwork::write_uniform_simplified_visible_curves(
  svg::SVG& svgWriter) const
{
  Color color(0, 0, 0, 1);
  CurveDiscretizationParameters curve_disc_params;

  // Get visible curves indices
  std::vector<NodeIndex> visible_curve_start_nodes;
  std::vector<NodeIndex> visible_curve_end_nodes;
  get_visible_curves(visible_curve_start_nodes, visible_curve_end_nodes);
  spdlog::info("{} start points and {} end points found",
               visible_curve_start_nodes.size(),
               visible_curve_end_nodes.size());

  // Get all curves
  std::vector<std::vector<PlanarPoint>> all_points;
  for (size_t i = 0; i < visible_curve_start_nodes.size(); ++i) {
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    std::vector<bool> is_cusp;
    discretize_curve(visible_curve_start_nodes[i],
                     visible_curve_end_nodes[i],
                     curve_disc_params,
                     points,
                     polyline,
                     is_cusp);
    all_points.push_back(points);
  }

  // Simplify curves
  std::vector<std::vector<PlanarPoint>> simplified_points;
  std::vector<std::vector<int>> simplified_polylines;
  simplify_curves(visible_curve_start_nodes,
                  visible_curve_end_nodes,
                  all_points,
                  simplified_points,
                  simplified_polylines);

  // Write curves
  for (size_t i = 0; i < simplified_polylines.size(); ++i) {
    add_curve_to_svg(simplified_points[i],
                     simplified_polylines[i],
                     svgWriter,
                     800,
                     400,
                     color);
  }
}

void
ProjectedCurveNetwork::write_uniform_closed_curves(svg::SVG& svgWriter) const
{
  Color color(0, 0, 0, 1);
  CurveDiscretizationParameters curve_disc_params;

  // Get closed curves indices
  std::vector<NodeIndex> closed_curve_start_nodes;
  std::vector<NodeIndex> closed_curve_end_nodes;
  get_closed_curves(closed_curve_start_nodes, closed_curve_end_nodes);
  spdlog::info("{} start points and {} end points found",
               closed_curve_start_nodes.size(),
               closed_curve_end_nodes.size());
  SPDLOG_INFO("Start nodes are {}",
              formatted_vector(closed_curve_start_nodes, ", "));
  SPDLOG_INFO("End nodes are {}",
              formatted_vector(closed_curve_end_nodes, ", "));

  // Write all curves
  for (size_t i = 0; i < closed_curve_start_nodes.size(); ++i) {
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    std::vector<bool> is_cusp;
    discretize_curve(closed_curve_start_nodes[i],
                     closed_curve_end_nodes[i],
                     curve_disc_params,
                     points,
                     polyline,
                     is_cusp);
    add_curve_to_svg(points, polyline, svgWriter, 800, 400, color);
  }
}

void
ProjectedCurveNetwork::serialize_closed_curves(
  const std::string& filename) const
{
  int prec = 17;
  std::ofstream out(filename, std::ios::out | std::ios::trunc);
  CurveDiscretizationParameters curve_disc_params;
  curve_disc_params.num_samples = 25;

  // Get closed curve nodes
  std::vector<NodeIndex> closed_curve_start_nodes;
  std::vector<NodeIndex> closed_curve_end_nodes;
  get_closed_curves(closed_curve_start_nodes, closed_curve_end_nodes);

  // Serialize closed curves
  for (size_t i = 0; i < closed_curve_start_nodes.size(); ++i) {
    // Discretize the curve
    std::vector<PlanarPoint> points;
    std::vector<int> polyline;
    std::vector<bool> is_cusp;
    discretize_curve(closed_curve_start_nodes[i],
                     closed_curve_end_nodes[i],
                     curve_disc_params,
                     points,
                     polyline,
                     is_cusp);
    // view_polygon(points, polyline, is_cusp);

    // Write the number of points in the curve
    out << "n " << points.size() << std::endl;

    // Write the curve points with cusp annotation to file
    for (size_t j = 0; j < points.size(); ++j) {
      out << (is_cusp[j] ? 1 : 0) << " ";
      out << std::setprecision(prec) << points[j][0] << " " << points[j][1]
          << std::endl;
    }
  }
}

void
ProjectedCurveNetwork::write(const std::string& output_path,
                             SVGOutputMode color_mode,
                             bool show_nodes) const
{
  Color color(0, 0, 0, 1);
  Eigen::Vector2i viewport = Eigen::Vector2i(800, 800);
  svg::SVG svgWriter(output_path, viewport);

  // Write curves
  if (color_mode == SVGOutputMode::uniform_segments) {
    write_uniform_segments(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_visible_segments) {
    write_uniform_visible_segments(svgWriter);
  } else if (color_mode == SVGOutputMode::contrast_invisible_segments) {
    write_contrast_invisible_segments(svgWriter);
  } else if (color_mode == SVGOutputMode::random_chains) {
    write_random_chains(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_chains) {
    write_uniform_chains(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_visible_chains) {
    write_uniform_visible_chains(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_closed_curves) {
    write_uniform_closed_curves(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_visible_curves) {
    write_uniform_visible_curves(svgWriter);
  } else if (color_mode == SVGOutputMode::uniform_simplified_visible_curves) {
    write_uniform_simplified_visible_curves(svgWriter);
  }

  // Write nodes
  if (show_nodes) {
    // Write nodes
    for (NodeIndex i = 0; i < num_nodes(); ++i) {
      PlanarPoint node_point = node_planar_point(i);
      if (is_boundary_cusp_node(i)) {
        spdlog::info("Writing boundary cusp node at {}", node_point);
        write_planar_point(
          node_point, svgWriter, 800, 400, Color(0, 0, 0.545, 1) // Blue
        );
      } else if (is_interior_cusp_node(i)) {
        spdlog::info("Writing interior cusp node at {}", node_point);
        write_planar_point(
          node_point, svgWriter, 800, 400, Color(0.537, 0.671, 0.890, 1) // Light blue
        );
      } else if (is_intersection_node(i)) {
        spdlog::info("Writing intersection node at {}", node_point);
        write_planar_point(
          node_point, svgWriter, 800, 400, Color(0.227, 0.420, 0.208, 1) // Green
        );
      }
    }
  }
}

std::vector<int> const&
ProjectedCurveNetwork::get_chain_start_nodes() const
{
  return m_chain_start_nodes;
}

void
ProjectedCurveNetwork::enumerate_quantitative_invisibility(
  std::vector<int>& quantitative_invisibility) const
{
  quantitative_invisibility.resize(num_segments());
  for (SegmentIndex i = 0; i < num_segments(); ++i) {
    quantitative_invisibility[i] =
      m_segments[i].get_quantitative_invisibility();
  }
}

// Add planar curves to the polyscope viewer
void
ProjectedCurveNetwork::add_planar_segments_to_viewer() const
{
  // Get planar curve list
  std::vector<RationalFunction<4, 2>> planar_curves;
  enumerate_planar_curves(planar_curves);
  CurveDiscretizationParameters curve_disc_params;

  // Discretize the planar curves
  std::vector<PlanarPoint> points;
  std::vector<std::vector<int>> polylines;
  discretize_curve_segments<4, 2>(
    planar_curves, curve_disc_params, points, polylines);
  MatrixXr points_mat = convert_nested_vector_to_matrix(points);
  std::vector<std::array<int, 2>> edges = convert_polylines_to_edges(polylines);
  polyscope::registerCurveNetwork2D(
    "planar_curve_network_segments", points_mat, edges);
  polyscope::getCurveNetwork("planar_curve_network_segments")
    ->setRadius(0.001);
  polyscope::getCurveNetwork("planar_curve_network_segments")
    ->setColor(glm::vec3(0.600, 0.000, 0.067));

  // Add QI values
  std::vector<int> quantitative_invisibility;
  enumerate_quantitative_invisibility(quantitative_invisibility);
  int num_segments = quantitative_invisibility.size();
  std::vector<int> QI_labels(curve_disc_params.num_samples * num_segments);
  for (int i = 0; i < num_segments; ++i) {
    for (int j = 0; j < curve_disc_params.num_samples; ++j) {
      QI_labels[curve_disc_params.num_samples * i + j] =
        quantitative_invisibility[i];
    }
  }
  MatrixXr colormap;
  generate_random_category_colormap(QI_labels, colormap);
  polyscope::getCurveNetwork("planar_curve_network_segments")
    ->addNodeColorQuantity("quantitative_invisibility", colormap);

  // Get visible curves separately
  std::vector<RationalFunction<4, 2>> visible_planar_curves;
  visible_planar_curves.reserve(planar_curves.size());
  for (size_t i = 0; i < planar_curves.size(); ++i) {
    if (quantitative_invisibility[i] == 0) {
      visible_planar_curves.push_back(planar_curves[i]);
    }
  }

  // Discretize the planar curves
  std::vector<PlanarPoint> visible_points;
  std::vector<std::vector<int>> visible_polylines;
  discretize_curve_segments<4, 2>(visible_planar_curves,
                                  curve_disc_params,
                                  visible_points,
                                  visible_polylines);
  MatrixXr visible_points_mat = convert_nested_vector_to_matrix(visible_points);
  std::vector<std::array<int, 2>> visible_edges =
    convert_polylines_to_edges(visible_polylines);
  polyscope::registerCurveNetwork2D(
    "visible_planar_curves", visible_points_mat, visible_edges);
  polyscope::getCurveNetwork("visible_planar_curves")->setRadius(0.0025);
  polyscope::getCurveNetwork("visible_planar_curves")
    ->setColor(glm::vec3(0.0, 0.0, 0.0));
}

// Add spatial curves to the polyscope viewer
void
ProjectedCurveNetwork::add_spatial_segments_to_viewer() const
{
  // Get spatial curve list
  std::vector<RationalFunction<4, 3>> spatial_curves;
  enumerate_spatial_curves(spatial_curves);
  CurveDiscretizationParameters curve_disc_params;

  // Discretize the spatial curves
  std::vector<SpatialVector> points;
  std::vector<std::vector<int>> polylines;
  discretize_curve_segments<4, 3>(
    spatial_curves, curve_disc_params, points, polylines);
  MatrixXr points_mat = convert_nested_vector_to_matrix(points);
  std::vector<std::array<int, 2>> edges = convert_polylines_to_edges(polylines);
  polyscope::registerCurveNetwork(
    "spatial_curve_network_segments", points_mat, edges);
  polyscope::getCurveNetwork("spatial_curve_network_segments")
    ->setRadius(0.001);
  polyscope::getCurveNetwork("spatial_curve_network_segments")
    ->setColor(glm::vec3(0.600, 0.000, 0.067));

  // Add QI values
  std::vector<int> quantitative_invisibility;
  enumerate_quantitative_invisibility(quantitative_invisibility);
  int num_segments = quantitative_invisibility.size();
  std::vector<int> QI_labels(curve_disc_params.num_samples * num_segments);
  for (int i = 0; i < num_segments; ++i) {
    for (int j = 0; j < curve_disc_params.num_samples; ++j) {
      QI_labels[curve_disc_params.num_samples * i + j] =
        quantitative_invisibility[i];
    }
  }
  MatrixXr colormap;
  generate_random_category_colormap(QI_labels, colormap);
  polyscope::getCurveNetwork("spatial_curve_network_segments")
    ->addNodeColorQuantity("quantitative_invisibility", colormap);

  // Get visible curves separately
  std::vector<RationalFunction<4, 3>> visible_spatial_curves;
  visible_spatial_curves.reserve(spatial_curves.size());
  for (size_t i = 0; i < spatial_curves.size(); ++i) {
    if (quantitative_invisibility[i] == 0) {
      visible_spatial_curves.push_back(spatial_curves[i]);
    }
  }

  // Discretize the spatial curves
  std::vector<SpatialVector> visible_points;
  std::vector<std::vector<int>> visible_polylines;
  discretize_curve_segments<4, 3>(visible_spatial_curves,
                                  curve_disc_params,
                                  visible_points,
                                  visible_polylines);
  MatrixXr visible_points_mat = convert_nested_vector_to_matrix(visible_points);
  std::vector<std::array<int, 2>> visible_edges =
    convert_polylines_to_edges(visible_polylines);
  polyscope::registerCurveNetwork(
    "visible_spatial_curves", visible_points_mat, visible_edges);
  polyscope::getCurveNetwork("visible_spatial_curves")->setRadius(0.0025);
  polyscope::getCurveNetwork("visible_spatial_curves")
    ->setColor(glm::vec3(0.0, 0.0, 0.0));
}

// Add spatial curves to the polyscope viewer
void
ProjectedCurveNetwork::add_spatial_nodes_to_viewer() const
{
  // Get all spatial nodes
  std::vector<SpatialVector> spatial_knot_nodes;
  std::vector<SpatialVector> spatial_marked_knot_nodes;
  std::vector<SpatialVector> spatial_intersection_nodes;
  std::vector<SpatialVector> spatial_interior_cusp_nodes;
  std::vector<SpatialVector> spatial_boundary_cusp_nodes;
  std::vector<SpatialVector> spatial_path_start_nodes;
  std::vector<SpatialVector> spatial_path_end_nodes;
  enumerate_spatial_nodes(spatial_knot_nodes,
                          spatial_marked_knot_nodes,
                          spatial_intersection_nodes,
                          spatial_interior_cusp_nodes,
                          spatial_boundary_cusp_nodes,
                          spatial_path_start_nodes,
                          spatial_path_end_nodes);

  // Get spatial tangents
  std::vector<SpatialVector> spatial_interior_cusp_in_tangents;
  std::vector<SpatialVector> spatial_interior_cusp_out_tangents;
  std::vector<SpatialVector> spatial_boundary_cusp_in_tangents;
  std::vector<SpatialVector> spatial_boundary_cusp_out_tangents;
  enumerate_cusp_spatial_tangents(spatial_interior_cusp_in_tangents,
                                  spatial_interior_cusp_out_tangents,
                                  spatial_boundary_cusp_in_tangents,
                                  spatial_boundary_cusp_out_tangents);

  // Register all spatial nodes
  polyscope::registerPointCloud("spatial_knot_nodes", spatial_knot_nodes);
  polyscope::getPointCloud("spatial_knot_nodes")->setEnabled(false);

  polyscope::registerPointCloud("spatial_marked_knot_nodes",
                                spatial_marked_knot_nodes);
  polyscope::getPointCloud("spatial_marked_knot_nodes")->setEnabled(false);

  polyscope::registerPointCloud("spatial_intersection_nodes",
                                spatial_intersection_nodes);
  polyscope::getPointCloud("spatial_intersection_nodes")
    ->setPointColor(glm::vec3(0.227, 0.420, 0.208));

  polyscope::registerPointCloud("spatial_interior_cusp_nodes",
                                spatial_interior_cusp_nodes);
  polyscope::getPointCloud("spatial_interior_cusp_nodes")
    ->setPointColor(glm::vec3(0.537, 0.671, 0.890));
  polyscope::getPointCloud("spatial_interior_cusp_nodes")
    ->addVectorQuantity("in_tangents", spatial_interior_cusp_in_tangents);
  polyscope::getPointCloud("spatial_interior_cusp_nodes")
    ->addVectorQuantity("out_tangents", spatial_interior_cusp_out_tangents);

  polyscope::registerPointCloud("spatial_boundary_cusp_nodes",
                                spatial_boundary_cusp_nodes);
  polyscope::getPointCloud("spatial_boundary_cusp_nodes")
    ->setPointColor(glm::vec3(0.0, 0.0, 0.545));
  polyscope::getPointCloud("spatial_boundary_cusp_nodes")
    ->addVectorQuantity("in_tangents", spatial_boundary_cusp_in_tangents);
  polyscope::getPointCloud("spatial_boundary_cusp_nodes")
    ->addVectorQuantity("out_tangents", spatial_boundary_cusp_out_tangents);

  polyscope::registerPointCloud("spatial_path_start_nodes",
                                spatial_path_start_nodes);
  polyscope::getPointCloud("spatial_path_start_nodes")->setEnabled(false);

  polyscope::registerPointCloud("spatial_path_end_nodes",
                                spatial_path_end_nodes);
  polyscope::getPointCloud("spatial_path_end_nodes")->setEnabled(false);
}

void
ProjectedCurveNetwork::clear()
{
  clear_geometry();
}

void
ProjectedCurveNetwork::clear_geometry()
{
  m_segments.clear();
  m_nodes.clear();
  m_chain_start_nodes.clear();
}

bool
ProjectedCurveNetwork::is_valid_projected_curve_network() const
{
  // Check that all segments and nodes are hit by iteration from the start node
  std::vector<bool> is_covered_node(num_nodes(), false);
  std::vector<bool> is_covered_segment(num_segments(), false);
  for (size_t i = 0; i < m_chain_start_nodes.size(); ++i) {
    NodeIndex ni = m_chain_start_nodes[i];
    is_covered_node[ni] = true;
    SegmentIndex start_si = out(ni);
    if (!is_valid_segment_index(start_si)) {
      spdlog::error("Start node is an end point");
      return false;
    }

    // Check chain from start node
    is_covered_segment[start_si] = true;
    SegmentChainIterator iter = get_segment_chain_iterator(start_si);
    for (; !iter.at_end_of_chain(); ++iter) {
      SegmentIndex si = *iter;
      is_covered_segment[si] = true;
      is_covered_node[to(si)] = true;
    }
  }

  size_t num_missed_nodes = 0;
  for (NodeIndex ni = 0; ni < num_nodes(); ++ni) {
    if (!is_covered_node[ni]) {
      num_missed_nodes++;
      spdlog::error("{} node {} is not covered by chain iteration",
                    m_nodes[ni].formatted_node(),
                    ni);
    }
  }

  size_t num_missed_segments = 0;
  for (SegmentIndex si = 0; si < num_segments(); ++si) {
    if (!is_covered_segment[si]) {
      num_missed_segments++;
      spdlog::error("Segment {} is not covered by chain iteration", si);
    }
  }

  if ((num_missed_segments > 0) || (num_missed_nodes > 0)) {
    spdlog::error(
      "Missed {} nodes and {} segments", num_missed_nodes, num_missed_segments);
  }

  return true;
}