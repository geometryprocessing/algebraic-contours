// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"

/// \file abstract_curve_network.h
///
/// Methods to build an curve network from minimal connectivity data.

/// An abstract curve network is a graph representing a finite set of possibly
/// intersecting directed curves. For simplicity, it is assumed that all
/// intersections are either transversal or T-nodes so that at most two curves
/// intersect at a node.
class AbstractCurveNetwork
{
public:
  // Indices
  typedef int NodeIndex;
  typedef int SegmentIndex;

  /// Default empty constructor.
  AbstractCurveNetwork();

  /// Construct the network from the basic topological information.
  ///
  /// @param[in] to_array: array mapping segments to their endpoints
  /// @param[in] out_array: array mapping nodes to their outgoing segment
  /// @param[in] intersection_array: list of intersection nodes
  AbstractCurveNetwork(const std::vector<NodeIndex>& to_array,
                       const std::vector<SegmentIndex>& out_array,
                       const std::vector<NodeIndex>& intersection_array);

  /// Update the basic topological information of the curve network.
  ///
  /// @param[in] to_array: array mapping segments to their endpoints
  /// @param[in] out_array: array mapping nodes to their outgoing segment
  /// @param[in] intersection_array: list of intersection nodes
  void update_topology(
    const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
    const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
    const std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array);

  /// Return the number of segments in the curve network.
  ///
  /// @return number of segments
  SegmentIndex get_num_segments() const { return m_to_array.size(); }

  /// Return the number of nodes in the curve network.
  ///
  /// @return number of nodes
  NodeIndex get_num_nodes() const { return m_out_array.size(); }

  /// Get the next segment after a given segment (or -1 if there is no next
  /// segment)
  ///
  /// @param[in] segment_index: query segment index
  /// @return next segment
  SegmentIndex next(SegmentIndex segment_index) const
  {
    if (!is_valid_segment_index(segment_index))
      return -1;
    return m_next_array[segment_index];
  }

  /// Get the previous segment after a given segment (or -1 if there is no
  /// previous segment)
  ///
  /// @param[in] segment_index: query segment index
  /// @return previous segment
  SegmentIndex prev(SegmentIndex segment_index) const
  {
    if (!is_valid_segment_index(segment_index))
      return -1;
    return m_prev_array[segment_index];
  }

  /// Get the node at the tip of the segment
  ///
  /// Note that this operation is valid for any valid segment
  ///
  /// @param[in] segment_index: query segment index
  /// @return to node of the segment
  NodeIndex to(SegmentIndex segment_index) const
  {
    if (!is_valid_segment_index(segment_index))
      return -1;
    return m_to_array[segment_index];
  }

  /// Get the node at the base of the segment
  ///
  /// Note that this operation is valid for any valid segment
  ///
  /// @param[in] segment_index: query segment index
  /// @return from node of the segment
  NodeIndex from(SegmentIndex segment_index) const
  {
    if (!is_valid_segment_index(segment_index))
      return -1;
    return m_from_array[segment_index];
  }

  /// Get the node that intersects the given node (or -1 if the node does not
  /// intersect another node)
  ///
  /// @param[in] node_index: query node index
  /// @return intersection node of the node
  NodeIndex intersection(NodeIndex node_index) const
  {
    if (!is_valid_node_index(node_index))
      return -1;
    return m_intersection_array[node_index];
  }

  /// Get the outgoing segment for the node (or -1 if none exists)
  ///
  /// @param[in] node_index: query node index
  /// @return out segment of the node
  SegmentIndex out(NodeIndex node_index) const
  {
    if (!is_valid_node_index(node_index))
      return -1;
    return m_out_array[node_index];
  }

  /// Get the incoming segment for the node (or -1 if none exists)
  ///
  /// @param[in] node_index: query node index
  /// @return in segment of the node
  SegmentIndex in(NodeIndex node_index) const
  {
    if (!is_valid_node_index(node_index))
      return -1;
    return m_in_array[node_index];
  }

  /// Determine if the node is on the boundary of a curve in the curve network.
  ///
  /// @param[in] node_index: query node index
  /// @return true iff the given node is a boundary node
  bool is_boundary_node(NodeIndex node_index) const;

  /// Determine if the node has an intersection.
  ///
  /// @param[in] node_index: query node index
  /// @return true iff the given node is an intersection node
  bool has_intersection_node(NodeIndex node_index) const;

  /// Determine if the node is a T-node, i.e., has an intersection and one of
  /// the two is on the boundary.
  ///
  /// Note that this is a weaker condition than having an intersection node and
  /// being and intersection node and is not simply a logical and of the two
  /// conditions.
  ///
  /// @param[in] node_index: query node index
  /// @return true iff the given node is a boundary intersection node
  bool is_tnode(NodeIndex node_index) const;

protected:
  // Initialization
  void init_abstract_curve_network();

  // Validity checks
  bool is_valid_segment_index(SegmentIndex segment_index) const;
  bool is_valid_node_index(NodeIndex node_index) const;
  bool is_valid_abstract_curve_network() const;

  // Clear
  void clear_topology();

private:
  // Segment topology
  std::vector<NodeIndex> m_to_array;
  std::vector<NodeIndex> m_from_array;
  std::vector<SegmentIndex> m_next_array;
  std::vector<SegmentIndex> m_prev_array;

  // Node topology
  std::vector<SegmentIndex> m_out_array;
  std::vector<NodeIndex> m_intersection_array;
  std::vector<SegmentIndex> m_in_array;
};

/// Build from map sending segments to their origin nodes.
///
/// @param[in] to_array: array mapping segments to their endpoints
/// @param[in] out_array: array mapping nodes to their outgoing segment
/// @param[out] out_array: array mapping segments to their origin points
void
build_from_array(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<AbstractCurveNetwork::NodeIndex>& from_array);

/// Check if input has valid indexing, meaning all to nodes are valid
///
/// Note that out may be invalid for some nodes if they are terminal
///
/// @param[in] to_array: array mapping segments to their endpoints
/// @param[in] out_array: array mapping nodes to their outgoing segment
/// @return true iff the curve data is valid
bool
is_valid_curve_data(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array);

/// Check if input describes a valid curve network
///
/// @param[in] to_array: array mapping segments to their endpoints
/// @param[in] out_array: array mapping nodes to their outgoing segment
/// @param[in] intersection_array: list of intersection nodes
/// @return true iff the curve network data is valid
bool
is_valid_minimal_curve_network_data(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  const std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array);
