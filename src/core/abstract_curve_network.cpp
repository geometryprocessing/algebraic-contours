// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "abstract_curve_network.h"

// ****************
// Helper Functions
// ****************

// Build next map from segments to the following segment or -1 if it is terminal
void
build_next_array(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& next_array)
{
  AbstractCurveNetwork::SegmentIndex num_segments = to_array.size();
  next_array.resize(num_segments);
  for (AbstractCurveNetwork::SegmentIndex si = 0; si < num_segments; ++si) {
    next_array[si] = out_array[to_array[si]];
  }
}

// Build prev map from segments to their previous segment or -1 if it is initial
void
build_prev_array(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<AbstractCurveNetwork::SegmentIndex>& prev_array)
{
  AbstractCurveNetwork::SegmentIndex num_segments = to_array.size();

  // Initialize prev array to -1
  prev_array.resize(num_segments);
  std::fill(prev_array.begin(), prev_array.end(), -1);

  // Find previous segments when they exist
  for (AbstractCurveNetwork::SegmentIndex si = 0; si < num_segments; ++si) {
    AbstractCurveNetwork::SegmentIndex next_segment = out_array[to_array[si]];
    if ((next_segment < 0) || (next_segment >= num_segments))
      continue;
    prev_array[next_segment] = si;
  }
}

// Build from map sending segments to their origin nodes.
void
build_from_array(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  std::vector<AbstractCurveNetwork::NodeIndex>& from_array)
{
  AbstractCurveNetwork::SegmentIndex num_segments = to_array.size();
  AbstractCurveNetwork::NodeIndex num_nodes = out_array.size();
  from_array.assign(num_segments, -1);
  for (AbstractCurveNetwork::NodeIndex ni = 0; ni < num_nodes; ++ni) {
    AbstractCurveNetwork::NodeIndex out_segment = out_array[ni];
    if ((out_segment < 0) || (out_segment >= num_segments))
      continue;
    from_array[out_segment] = ni;
  }
}

// Build in map from nodes to incoming segments or -1 if they are initial
void
build_in_array(const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
               const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
               std::vector<AbstractCurveNetwork::SegmentIndex>& in_array)
{
  AbstractCurveNetwork::SegmentIndex num_segments = to_array.size();
  AbstractCurveNetwork::NodeIndex num_nodes = out_array.size();
  in_array.resize(num_nodes);
  std::fill(in_array.begin(), in_array.end(), -1);
  for (AbstractCurveNetwork::SegmentIndex si = 0; si < num_segments; ++si) {
    in_array[to_array[si]] = si;
  }
}

// Check if input has valid indexing
bool
is_valid_curve_data(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array)
{
  AbstractCurveNetwork::SegmentIndex num_segments = to_array.size();
  AbstractCurveNetwork::NodeIndex num_nodes = out_array.size();

  // Check all to nodes are valid
  for (AbstractCurveNetwork::SegmentIndex si = 0; si < num_segments; ++si) {
    if ((to_array[si] < 0) || (to_array[si] >= num_nodes)) {
      spdlog::error("Segment {} is invalid with to node {}", si, to_array[si]);
      return false;
    }
  }

  return true;
}

// Check if input describes a valid curve network
bool
is_valid_minimal_curve_network_data(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  const std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array)
{
  size_t num_segments = to_array.size();
  size_t num_nodes = out_array.size();

  // Check sizes
  if (to_array.size() != num_segments) {
    spdlog::error("to domain not in bijection with number of segments");
    return false;
  }
  if (out_array.size() != num_nodes) {
    spdlog::error("out domain not in bijection with number of nodes");
    return false;
  }
  if (intersection_array.size() != num_nodes) {
    spdlog::error("out domain not in bijection with number of nodes");
    return false;
  }

  // Check all out nodes are valid (to and intersection array can have invalid
  // nodes)
  if (!is_valid_curve_data(to_array, out_array)) {
    return false;
  }

  return true;
}

// **********************
// Abstract Curve Network
// **********************

// **************
// Public methods
// **************

AbstractCurveNetwork::AbstractCurveNetwork()
{
  // Clear all information
  clear_topology();
}

AbstractCurveNetwork::AbstractCurveNetwork(
  const std::vector<NodeIndex>& to_array,
  const std::vector<SegmentIndex>& out_array,
  const std::vector<NodeIndex>& intersection_array)
  : m_to_array(to_array)
  , m_out_array(out_array)
  , m_intersection_array(intersection_array)
{
  // Check input validity
#if CHECK_VALIDITY
  if (!is_valid_minimal_curve_network_data(
        to_array, out_array, intersection_array)) {
    spdlog::error("Could not build abstract curve network");
    clear_topology();
    return;
  }
#endif

  // Build curve network
  init_abstract_curve_network();

  // Check validity
#if CHECK_VALIDITY
  if (!is_valid_abstract_curve_network()) {
    spdlog::error("Inconsistent abstract curve network built");
    clear_topology();
    return;
  }
#endif
}

void
AbstractCurveNetwork::update_topology(
  const std::vector<AbstractCurveNetwork::NodeIndex>& to_array,
  const std::vector<AbstractCurveNetwork::SegmentIndex>& out_array,
  const std::vector<AbstractCurveNetwork::NodeIndex>& intersection_array)
{
  // Check input validity
  assert(is_valid_minimal_curve_network_data(
    to_array, out_array, intersection_array));
  if (!is_valid_minimal_curve_network_data(
        to_array, out_array, intersection_array)) {
    spdlog::error("Could not build abstract curve network");
    clear_topology();
    return;
  }

  // Set input arrays
  m_to_array = to_array;
  m_out_array = out_array;
  m_intersection_array = intersection_array;

  // Build curve network
  init_abstract_curve_network();

  // Check validity
  if (!is_valid_abstract_curve_network()) {
    spdlog::error("Inconsistent abstract curve network built");
    clear_topology();
    return;
  }
}

bool
AbstractCurveNetwork::is_boundary_node(NodeIndex node_index) const
{
  // Invalid nodes are not boundary nodes
  if (!is_valid_node_index(node_index))
    return false;

  // Nodes without in or out segments are on the boundary
  if (!is_valid_segment_index(in(node_index)))
    return true;
  if (!is_valid_segment_index(out(node_index)))
    return true;

  // Otherwise, interior node
  return false;
}

bool
AbstractCurveNetwork::has_intersection_node(NodeIndex node_index) const
{
  // Invalid nodes do not have intersection nodes
  if (!is_valid_node_index(node_index))
    return false;

  // Check if intersection node exists
  return (is_valid_node_index(intersection(node_index)));
}

bool
AbstractCurveNetwork::is_tnode(NodeIndex node_index) const
{
  // Invalid nodes are not T-nodes
  if (!is_valid_node_index(node_index))
    return false;

  // Nodes without intersections are not T-nodes
  if (!has_intersection_node(node_index))
    return false;

  // T-nodes are on the boundary or intersect a boundary node
  if (is_boundary_node(node_index))
    return true;
  if (is_boundary_node(intersection(node_index)))
    return true;

  // Otherwise, interior intersection node
  return false;
}

// *****************
// Protected methods
// *****************

// Implementation of the main constructor
void
AbstractCurveNetwork::init_abstract_curve_network()
{
  // Build next arrays
  build_next_array(m_to_array, m_out_array, m_next_array);
  build_prev_array(m_to_array, m_out_array, m_prev_array);
  build_from_array(m_to_array, m_out_array, m_from_array);
  build_in_array(m_to_array, m_out_array, m_in_array);
}

// Determine if the index describes a segment of the curve network
bool
AbstractCurveNetwork::is_valid_segment_index(
  AbstractCurveNetwork::SegmentIndex segment_index) const
{
  // Ensure in bounds for segment list
  if (segment_index < 0)
    return false;
  if (segment_index >= get_num_segments())
    return false;

  return true;
}

// Determine if the index describes a node of the curve network
bool
AbstractCurveNetwork::is_valid_node_index(
  AbstractCurveNetwork::NodeIndex node_index) const
{
  // Ensure in bounds for node list
  if (node_index < 0)
    return false;
  if (node_index >= get_num_nodes())
    return false;

  return true;
}

// Clear all member data
void
AbstractCurveNetwork::clear_topology()
{
  m_next_array.clear();
  m_prev_array.clear();
  m_to_array.clear();
  m_from_array.clear();
  m_intersection_array.clear();
  m_out_array.clear();
  m_in_array.clear();
}

// ***************
// Private methods
// ***************

// General validity checker for the network topology
bool
AbstractCurveNetwork::is_valid_abstract_curve_network() const
{
  size_t num_segments = get_num_segments();
  size_t num_nodes = get_num_nodes();

  // Array size checks
  if (m_next_array.size() != num_segments) {
    spdlog::error("Inconsistent next array");
    return false;
  }
  if (m_prev_array.size() != num_segments) {
    spdlog::error("Inconsistent prev array");
    return false;
  }
  if (m_to_array.size() != num_segments) {
    spdlog::error("Inconsistent to array");
    return false;
  }
  if (m_from_array.size() != num_segments) {
    spdlog::error("Inconsistent from array");
    return false;
  }
  if (m_intersection_array.size() != num_nodes) {
    spdlog::error("Inconsistent intersection array");
    return false;
  }
  if (m_out_array.size() != num_nodes) {
    spdlog::error("Inconsistent out array");
    return false;
  }
  if (m_in_array.size() != num_nodes) {
    spdlog::error("Inconsistent in array");
    return false;
  }

  // Check segment topology
  for (SegmentIndex si = 0; si < get_num_segments(); ++si) {
    // Check to node
    if (!is_valid_node_index(to(si))) {
      spdlog::error("To does not have a valid endpoint for segment {}", si);
      return false;
    }
    if (in(to(si)) != si) {
      spdlog::error("in(to(s)) is not the identity for segment {}", si);
      return false;
    }

    // Check from node
    if (!is_valid_node_index(from(si))) {
      spdlog::error("From does not have a valid endpoint for segment {}", si);
      return false;
    }
    if (out(from(si)) != si) {
      spdlog::error("out(from(s)) is not the identity for segment {}", si);
      return false;
    }

    // Check next segment is consistent if it exists
    if (is_valid_segment_index(next(si))) {
      if (prev(next(si)) != si) {
        spdlog::error(
          "prev(next(s)) is not the identity for nonterminal segment {}", si);
        spdlog::error("next(s) is {}", next(si));
        return false;
      }
    }
    // Check to node is an endpoint if the next segment does not exist
    else {
      if (is_valid_segment_index(out(to(si)))) {
        spdlog::error("Terminal segment {} does not have a terminal endpoint",
                      si);
        return false;
      }
    }

    // Check prev segment is consistent if it exists
    if (is_valid_segment_index(prev(si))) {
      if (next(prev(si)) != si) {
        spdlog::error(
          "next(prev(s)) is not the identity for noninitial segment {}", si);
        return false;
      }
    }
    // Check to node is an endpoint if the next segment does not exist
    else {
      if (is_valid_segment_index(in(from(si)))) {
        spdlog::error("Initial segment {} does not have a initial start point",
                      si);
        return false;
      }
    }
  }

  // Check node topology
  std::vector<bool> is_out_segment(get_num_segments(), false);
  std::vector<bool> is_in_segment(get_num_segments(), false);
  for (NodeIndex ni = 0; ni < get_num_nodes(); ++ni) {
    // Check the outgoing segment comes from the node if it exists
    if (is_valid_segment_index(out(ni))) {
      if (from(out(ni)) != ni) {
        spdlog::error(
          "from(out(n)) is not the identity for nonterminal node {}", ni);
        return false;
      }
      is_out_segment[out(ni)] = true;
    }

    // Check the incoming segment goes to the node if it exists
    if (is_valid_segment_index(in(ni))) {
      if (to(in(ni)) != ni) {
        spdlog::error("to(in(n)) is not the identity for non ininital node {}",
                      ni);
        return false;
      }
      is_in_segment[in(ni)] = true;
    }

    // Check the intersection is a closed order 2 loop if it exists
    if (is_valid_node_index(intersection(ni))) {
      if (intersection(intersection(ni)) != ni) {
        spdlog::error("Intersection is order 2 for intersection node {}", ni);
        return false;
      }
    }
  }

  // Check all segments originate from some node
  if (vector_contains<bool>(is_out_segment, false)) {
    spdlog::error("Segment does not have a starting node");
    return false;
  }

  // Check all segments go into some node
  if (vector_contains<bool>(is_in_segment, false)) {
    spdlog::error("Segment does not have a terminal node");
    return false;
  }

  return true;
}
