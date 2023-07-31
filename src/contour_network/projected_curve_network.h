// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "abstract_curve_network.h"
#include "common.h"
#include "conic.h"
#include "intersection_data.h"
#include "rational_function.h"
#include "svg.h"

/// \file projected_curve_network.h
///
/// Methods to compute a simple planar curve network from annotated plane curve
/// soup.

enum class SVGOutputMode
{
  uniform_segments,            // All contours in uniform color
  uniform_visible_segments,    // Visible contours in uniform color
  contrast_invisible_segments, // Visible and invisible segments in a different
                               // color
  random_chains,               // Chains in random colors
  uniform_chains,              // Chains in uniform color
  uniform_visible_chains,      // Visible chains in uniform color
  uniform_visible_curves, // Visible curves with no breaks at special points
  uniform_closed_curves,  // All closed curves with no breaks at special points
  uniform_simplified_visible_curves // All visible closed curves with
                                    // simplification

};

class SegmentGeometry
{
public:
  SegmentGeometry() {}

  SegmentGeometry(const Conic& parameter_curve,
                  const RationalFunction<4, 3>& spatial_curve,
                  const RationalFunction<4, 2>& planar_curve,
                  const std::map<std::string, int> segment_labels)
  {
    m_planar_curve = planar_curve;
    m_spatial_curve = spatial_curve;
    m_parameter_curve = parameter_curve;
    m_segment_labels = segment_labels;
    m_quantitative_invisibility = -1;

    if (!is_valid_segment_geometry()) {
      spdlog::error("Invalid segment made");
    }
  }

  // Geometric data
  RationalFunction<4, 2> const& get_planar_curve() const
  {
    return m_planar_curve;
  }

  RationalFunction<4, 3> const& get_spatial_curve() const
  {
    return m_spatial_curve;
  }

  Conic const& get_parameter_curve() const { return m_parameter_curve; }

  void set_segment_label(const std::string& label_name, int new_segment_label)
  {
    m_segment_labels[label_name] = new_segment_label;
  }

  int get_segment_label(const std::string& label_name) const
  {
    return m_segment_labels.at(label_name);
  }

  void set_quantitative_invisibility(int new_quantitative_invisibility)
  {
    if (new_quantitative_invisibility < 0) {
      spdlog::error("Cannot set negative segment QI");
      return;
    }
    m_quantitative_invisibility = new_quantitative_invisibility;
  }

  int get_quantitative_invisibility() const
  {
    return m_quantitative_invisibility;
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const SegmentGeometry& segment_geometry)
  {
    out << "Parameter curve segment: " << segment_geometry.get_parameter_curve()
        << std::endl;
    out << "Spatial curve segment: " << segment_geometry.get_spatial_curve()
        << std::endl;
    out << "Planar curve segment: " << segment_geometry.get_planar_curve()
        << std::endl;

    return out;
  }

  // Split segment into two new segments at a given knot value
  void split_at_knot(double knot,
                     SegmentGeometry& lower_segment,
                     SegmentGeometry& upper_segment)
  {
    // Split all three segment curves at the knot
    Conic lower_parameter_curve;
    RationalFunction<4, 3> lower_spatial_curve;
    RationalFunction<4, 2> lower_planar_curve;
    Conic upper_parameter_curve;
    RationalFunction<4, 3> upper_spatial_curve;
    RationalFunction<4, 2> upper_planar_curve;
    m_parameter_curve.split_at_knot(
      knot, lower_parameter_curve, upper_parameter_curve);
    m_planar_curve.split_at_knot(knot, lower_planar_curve, upper_planar_curve);
    m_spatial_curve.split_at_knot(
      knot, lower_spatial_curve, upper_spatial_curve);

    // Build new segments from the split curves
    lower_segment = SegmentGeometry(lower_parameter_curve,
                                    lower_spatial_curve,
                                    lower_planar_curve,
                                    m_segment_labels);
    upper_segment = SegmentGeometry(upper_parameter_curve,
                                    upper_spatial_curve,
                                    upper_planar_curve,
                                    m_segment_labels);
  }

private:
  RationalFunction<4, 2> m_planar_curve;
  RationalFunction<4, 3> m_spatial_curve;
  Conic m_parameter_curve;
  std::map<std::string, int> m_segment_labels;
  int m_quantitative_invisibility;

  bool is_valid_segment_geometry() const
  {
    double t0 = m_planar_curve.domain().get_lower_bound();
    double t1 = m_planar_curve.domain().get_upper_bound();
    if (m_spatial_curve.domain().get_lower_bound() != t0) {
      spdlog::error("Lower bound error");
      return false;
    }
    if (m_parameter_curve.domain().get_lower_bound() != t0) {
      spdlog::error("Lower bound error");
      return false;
    }
    if (m_spatial_curve.domain().get_upper_bound() != t1) {
      spdlog::error("Upper bound error");
      return false;
    }
    if (m_parameter_curve.domain().get_upper_bound() != t1) {
      spdlog::error("Upper bound error");
      return false;
    }

    return true;
  }
};

class NodeGeometry
{
public:
  NodeGeometry()
  {
    m_node_type = knot;
    m_quantitative_invisibility = -1;
  }

  void mark_as_knot() { m_node_type = knot; }

  void mark_as_marked_knot() { m_node_type = marked_knot; }

  void mark_as_boundary_cusp() { m_node_type = boundary_cusp; }

  void mark_as_interior_cusp() { m_node_type = interior_cusp; }

  void mark_as_intersection() { m_node_type = intersection; }

  void mark_as_path_start_node() { m_node_type = path_start_node; }

  void mark_as_path_end_node() { m_node_type = path_end_node; }

  bool is_knot() const { return (m_node_type == knot); }

  bool is_marked_knot() const { return (m_node_type == marked_knot); }

  bool is_intersection() const { return (m_node_type == intersection); }

  bool is_interior_cusp() const { return (m_node_type == interior_cusp); }

  bool is_boundary_cusp() const { return (m_node_type == boundary_cusp); }

  bool is_path_start_node() const { return (m_node_type == path_start_node); }

  bool is_path_end_node() const { return (m_node_type == path_end_node); }

  void set_quantitative_invisibility(int new_quantitative_invisibility)
  {
    m_qi_set = true;
    if (new_quantitative_invisibility < 0) {
      spdlog::error("Cannot set negative node QI");
      return;
    }
    m_quantitative_invisibility = new_quantitative_invisibility;
  }

  int get_quantitative_invisibility() const
  {
    return m_quantitative_invisibility;
  }

  bool quantitative_invisibility_is_set() const { return m_qi_set; }
  void mark_quantitative_invisibility_as_set() { m_qi_set = true; }

  std::string formatted_node() const
  {
    if (m_node_type == knot)
      return "knot";
    else if (m_node_type == marked_knot)
      return "marked_knot";
    else if (m_node_type == boundary_cusp)
      return "boundary_cusp";
    else if (m_node_type == interior_cusp)
      return "interior_cusp";
    else if (m_node_type == intersection)
      return "intersection";
    else if (m_node_type == path_end_node)
      return "path_end_node";
    else if (m_node_type == path_start_node)
      return "path_start_node";
    else
      return "unknown";
  }

private:
  // Geometric data
  // TODO boundary and interior cusps are no longer clear in the general
  // context. Should rename to something more general like C1 cusp and C0
  // cusp, etc.
  enum
  {
    knot,
    marked_knot,
    boundary_cusp,
    interior_cusp,
    intersection,
    path_end_node,
    path_start_node
  } m_node_type;

  int m_quantitative_invisibility;
  bool m_qi_set = false;
};

/// A projected curve network is a curve network of intersecting planar curves
/// arising from the projection of spatial curves to the xy plane.
class ProjectedCurveNetwork : public AbstractCurveNetwork
{
public:
  /// Default constructor
  ProjectedCurveNetwork();

  /// Construct the curve network from the relevant annotated geometric
  /// information.
  ///
  /// @param[in] parameter_segments: uv domain quadratic curves parametrizing
  ///     the other curves
  /// @param[in] spatial_segments: spatial rational curves before projection
  /// @param[in] planar_segments: planar rational curves after projection
  /// @param[in] chain_labels: list of maps of labels for each segment (e.g.,
  /// patch label)
  /// @param[in] chains: list of lists of chained curve indices
  /// @param[in] chain_labels: list of chain labels for each segment
  /// @param[in] interior_cusps: list of lists of cusp points per segment
  /// @param[in] has_cusp_at_base: list of bools per segment indicating if the
  ///     segment base node is a cusp
  /// @param[in] intersections: list of lists of intersection points per segment
  /// @param[in] intersection_indices: list of lists of indices of curves
  /// corresponding
  ///     to intersection points per segment
  ProjectedCurveNetwork(
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
    int num_intersections);

  // Counts
  AbstractCurveNetwork::SegmentIndex num_segments() const;
  AbstractCurveNetwork::NodeIndex num_nodes() const;

  // Segment geometry
  Conic const& segment_parameter_curve(SegmentIndex segment_index) const;
  RationalFunction<4, 2> const& segment_planar_curve(
    SegmentIndex segment_index) const;
  RationalFunction<4, 3> const& segment_spatial_curve(
    SegmentIndex segment_index) const;
  int get_segment_label(SegmentIndex segment_index,
                        const std::string& label_name) const;
  int get_segment_quantitative_invisibility(SegmentIndex segment_index) const;
  void set_segment_quantitative_invisibility(SegmentIndex segment_index,
                                             int new_quantitative_invisibility);

  /// Get list of all parameter curves in the network
  ///
  /// @param[out] parameter_curves: list of all parameter curves
  void enumerate_parameter_curves(std::vector<Conic>& parameter_curves) const;

  /// Get list of all spatial curves in the network
  ///
  /// @param[out] spatial_curves: list of all spatial curves
  void enumerate_spatial_curves(
    std::vector<RationalFunction<4, 3>>& spatial_curves) const;

  /// Get list of all planar curves in the network
  ///
  /// @param[out] planar_curves: list of all planar curves
  void enumerate_planar_curves(
    std::vector<RationalFunction<4, 2>>& planar_curves) const;

  /// Get list of all annotated spatial nodes in the network
  ///
  /// @param[out] spatial_knot_nodes: list of all spatial knot nodes
  /// @param[out] spatial_marked_knot_nodes: list of all spatial marked knot
  /// nodes
  /// @param[out] spatial_intersection_nodes: list of all spatial intersection
  /// nodes
  /// @param[out] spatial_interior_cusp_nodes: list of all spatial interior cusp
  /// nodes
  /// @param[out] spatial_boundary_cusp_nodes: list of all spatial boundary cusp
  /// nodes
  /// @param[out] spatial_path_start_nodes: list of all spatial path start nodes
  /// @param[out] spatial_path_end_nodes: list of all spatial path end nodes
  void enumerate_spatial_nodes(
    std::vector<SpatialVector>& spatial_knot_nodes,
    std::vector<SpatialVector>& spatial_marked_knot_nodes,
    std::vector<SpatialVector>& spatial_intersection_nodes,
    std::vector<SpatialVector>& spatial_interior_cusp_nodes,
    std::vector<SpatialVector>& spatial_boundary_cusp_nodes,
    std::vector<SpatialVector>& spatial_path_start_nodes,
    std::vector<SpatialVector>& spatial_path_end_nodes) const;

  /// Get list of all annotated planar nodes in the network
  ///
  /// @param[out] planar_knot_nodes: list of all planar knot nodes
  /// @param[out] planar_marked_knot_nodes: list of all planar marked knot nodes
  /// @param[out] planar_intersection_nodes: list of all planar intersection
  /// nodes
  /// @param[out] planar_interior_cusp_nodes: list of all planar interior cusp
  /// nodes
  /// @param[out] planar_boundary_cusp_nodes: list of all planar boundary cusp
  /// nodes
  /// @param[out] planar_path_start_nodes: list of all planar path start nodes
  /// @param[out] planar_path_end_nodes: list of all planar path end nodes
  void enumerate_planar_nodes(
    std::vector<PlanarPoint>& planar_knot_nodes,
    std::vector<PlanarPoint>& planar_marked_knot_nodes,
    std::vector<PlanarPoint>& planar_intersection_nodes,
    std::vector<PlanarPoint>& planar_interior_cusp_nodes,
    std::vector<PlanarPoint>& planar_boundary_cusp_nodes,
    std::vector<PlanarPoint>& planar_path_start_nodes,
    std::vector<PlanarPoint>& planar_path_end_nodes) const;

  void enumerate_cusp_spatial_tangents(
    std::vector<SpatialVector>& spatial_interior_cusp_in_tangents,
    std::vector<SpatialVector>& spatial_interior_cusp_out_tangents,
    std::vector<SpatialVector>& spatial_boundary_cusp_in_tangents,
    std::vector<SpatialVector>& spatial_boundary_cusp_out_tangents) const;

  /// Get list of all quantitative invisibility values in the network
  ///
  /// @param[out] quantitative_invisibility: list of all QI values
  void enumerate_quantitative_invisibility(
    std::vector<int>& quantitative_invisibility) const;

  /// Get list of all start nodes in the network
  ///
  /// @return list of all chain start nodes
  std::vector<int> const& get_chain_start_nodes() const;

  // Node geometry
  bool is_knot_node(NodeIndex node_index) const;
  bool is_marked_knot_node(NodeIndex node_index) const;
  bool is_intersection_node(NodeIndex node_index) const;
  bool is_interior_cusp_node(NodeIndex node_index) const;
  bool is_boundary_cusp_node(NodeIndex node_index) const;
  bool is_path_start_node(NodeIndex node_index) const;
  bool is_path_end_node(NodeIndex node_index) const;
  PlanarPoint node_planar_point(NodeIndex node_index) const;
  SpatialVector node_spatial_point(NodeIndex node_index) const;
  SpatialVector node_spatial_in_tangent(NodeIndex node_index) const;
  SpatialVector node_spatial_out_tangent(NodeIndex node_index) const;
  SpatialVector node_spatial_tangent(NodeIndex node_index) const;
  PlanarPoint node_planar_in_tangent(NodeIndex node_index) const;
  PlanarPoint node_planar_out_tangent(NodeIndex node_index) const;
  PlanarPoint node_planar_tangent(NodeIndex node_index) const;
  int get_node_quantitative_invisibility(NodeIndex node_index) const;
  void set_node_quantitative_invisibility(NodeIndex node_index,
                                          int new_quantitative_invisibility);
  bool node_quantitative_invisibility_is_set(NodeIndex node_index) const
  {
    return m_nodes[node_index].quantitative_invisibility_is_set();
  }
  void mark_node_quantitative_invisibility_as_set(NodeIndex node_index)
  {
    m_nodes[node_index].mark_quantitative_invisibility_as_set();
  }

  /// @brief Iterator over contour segment chains until FIXME is found.
  class SegmentChainIterator
  {
  public:
    SegmentChainIterator(const ProjectedCurveNetwork& parent,
                         SegmentIndex segment_index)
      : m_parent(parent)
    {
      m_current_segment_index = segment_index;
      m_is_end_of_chain = false;
      m_is_reverse_end_of_chain = false;
      spdlog::trace("Iterator initialized to segment {}",
                    m_current_segment_index);
    }

    SegmentChainIterator& operator++()
    {
      // Check if off end of chain
      if (m_is_end_of_chain || (m_is_reverse_end_of_chain))
        return *this;

      // Check if next is end of chain (i.e. the next node is not a knot node)
      NodeIndex to_node = m_parent.to(m_current_segment_index);
      if (!m_parent.is_knot_node(to_node))
      // if
      // (!m_parent.is_valid_segment_index(m_parent.next(m_current_segment_index)))
      {
        m_current_segment_index = -1;
        m_is_end_of_chain = true;
      } else {
        m_current_segment_index = m_parent.next(m_current_segment_index);
      }

      spdlog::trace("Iterator moved to segment {}", m_current_segment_index);

      return *this;
    }

    SegmentChainIterator operator++(int)
    {
      SegmentChainIterator temp = *this;
      ++*this;
      return temp;
    }

    SegmentChainIterator& operator--()
    {
      // Check if off end of chain
      if (m_is_end_of_chain || (m_is_reverse_end_of_chain))
        return *this;

      // Check if prev is reverse end of chain (i.e. the prev node is not a knot
      // node)
      NodeIndex from_node = m_parent.from(m_current_segment_index);
      if (!m_parent.is_knot_node(from_node))
      // if
      // (!m_parent.is_valid_segment_index(m_parent.prev(m_current_segment_index)))
      {
        m_current_segment_index = -1;
        m_is_reverse_end_of_chain = true;
      } else {
        m_current_segment_index = m_parent.prev(m_current_segment_index);
      }

      spdlog::trace("Iterator moved to segment {}", m_current_segment_index);

      return *this;
    }

    SegmentChainIterator operator--(int)
    {
      SegmentChainIterator temp = *this;
      --*this;
      return temp;
    }

    bool at_end_of_chain() { return m_is_end_of_chain; }

    bool at_reverse_end_of_chain() { return m_is_reverse_end_of_chain; }

    SegmentIndex operator*() { return m_current_segment_index; }

  private:
    const ProjectedCurveNetwork& m_parent;
    SegmentIndex m_current_segment_index;
    bool m_is_end_of_chain;
    bool m_is_reverse_end_of_chain;
  };

  /// Build a segment chain iterator from some starting segment index.
  ///
  /// @param segment_index: index to start the chain iterator at
  /// @return chain iterator for the given starting segment
  SegmentChainIterator get_segment_chain_iterator(
    SegmentIndex segment_index) const;

  /// Output visible planar curves as an SVG file, with visibility determined by
  /// having a non-positive quantitative invisibility.
  ///
  /// @param[in] output_path: filepath for the output SVG file
  /// @param[in] color_mode: choice of segment coloring
  /// @param[in] show_nodes: show nodes iff true
  void write(const std::string& output_path,
             SVGOutputMode color_mode = SVGOutputMode::uniform_visible_chains,
             bool show_nodes = false) const;

  /// Write closed curves to file as polygons with cusps indicated in a parallel
  /// list as 1 for cusp and 0 for no cusp.
  ///
  /// @param[in] filename: file to write the serialized contours to
  void serialize_closed_curves(const std::string& filename) const;

  /// Clear all topology and geometry data
  void clear();

  /// Add spatial curve network segments and annotated nodes to polyscope
  /// viewer.
  void add_spatial_network_to_viewer() const;

  /// View spatial curve network segments and annotated nodes.
  void spatial_network_viewer() const;

  /// Get start and end nodes for the closed curves of the surface. For open
  /// curves, this is the start node of the curve.
  ///
  /// @param[out] closed_curve_start_nodes: closed curve start node indices
  /// @param[out] closed_curve_end_nodes: closed curve end node indices
  void get_closed_curves(std::vector<NodeIndex>& closed_curve_start_nodes,
                         std::vector<NodeIndex>& closed_curve_end_nodes) const;

  /// Get start and end nodes for the visible curves of the surface. For open
  /// curves, this is the start node of the curve.
  ///
  /// @param[out] visible_curve_start_nodes: visible curve start node indices
  /// @param[out] visible_curve_end_nodes: visible curve end node indices
  void get_visible_curves(
    std::vector<NodeIndex>& visible_curve_start_nodes,
    std::vector<NodeIndex>& visible_curve_end_nodes) const;

protected:
  void init_projected_curve_network(
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
    int num_intersections);

  bool is_valid_projected_curve_network() const;

private:
  void init_chain_start_nodes();

  void discretize_segment_chain(SegmentChainIterator& iter,
                                std::vector<PlanarPoint>& points,
                                std::vector<int>& polyline) const;
  void discretize_curve(NodeIndex start_node_index,
                        NodeIndex end_node_index,
                        const CurveDiscretizationParameters& curve_disc_params,
                        std::vector<PlanarPoint>& points,
                        std::vector<int>& polyline,
                        std::vector<bool>& is_cusp) const;
  void simplify_curves(
    const std::vector<NodeIndex>& visible_curve_start_nodes,
    const std::vector<NodeIndex>& visible_curve_end_nodes,
    const std::vector<std::vector<PlanarPoint>>& all_points,
    std::vector<std::vector<PlanarPoint>>& simplified_points,
    std::vector<std::vector<int>>& simplified_polylines) const;
  void write_uniform_segments(svg::SVG& svgWriter) const;
  void write_uniform_visible_segments(svg::SVG& svgWriter) const;
  void write_contrast_invisible_segments(svg::SVG& svgWriter) const;
  void write_random_chains(svg::SVG& svgWriter) const;
  void write_uniform_chains(svg::SVG& svgWriter) const;
  void write_uniform_visible_chains(svg::SVG& svgWriter) const;
  void write_uniform_visible_curves(svg::SVG& svgWriter) const;
  void write_uniform_closed_curves(svg::SVG& svgWriter) const;
  void write_uniform_simplified_visible_curves(svg::SVG& svgWriter) const;

  void add_planar_segments_to_viewer() const;

  void add_spatial_segments_to_viewer() const;

  void add_spatial_nodes_to_viewer() const;

  void clear_geometry();

  // Geometry
  std::vector<SegmentGeometry> m_segments;
  std::vector<NodeGeometry> m_nodes;
  std::vector<NodeIndex> m_chain_start_nodes;
};
