// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "convex_polygon.h"
#include "evaluate_surface_normal.h"
#include "polynomial_function.h"
#include "polyscope/point_cloud.h"
#include <igl/is_vertex_manifold.h>


/// \file quadratic_spline_surface.h
///
/// Representation for quadratic surface patches with convex domains.

/// A quadratic surface patch with convex polygonal domain.
///
/// Supports:
/// - evaluation and sampling of points and normals on the surface
/// - triangulation
/// - conversion to Bezier form
/// - bounding box computation
/// - boundary curve parameterization
/// - cone point annotation
/// - (de)serialization
class QuadraticSplineSurfacePatch
{
public:
  /// @brief Default constructor
  QuadraticSplineSurfacePatch();

  /// Construct the patch from a domain polygon surface coefficients in order:
  ///   a_0 a_u a_v a_uv a_uu a_vv
  ///
  /// @param[in] surface_mapping_coeffs: surface mapping coefficients
  /// @param[in] domain: uv domain boundaries
  QuadraticSplineSurfacePatch(const Matrix6x3r& surface_mapping_coeffs,
                              const ConvexPolygon& domain);

  /// Get the dimension of the surface patch ambient space
  ///
  /// @return dimension of the ambient space
  int dimension() const;

  /// Mark one of the vertices as a cone
  ///
  /// @param[in] cone_index: index of the cone in the triangle
  void mark_cone(int cone_index);

  /// Determine if the patch has a cone
  ///
  /// @return true iff the patch has a cone
  bool has_cone() const;

  /// Get the cone index, or -1 if none exists
  ///
  /// @return true iff the patch has a cone
  int get_cone() const;

  /// Get the surface mapping coefficients
  ///
  /// @return reference to the surface mapping
  Matrix6x3r const& get_surface_mapping() const;

  /// Get the surface normal mapping coefficients
  ///
  /// @return reference to the surface normal mapping
  Matrix6x3r const& get_normal_mapping() const;

  /// Get the surface mapping coefficients with normalized domain.
  ///
  /// Compute them if they haven't been computed yet
  ///
  /// @return reference to the normalized surface mapping
  Matrix6x3r const& get_normalized_surface_mapping() const;

  /// Get the surface mapping coefficients with normalized domain.
  ///
  /// Compute them if they haven't been computed yet
  ///
  /// @return reference to the bezier points
  Matrix6x3r const& get_bezier_points() const;

  /// @brief Compute the bounding box for the surface patch
  ///
  /// @param[out] min_point: minimum coordinates bounding box point
  /// @param[out] max_point: maximum coordinates bounding box point
  void get_bounding_box(SpatialVector& min_point,
                        SpatialVector& max_point) const
  {
    min_point = m_min_point;
    max_point = m_max_point;
  }

  /// @brief Compute the minimum point of the bounding box for the surface patch
  ///
  /// @return min_point: minimum coordinates bounding box point
  SpatialVector const& get_bounding_box_min_point() const
  {
    return m_min_point;
  }
  double get_bbox_x_min() const { return m_min_point[0]; }
  double get_bbox_y_min() const { return m_min_point[1]; }

  /// @brief Compute the maximum point of the bounding box for the surface patch
  ///
  /// @return max_point: maximum coordinates bounding box point
  SpatialVector const& get_bounding_box_max_point() const
  {
    return m_max_point;
  }
  double get_bbox_x_max() const { return m_max_point[0]; }
  double get_bbox_y_max() const { return m_max_point[1]; }

  /// Get the convex domain of the patch.
  ///
  /// @return reference to the convex domain
  ConvexPolygon const& get_domain() const;

  /// Get the patch boundaries as spatial curves.
  ///
  /// @param[out] patch_boundaries: patch boundary spatial curves
  void get_patch_boundaries(
    std::array<RationalFunction<4, 3>, 3>& patch_boundaries) const;

  /// Construct a spline surface patch with the same image but where the domain
  /// is normalized to the triangle u + v <= 1 in the positive quadrant.
  ///
  /// @param[out] normalized_spline_surface_patch: normalized patch
  void normalize_patch_domain(
    QuadraticSplineSurfacePatch& normalized_spline_surface_patch) const;

  /// Given a normalized domain point in the triangle u + v <= 1, map it to the
  /// corresponding point in the patch domain
  ///
  /// @param[in] normalized_domain_point: normalized (barycentric) domain point
  /// @return corresponding point in the domain triangle
  PlanarPoint denormalize_domain_point(
    PlanarPoint& normalized_domain_point) const;

  /// Evaluate the surface at a given domain point.
  ///
  /// @param[in] domain_point: domain evaluation point
  /// @param[in] surface_point: image of the domain point on the surface
  void evaluate(const PlanarPoint& domain_point,
                SpatialVector& surface_point) const;

  /// Evaluate the surface normal at a given domain point.
  ///
  /// @param[in] domain_point: domain evaluation point
  /// @param[in] surface_point: surface normal at the image of the domain point
  void evaluate_normal(const PlanarPoint& domain_point,
                       SpatialVector& surface_normal) const;

  /// Sample points on the surface.
  ///
  /// @param[in] sampling_density: sampling density parameter
  /// @param[out] spline_surface_patch_points: sampled points on the surface
  void sample(size_t sampling_density,
              std::vector<SpatialVector>& spline_surface_patch_points) const;

  /// Triangulate the surface patch.
  ///
  /// @param[in] num_refinements: number of refinements of the domain to perform
  /// @param[out] V: triangulated patch vertex positions
  /// @param[out] F: triangulated patch faces
  /// @param[out] N: triangulated patch vertex normals
  void triangulate(size_t num_refinements,
                   Eigen::MatrixXd& V,
                   Eigen::MatrixXi& F,
                   Eigen::MatrixXd& N) const;

  /// Add triangulated patch to the polyscope viewer.
  ///
  /// @param[in] patch_name: name to assign the patch in the viewer
  void add_patch_to_viewer(std::string patch_name = "surface_patch") const;

  /// Overloaded insertion operator.
  friend std::ostream& operator<<(
    std::ostream& out,
    const QuadraticSplineSurfacePatch& spline_surface_patch);

  /// Write the patch information to the output stream in the format
  ///   c a_0 a_u a_v a_uv a_uu a_vv
  ///   p1 p1_u p1_v
  ///   p2 p2_u p2_v
  ///   p3 p3_u p3_v
  ///
  /// @param[out] out: stream to write serialization to
  void serialize(std::ostream& out) const;

  /// Write patch to file.
  ///
  /// @param[in] filename: file to write serialized patch to
  void write_patch(const std::string& filename) const;

private:
  std::string formatted_patch() const;

  // Core independent data
  Matrix6x3r m_surface_mapping_coeffs;
  ConvexPolygon m_domain;

  // Inferred dependent data
  Matrix6x3r m_normal_mapping_coeffs;
  Matrix6x3r m_normalized_surface_mapping_coeffs;
  Matrix6x3r m_bezier_points;
  SpatialVector m_min_point;
  SpatialVector m_max_point;

  // Additional cone marker to handle degenerate configurations
  int m_cone_index;
};
