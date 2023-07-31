// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_contours.h"

#include "compute_boundaries.h"
#include "intersect_conic.h"
#include "parametrize_conic.h"
#include "validity.h"

// Return false if the contour domain curve segments are not contained in the
// patch domain based on a simple sampling heuristic
bool
are_contained_in_patch_heuristic(
  const std::vector<Conic>& contour_domain_curve_segments,
  const QuadraticSplineSurfacePatch& spline_surface_patch)
{
  for (size_t i = 0; i < contour_domain_curve_segments.size(); ++i) {
    // Get curve sample points at the midpoint and near endpoints
    PlanarPoint start_point, mid_point, end_point;
    contour_domain_curve_segments[i].evaluate_normalized_coordinate(
      0.01, start_point);
    contour_domain_curve_segments[i].evaluate_normalized_coordinate(0.5,
                                                                    mid_point);
    contour_domain_curve_segments[i].evaluate_normalized_coordinate(0.99,
                                                                    end_point);

    // Check domain containment
    const ConvexPolygon& domain = spline_surface_patch.get_domain();
    if (!domain.contains(start_point))
      return false;
    if (!domain.contains(mid_point))
      return false;
    if (!domain.contains(end_point))
      return false;
  }

  return true;
}

void
compute_contour_equation(const Matrix6x3r& normal_mapping_coeffs,
                         const Matrix3x3r& frame,
                         Eigen::Matrix<double, 6, 1>& contour_equation_coeffs)
{
  assert(is_valid_spatial_mapping(normal_mapping_coeffs));
  assert(is_valid_frame(frame));

  Eigen::Matrix<double, 3, 1> tau = frame.col(2);
  contour_equation_coeffs = normal_mapping_coeffs * tau;

  // Normalize equation coefficients for stability
  contour_equation_coeffs /= contour_equation_coeffs.norm();
}

void
compute_quadratic_surface_contours(
  const Matrix6x3r& surface_mapping_coeffs,
  const Matrix6x3r& normal_mapping_coeffs,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments)
{
  assert(is_valid_spatial_mapping(surface_mapping_coeffs));
  assert(is_valid_spatial_mapping(normal_mapping_coeffs));
  assert(is_valid_frame(frame));

  SPDLOG_TRACE("Computing contours for quadratic surface {} with normals {}",
               formatted_bivariate_quadratic_mapping(surface_mapping_coeffs),
               formatted_bivariate_quadratic_mapping(normal_mapping_coeffs));

  // Compute the coefficients of the implicit contour function
  Eigen::Matrix<double, 6, 1> contour_equation_coeffs;
  compute_contour_equation(
    normal_mapping_coeffs, frame, contour_equation_coeffs);
  SPDLOG_TRACE("Contour function: {}",
               formatted_bivariate_quadratic_mapping(contour_equation_coeffs));

  // Parametrize the contour in the parametric domain
  parametrize_conic(contour_equation_coeffs, contour_domain_curve_segments);
  SPDLOG_TRACE("Contour domain curve segments:\n{}",
               formatted_vector(contour_domain_curve_segments));

  // Lift the contour to the surface
  contour_segments.resize(contour_domain_curve_segments.size());
  for (size_t k = 0; k < contour_domain_curve_segments.size(); ++k) {
    contour_domain_curve_segments[k].pullback_quadratic_function<3>(
      surface_mapping_coeffs, contour_segments[k]);
  }
  SPDLOG_TRACE("Contour segments:\n{}", formatted_vector(contour_segments));
}

// Get all contour segments for a spline surface patch
void
compute_spline_surface_patch_contours(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<std::pair<int, int>>& line_intersection_indices)
{
  spdlog::trace("Computing contours for spline surface patch");
  contour_domain_curve_segments.clear();
  contour_segments.clear();
  line_intersection_indices.clear();

  // Get surface mapping
  Matrix6x3r surface_mapping_coeffs =
    spline_surface_patch.get_surface_mapping();
  SPDLOG_TRACE("Patch surface mapping coefficients: {}",
               formatted_bivariate_quadratic_mapping(surface_mapping_coeffs));

  // Get surface normal mapping
  Matrix6x3r normal_mapping_coeffs = spline_surface_patch.get_normal_mapping();
  SPDLOG_TRACE("Patch normal mapping coefficients: {}",
               formatted_bivariate_quadratic_mapping(normal_mapping_coeffs));

  // Get implicit contour equation
  Eigen::Matrix<double, 6, 1> contour_equation_coeffs;
  compute_contour_equation(
    normal_mapping_coeffs, frame, contour_equation_coeffs);
  SPDLOG_TRACE("Patch contour equation coefficients: {}",
               formatted_bivariate_quadratic_mapping(contour_equation_coeffs));

  // Get full quadratic contours
  SPDLOG_TRACE("Parametrizing patch contour domain curves");
  std::vector<Conic> contour_domain_curves;
  parametrize_conic(contour_equation_coeffs, contour_domain_curves);
  SPDLOG_TRACE("Domain curves: {}", formatted_vector(contour_domain_curves));

  // Intersect contour domain curves with patch boundaries
  SPDLOG_TRACE("Intersecting domain curves with patch domain boundary");
  const ConvexPolygon& domain = spline_surface_patch.get_domain();
  for (size_t i = 0; i < contour_domain_curves.size(); ++i) {
    Conic current_contour_domain_curve = contour_domain_curves[i];
    std::vector<Conic> current_contour_domain_curve_segments;
    std::vector<std::pair<int, int>> current_contour_line_intersection_indices;
    intersect_conic_with_convex_polygon(
      current_contour_domain_curve,
      domain,
      current_contour_domain_curve_segments,
      current_contour_line_intersection_indices);
    append(contour_domain_curve_segments,
           current_contour_domain_curve_segments);
    append(line_intersection_indices,
           current_contour_line_intersection_indices);
    spdlog::trace("Contour domain curve split into {} segments",
                  current_contour_domain_curve_segments.size());
  }

  // Lift contour domain curves to the surface
  spdlog::trace("Lifting domain curves to the surface");
  for (size_t i = 0; i < contour_domain_curve_segments.size(); ++i) {
    Conic contour_domain_curve_segment = contour_domain_curve_segments[i];
    RationalFunction<4, 3> contour_segment;
    contour_domain_curve_segment.pullback_quadratic_function<3>(
      surface_mapping_coeffs, contour_segment);
    contour_segments.push_back(contour_segment);
    spdlog::trace("Domain curve {} lifted to {}",
                  contour_domain_curve_segment,
                  contour_segment);
  }

  assert(are_contained_in_patch_heuristic(contour_domain_curve_segments,
                                          spline_surface_patch));
}

// Get all contour segments for a spline surface patch
void
compute_spline_surface_cone_patch_contours(
  const QuadraticSplineSurfacePatch& spline_surface_patch,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<std::pair<int, int>>& line_intersection_indices)
{
  spdlog::trace("Computing contours for spline surface patch");
  contour_domain_curve_segments.clear();
  contour_segments.clear();
  line_intersection_indices.clear();

  // Check for cone
  if (!spline_surface_patch.has_cone()) {
    spdlog::error(
      "Tried to compute cone patch contours in a patch without a cone");
    return;
  }
  int cone_corner_index = spline_surface_patch.get_cone();

  // Get surface mapping
  Matrix6x3r surface_mapping_coeffs =
    spline_surface_patch.get_surface_mapping();
  SPDLOG_TRACE("Patch surface mapping coefficients: {}",
               formatted_bivariate_quadratic_mapping(surface_mapping_coeffs));

  // Get surface normal mapping
  Matrix6x3r normal_mapping_coeffs = spline_surface_patch.get_normal_mapping();
  SPDLOG_TRACE("Patch normal mapping coefficients: {}",
               formatted_bivariate_quadratic_mapping(normal_mapping_coeffs));

  // Get implicit contour equation
  Vector6r contour_equation_coeffs;
  compute_contour_equation(
    normal_mapping_coeffs, frame, contour_equation_coeffs);
  SPDLOG_TRACE(
    "Patch contour equation coefficients: {}",
    formatted_bivariate_quadratic_mapping(contour_equation_coeffs, 17));

  // Get full quadratic contours
  SPDLOG_TRACE("Parametrizing cone patch contour domain curves");
  std::vector<Conic> contour_domain_curves;
  parametrize_cone_patch_conic(contour_equation_coeffs, contour_domain_curves);
  SPDLOG_TRACE("Domain curves: {}", formatted_vector(contour_domain_curves));

  // Intersect contour domain curves with patch boundaries
  const ConvexPolygon& domain = spline_surface_patch.get_domain();
  spdlog::trace("Intersecting domain curves with patch domain {}",
                domain.get_vertices());
  for (size_t i = 0; i < contour_domain_curves.size(); ++i) {
    Conic current_contour_domain_curve = contour_domain_curves[i];

    // Check if the contour is stable
    bool conic_intersects = check_if_conic_intersects_cone_patch_domain(
      current_contour_domain_curve, domain, cone_corner_index);
    if (!conic_intersects) {
      spdlog::trace("Line through conic does not intersect the domain");
      continue;
    }

    // Intersect the contour domain with the domain robustly
    Conic current_contour_domain_curve_segment;
    std::pair<int, int> current_contour_line_intersection_indices;
    conic_intersects =
      intersect_conic_in_cone_patch(current_contour_domain_curve,
                                    domain,
                                    cone_corner_index,
                                    current_contour_domain_curve_segment,
                                    current_contour_line_intersection_indices);
    if (!conic_intersects) {
      spdlog::trace("Conic ray does not intersect the domain");
      continue;
    }

    // Add the contour if it exists
    contour_domain_curve_segments.push_back(
      current_contour_domain_curve_segment);
    line_intersection_indices.push_back(
      current_contour_line_intersection_indices);
  }

  // Lift contour domain curves to the surface
  spdlog::trace("Lifting cone patch domain curves to the surface");
  for (size_t i = 0; i < contour_domain_curve_segments.size(); ++i) {
    const Conic& contour_domain_curve_segment =
      contour_domain_curve_segments[i];
    RationalFunction<4, 3> contour_segment;
    contour_domain_curve_segment.pullback_quadratic_function<3>(
      surface_mapping_coeffs, contour_segment);
    contour_segments.push_back(contour_segment);
    spdlog::debug("Cone patch domain curve {} lifted to {}",
                  contour_domain_curve_segment,
                  contour_segment);
  }

  assert(are_contained_in_patch_heuristic(contour_domain_curve_segments,
                                          spline_surface_patch));
}

void
compute_spline_surface_contours(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  std::vector<std::pair<int, int>>& line_intersection_indices)
{
  assert(is_valid_frame(frame));

  contour_domain_curve_segments.clear();
  contour_segments.clear();
  contour_patch_indices.clear();
  line_intersection_indices.clear();

  // Compute contours for each quadratic patch
  for (QuadraticSplineSurface::PatchIndex patch_index = 0;
       patch_index < spline_surface.num_patches();
       ++patch_index) {
    QuadraticSplineSurfacePatch const& spline_surface_patch =
      spline_surface.get_patch(patch_index);

    std::vector<Conic> patch_contour_domain_curve_segments;
    std::vector<RationalFunction<4, 3>> patch_contour_segments;
    std::vector<std::pair<int, int>> patch_line_intersection_indices;
    if (spline_surface_patch.has_cone()) {
      spdlog::trace("Parametrizing cone patch {}", patch_index);
      compute_spline_surface_cone_patch_contours(
        spline_surface_patch,
        frame,
        patch_contour_domain_curve_segments,
        patch_contour_segments,
        patch_line_intersection_indices);
    } else {
      compute_spline_surface_patch_contours(spline_surface_patch,
                                            frame,
                                            patch_contour_domain_curve_segments,
                                            patch_contour_segments,
                                            patch_line_intersection_indices);
    }

    // Append patch indices
    size_t num_patch_contours = patch_contour_segments.size();
    std::vector<QuadraticSplineSurface::PatchIndex> patch_index_list(
      num_patch_contours, patch_index);
    append(contour_domain_curve_segments, patch_contour_domain_curve_segments);
    append(contour_segments, patch_contour_segments);
    append(contour_patch_indices, patch_index_list);
    append(line_intersection_indices, patch_line_intersection_indices);
    assert(contour_domain_curve_segments.size() == contour_segments.size());
    assert(contour_patch_indices.size() == contour_segments.size());
  }

  spdlog::info("Found {} contour segments", contour_segments.size());
}

void
compute_spline_surface_boundaries(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<std::pair<int, int>>& patch_boundary_edges,
  std::vector<Conic>& boundary_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& boundary_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& boundary_patch_indices)
{
  size_t num_boundary_edges = patch_boundary_edges.size();
  boundary_domain_curve_segments.reserve(num_boundary_edges);
  boundary_segments.reserve(num_boundary_edges);
  boundary_patch_indices.reserve(num_boundary_edges);

  // Get patch boundary
  for (size_t i = 0; i < num_boundary_edges; ++i) {
    // Get spline surface patch
    int patch_index = patch_boundary_edges[i].first;
    int patch_edge_index = patch_boundary_edges[i].second;
    QuadraticSplineSurfacePatch const& spline_surface_patch =
      spline_surface.get_patch(patch_index);

    // Get patch domain boundaries
    std::array<LineSegment, 3> patch_domain_boundaries;
    spline_surface_patch.get_domain().parametrize_patch_boundaries(
      patch_domain_boundaries);
    boundary_domain_curve_segments.push_back(
      patch_domain_boundaries[patch_edge_index]);

    // Get patch boundaries
    std::array<RationalFunction<4, 3>, 3> patch_boundaries;
    spline_surface_patch.get_patch_boundaries(patch_boundaries);
    boundary_segments.push_back(patch_boundaries[patch_edge_index]);

    // Record patch index
    boundary_patch_indices.push_back(patch_index);
  }
}

double
compute_boundary_intersection_parameter(
  const Conic& boundary_domain_curve_segment,
  const PlanarPoint& intersection_point)
{
  // The boundary domain curve is always a line
  Matrix3x2r P_coeffs = boundary_domain_curve_segment.get_numerators();
  PlanarPoint x0 = P_coeffs.row(0) - intersection_point;
  PlanarPoint d = P_coeffs.row(1);
  return -x0.dot(d) / d.dot(d);
}

void
compute_spline_surface_boundary_intersections(
  const QuadraticSplineSurface& spline_surface,
  const std::vector<Conic>& contour_domain_curve_segments,
  const std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  const std::vector<std::pair<int, int>>& line_intersection_indices,
  const std::vector<std::pair<int, int>>& patch_boundary_edges,
  const std::vector<Conic>& boundary_domain_curve_segments,
  std::vector<std::vector<IntersectionData>>& contour_intersections,
  int& num_intersections)
{
  size_t num_patches = spline_surface.num_patches();
  size_t num_interior_contours = contour_domain_curve_segments.size();

  // Build map from spine surface edges to boundary edge indices (or -1)
  std::vector<std::array<int, 3>> patch_boundary_contour_map(num_patches,
                                                             { -1, -1, -1 });
  std::vector<bool> is_boundary_patch(num_patches, false);
  spdlog::debug("Intersecting with {} patch boundaries",
                patch_boundary_edges.size());
  for (size_t i = 0; i < patch_boundary_edges.size(); ++i) {
    int patch_index = patch_boundary_edges[i].first;
    int patch_edge_index = patch_boundary_edges[i].second;
    spdlog::debug("Boundary {}: {}, {}", i, patch_index, patch_edge_index);
    patch_boundary_contour_map[patch_index][patch_edge_index] = i;
    is_boundary_patch[patch_index] = true;
  }

  // Compute intersections for contours with intersections on the boundary
  for (size_t i = 0; i < contour_domain_curve_segments.size(); ++i) {
    size_t patch_index = contour_patch_indices[i];
    // FIXME This +1 is from an indexing inconsistency somewhere that should be
    // fixed
    int start_line_intersection_index =
      (line_intersection_indices[i].first + 1) % 3;
    int end_line_intersection_index =
      (line_intersection_indices[i].second + 1) % 3;
    if (is_boundary_patch[patch_index]) {
      spdlog::debug("Contour {} in boundary patch {}", i, patch_index);
      spdlog::debug("Contour {} has line intersection indices {} and {}",
                    i,
                    start_line_intersection_index,
                    end_line_intersection_index);
      if ((start_line_intersection_index >= 0) &&
          (end_line_intersection_index >= 0)) {
        spdlog::debug(
          "Contour {} has boundary contours {} and {}",
          i,
          patch_boundary_contour_map[patch_index]
                                    [start_line_intersection_index],
          patch_boundary_contour_map[patch_index][end_line_intersection_index]);
      }
    }

    // Handle start point
    if (start_line_intersection_index >= 0) {
      int start_boundary_contour =
        patch_boundary_contour_map[patch_index][start_line_intersection_index];
      if (start_boundary_contour >= 0) {
        spdlog::debug("Start point intersection for interior contour {} in "
                      "patch {} with boundary contour {} on patch line {}",
                      i,
                      patch_index,
                      start_boundary_contour,
                      start_line_intersection_index);
        spdlog::debug("Interior domain contour: {}",
                      contour_domain_curve_segments[i]);
        spdlog::debug("Boundary domain contour: {}",
                      boundary_domain_curve_segments[start_boundary_contour]);
        // Build interior intersection data
        IntersectionData interior_intersection_data;
        interior_intersection_data.knot =
          contour_domain_curve_segments[i].domain().get_lower_bound();
        // interior_intersection_data.knot =
        // contour_domain_curve_segments[i].domain().get_upper_bound();
        interior_intersection_data.intersection_index =
          num_interior_contours + start_boundary_contour;
        interior_intersection_data.intersection_knot =
          compute_boundary_intersection_parameter(
            boundary_domain_curve_segments[start_boundary_contour],
            contour_domain_curve_segments[i](interior_intersection_data.knot));
        interior_intersection_data.is_base = true;
        interior_intersection_data.id = num_intersections;
        contour_intersections[i].push_back(interior_intersection_data);

        // Build complementary boundary intersection data
        IntersectionData boundary_intersection_data;
        boundary_intersection_data.knot =
          interior_intersection_data.intersection_knot;
        boundary_intersection_data.intersection_index = i;
        boundary_intersection_data.intersection_knot =
          interior_intersection_data.knot;
        boundary_intersection_data.id = num_intersections;
        contour_intersections[num_interior_contours + start_boundary_contour]
          .push_back(boundary_intersection_data);
        num_intersections++;
      }
    } else {
      spdlog::debug("Contour {} start in patch {} does not intersect an edge",
                    i,
                    patch_index);
    }

    // Handle endpoint
    if (end_line_intersection_index >= 0) {
      int end_boundary_contour =
        patch_boundary_contour_map[patch_index][end_line_intersection_index];
      if (end_boundary_contour >= 0) {
        spdlog::debug("End point intersection for interior contour {} in patch "
                      "{} with boundary contour {} on patch line {}",
                      i,
                      patch_index,
                      end_boundary_contour,
                      end_line_intersection_index);
        spdlog::debug("Interior domain contour: {}",
                      contour_domain_curve_segments[i]);
        spdlog::debug("Boundary domain contour: {}",
                      boundary_domain_curve_segments[end_boundary_contour]);
        // Build interior intersection data
        IntersectionData interior_intersection_data;
        interior_intersection_data.knot =
          contour_domain_curve_segments[i].domain().get_upper_bound();
        // interior_intersection_data.knot =
        // contour_domain_curve_segments[i].domain().get_lower_bound();
        interior_intersection_data.intersection_index =
          num_interior_contours + end_boundary_contour;
        interior_intersection_data.intersection_knot =
          compute_boundary_intersection_parameter(
            boundary_domain_curve_segments[end_boundary_contour],
            contour_domain_curve_segments[i](interior_intersection_data.knot));
        interior_intersection_data.is_tip = true;
        interior_intersection_data.id = num_intersections;
        contour_intersections[i].push_back(interior_intersection_data);

        // Build complementary boundary intersection data
        IntersectionData boundary_intersection_data;
        boundary_intersection_data.knot =
          interior_intersection_data.intersection_knot;
        boundary_intersection_data.intersection_index = i;
        boundary_intersection_data.intersection_knot =
          interior_intersection_data.knot;
        boundary_intersection_data.id = num_intersections;
        contour_intersections[num_interior_contours + end_boundary_contour]
          .push_back(boundary_intersection_data);
        num_intersections++;
      }
    } else {
      spdlog::debug("Contour {} end in patch {} does not intersect an edge",
                    i,
                    patch_index);
    }
  }
  spdlog::info("{} interior and boundary intersections found",
               num_intersections);
}

void
compute_spline_surface_contours_and_boundaries(
  const QuadraticSplineSurface& spline_surface,
  const Matrix3x3r& frame,
  const std::vector<std::pair<int, int>>& patch_boundary_edges,
  std::vector<Conic>& contour_domain_curve_segments,
  std::vector<RationalFunction<4, 3>>& contour_segments,
  std::vector<QuadraticSplineSurface::PatchIndex>& contour_patch_indices,
  std::vector<bool>& contour_is_boundary,
  std::vector<std::vector<IntersectionData>>& contour_intersections,
  int& num_intersections)
{
  num_intersections = 0;

  // Compute contours
  std::vector<std::pair<int, int>> line_intersection_indices;
  compute_spline_surface_contours(spline_surface,
                                  frame,
                                  contour_domain_curve_segments,
                                  contour_segments,
                                  contour_patch_indices,
                                  line_intersection_indices);
  contour_is_boundary = std::vector<bool>(contour_segments.size(), false);

  // Compute boundaries
  std::vector<Conic> boundary_domain_curve_segments;
  std::vector<RationalFunction<4, 3>> boundary_segments;
  std::vector<QuadraticSplineSurface::PatchIndex> boundary_patch_indices;
  compute_spline_surface_boundaries(spline_surface,
                                    patch_boundary_edges,
                                    boundary_domain_curve_segments,
                                    boundary_segments,
                                    boundary_patch_indices);
  std::vector<bool> boundary_is_boundary(boundary_segments.size(), true);

  // Compute intersections of the contours and the boundaries
  contour_intersections.resize(contour_domain_curve_segments.size() +
                               boundary_domain_curve_segments.size());
  for (size_t i = 0; i < contour_intersections.size(); ++i) {
    contour_intersections[i].clear();
  }
  compute_spline_surface_boundary_intersections(spline_surface,
                                                contour_domain_curve_segments,
                                                contour_patch_indices,
                                                line_intersection_indices,
                                                patch_boundary_edges,
                                                boundary_domain_curve_segments,
                                                contour_intersections,
                                                num_intersections);

  append(contour_domain_curve_segments, boundary_domain_curve_segments);
  append(contour_segments, boundary_segments);
  append(contour_patch_indices, boundary_patch_indices);
  append(contour_is_boundary, boundary_is_boundary);
}

void
pad_contours(std::vector<Conic>& contour_domain_curve_segments,
             std::vector<RationalFunction<4, 3>>& contour_segments,
             std::vector<RationalFunction<4, 2>>& planar_contour_segments,
             double pad_amount)
{
  if (pad_amount <= 0.0)
    return;

  for (size_t i = 0; i < contour_domain_curve_segments.size(); ++i) {
    contour_domain_curve_segments[i].domain().pad_lower_bound(pad_amount);
    contour_domain_curve_segments[i].domain().pad_upper_bound(pad_amount);
    contour_segments[i].domain().pad_lower_bound(pad_amount);
    contour_segments[i].domain().pad_upper_bound(pad_amount);
    planar_contour_segments[i].domain().pad_lower_bound(pad_amount);
    planar_contour_segments[i].domain().pad_upper_bound(pad_amount);
  }
}