// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "twelve_split_spline.h"
#include "compute_boundaries.h"

TwelveSplitSplineSurface::TwelveSplitSplineSurface(
  const Eigen::MatrixXd& V,
  const AffineManifold& affine_manifold,
  const OptimizationParameters& optimization_params,
  std::vector<std::vector<int>>& face_to_patch_indices,
  std::vector<int>& patch_to_face_indices,
  Eigen::SparseMatrix<double>& fit_matrix,
  Eigen::SparseMatrix<double>& energy_hessian,
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>&
    energy_hessian_inverse)
  : m_affine_manifold(affine_manifold)
{
  // Generate normals
  MatrixXr N;
  generate_face_normals(V, affine_manifold, N);

  // Generate fit matrix by setting the parametrized quadratic surface mapping
  // factor to zero
  double fit_energy;
  VectorXr fit_derivatives;
  OptimizationParameters optimization_params_fit = optimization_params;
  optimization_params_fit.parametrized_quadratic_surface_mapping_factor = 0.0;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> fit_matrix_inverse;
  build_twelve_split_spline_energy_system(V,
                                          N,
                                          affine_manifold,
                                          optimization_params_fit,
                                          fit_energy,
                                          fit_derivatives,
                                          fit_matrix,
                                          fit_matrix_inverse);

  // Build full energy hessian system
  double energy;
  VectorXr derivatives;
  build_twelve_split_spline_energy_system(V,
                                          N,
                                          affine_manifold,
                                          optimization_params,
                                          energy,
                                          derivatives,
                                          energy_hessian,
                                          energy_hessian_inverse);

  // Build optimized corner and midpoint data
  generate_optimized_twelve_split_position_data(V,
                                                affine_manifold,
                                                fit_matrix,
                                                energy_hessian_inverse,
                                                m_corner_data,
                                                m_midpoint_data);

  // Get cone corners
  std::vector<std::array<bool, 3>> is_cone_corner;
  affine_manifold.compute_cones_corners(is_cone_corner);

  // Initialize position data and patches
  init_twelve_split_patches(m_corner_data,
                            m_midpoint_data,
                            is_cone_corner,
                            face_to_patch_indices,
                            patch_to_face_indices);
}

TwelveSplitSplineSurface::TwelveSplitSplineSurface(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  const std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data,
  std::vector<std::vector<int>>& face_to_patch_indices,
  std::vector<int>& patch_to_face_indices)
  : m_corner_data(corner_data)
  , m_midpoint_data(midpoint_data)
{
  int num_faces = corner_data.size();
  spdlog::info("Building surface directly from position data");

  // Assume no cones
  std::vector<std::array<bool, 3>> is_cone_corner(num_faces,
                                                  { false, false, false });

  // Initialize position data and patches
  init_twelve_split_patches(corner_data,
                            midpoint_data,
                            is_cone_corner,
                            face_to_patch_indices,
                            patch_to_face_indices);
}

void
TwelveSplitSplineSurface::update_positions(
  const Eigen::MatrixXd& V,
  const Eigen::SparseMatrix<double>& fit_matrix,
  const Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>&
    energy_hessian_inverse)
{
  AffineManifold const& affine_manifold = get_affine_manifold();

  // Generate normals
  MatrixXr N;
  generate_face_normals(V, affine_manifold, N);

  // Build optimized corner and midpoint data
  generate_optimized_twelve_split_position_data(V,
                                                affine_manifold,
                                                fit_matrix,
                                                energy_hessian_inverse,
                                                m_corner_data,
                                                m_midpoint_data);

  // Get cone corners
  std::vector<std::array<bool, 3>> is_cone_corner;
  affine_manifold.compute_cones_corners(is_cone_corner);

  // Initialize position data and patches
  std::vector<std::vector<int>> face_to_patch_indices;
  std::vector<int> patch_to_face_indices;
  init_twelve_split_patches(m_corner_data,
                            m_midpoint_data,
                            is_cone_corner,
                            face_to_patch_indices,
                            patch_to_face_indices);
}

void
TwelveSplitSplineSurface::add_position_data_to_viewer() const
{
  // Add corner position data if it exists
  if (!m_corner_data.empty()) {
    MatrixXr position_matrix;
    MatrixXr first_derivative_matrix;
    MatrixXr second_derivative_matrix;
    generate_corner_data_matrices(m_corner_data,
                                  position_matrix,
                                  first_derivative_matrix,
                                  second_derivative_matrix);
    polyscope::registerPointCloud("corner data", position_matrix);
    polyscope::getPointCloud("corner data")
      ->addVectorQuantity("first derivatives", first_derivative_matrix);
    polyscope::getPointCloud("corner data")
      ->addVectorQuantity("second derivatives", second_derivative_matrix);
  }

  // Add midpoint position data if it (and the corner data) exists
  if ((!m_corner_data.empty()) && (!m_midpoint_data.empty())) {
    MatrixXr position_matrix;
    MatrixXr tangent_derivative_matrix;
    MatrixXr normal_derivative_matrix;
    generate_midpoint_data_matrices(m_corner_data,
                                    m_midpoint_data,
                                    position_matrix,
                                    tangent_derivative_matrix,
                                    normal_derivative_matrix);
    polyscope::registerPointCloud("midpoint data", position_matrix);
    polyscope::getPointCloud("midpoint data")
      ->addVectorQuantity("tangent derivatives", tangent_derivative_matrix);
    polyscope::getPointCloud("midpoint data")
      ->addVectorQuantity("normal derivatives", normal_derivative_matrix);
  }
}

void
TwelveSplitSplineSurface::view(Eigen::Matrix<double, 3, 1> color,
                             int num_subdivisions) const
{
  add_surface_to_viewer(color, num_subdivisions);
  add_position_data_to_viewer();
  polyscope::show();
}

void
TwelveSplitSplineSurface::clear()
{
  m_affine_manifold.clear();
  m_corner_data.clear();
  m_midpoint_data.clear();
  m_patches.clear();
}

// Initialize patches for a twelve split Powell-Sabin type surface
void
TwelveSplitSplineSurface::init_twelve_split_patches(
  const std::vector<std::array<TriangleCornerFunctionData, 3>>& corner_data,
  const std::vector<std::array<TriangleMidpointFunctionData, 3>>& midpoint_data,
  const std::vector<std::array<bool, 3>>& is_cone_corner,
  std::vector<std::vector<int>>& face_to_patch_indices,
  std::vector<int>& patch_to_face_indices)
{
  int num_faces = corner_data.size();

  // Get number of patches per face
  int patches_per_face = 12;
  int num_patches = patches_per_face * num_faces;

  // Get general patch domains to use for all faces
  std::array<std::array<Eigen::Matrix<double, 3, 1>, 3>, 12> patch_boundaries;
  generate_twelve_split_spline_patch_patch_boundaries(patch_boundaries);
  std::vector<ConvexPolygon> domains(patches_per_face);
  for (int i = 0; i < patches_per_face; ++i) {
    domains[i] = ConvexPolygon(patch_boundaries[i]);
  }

  // Generate map from patches to input mesh corners
  std::array<std::pair<int, int>, 12> patch_to_corner_map;
  generate_twelve_split_spline_patch_patch_to_corner_map(patch_to_corner_map);

  // Clear face to patch mappings
  face_to_patch_indices.resize(num_faces);
  for (int i = 0; i < num_faces; ++i) {
    face_to_patch_indices[i].clear();
  }
  patch_to_face_indices.clear();
  patch_to_face_indices.reserve(num_patches);

  // Iterate over face position data
  m_patches.clear();
  m_patches.reserve(num_patches);
  for (int face_index = 0; face_index < num_faces; ++face_index) {
    // Get surface mappings
    std::array<Eigen::Matrix<double, 6, 3>, 12> surface_mappings;
    generate_twelve_split_spline_patch_surface_mapping<double>(
      corner_data[face_index], midpoint_data[face_index], surface_mappings);

    // Add patches
    face_to_patch_indices[face_index].clear();
    for (int j = 0; j < patches_per_face; ++j) {
      // Add patch to surface
      m_patches.push_back(
        QuadraticSplineSurfacePatch(surface_mappings[j], domains[j]));

      // Mark cones
      int corner_index = patch_to_corner_map[j].first;
      if ((corner_index >= 0) && (is_cone_corner[face_index][corner_index])) {
        int patch_cone_index = patch_to_corner_map[j].second;
        m_patches.back().mark_cone(patch_cone_index);
      }

      // Update indices
      int new_patch_index = m_patches.size() - 1;
      patch_to_face_indices.push_back(face_index);
      face_to_patch_indices[face_index].push_back(new_patch_index);
      assert(patch_to_face_indices[face_to_patch_indices[face_index].back()] ==
             face_index);
      assert(face_to_patch_indices[patch_to_face_indices.back()].back() ==
             new_patch_index);
    }
  }

  // Initialize hash tables
  compute_patch_hash_tables();
}

void
TwelveSplitSplineSurface::generate_face_normals(
  const Eigen::MatrixXd& V,
  const AffineManifold& affine_manifold,
  Eigen::MatrixXd& N)
{
  Eigen::MatrixXi const& F = affine_manifold.get_faces();

  // Compute the cones of the affine manifold
  std::vector<AffineManifold::Index> cones;
  affine_manifold.compute_cones(cones);

  // Get vertex normals
  Eigen::MatrixXd N_vertices;
  igl::per_vertex_normals(V, F, N_vertices);

  // Set the face one ring normals of the cone vertices to the cone vertex
  // normal
  N.setZero(F.rows(), 3);
  for (size_t i = 0; i < cones.size(); ++i) {
    int ci = cones[i];
    VertexManifoldChart const& chart = affine_manifold.get_vertex_chart(ci);
    for (size_t j = 0; j < chart.face_one_ring.size(); ++j) {
      int fj = chart.face_one_ring[j];
      N.row(fj) = N_vertices.row(ci);
    }
  }
}

void
generate_twelve_split_spline_patch_patch_boundaries(
  std::array<std::array<Eigen::Matrix<double, 3, 1>, 3>, 12>& patch_boundaries)
{
  size_t num_patches = 12;
  size_t num_boundaries = 3;
  size_t num_coeffs = 3;

  // Get boundary coefficients
  double bound_coeffs[12][3][3];
  PS12tri_bounds_coeffs(bound_coeffs);

  // Reorganize boundary coefficients
  for (size_t i = 0; i < num_patches; ++i) {
    for (size_t j = 0; j < num_boundaries; ++j) {
      for (size_t k = 0; k < num_coeffs; ++k) {
        patch_boundaries[i][j][k] = bound_coeffs[i][j][k];
      }
    }
  }
}

void
generate_twelve_split_spline_patch_patch_to_corner_map(
  std::array<std::pair<int, int>, 12>& patch_to_corner_map)
{
  // First six patches are interior
  for (size_t i = 0; i < 6; ++i) {
    patch_to_corner_map[i] = std::make_pair(-1, -1);
  }

  // Hand code the rest
  patch_to_corner_map[6] = std::make_pair(0, 1);
  patch_to_corner_map[7] = std::make_pair(1, 1);
  patch_to_corner_map[8] = std::make_pair(1, 1);
  patch_to_corner_map[9] = std::make_pair(2, 1);
  patch_to_corner_map[10] = std::make_pair(2, 1);
  patch_to_corner_map[11] = std::make_pair(0, 1);
}

void
generate_twelve_split_domain_areas(const PlanarPoint& v0,
                                   const PlanarPoint& v1,
                                   const PlanarPoint& v2,
                                   std::array<double, 12>& patch_areas)
{
  PlanarPoint m = (v0 + v1 + v2) / 3.0;
  PlanarPoint e01 = (v0 + v1) / 2.0;
  PlanarPoint e12 = (v1 + v2) / 2.0;
  PlanarPoint e20 = (v2 + v0) / 2.0;
  PlanarPoint f0 = (e01 + e20) / 2.0;
  PlanarPoint f1 = (e12 + e01) / 2.0;
  PlanarPoint f2 = (e20 + e12) / 2.0;

  patch_areas[0] = area_from_positions(e20, f0, m);
  patch_areas[1] = area_from_positions(f0, e01, m);
  patch_areas[2] = area_from_positions(e01, f1, m);
  patch_areas[3] = area_from_positions(f1, e12, m);
  patch_areas[4] = area_from_positions(e12, f2, m);
  patch_areas[5] = area_from_positions(f2, e20, m);
  patch_areas[6] = area_from_positions(v1, e12, f1);
  patch_areas[7] = area_from_positions(v2, f2, e12);
  patch_areas[8] = area_from_positions(v2, e20, f2);
  patch_areas[9] = area_from_positions(v0, f0, e20);
  patch_areas[10] = area_from_positions(v0, e01, f0);
  patch_areas[11] = area_from_positions(v1, f1, e01);
}

void
compute_twelve_split_spline_patch_boundary_edges(
  const Eigen::MatrixXi& F,
  const std::vector<std::vector<int>>& face_to_patch_indices,
  std::vector<std::pair<int, int>>& patch_boundary_edges)
{
  patch_boundary_edges.clear();
  spdlog::info("Computing patch boundary edges for mesh with {} faces",
               F.rows());

  // Validate input
  if (face_to_patch_indices.size() != static_cast<size_t>(F.rows())) {
    spdlog::error(
      "Incompatible number of mesh faces ({}) and face to patch mappings ({})",
      F.rows(),
      face_to_patch_indices.size());
    return;
  }

  // Get face boundary edges
  std::vector<std::pair<int, int>> face_boundary_edges;
  compute_face_boundary_edges(F, face_boundary_edges);

  // Get boundary patch corners
  patch_boundary_edges.reserve(2 * face_boundary_edges.size());
  for (size_t i = 0; i < face_boundary_edges.size(); ++i) {
    // Get two patch edge corners corresponding to the face edge
    // WARNING: There are some magic numbers here from the construction
    // of the twelve split
    // FIXME Double check these numbers
    int face_index = face_boundary_edges[i].first;
    int face_vertex_index = (face_boundary_edges[i].second + 1) % 3;

    int first_patch_index =
      face_to_patch_indices[face_index][6 + (2 * face_vertex_index)];
    int first_patch_vertex_index = 1;
    int second_patch_index =
      face_to_patch_indices[face_index][7 + (2 * face_vertex_index)];
    int second_patch_vertex_index = 0;

    // Skip faces without a patch
    if ((first_patch_index < 0) || (second_patch_index < 0)) {
      continue;
    }

    // Add patch boundary edges
    patch_boundary_edges.push_back(
      std::make_pair(first_patch_index, first_patch_vertex_index));
    patch_boundary_edges.push_back(
      std::make_pair(second_patch_index, second_patch_vertex_index));
  }
}