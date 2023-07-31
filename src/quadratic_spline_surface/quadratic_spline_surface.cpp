// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "quadratic_spline_surface.h"
#include "convex_polygon.h"
#include "polynomial_function.h"
#include "polyscope/point_cloud.h"
#include "twelve_split_spline.h"
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>

// ************************
// Quadratic Spline Surface
// ************************

// **************
// Public Methods
// **************

QuadraticSplineSurface::QuadraticSplineSurface()
{
  clear();
}

QuadraticSplineSurface::QuadraticSplineSurface(
  std::vector<QuadraticSplineSurfacePatch>& patches)
{
  clear();
  m_patches = patches;
  compute_patch_hash_tables();
}

QuadraticSplineSurface
QuadraticSplineSurface::subsurface(
  const std::vector<PatchIndex>& patch_indices) const
{
  std::vector<QuadraticSplineSurfacePatch> sub_patches;
  sub_patches.clear();
  sub_patches.reserve(patch_indices.size());
  for (size_t i = 0; i < patch_indices.size(); ++i) {
    sub_patches.push_back(m_patches[patch_indices[i]]);
  }

  QuadraticSplineSurface subsurface_spline(sub_patches);
  return subsurface_spline;
}

void
QuadraticSplineSurface::clear()
{
  m_patches.clear();
}

void
QuadraticSplineSurface::triangulate_patch(const PatchIndex& patch_index,
                                          int num_refinements,
                                          Eigen::MatrixXd& V,
                                          Eigen::MatrixXi& F,
                                          Eigen::MatrixXd& N) const
{
  get_patch(patch_index).triangulate(num_refinements, V, F, N);
}

void
QuadraticSplineSurface::discretize(
  const SurfaceDiscretizationParameters& surface_disc_params,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& N) const
{
  V.resize(0, 0);
  F.resize(0, 0);
  N.resize(0, 0);
  std::vector<Eigen::MatrixXd> V_vec(0);
  std::vector<Eigen::MatrixXi> F_vec(0);
  std::vector<Eigen::MatrixXd> N_vec(0);
  int num_subdivisions = surface_disc_params.num_subdivisions;
  if (empty())
    return;

  // Build triangulated surface in place
  PatchIndex patch_index = 0;
  triangulate_patch(patch_index, num_subdivisions, V, F, N);
  int num_patch_vertices = V.rows();
  int num_patch_faces = F.rows();
  V.conservativeResize(num_patch_vertices * num_patches(), 3);
  F.conservativeResize(num_patch_faces * num_patches(), 3);
  N.conservativeResize(num_patch_vertices * num_patches(), 3);
  ++patch_index;

  for (; patch_index < num_patches(); ++patch_index) {
    Eigen::MatrixXd V_patch;
    Eigen::MatrixXi F_patch;
    Eigen::MatrixXd N_patch;
    triangulate_patch(patch_index, num_subdivisions, V_patch, F_patch, N_patch);
    V.block(num_patch_vertices * patch_index, 0, num_patch_vertices, V.cols()) =
      V_patch;
    F.block(num_patch_faces * patch_index, 0, num_patch_faces, F.cols()) =
      F_patch + Eigen::MatrixXi::Constant(
                  num_patch_faces, F.cols(), num_patch_vertices * patch_index);
    N.block(num_patch_vertices * patch_index, 0, num_patch_vertices, N.cols()) =
      N_patch;
  }

  spdlog::info("{} surface vertices", V.rows());
  spdlog::info("{} surface faces", F.rows());
  spdlog::info("{} surface normals", N.rows());
}

std::tuple<Eigen::MatrixXd, // V
           Eigen::MatrixXi, // F
           Eigen::MatrixXd  // N
           >
QuadraticSplineSurface::discretize(
  const SurfaceDiscretizationParameters& surface_disc_params) const
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;
  discretize(surface_disc_params, V, F, N);
  return std::make_tuple(V, F, N);
}

void
QuadraticSplineSurface::discretize_patch_boundaries(
  std::vector<SpatialVector>& points,
  std::vector<std::vector<int>>& polylines) const
{
  points.clear();
  polylines.clear();

  for (PatchIndex patch_index = 0; patch_index < num_patches(); ++patch_index) {
    std::array<LineSegment, 3> patch_boundaries;
    auto& spline_surface_patch = get_patch(patch_index);
    spline_surface_patch.get_domain().parametrize_patch_boundaries(
      patch_boundaries);
    for (size_t k = 0; k < patch_boundaries.size(); ++k) {
      // Get points on the boundary curve
      std::vector<PlanarPoint> parameter_points_k;
      patch_boundaries[k].sample_points(5, parameter_points_k);
      std::vector<SpatialVector> points_k(parameter_points_k.size());
      for (size_t i = 0; i < parameter_points_k.size(); ++i) {
        spline_surface_patch.evaluate(parameter_points_k[i], points_k[i]);
      }

      // Build polyline for the given curve
      std::vector<int> polyline;
      polyline.resize(points_k.size());
      for (size_t l = 0; l < points_k.size(); ++l) {
        polyline[l] = points.size() + l;
      }

      append(points, points_k);
      polylines.push_back(polyline);
    }
  }
}

void
QuadraticSplineSurface::save_obj(const std::string& filename) const
{
  // Generate mesh discretization
  Eigen::MatrixXd V, TC;
  Eigen::MatrixXi F, FTC;
  Eigen::MatrixXd N;
  SurfaceDiscretizationParameters surface_disc_params;
  discretize(surface_disc_params, V, F, N);

  // Write mesh to file
  igl::writeOBJ(filename, V, F, N, F, TC, FTC);
}

void
QuadraticSplineSurface::add_surface_to_viewer(Eigen::Matrix<double, 3, 1> color,
                                              int num_subdivisions) const
{
  // Generate mesh discretization
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;
  SurfaceDiscretizationParameters surface_disc_params;
  surface_disc_params.num_subdivisions = num_subdivisions;
  discretize(surface_disc_params, V, F, N);

  // Add surface mesh
  polyscope::init();
  polyscope::registerSurfaceMesh("surface", V, F)
    ->setEdgeWidth(0);
  polyscope::getSurfaceMesh("surface")->setSurfaceColor(
    glm::vec3(color[0], color[1], color[2]));

  // Discretize patch boundaries
  std::vector<SpatialVector> boundary_points;
  std::vector<std::vector<int>> boundary_polylines;
  discretize_patch_boundaries(boundary_points, boundary_polylines);

  // View contour curve network
  MatrixXr boundary_points_mat =
    convert_nested_vector_to_matrix(boundary_points);
  std::vector<std::array<int, 2>> boundary_edges =
    convert_polylines_to_edges(boundary_polylines);
  polyscope::registerCurveNetwork(
    "patch_boundaries", boundary_points_mat, boundary_edges);
  polyscope::getCurveNetwork("patch_boundaries")
    ->setColor(glm::vec3(0.670, 0.673, 0.292));
  polyscope::getCurveNetwork("patch_boundaries")->setRadius(0.0005);
  polyscope::getCurveNetwork("patch_boundaries")->setRadius(0.0005);
  polyscope::getCurveNetwork("patch_boundaries")->setEnabled(false);
}

void
QuadraticSplineSurface::view(Eigen::Matrix<double, 3, 1> color,
                             int num_subdivisions) const
{
  add_surface_to_viewer(color, num_subdivisions);
  polyscope::show();
}

void
QuadraticSplineSurface::screenshot(const std::string& filename,
                                   SpatialVector camera_position,
                                   SpatialVector camera_target,
                                   bool use_orthographic) const
{
  add_surface_to_viewer();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  glm::vec3 glm_camera_position = { camera_position[0],
                                    camera_position[1],
                                    camera_position[2] };
  glm::vec3 glm_camera_target = { camera_target[0],
                                  camera_target[1],
                                  camera_target[2] };
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  if (use_orthographic) {
    polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  }
  else {
    polyscope::view::projectionMode = polyscope::ProjectionMode::Perspective;
  }
  polyscope::screenshot(filename);
  spdlog::info("Screenshot saved to {}", filename);
}

void
QuadraticSplineSurface::serialize(std::ostream& out) const
{
  for (size_t i = 0; i < m_patches.size(); ++i) {
    m_patches[i].serialize(out);
  }
}

void
QuadraticSplineSurface::deserialize(std::istream& in)
{
  m_patches.clear();

  std::string line;
  while (std::getline(in, line)) {
    // TODO Add better checking and optional comments
    // if (line != "patch")
    //{
    //  spdlog::error("Could not deserialize spline");
    //  m_patches.clear();
    //  return;
    //}

    // Read coordinates
    Matrix6x3r surface_mapping_coeffs;
    for (int coord = 0; coord < 3; ++coord) {
      std::getline(in, line);
      std::istringstream iss(line);

      // Get line label
      std::string label;
      if (!(iss >> label)) {
        spdlog::error("Could not deserialize spline");
        m_patches.clear();
        return;
      }

      // Read coefficients
      for (int i = 0; i < 6; ++i) {
        if (!(iss >> surface_mapping_coeffs(i, coord))) {
          spdlog::error("Could not deserialize spline");
          m_patches.clear();
          return;
        }
      }
    }

    // Read convex domain vertices
    Matrix3x2r vertices;
    for (int i = 0; i < 3; ++i) {
      std::getline(in, line);
      std::istringstream iss(line);

      // Get line label
      std::string label;
      if (!(iss >> label)) {
        spdlog::error("Could not deserialize spline");
        m_patches.clear();
        return;
      }

      for (int j = 0; j < 2; ++j) {
        if (!(iss >> vertices(i, j))) {
          spdlog::error("Could not deserialize spline");
          m_patches.clear();
          return;
        }
      }
    }
    ConvexPolygon domain(vertices);

    // Add patch to the spline surface
    m_patches.push_back(
      QuadraticSplineSurfacePatch(surface_mapping_coeffs, domain));
  }
}

void
QuadraticSplineSurface::write_spline(const std::string& filename) const
{
  spdlog::info("Writing spline to {}", filename);
  std::ofstream output_file(filename, std::ios::out | std::ios::trunc);
  serialize(output_file);
  output_file.close();
}

void
QuadraticSplineSurface::read_spline(const std::string& filename)
{
  std::ifstream input_file;
  input_file.open(filename);
  deserialize(input_file);
  input_file.close();
}

std::pair<int, int>
QuadraticSplineSurface::compute_hash_indices(const PlanarPoint& point) const
{
  int hash_x = (point[0] - patches_bbox_x_min) / hash_x_interval;
  int hash_y = (point[1] - patches_bbox_y_min) / hash_y_interval;
  if ((hash_x < 0) || (hash_x >= HASH_TABLE_SIZE)) {
    spdlog::error("x hash index out of bounds");
    hash_x = std::max(std::min(hash_x, HASH_TABLE_SIZE - 1), 0);
  }
  if ((hash_y < 0) || (hash_y >= HASH_TABLE_SIZE)) {
    spdlog::error("y hash index out of bounds");
    hash_y = std::max(std::min(hash_y, HASH_TABLE_SIZE - 1), 0);
  }
  return std::make_pair(hash_x, hash_y);
}

void
QuadraticSplineSurface::compute_patch_hash_tables()
{
  int num_patch = num_patches();
  int hash_size_x = HASH_TABLE_SIZE;
  int hash_size_y = HASH_TABLE_SIZE;

  // clear the hash table;
  for (int i = 0; i < hash_size_x; i++) {
    for (int j = 0; j < hash_size_y; j++) {
      hash_table[i][j].clear();
    }
  }

  // compute bounding box for all the patches
  compute_patches_bbox();
  double x_min, x_max, y_min, y_max;
  x_min = patches_bbox_x_min;
  x_max = patches_bbox_x_max;
  y_min = patches_bbox_y_min;
  y_max = patches_bbox_y_max;

  for (int i = 1; i < num_patch; i++) {
    if (x_min > m_patches[i].get_bbox_x_min())
      x_min = m_patches[i].get_bbox_x_min();
    if (x_max < m_patches[i].get_bbox_x_max())
      x_max = m_patches[i].get_bbox_x_max();
    if (y_min > m_patches[i].get_bbox_y_min())
      y_min = m_patches[i].get_bbox_y_min();
    if (y_max < m_patches[i].get_bbox_y_max())
      y_max = m_patches[i].get_bbox_y_max();
  }

  double x_interval = (x_max - x_min) / hash_size_x;
  double y_interval = (y_max - y_min) / hash_size_y;

  hash_x_interval = x_interval;
  hash_y_interval = y_interval;

  // reverse_hash_table.reserve(num_patch);

  double eps = 1e-10;

  // hash into each box
  for (int i = 0; i < num_patch; i++) {
    int left_x = (m_patches[i].get_bbox_x_min() - eps - x_min) / x_interval;
    int right_x =
      hash_size_x -
      int((x_max - m_patches[i].get_bbox_x_max() - eps) / x_interval) - 1;
    int left_y = (m_patches[i].get_bbox_y_min() - eps - y_min) / y_interval;
    int right_y =
      hash_size_y -
      int((y_max - m_patches[i].get_bbox_y_max() - eps) / y_interval) - 1;

    for (int j = left_x; j <= right_x; j++) {
      for (int k = left_y; k <= right_y; k++) {
        hash_table[j][k].push_back(i);
        // reverse_hash_table[i].push_back(std::make_pair(j, k));
      }
    }
  }
}

// ***************
// Private Methods
// ***************

// Determine if a patch index is valid
bool
QuadraticSplineSurface::is_valid_patch_index(
  QuadraticSplineSurface::PatchIndex patch_index) const
{
  if (patch_index < 0)
    return false;
  if (patch_index >= num_patches())
    return false;

  return true;
}

// Compute bounding boxes for the patches
void
QuadraticSplineSurface::compute_patches_bbox()
{
  double x_min, x_max, y_min, y_max;
  x_min = m_patches[0].get_bbox_x_min();
  x_max = m_patches[0].get_bbox_x_max();
  y_min = m_patches[0].get_bbox_y_min();
  y_max = m_patches[0].get_bbox_y_max();

  for (PatchIndex i = 1; i < num_patches(); i++) {
    if (x_min > m_patches[i].get_bbox_x_min())
      x_min = m_patches[i].get_bbox_x_min();
    if (x_max < m_patches[i].get_bbox_x_max())
      x_max = m_patches[i].get_bbox_x_max();
    if (y_min > m_patches[i].get_bbox_y_min())
      y_min = m_patches[i].get_bbox_y_min();
    if (y_max < m_patches[i].get_bbox_y_max())
      y_max = m_patches[i].get_bbox_y_max();
  }

  patches_bbox_x_min = x_min;
  patches_bbox_x_max = x_max;
  patches_bbox_y_min = y_min;
  patches_bbox_y_max = y_max;
}
