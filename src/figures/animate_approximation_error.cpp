// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "common.h"
#include "apply_transformation.h"
#include "compute_boundaries.h"
#include "contour_network.h"
#include "generate_transformation.h"
#include "globals.cpp"
#include "twelve_split_spline.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <CLI/CLI.hpp>

/// Executable to animate the interpolation of a rotating figure by projecting, interpolating, and
/// inverting the projection to give a fixed camera perspective of the interpolation.

void normalize_and_rotate_vertices(
  const MatrixXr &input_V,
  const Matrix3x3r &frame,
  Eigen::MatrixXd &output_V
) {
  // Compute mesh midpoint and bounding box diagonal
  SpatialVector min_point;
  SpatialVector max_point;
  compute_point_cloud_bounding_box(input_V, min_point, max_point);
  SpatialVector mesh_midpoint = 0.5 * (max_point + min_point);
  SpatialVector bounding_box_diagonal = max_point - min_point;
  spdlog::info("Initial mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Initial mesh midpoint: {}", mesh_midpoint);

  // Normalize the vertices
  double scale_factor =  bounding_box_diagonal.norm();
  size_t num_vertices = input_V.rows();
  output_V.resize(num_vertices, 3);
  for (size_t i = 0; i < num_vertices; ++i)
  {
    output_V.row(i) = (input_V.row(i) - mesh_midpoint) / scale_factor;
  }

  // Generate rotation matrix
  spdlog::info("Projecting onto frame:\n{}", frame);
  Eigen::Matrix<double, 4, 4> frame_rotation_matrix = rotate_frame_projective_matrix(frame);
  
  // Generate translation matrix
  double z_distance = 1.1;
  SpatialVector translation(3);
  translation << 0.0,
                 0.0,
                 z_distance;
  Eigen::Matrix<double, 4, 4> translation_matrix = translation_projective_matrix(translation);

  // Apply the transformations
  Eigen::Matrix<double, 4, 4> projective_transformation = translation_matrix * frame_rotation_matrix;
  spdlog::info("Apply transformation:\n{}", projective_transformation);
  apply_transformation_to_vertices_in_place(output_V, projective_transformation);

  // Check final midpoint location
  compute_point_cloud_bounding_box(output_V, min_point, max_point);
  mesh_midpoint = 0.5 * (max_point + min_point);
  spdlog::info("Final mesh bounding box: {}, {}", min_point, max_point);
  spdlog::info("Final mesh midpoint: {}", mesh_midpoint);
}

void project_vertices(
  const MatrixXr &input_V,
  Eigen::MatrixXd &output_V
) {
  size_t num_vertices = input_V.rows();

  // Generate matrix to send the origin to infinity
  double camera_to_plane_distance = 1.0;
  Eigen::Matrix<double, 4, 4> projection_matrix = origin_to_infinity_projective_matrix(
    camera_to_plane_distance
  );

  // Apply the transformation
  output_V = input_V;
  Eigen::Matrix<double, 4, 4> projective_transformation = projection_matrix;
  spdlog::info("Apply transformation:\n{}", projective_transformation);
  apply_transformation_to_vertices_in_place(output_V, projective_transformation);

  // Reflect z coordinate
  for (size_t i = 0; i < num_vertices; ++i)
  {
    //output_V(i, 2) = -output_V(i, 2); // Reflect z coordinate
  }
}

void unproject_vertices(
  const MatrixXr &input_V,
  Eigen::MatrixXd &output_V
) {
  size_t num_vertices = input_V.rows();

  // Generate matrix to send the origin to infinity
  double camera_to_plane_distance = 1.0;
  Eigen::Matrix<double, 4, 4> projection_matrix = origin_to_infinity_projective_matrix(
    camera_to_plane_distance
  );

  // Reflect z coordinate
  output_V = input_V;
  for (size_t i = 0; i < num_vertices; ++i)
  {
    //output_V(i, 2) = -output_V(i, 2); // Reflect z coordinate
  }

  // Apply the inverse transformation
  Eigen::Matrix<double, 4, 4> projective_transformation = projection_matrix.inverse();
  spdlog::info("Apply transformation:\n{}", projective_transformation);
  apply_transformation_to_vertices_in_place(output_V, projective_transformation);
}

// Function to compute the pointwise distance between two vertex sets
void
compute_vertex_distance(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::VectorXd& distances
) {
  int num_vertices = V.rows();
  if (W.rows() != num_vertices) return;

  // Compute distances by row
  distances.resize(num_vertices);
  for (int i = 0; i < num_vertices; ++i)
  {
    distances[i] = (V.row(i) - W.row(i)).norm();
  }
}

// Helper function to compute the distance between two planar projected points
// with orthographic projection in the z direction
double
compute_vector_planar_distance(
  const Eigen::VectorXd& v,
  const Eigen::VectorXd& w
) {
  Eigen::VectorXd d = v - w;
  return std::sqrt(d[0] * d[0] + d[1] * d[1]);
}

// Compute the pointwise distance between two planar projected vertex sets
// with orthographic projection in the z direction
void
compute_vertex_planar_distance(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::VectorXd& planar_distances
) {
  int num_vertices = V.rows();
  if (W.rows() != num_vertices) return;

  // Compute distances by row
  planar_distances.resize(num_vertices);
  for (int i = 0; i < num_vertices; ++i)
  {
    planar_distances[i] = compute_vector_planar_distance(V.row(i), W.row(i));
  }
}

// Screenshot a surface with a distance scalar field
void
screenshot_surface(
  const std::string& filename,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXd& distances
) {
  polyscope::registerSurfaceMesh("surface", V, F);
  polyscope::getSurfaceMesh("surface")
    ->addVertexScalarQuantity("distances", distances)
    ->setColorMap("reds")
    ->setMapRange(std::make_pair(0, 1e-4))
    ->setEnabled(true);
  glm::vec3 glm_camera_position = { 0, 0, 0 };
  glm::vec3 glm_camera_target = { 0, 0, 1 };
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::screenshot(filename);
  polyscope::removeAllStructures();
}

// Screenshot a surface with flat shading
void
screenshot_planar_surface(
  const std::string& filename,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
) {
  polyscope::registerSurfaceMesh("planar surface", V, F)
    ->setMaterial("flat")
    ->setSurfaceColor(glm::vec3(0, 0, 0));
  glm::vec3 glm_camera_position = { 0, 0, 0 };
  glm::vec3 glm_camera_target = { 0, 0, 1 };
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::screenshot(filename);
  polyscope::removeAllStructures();
}

int main(int argc, char *argv[]) {
  CLI::App app{"Generate approximation error animation for a given mesh."};
  std::string input_filename = "";
  std::string output_dir = "./";
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  app.add_option("-o,--output", output_dir, "Output directory")
    ->check(CLI::ExistingDirectory);
  CLI11_PARSE(app, argc, argv);

  // Set logging and discretization level
  DISCRETIZATION_LEVEL = 2;
  spdlog::set_level(spdlog::level::off);

	// Initialize polyscope
  polyscope::init();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

  // Get input mesh
  Eigen::MatrixXd V, uv, N;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);
  MatrixXr initial_V = V;

  // Generate quadratic spline
  OptimizationParameters optimization_params;
  std::vector<std::vector<int>> face_to_patch_indices;
  std::vector<int> patch_to_face_indices;
  Eigen::SparseMatrix<double> fit_matrix;
  Eigen::SparseMatrix<double> energy_hessian;
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>
      energy_hessian_inverse;
  AffineManifold affine_manifold(F, uv, FT);
  TwelveSplitSplineSurface spline_surface(
      initial_V, affine_manifold,
      optimization_params, face_to_patch_indices, patch_to_face_indices,
      fit_matrix, energy_hessian, energy_hessian_inverse);

  // Get perturbation frame
  MatrixXr pert_frame(3, 3);
  pert_frame << 0.9218127938, -0.3785297887, 0.08352467839, 0.3815369887,
       0.9240681946, -0.0229673378, -0.06848867707, 0.0530393400,
       0.9962410000;

  // Iterate through rotation frames
  SurfaceDiscretizationParameters surface_disc_params;
  surface_disc_params.num_subdivisions = 2;
  double pert = 1e-2;
  double theta = 0;
  double dx = 1;
  int num_frames = 360;
  for (int i = 0; i < num_frames; ++i)
  {
    // Rotate the mesh
    theta = dx * i + pert;
    MatrixXr rotation_matrix = y_axis_rotation_projective_matrix(theta);
    Matrix3x3r frame = rotation_matrix.block(0, 0, 3, 3) * pert_frame;
    Eigen::MatrixXd V_rotated;
    normalize_and_rotate_vertices(initial_V, frame, V_rotated);

    // Get the unprojected spline surface
    spline_surface.update_positions(
        V_rotated, fit_matrix, energy_hessian_inverse);

    // Triangulate the unprojected surface (in perspective space)
    Eigen::MatrixXd V_unprojected_perspective;
    Eigen::MatrixXi F_unprojected;
    Eigen::MatrixXd N_unprojected;
    spline_surface.discretize(
      surface_disc_params,
      V_unprojected_perspective,
      F_unprojected,
      N_unprojected
    );

    // Project the input vertices to orthographic space
    Eigen::MatrixXd V_orthogonal;
    project_vertices(V_rotated, V_orthogonal);

    // Get the projected vertex spline surface
    spline_surface.update_positions(
        V_orthogonal, fit_matrix, energy_hessian_inverse);

    // Triangulate the projected surface (in orthographic space)
    Eigen::MatrixXd V_projected_orthogonal;
    Eigen::MatrixXi F_projected;
    Eigen::MatrixXd N_projected;
    spline_surface.discretize(
      surface_disc_params,
      V_projected_orthogonal,
      F_projected,
      N_projected
    );

    // Unproject the projected surface triangulation to perspective space
    Eigen::MatrixXd V_projected_perspective;
    unproject_vertices(V_projected_orthogonal, V_projected_perspective);

    // Project the unprojected surface triangulation to orthographic space
    Eigen::MatrixXd V_unprojected_orthogonal;
    project_vertices(V_unprojected_perspective, V_unprojected_orthogonal);

    // Get the distances between the mesh vertices in perspective space
    Eigen::VectorXd distances;
    compute_vertex_distance(V_unprojected_perspective, V_projected_perspective, distances);

    // Write unprojected surface output
    std::string output_path;
    output_path = join_path(output_dir, "unprojected_frame_" + std::to_string(i) + ".png");
    screenshot_surface(output_path, V_unprojected_perspective, F_unprojected, distances);
    output_path = join_path(output_dir, "unprojected_planar_frame_" + std::to_string(i) + ".png");
    screenshot_planar_surface(output_path, V_unprojected_perspective, F_unprojected);

    // Write projected surface output
    output_path = join_path(output_dir, "projected_frame_" + std::to_string(i) + ".png");
    screenshot_surface(output_path, V_projected_perspective, F_projected, distances);
    output_path = join_path(output_dir, "projected_planar_frame_" + std::to_string(i) + ".png");
    screenshot_planar_surface(output_path, V_projected_perspective, F_projected);

  }
}
