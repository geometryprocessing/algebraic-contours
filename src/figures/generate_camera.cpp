// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "apply_transformation.h"
#include "common.h"
#include "generate_transformation.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include <igl/readOBJ.h>
#include <CLI/CLI.hpp>


int main(int argc, char *argv[]) {
  // Get command line arguments
  CLI::App app{"Generate perspective figure example for a given mesh."};
  std::string input_filename = "";
  app.add_option("-i,--input", input_filename, "Mesh filepath")
    ->check(CLI::ExistingFile)
    ->required();
  CLI11_PARSE(app, argc, argv);

  // Set logger level
  spdlog::set_level(spdlog::level::off);

  // Get input mesh
  Eigen::MatrixXd V, uv, N, orthographic_V;
  Eigen::MatrixXi F, FT, FN;
  igl::readOBJ(input_filename, V, uv, N, F, FT, FN);

  // Build perspective viewer
  polyscope::init();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

  // Look at viewer
  glm::vec3 glm_camera_position = {0.0, 0.0, 0.0};
  glm::vec3 glm_camera_target = {0.0, 0.0, 4.0};
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::show();
  polyscope::removeAllStructures();

  // Project mesh and frustum to send the camera to infinity while fixing the
  // plane z = 1
  MatrixXr projection_matrix = origin_to_infinity_projective_matrix(1.0);

  // Build orthographic viewer
  polyscope::registerSurfaceMesh("original_mesh", V, F);
  polyscope::getSurfaceMesh("original_mesh")
      ->setSurfaceColor(glm::vec3(0.670, 0.673, 0.292));

  float rx = 162.398;
  float ry = -47.3441;
  float rz = 181.622;
  float tx = -2.20568;
  float ty = 1.21006;
  float tz = -1.91432;
  auto viewer_callback = [&]() {
    ImGui::InputFloat("rx", &rx);
    ImGui::InputFloat("ry", &ry);
    ImGui::InputFloat("rz", &rz);
    ImGui::InputFloat("tx", &tx);
    ImGui::InputFloat("ty", &ty);
    ImGui::InputFloat("tz", &tz);
    if (ImGui::Button("Project")) {
      MatrixXr projective_matrix = blender_rotation_projective_matrix(rx, ry, rz, tx, ty, tz);
      apply_transformation_to_vertices(V, projective_matrix, orthographic_V);
      polyscope::registerSurfaceMesh("orthographic_mesh", orthographic_V, F);
    }
    if (ImGui::Button("Write Matrix")) {
      MatrixXr projective_matrix = blender_rotation_projective_matrix(rx, ry, rz, tx, ty, tz);
      write_camera_matrix("camera.txt", projective_matrix);
    }
  };
  polyscope::state::userCallback = viewer_callback;

  // Look at viewer
  glm_camera_position = { 0.0, 0.0, 0.0 };
  glm_camera_target = { 0.0, 0.0, 1.0 };
  //polyscope::view::projectionMode = polyscope::ProjectionMode::Orthographic;
  polyscope::view::lookAt(glm_camera_position, glm_camera_target);
  polyscope::show();
  polyscope::removeAllStructures();
}
