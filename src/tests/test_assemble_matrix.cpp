#include "affine_manifold.h"
#include "common.h"
#include "compute_local_twelve_split_hessian.h"
#include "optimize_spline_surface.h"
#include "twelve_split_spline.h"
#include "generate_position_data.h"
#include "generate_shapes.h"
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <iostream>
#include <vector>

/// Purpose: We currently have multiple methods of computing optimized quadratic
/// surfaces. In order to validate each and more easily resolve inconsistencies,
/// we run the current default method using autodiff variables on a quadratic
/// surface with extensive logging.

// Display quadratic surface mapping output to console
//void display_quadratic_surface_mapping(const Eigen::MatrixXd &uv,
//                                       const Eigen::MatrixXi &F,
//                                       double uv_coeff, double uu_coeff,
//                                       double vv_coeff) {
//  spdlog::info("Using quadratic {} * uv + {} * u^2 + {} * v^2", uv_coeff,
//               uu_coeff, vv_coeff);
//
//  // Build parametric domain manifold
//  ParametricAffineManifold parametric_affine_manifold(F, uv);
//  spdlog::info("Using uv domain triangle\n{}", uv);
//
//  // Build halfedge for edge index management
//  std::vector<std::vector<size_t>> corner_to_he;
//  std::vector<std::pair<size_t, size_t>> he_to_corner;
//  Halfedge halfedge(F, corner_to_he, he_to_corner);
//
//  // Generate functors to evaluate the quadratic and its gradients
//  QuadraticPositionFunction position_func(uv_coeff, uu_coeff, vv_coeff);
//  QuadraticGradientFunction gradient_func(uv_coeff, uu_coeff, vv_coeff);
//
//  // Display function corner and midpoint data
//  std::vector<std::array<TriangleCornerFunctionData, 3>> corner_data;
//  std::vector<std::array<TriangleMidpointFunctionData, 3>> midpoint_data;
//  generate_parametric_affine_manifold_corner_data(
//      position_func, gradient_func, parametric_affine_manifold, corner_data);
//  generate_parametric_affine_manifold_midpoint_data(
//      gradient_func, parametric_affine_manifold, midpoint_data);
//  for (size_t i = 0; i < 3; ++i) {
//    spdlog::info("Data for corner {}", i);
//    spdlog::info("Fuction value: {}",
//                 corner_data[0][i].function_value.transpose());
//    spdlog::info("First edge derivative: {}",
//                 corner_data[0][i].first_edge_derivative.transpose());
//    spdlog::info("Second edge derivative: {}",
//                 corner_data[0][i].second_edge_derivative.transpose());
//  }
//  for (size_t i = 0; i < 3; ++i) {
//    spdlog::info("Data for midpoint opposite corner {}", i);
//    spdlog::info("Midpoint derivative: {}",
//                 midpoint_data[0][i].normal_derivative.transpose());
//  }
//
//  // Generate the surface mapping
//  std::array<Eigen::Matrix<double, 6, 3>, 12> surface_mappings;
//  generate_twelve_split_spline_patch_surface_mapping<double>(
//      corner_data[0], midpoint_data[0], surface_mappings);
//  for (size_t i = 0; i < 12; ++i) {
//    spdlog::info("Surface mapping for patch {}:\n{}", i, surface_mappings[i]);
//  }
//
//  // Reorganize surface mapping into x, y, z coordinates
//  Eigen::Matrix<double, 12, 6> x_mat, y_mat, z_mat;
//  for (size_t i = 0; i < 12; ++i) {
//    for (size_t j = 0; j < 6; ++j) {
//      x_mat(i, j) = surface_mappings[i](j, 0);
//      y_mat(i, j) = surface_mappings[i](j, 1);
//      z_mat(i, j) = surface_mappings[i](j, 2);
//    }
//  }
//  spdlog::info("x coefficient matrix\n{}", x_mat);
//  spdlog::info("y coefficient matrix\n{}", y_mat);
//  spdlog::info("z coefficient matrix\n{}", z_mat);
//}
//
//bool generate_twelve_split_quadratic_surface_mapping(
//    const Eigen::MatrixXd &uv, const Eigen::MatrixXi &F, double uv_coeff,
//    double uu_coeff, double vv_coeff,
//    std::array<SpatialVector, 3> &vertex_positions_T,
//    std::array<Matrix2x3r, 3> &vertex_gradients_T,
//    std::array<SpatialVector, 3> &edge_gradients_T,
//    std::array<SpatialVector, 3> &initial_vertex_positions_T,
//    std::array<PlanarPoint, 3> &face_vertex_uv_positions,
//    std::array<Matrix2x2r, 3> &corner_to_corner_uv_positions,
//    std::array<PlanarPoint, 3> &midpoint_to_corner_uv_positions,
//    std::array<bool, 3> &reverse_edge_orientations) {
//  spdlog::info("Using quadratic {} * uv + {} * u^2 + {} * v^2", uv_coeff,
//               uu_coeff, vv_coeff);
//
//  // Build parametric domain manifold
//  ParametricAffineManifold parametric_affine_manifold(F, uv);
//  spdlog::info("Using uv domain triangle\n{}", uv);
//
//  // Build halfedge for edge index management
//  std::vector<std::vector<size_t>> corner_to_he;
//  std::vector<std::pair<size_t, size_t>> he_to_corner;
//  Halfedge halfedge(F, corner_to_he, he_to_corner);
//
//  // Generate functors to evaluate the quadratic and its gradients
//  QuadraticPositionFunction position_func(uv_coeff, uu_coeff, vv_coeff);
//  QuadraticGradientFunction gradient_func(uv_coeff, uu_coeff, vv_coeff);
//
//  // Get position and gradient values
//  MatrixXr V;
//  std::vector<Matrix2x3r> vertex_gradients;
//  std::vector<std::array<SpatialVector, 3>> edge_gradients;
//  generate_parametric_affine_manifold_vertex_positions(
//      position_func, parametric_affine_manifold, V);
//  generate_parametric_affine_manifold_vertex_gradients(
//      gradient_func, parametric_affine_manifold, vertex_gradients);
//  generate_parametric_affine_manifold_edge_gradients(
//      gradient_func, parametric_affine_manifold, halfedge, edge_gradients);
//
//  // Display values
//  spdlog::info("Vertex positions:\n{}", V);
//  for (size_t i = 0; i < vertex_gradients.size(); ++i) {
//    spdlog::info("Vertex gradients for vertex {}:\n{}", i, vertex_gradients[i]);
//  }
//  for (size_t i = 0; i < 3; ++i) {
//    spdlog::info("Edge gradient opposite vertex {}: {}", i,
//                 edge_gradients[0][i].transpose());
//  }
//
//  // Display uv coordinates
//  // spdlog::info("Affine manifold uv positions:");
//  // PlanarPoint uv_coords;
//  // parametric_affine_manifold.get_vertex_global_uv(0, uv_coords);
//  // spdlog::info(uv_coords.transpose());
//  // parametric_affine_manifold.get_vertex_global_uv(1, uv_coords);
//  // spdlog::info(uv_coords.transpose());
//  // parametric_affine_manifold.get_vertex_global_uv(2, uv_coords);
//  // spdlog::info(uv_coords.transpose());
//
//  // Compute the hessian for this configuration
//  OptimizationParameters optimization_params;
//  optimization_params.parametrized_quadratic_surface_mapping_factor = 1;
//  TriangleEnergy<TwelveSplitVar> triangle_energy(optimization_params);
//
//  // Build vertex positions and gradients
//  int i = 0;
//  int j = 1;
//  int k = 2;
//  int face_index = 0;
//  size_t num_vertices = V.rows();
//  std::vector<SpatialVector> vertex_positions(num_vertices);
//  std::vector<SpatialVector> initial_vertex_positions(num_vertices);
//  for (size_t i = 0; i < num_vertices; ++i) {
//    vertex_positions[i] = V.row(i);
//    initial_vertex_positions[i] = V.row(i);
//  }
//
//  build_face_variable_vector(vertex_positions, i, j, k, vertex_positions_T);
//  build_face_variable_vector(vertex_gradients, i, j, k, vertex_gradients_T);
//  build_face_variable_vector(initial_vertex_positions, i, j, k,
//                             initial_vertex_positions_T);
//  edge_gradients_T = edge_gradients[face_index];
//
//  // Get the global uv values for the face vertices
//  parametric_affine_manifold.get_face_global_uv(face_index,
//                                                face_vertex_uv_positions);
//
//  // Get corner uv positions for the given face corners
//  parametric_affine_manifold.get_face_corner_charts(
//      face_index, corner_to_corner_uv_positions);
//
//  // Get corner uv positions for the given face edges
//  std::array<MatrixXr, 3> face_edge_uv_positions;
//  parametric_affine_manifold.get_face_edge_charts(face_index,
//                                                  face_edge_uv_positions);
//
//  // Get only the uv position from the edge midpoint to the opposite corner
//  for (size_t i = 0; i < 3; ++i) {
//    midpoint_to_corner_uv_positions[i] = face_edge_uv_positions[i].row(1);
//    spdlog::info("Midpoint to corner uv position {}: {}", i,
//                 midpoint_to_corner_uv_positions[i]);
//  }
//
//  // Get edge orientations
//  for (size_t i = 0; i < 3; ++i) {
//    EdgeManifoldChart const &chart =
//        parametric_affine_manifold.get_edge_chart(face_index, i);
//    reverse_edge_orientations[i] = (chart.top_face_index != face_index);
//  }
//
//  //   double local_energy;
//  //   TwelveSplitGradient local_derivatives;
//  //   TwelveSplitHessian local_hessian;
//  //   compute_local_twelve_split_energy_quadratic(
//  //       vertex_positions_T, vertex_gradients_T, edge_gradients_T,
//  //       initial_vertex_positions_T, face_vertex_uv_positions,
//  //       corner_to_corner_uv_positions, midpoint_to_corner_uv_positions,
//  //       reverse_edge_orientations, optimization_params, local_energy,
//  //       local_derivatives, local_hessian);
//
//  // Convert values to differentiable variables
//  // WARNING Old method
//  //	std::array<TwelveSplitVector3v, 3> face_vertex_positions,
//  // face_edge_gradients; 	std::array<TwelveSplitMatrix2x3v, 3>
//  // face_vertex_gradients;
//  //	initialize_local_twelve_split_variables<TwelveSplitVar,
//  // TwelveSplitVector3v, TwelveSplitMatrix2x3v>(
//  // vertex_positions_T, 		vertex_gradients_T,
//  // edge_gradients_T,
//  //		face_vertex_positions,
//  //		face_vertex_gradients,
//  //		face_edge_gradients
//  //	);
//  //
//  //	// Add energy for the current face to the total energy
//  //	TwelveSplitVar face_energy = triangle_energy(
//  //		face_vertex_positions,
//  //		face_vertex_gradients,
//  //		face_edge_gradients,
//  //		initial_vertex_positions_T,
//  //		face_vertex_uv_positions,
//  //		corner_to_corner_uv_positions,
//  //		midpoint_to_corner_uv_positions,
//  //		reverse_edge_orientations
//  //	);
//  //  double local_energy;
//  //  TwelveSplitGradient local_derivatives;
//  //  TwelveSplitHessian local_hessian;
//  //  local_energy = compute_variable_value<TwelveSplitVar>(face_energy);
//  //  compute_variable_gradient<TwelveSplitGradient,
//  //  TwelveSplitVar>(face_energy, local_derivatives);
//  //  compute_variable_hessian<TwelveSplitHessian, TwelveSplitVar>(face_energy,
//  //  local_hessian);
//
//  // Display output
//  //   spdlog::info("Local face energy {}", local_energy);
//  //   spdlog::info("Local face derivative {}", local_derivatives);
//  //   spdlog::info("Local face hessian {}", local_hessian);
//
//  //   std::array<double, 12> patch_areas;
//  //   generate_twelve_split_domain_areas(uv.row(0), uv.row(1), uv.row(2),
//  //                                      patch_areas);
//
//  //   for (size_t i = 0; i < 12; ++i) {
//  //     spdlog::info("Face patch areas: {}", patch_areas[i]);
//  //   }
//}
//
//TEST_CASE("generate_input_data", "[assemble_hessian]") {
//  spdlog::set_level(spdlog::level::info);
//
//  // Build affine manifold domain
//  Eigen::MatrixXd uv(3, 2);
//  Eigen::MatrixXi F(1, 3);
//  double in1, in2, in3, in4, in5, in6;
//  double c_uv, c_uu, c_vv;
//
//  // uv <<
//  //   1.0,  0.0,
//  //   0.0,  1.0,
//  //   0.0,  0.0;
//  //   uv << 2.0, 0.0, 0.0, 2.0, 0.0, 0.0;
//
//  std::cin >> in1 >> in2 >> in3 >> in4 >> in5 >> in6;
//  std::cin >> c_uv >> c_uu >> c_vv;
//  uv << in1, in2, in3, in4, in5, in6;
//  std::cout << uv << std::endl;
//  //   uv << 1.0, 0.0, 0.0, 2.0, 0.0, 0.0;
//  //   uv << 0.8, 0.02, 0.1, 0.85, 0.5, -0.05;
//  F << 0, 1, 2;
//
//  std::array<SpatialVector, 3> vertex_positions_T;
//  std::array<Matrix2x3r, 3> vertex_gradients_T;
//  std::array<SpatialVector, 3> edge_gradients_T;
//  std::array<SpatialVector, 3> initial_vertex_positions_T;
//  std::array<PlanarPoint, 3> face_vertex_uv_positions;
//  std::array<Matrix2x2r, 3> corner_to_corner_uv_positions;
//  std::array<PlanarPoint, 3> midpoint_to_corner_uv_positions;
//  std::array<bool, 3> reverse_edge_orientations;
//  OptimizationParameters optimization_params;
//  optimization_params.parametrized_quadratic_surface_mapping_factor = 1;
//  optimization_params.position_difference_factor = 1000;
//
//  //   generate_twelve_split_quadratic_surface_mapping(uv, F, 0.0, 0.0,
//  //   0.0);
//  generate_twelve_split_quadratic_surface_mapping(
//      uv, F, c_uv, c_uu, c_vv, vertex_positions_T, vertex_gradients_T,
//      edge_gradients_T, initial_vertex_positions_T, face_vertex_uv_positions,
//      corner_to_corner_uv_positions, midpoint_to_corner_uv_positions,
//      reverse_edge_orientations);
//  //   generate_twelve_split_quadratic_surface_mapping(uv, F, 0.0, 1.0, 0.0);
//  //   generate_twelve_split_quadratic_surface_mapping(uv, F, 0.0, 0.0, 1.0);
//  //   generate_twelve_split_quadratic_surface_mapping(uv, F, 1.0, 2.0, 3.0);
//
//  Eigen::Matrix<double, 12, 3> r;
//  r.row(0) = vertex_positions_T[0];
//  r.row(1) = vertex_positions_T[1];
//  r.row(2) = vertex_positions_T[2];
//  r.row(3) = vertex_gradients_T[0].row(0);
//  r.row(4) = vertex_gradients_T[0].row(1);
//  r.row(5) = vertex_gradients_T[1].row(0);
//  r.row(6) = vertex_gradients_T[1].row(1);
//  r.row(7) = vertex_gradients_T[2].row(0);
//  r.row(8) = vertex_gradients_T[2].row(1);
//  r.row(9) = edge_gradients_T[0];
//  r.row(10) = edge_gradients_T[1];
//  r.row(11) = edge_gradients_T[2];
//
//  Eigen::Matrix<double, 3, 2> uv_p;
//  uv_p.row(0) = face_vertex_uv_positions[0];
//  uv_p.row(1) = face_vertex_uv_positions[1];
//  uv_p.row(2) = face_vertex_uv_positions[2];
//
//  std::cout << "uv_p:" << std::endl;
//  std::cout << uv_p << std::endl;
//
//  auto C_gl = get_C_gl(uv_p, corner_to_corner_uv_positions, reverse_edge_orientations);
//  std::cout << "C_gl:" << std::endl;
//  std::cout << C_gl << std::endl << std::endl;
//
//  std::cout << "C_gl * r:" << std::endl;
//  std::cout << C_gl * r << std::endl << std::endl;
//
//  auto Tquad_Cder_Csub = get_Tquad_Cder_Csub();
//
//  std::cout << "Tquad_Cder_Csub * C_gl * r:" << std::endl;
//  std::cout << Tquad_Cder_Csub * C_gl * r << std::endl << std::endl;
//
//  auto R_quad = get_R_quad(uv_p);
//
//  std::cout << "R_quad:" << std::endl;
//  std::cout << R_quad << std::endl << std::endl;
//
//  auto S = get_S(uv_p);
//
//  std::cout << "S:" << std::endl;
//  std::cout << S << std::endl << std::endl;
//
//  std::cout << "G * r:" << std::endl;
//  std::cout << R_quad * Tquad_Cder_Csub * C_gl * r << std::endl << std::endl;
//
//  // autodiff
//  TwelveSplitGradient local_derivatives_autodiff;
//  TwelveSplitHessian local_hessian_autodiff;
//  double local_energy_autodiff;
//
//  compute_local_twelve_split_energy_quadratic(
//      vertex_positions_T, vertex_gradients_T, edge_gradients_T,
//      initial_vertex_positions_T, face_vertex_uv_positions,
//      corner_to_corner_uv_positions, midpoint_to_corner_uv_positions,
//      reverse_edge_orientations, optimization_params, local_energy_autodiff,
//      local_derivatives_autodiff, local_hessian_autodiff);
//
//  std::cout << "energy autodiff: " << local_energy_autodiff << std::endl
//            << std::endl;
//
//  std::cout << "derivative autodiff: " << std::endl
//            << local_derivatives_autodiff << std::endl
//            << std::endl;
//
//  std::cout << "hessian autodiff: " << std::endl
//            << local_hessian_autodiff << std::endl
//            << std::endl;
//
//  // assemble
//  TwelveSplitGradient local_derivatives_assemble;
//  TwelveSplitHessian local_hessian_assemble;
//  double local_energy_assemble;
//
//  std::array<bool, 3> is_cone = {false, false, false};
//  std::array<bool, 3> is_cone_adjacent = {false, false, false};
//  SpatialVector normal;
//  normal.setZero();
//  compute_local_twelve_split_energy_quadratic_manually(
//      vertex_positions_T, vertex_gradients_T, edge_gradients_T,
//      initial_vertex_positions_T, face_vertex_uv_positions, corner_to_corner_uv_positions,
//      reverse_edge_orientations, is_cone, is_cone_adjacent, normal, optimization_params, local_energy_assemble,
//      local_derivatives_assemble, local_hessian_assemble);
//
//  std::cout << "energy mine: " << local_energy_assemble << std::endl
//            << std::endl;
//
//  std::cout << "derivative mine: " << std::endl
//            << local_derivatives_assemble << std::endl
//            << std::endl;
//
//  std::cout << "hessian mine: " << std::endl
//            << local_hessian_assemble << std::endl
//            << std::endl;
//
//  std::cout << "energy diff: "
//            << abs(local_energy_autodiff - local_energy_assemble) << std::endl;
//
//  std::cout << "derivatives norm diff: "
//            << (local_derivatives_autodiff - local_derivatives_assemble).norm()
//            << std::endl;
//
//  std::cout << "hessian norm diff: "
//            << (local_hessian_autodiff - local_hessian_assemble).norm()
//            << std::endl;
//}