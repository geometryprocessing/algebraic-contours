#include <catch2/catch_test_macros.hpp>

#include "common.h"
#include "generate_shapes.h"
#include "optimize_spline_surface.h"
#include <fenv.h>


//void build_corner_variables(
//  const MatrixXr &P,
//  const MatrixXr &G,
//  DynamicVector3v &pi,
//  DynamicVector3v &pj,
//  DynamicVector3v &pk,
//  DynamicMatrix2x3v &Gi
//) {
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(0), pi, 0, 15);
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(1), pj, 3, 15);
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(2), pk, 6, 15);
//  build_independent_variable_matrix<DynamicVar, DynamicMatrix2x3v>(G, Gi, 9, 15);
//}
//
//
//void build_triangle_variables(
//  const MatrixXr &P,
//  const MatrixXr &G,
//  DynamicVector3v &pi,
//  DynamicVector3v &pj,
//  DynamicVector3v &pk,
//  DynamicMatrix2x3v &Gi,
//  DynamicMatrix2x3v &Gj,
//  DynamicMatrix2x3v &Gk
//) {
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(0), pi, 0, 27);
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(1), pj, 3, 27);
//  build_independent_variable_vector<DynamicVar, DynamicVector3v>(P.row(2), pk, 6, 27);
//  build_independent_variable_matrix<DynamicVar, DynamicMatrix2x3v>(G.block(0, 0, 2, 3), Gi, 9, 27);
//  build_independent_variable_matrix<DynamicVar, DynamicMatrix2x3v>(G.block(2, 0, 2, 3), Gj, 15, 27);
//  build_independent_variable_matrix<DynamicVar, DynamicMatrix2x3v>(G.block(4, 0, 2, 3), Gk, 21, 27);
//}

// 2u^2 + v^2
void build_quadratic_biharmonic_square(
  MatrixXr &V,
  std::vector<Matrix2x3r> &gradients,
  std::vector<std::vector<MatrixXr>> &layout
) {
  V.resize(5, 3);
  V <<
    0, 0, 0,
    2, 0, 8,
    2, 2, 12,
    0, 2, 4,
    1, 1, 3;

  Eigen::MatrixXi F(4, 3);
  F <<
    0, 1, 4,
    1, 2, 4,
    2, 3, 4,
    3, 0, 4;

  layout.resize(4);
  for (size_t face_index = 0; face_index < 4; ++face_index)
  {
    layout[face_index].resize(3);
    for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
    {
      size_t i = F(face_index, vertex_index);
      size_t j = F(face_index, (vertex_index + 1) % 3);
      size_t k = F(face_index, (vertex_index + 2) % 3);

      VectorXr vi = V.row(i);
      VectorXr vj = V.row(j);
      VectorXr vk = V.row(k);

      layout[face_index][vertex_index].resize(2, 2);
      layout[face_index][vertex_index](0,0) = vj(0) - vi(0);
      layout[face_index][vertex_index](0,1) = vj(1) - vi(1);
      layout[face_index][vertex_index](1,0) = vk(0) - vi(0);
      layout[face_index][vertex_index](1,1) = vk(1) - vi(1);
    }
  }

  gradients.resize(5);
  gradients[0] <<
    1, 0, 0,
    0, 1, 0;
  gradients[1] <<
    1, 0, 8,
    0, 1, 0;
  gradients[2] <<
    1, 0, 8,
    0, 1, 4;
  gradients[3] <<
    1, 0, 0,
    0, 1, 4;
  gradients[4] <<
    1, 0, 4,
    0, 1, 2;
}

SpatialVector quadratic_torus_point(
  double u,
  double v,
  double R1=5.0,
  double R2=1.0
) {
  SpatialVector point(3);
  point <<
    R1 * (1 - (u * u / 2)) + R2 * (1 - (v * v / 2) - (u * u / 2)),
    (R1 + R2) * u,
    R2 * v;
  return point;
}

Matrix2x3r quadratic_torus_gradient(
  double u,
  double v,
  double R1=5.0,
  double R2=1.0
) {
  Matrix2x3r gradient(2, 3);
  gradient <<
    -R1 * u - R2 * u, (R1 + R2),  0,
             -R2 * v,         0, R2;

  return gradient;
}


void build_quadratic_torus_square(
  MatrixXr &V,
  std::vector<Matrix2x3r> &gradients,
  std::vector<std::vector<MatrixXr>> &layout
) {
  double pi = M_PI;
  V.resize(5, 3);
  V.row(0) = quadratic_torus_point(-pi/4.0, -pi/2.0);
  V.row(1) = quadratic_torus_point( pi/4.0, -pi/2.0);
  V.row(2) = quadratic_torus_point( pi/4.0,  pi/2.0);
  V.row(3) = quadratic_torus_point(-pi/4.0,  pi/2.0);
  V.row(4) = quadratic_torus_point(0.0, 0.0);

  Eigen::MatrixXi F(4, 3);
  F <<
    0, 1, 4,
    1, 2, 4,
    2, 3, 4,
    3, 0, 4;

  MatrixXr U(5, 2);
  U <<
    -pi/4.0, -pi/2.0,
     pi/4.0, -pi/2.0,
     pi/4.0,  pi/2.0,
    -pi/4.0,  pi/2.0,
        0.0,     0.0;

  layout.resize(4);
  for (size_t face_index = 0; face_index < 4; ++face_index)
  {
    layout[face_index].resize(3);
    for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
    {
      size_t i = F(face_index, vertex_index);
      size_t j = F(face_index, (vertex_index + 1) % 3);
      size_t k = F(face_index, (vertex_index + 2) % 3);

      PlanarPoint ui = U.row(i);
      PlanarPoint uj = U.row(j);
      PlanarPoint uk = U.row(k);

      layout[face_index][vertex_index].resize(2, 2);
      layout[face_index][vertex_index](0,0) = uj(0) - ui(0);
      layout[face_index][vertex_index](0,1) = uj(1) - ui(1);
      layout[face_index][vertex_index](1,0) = uk(0) - ui(0);
      layout[face_index][vertex_index](1,1) = uk(1) - ui(1);
    }
  }

  gradients.resize(5);
  gradients[0].resize(3, 2);
  gradients[0] = quadratic_torus_gradient(-pi/4.0, -pi/2.0);
  gradients[1] = quadratic_torus_gradient( pi/4.0, -pi/2.0);
  gradients[2] = quadratic_torus_gradient( pi/4.0,  pi/2.0);
  gradients[3] = quadratic_torus_gradient(-pi/4.0,  pi/2.0);
  gradients[4] = quadratic_torus_gradient(0.0, 0.0);
}

// u^3 + v^3
void build_cubic_biharmonic_square(
  MatrixXr &V,
  std::vector<Matrix2x3r> &gradients,
  std::vector<std::vector<MatrixXr>> &layout
) {
  V.resize(5, 3);
  V <<
    0, 0, 0,
    2, 0, 8,
    2, 2, 16,
    0, 2, 8,
    1, 1, 2;

  Eigen::VectorXi F(4, 3);
  F <<
    0, 1, 4,
    1, 2, 4,
    2, 3, 4,
    3, 0, 4;

  layout.resize(4);
  for (size_t face_index = 0; face_index < 4; ++face_index)
  {
    layout[face_index].resize(3);
    for (size_t i = 0; i < 3; ++i)
    {
      size_t j = (i + 1) % 3;
      size_t k = (i + 2) % 3;

      SpatialVector vi = V.row(i);
      SpatialVector vj = V.row(j);
      SpatialVector vk = V.row(k);

      layout[face_index][i].resize(2, 2);
      layout[face_index][i](0,0) = vj(0) - vi(0);
      layout[face_index][i](0,1) = vj(1) - vi(1);
      layout[face_index][i](1,0) = vk(0) - vi(0);
      layout[face_index][i](1,1) = vk(1) - vi(1);
    }
  }

  gradients.resize(5);
  gradients[0] <<
    1, 0, 0,
    0, 1, 0;
  gradients[1] <<
    1, 0, 12,
    0, 1, 0;
  gradients[2] <<
    1, 0, 12,
    0, 1, 12;
  gradients[3] <<
    1, 0, 0,
    0, 1, 12;
  gradients[4] <<
    1, 0, 3,
    0, 1, 3;
}


TEST_CASE ( "Single triangle affine manifold energy can be computed", "[optimize_spline_surface]")
{
  MatrixXr initial_V(3, 3);
  std::vector<Matrix2x3r> initial_gradients(3);
  
  MatrixXr optimized_V;
  std::vector<Matrix2x3r> optimized_gradients;
  std::vector<size_t> variable_vertices;
  arange(initial_V.rows(), variable_vertices);
  OptimizationParameters optimization_params;
  //optimization_params.least_squares_factor = 1.0;
  optimization_params.position_difference_factor = 1.0;

  // Initialize equilateral triangle affine manifold
  Eigen::MatrixXi F(1, 3);
  std::vector<std::vector<double>> l(1);
  F << 0, 1, 2;
  //l[0] = { sqrt(2.0), 1.0, 1.0 };
  //AffineManifold affine_manifold(F, l);

  SECTION ( "Zero" )
  {
//    initial_V << 
//      0, 0, 0,
//      0, 0, 0,
//      0, 0, 0;
//    initial_gradients[0] <<
//      0, 0, 0,
//      0, 0, 0;
//    initial_gradients[1] <<
//      0, 0, 0,
//      0, 0, 0;
//    initial_gradients[2] <<
//      0, 0, 0,
//      0, 0, 0;
//
//    double energy = optimize_six_split_spline_surface(
//      initial_V,
//      initial_gradients,
//      affine_manifold,
//      variable_vertices,
//      optimization_params,
//      optimized_V,
//      optimized_gradients
//    );
//
//    spdlog::trace("Zero edge energy: {}", energy);
//    REQUIRE ( float_equal(energy, 0.0) );
//    REQUIRE ( matrix_equal(initial_V, optimized_V) );
//    REQUIRE ( matrix_equal(initial_gradients[0], optimized_gradients[0]) );
//    REQUIRE ( matrix_equal(initial_gradients[1], optimized_gradients[1]) );
//    REQUIRE ( matrix_equal(initial_gradients[2], optimized_gradients[2]) );
  }

  SECTION ( "Zero with random gradients" )
  {
//    initial_V << 
//      0, 0, 0,
//      0, 0, 0,
//      0, 0, 0;
//    initial_gradients[0] <<
//      0, 1, 3,
//      9, 2, 1;
//    initial_gradients[1] <<
//      1, 1, 6,
//      1, 2, 1;
//    initial_gradients[2] <<
//      1, 2, 1,
//      3, 1, 3;
//    Matrix2x3r expected_gradient(2, 3);
//    expected_gradient <<
//      0, 0, 0,
//      0, 0, 0;
//
//    double energy = optimize_six_split_spline_surface(
//      initial_V,
//      initial_gradients,
//      affine_manifold,
//      variable_vertices,
//      optimization_params,
//      optimized_V,
//      optimized_gradients
//    );
//
//    spdlog::trace("Zero edge energy: {}", energy);
//    spdlog::trace("Optimized gradients:\n{}", formatted_vector(optimized_gradients));
//    REQUIRE ( float_equal(energy, 0.0) );
//    REQUIRE ( matrix_equal(initial_V, optimized_V) );
//    REQUIRE ( matrix_equal(expected_gradient, optimized_gradients[0]) );
//    REQUIRE ( matrix_equal(expected_gradient, optimized_gradients[1]) );
//    REQUIRE ( matrix_equal(expected_gradient, optimized_gradients[2]) );
  }

  SECTION ( "Planar triangle" )
  {
//    initial_V << 
//      0, 0, 0,
//      1, 0, 0,
//      0, 1, 0;
//    initial_gradients[0] <<
//      0, 1, 3,
//      9, 2, 1;
//    initial_gradients[1] <<
//      1, 1, 6,
//      1, 2, 1;
//    initial_gradients[2] <<
//      1, 2, 1,
//      3, 1, 3;
//    Matrix2x3r expected_gradient(2, 3);
//    expected_gradient <<
//      1, 0, 0,
//      0, 1, 0;
//
//    double energy = optimize_six_split_spline_surface(
//      initial_V,
//      initial_gradients,
//      affine_manifold,
//      variable_vertices,
//      optimization_params,
//      optimized_V,
//      optimized_gradients
//    );
//
//    spdlog::trace("Zero edge energy: {}", energy);
//    spdlog::trace("Optimized gradients:\n{}", formatted_vector(optimized_gradients));
//    REQUIRE ( float_equal(energy, 0.0) );
//    REQUIRE ( matrix_equal(initial_V, optimized_V) );
//    REQUIRE ( matrix_equal(expected_gradient, optimized_gradients[0]) );
  }
}


struct QuadraticTestParameters
{
  size_t resolution = 3;
  double u_curvature = 0;
  double v_curvature = 0;
  double uv_curvature = 0;
  double dx = 0;
  double dy = 0;
  double dz = 0;
  bool tile_perturbation = false;
};

void generate_quadratic_test(
  const QuadraticTestParameters &quadratic_test_params,
  MatrixXr &V,
  std::vector<Matrix2x3r> &gradients,
  std::vector<int> &variable_vertices,
  AffineManifold &affine_manifold
) {
  size_t resolution = quadratic_test_params.resolution;
  double u_curvature = quadratic_test_params.u_curvature;
  double v_curvature = quadratic_test_params.v_curvature;
  double uv_curvature = quadratic_test_params.uv_curvature;
  double dx = quadratic_test_params.dx;
  double dy = quadratic_test_params.dy;
  double dz = quadratic_test_params.dz;
  bool tile_perturbation = quadratic_test_params.tile_perturbation;

  // Build uv layout 
  std::vector< std::vector<PlanarPoint> > layout_point_grid;
  generate_global_layout_grid(resolution, layout_point_grid);

  size_t perturbed_index = (resolution / 2);
  size_t perturbed_u_index = perturbed_index;
  size_t perturbed_v_index = perturbed_index;
  if (tile_perturbation)
  {
    for (size_t i = 1; i < resolution; i += 2)
    {
      for (size_t j = 1; j < resolution; j += 2)
      {
      spdlog::info("Perturbing index {}, {}", i, j);
      layout_point_grid[i][j][0] += dx;
      layout_point_grid[i][j][1] += dy;
      }
    }
  } else {
    spdlog::info("Perturbing index {}, {}", perturbed_u_index, perturbed_v_index);
    layout_point_grid[perturbed_u_index][perturbed_v_index][0] += dx;
    layout_point_grid[perturbed_u_index][perturbed_v_index][1] += dy;
  }

  // Generate quadratic function grid points
  std::vector< std::vector<SpatialVector> > point_grid;
  generate_quadratic_grid(
    layout_point_grid,
    u_curvature,
    v_curvature,
    uv_curvature,
    point_grid
  );

  // Build mesh from the grid points
  Eigen::MatrixXi F;
  std::vector<std::vector<double>> l;
  bool closed_surface = false;
  generate_mesh_from_grid(
    point_grid,
    layout_point_grid,
    V,
    F,
    l,
    closed_surface
  );
  // FIXME V(perturbed_u_index, perturbed_v_index) -= quadratic_test_params.dz;
  spdlog::debug("Initial spatial points:\n{}", V);

  // Build an affine manifold with global layout
  MatrixXr uv;
  Eigen::MatrixXi layout_F;
  std::vector<std::vector<double>> layout_l;
  generate_mesh_from_grid(layout_point_grid, uv, layout_F, layout_l, closed_surface);
  affine_manifold = ParametricAffineManifold(F, uv);
  spdlog::debug("Global UV points:\n{}", uv);

  // Generate gradients for the quadratic at the grid points
  std::vector< std::vector<Matrix2x3r> > gradient_grid;
  generate_quadratic_gradients_grid(
    layout_point_grid,
    u_curvature,
    v_curvature,
    uv_curvature,
    gradient_grid
  );

  // Generate gradients and fixed degrees of freedom from the grid
  generate_gradients_from_grid(gradient_grid, gradients, variable_vertices);
}

//double compute_quadratic_test_energy(
//  const QuadraticTestParameters &quadratic_test_params,
//  bool optimize = false
//) {
//  MatrixXr V;
//  std::vector<Matrix2x3r> gradients;
//  std::vector<int> variable_vertices;
//  AffineManifold affine_manifold;
//  generate_quadratic_test(
//    quadratic_test_params,
//    V,
//    gradients,
//    variable_vertices,
//    affine_manifold
//  );
//
//  // Set the optimization parameters
//  OptimizationParameters optimization_params;
//  optimization_params.parametrized_quadratic_surface_mapping_factor = 1.0;
//  optimization_params.position_difference_factor = 1.0;
//
//  // Optimize the analytic initial values
//  MatrixXr optimized_V;
//  std::vector<Matrix2x3r> optimized_gradients;
//  double energy = optimize_six_split_spline_surface(
//    V,
//    gradients,
//    affine_manifold,
//    variable_vertices,
//    optimization_params,
//    optimized_V,
//    optimized_gradients
//  );
//
//  return energy;
//}

bool quadratic_refinement_test(
  const QuadraticTestParameters &quadratic_test_params
) {
  spdlog::set_level(spdlog::level::info);
  spdlog::info("Testing quadratic refinement");
  bool success = true;

  // Build test parameters for the quadratic refinement
  QuadraticTestParameters quadratic_refinement_test_params;
  quadratic_refinement_test_params = quadratic_test_params;

  // Scale proportionally to the resolution
  int resolution = quadratic_refinement_test_params.resolution;
  double k2 = (resolution - 1);
  quadratic_refinement_test_params.u_curvature /= k2;
  quadratic_refinement_test_params.v_curvature /= k2;
  quadratic_refinement_test_params.uv_curvature /= k2;

  // Compute the energy
  //double energy = compute_quadratic_test_energy(quadratic_refinement_test_params, true);
  ///spdlog::info("{}x{} energy: {}", resolution, resolution, energy);

  return success;
}



//void test()
//{
//  std::vector<std::vector<int>> face_to_patch_indices;
//  std::vector<int> patch_to_face_indices;
//  QuadraticSplineSurface initial_spline_surface(
//    initial_V,
//    initial_gradients,
//    affine_manifold,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//  //initial_spline_surface.view();
//  double pert = 1e-6;
//  MatrixXr V_pert(10, 3);
//  V_pert.setZero(10, 3);
//  V_pert(0, 0) = pert;
//  V_pert(1, 1) = pert;
//  V_pert(2, 2) = pert;
//  MatrixXr gradients_pert(20, 3);
//  gradients_pert.setZero(20, 3);
//  gradients_pert(6, 0) = pert;
//  gradients_pert(8, 1) = pert;
//  gradients_pert(10, 2) = pert;
//  gradients_pert(13, 0) = pert;
//  gradients_pert(15, 1) = pert;
//  gradients_pert(17, 2) = pert;
//  MatrixXr analytic_hessian;
//  VectorXr analytic_derivatives;
//  double analytic_energy;
//  for (size_t i = 9; i < 10; ++i)
//  {
//    spdlog::info("Perturbation {} with\n{}\nand\n{}", i, V_pert.row(i), gradients_pert.block(2 * i, 0, 2, 3));
//    MatrixXr initial_V_pert = initial_V;
//    initial_V_pert.row(4) += V_pert.row(i);
//    std::vector<MatrixXr> initial_gradients_pert = initial_gradients;
//    initial_gradients_pert[4] += gradients_pert.block(2 * i, 0, 2, 3) ;
//    initialize_six_split_variables(
//      initial_V_pert,
//      initial_gradients_pert,
//      variable_vertices,
//      positions,
//      gradients,
//      false
//    );
//    generate_six_split_energy_quadratic(
//      positions,
//      gradients,
//      initial_V,
//      affine_manifold,
//      triangle_energy,
//      analytic_hessian,
//      analytic_derivatives,
//      analytic_energy
//    );
//    assert( matrix_equal(analytic_hessian, hessian) );
//  }
//  spdlog::info("Initial variable values: {}", initial_variable_values);
//  spdlog::info("Analytic variable values: {}", variable_values);
//
//  optimized_V = initial_V;
//  optimized_gradients = initial_gradients;
//  spdlog::info("z coord derivative {}", analytic_derivatives[2]);
//  optimized_V(4, 2) -= 10 * analytic_derivatives[2];
//  initialize_six_split_variables(
//    optimized_V,
//    optimized_gradients,
//    variable_vertices,
//    positions,
//    gradients,
//    false
//  );
//  generate_six_split_energy_quadratic(
//    positions,
//    gradients,
//    initial_V,
//    affine_manifold,
//    triangle_energy,
//    hessian,
//    derivatives,
//    energy
//  );
//  QuadraticSplineSurface first_optimized_spline_surface(
//    optimized_V,
//    optimized_gradients,
//    affine_manifold,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//  first_optimized_spline_surface.view();
//
//  if ( !float_equal(residual.norm(), 0.0) ) return false;
//  return true;
//
//
//  spdlog::debug("Optimized energy: {}", energy);
//  spdlog::debug("Optimized positions:\n{}", optimized_V);
//  spdlog::debug("Optimized gradients:\n{}", formatted_vector(optimized_gradients));
//  QuadraticSplineSurface second_optimized_spline_surface(
//    optimized_V,
//    optimized_gradients,
//    affine_manifold,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//  first_optimized_spline_surface.view();
//
//  // Check for errors 
//  bool success = true;
//  for (size_t i = 0; i < optimized_gradients.size(); ++i)
//  {
//    if( !matrix_equal( analytic_gradients[i], optimized_gradients[i] ) )
//    {
//      spdlog::error(
//        "Optimized gradient\n{}\nnot equal to expected value\n{}",
//        optimized_gradients[i],
//        analytic_gradients[i]
//      );
//      spdlog::error(
//        "Corresponding optimized vertex is \n{}\nwith expected value{}",
//        optimized_V.row(i),
//        initial_V.row(i)
//      );
//      success = false;
//    }
//  }
//  // Get initial residual
//  generate_six_split_variable_value_vector(
//    positions,
//    gradients,
//    variable_vertices,
//    variable_values
//  );
//  compute_quadratic_residual(
//    variable_values,
//    derivatives,
//    hessian,
//    residual
//  );
//  spdlog::info("Initial residual: {}", residual);
//  if (!success) return false;
//  return true;
//
//  spdlog::info("Testing with one gradient set to 0");
//  initial_gradients[resolution + 1] <<
//    0, 0, 0,
//    0, 0, 0;
//  energy = optimize_six_split_spline_surface(
//    initial_V,
//    initial_gradients,
//    affine_manifold,
//    variable_vertices,
//    optimization_params,
//    optimized_V,
//    optimized_gradients
//  );
//
//  QuadraticSplineSurface optimized_spline_surface(
//    optimized_V,
//    optimized_gradients,
//    affine_manifold,
//    face_to_patch_indices,
//    patch_to_face_indices
//  );
//  optimized_spline_surface.view();
//
//  spdlog::debug("Perturbed start energy: {}", energy);
//  spdlog::debug("Optimized positions:\n{}", optimized_V);
//  spdlog::debug("Optimized gradients:\n{}", formatted_vector(optimized_gradients));
//
//  for (size_t i = 0; i < optimized_gradients.size(); ++i)
//  {
//    if( !matrix_equal( analytic_gradients[i], optimized_gradients[i] ) ) return false;
//  }
//
//  return true;
//}

TEST_CASE ( "Quadratic affine manifold energy can be computed", "[optimize_spline_surface]")
{
  MatrixXr initial_V;
  std::vector<Matrix2x3r> initial_gradients;
  MatrixXr optimized_V;
  std::vector<Matrix2x3r> optimized_gradients;

  SECTION ( "Small planar" )
  {
    QuadraticTestParameters quadratic_test_params;
    //bool success = quadratic_residual_test(quadratic_test_params);
    //REQUIRE ( success );
  }

//  SECTION ( "Large planar" )
//  {
//    size_t resolution = 4;
//    double u_curvature = 0;
//    double v_curvature = 0;
//    bool success = quadratic_surface_test(
//      resolution,
//      u_curvature,
//      v_curvature,
//      initial_V,
//      initial_gradients,
//      optimized_V,
//      optimized_gradients
//    );
//
//    REQUIRE ( success );
//  }

  SECTION ( "Small quadratic" )
  {
    QuadraticTestParameters quadratic_test_params;
    quadratic_test_params.u_curvature = 2;
    quadratic_test_params.dx = 0.1;
    quadratic_test_params.tile_perturbation = true;
    bool success = true;
    bool test_success;
    //test_success = quadratic_residual_test(quadratic_test_params);
    success = success && test_success;

    test_success = quadratic_refinement_test(quadratic_test_params);
    success = success && test_success;

    //quadratic_test_params.resolution = 5;
    //test_success = quadratic_refinement_test(quadratic_test_params);
    //success = success && test_success;
    
    //quadratic_test_params.resolution = 7;
    //test_success = quadratic_refinement_test(quadratic_test_params);
    //success = success && test_success;

    //quadratic_test_params.resolution = 9;
    //test_success = quadratic_refinement_test(quadratic_test_params);
    //success = success && test_success;
    //REQUIRE ( success );
  }

//  SECTION ( "Large quadratic" )
//  {
//    size_t resolution = 4;
//    double u_curvature = 0;
//    double v_curvature = 1;
//    bool success = quadratic_surface_test(
//      resolution,
//      u_curvature,
//      v_curvature,
//      initial_V,
//      initial_gradients,
//      optimized_V,
//      optimized_gradients
//    );
//
//    REQUIRE ( success );
//
//  }
}

/*
TEST_CASE ( "Triangle quadratic surface mapping energy can be computed", "[optimize_spline_surface]")
{
  MatrixXr P(3, 3);
  MatrixXr G(6, 3);
  DynamicVectorXv pi, pj, pk;
  DynamicMatrixXv Gi, Gj, Gk;
  MatrixXr real_Gi, real_Gj, real_Gk;
  VectorXr gradient;
  MatrixXr hessian;

  SECTION ( "Planar triangle" )
  {
    // Setup planar triangle
    P << 
      0, 0, 0,
      1, 0, 0,
      0, 1, 0;
    G <<
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0;
    MatrixXr Ui(2,2);
    MatrixXr Uj(2,2);
    MatrixXr Uk(2,2);
    Ui <<
      1, 0,
      0, 1;
    Uj <<
      -1, 1,
      -1, 0;
    Uk <<
      0, -1,
      1, -1;

    // Build independent variables for the triangle
    build_triangle_variables(P, G, pi, pj, pk, Gi, Gj, Gk);

    // Compute initial differentiable energy 
    TriangleEnergy triangle_energy;
    DynamicVar energy = triangle_energy.triangle_quadratic_surface_mapping_energy(pi, pj, pk, Gi, Gj, Gk, Ui, Uj, Uk);
    double energy_value = compute_variable_value(energy);
    compute_variable_gradient(energy, gradient);
    compute_variable_hessian(energy, hessian);
    spdlog::trace("Energy: {}", energy_value);
    spdlog::trace("Gradient:\n{}", gradient);

    // Get Newton direction for fixed positions
    VectorXr newton_direction;
    compute_fixed_position_newton_direction(
      gradient,
      hessian,
      9,
      18,
      newton_direction
    );

    // Update variables with Newton direction
    spdlog::trace("Newton direction:\n{}", newton_direction);
    update_independent_variable_matrix(newton_direction, Gi, 9);
    update_independent_variable_matrix(newton_direction, Gj, 15);
    update_independent_variable_matrix(newton_direction, Gk, 21);

    // Check new gradient values
    extract_variable_matrix_values(Gi, real_Gi);
    extract_variable_matrix_values(Gj, real_Gj);
    extract_variable_matrix_values(Gk, real_Gk);
    spdlog::trace("New gradients:\n{}\n{}\n{}", real_Gi, real_Gj, real_Gk);

    // Compute new energy with new gradients (should be 0)
    energy = triangle_energy.triangle_quadratic_surface_mapping_energy(pi, pj, pk, Gi, Gj, Gk, Ui, Uj, Uk);
    energy_value = compute_variable_value(energy);
    compute_variable_gradient(energy, gradient);
    spdlog::trace("New energy: {}", energy_value);
    spdlog::trace("New gradient:\n{}", gradient);
    REQUIRE ( float_equal(compute_variable_value(energy), 0.0) );
  }


  SECTION ( "Quadratic biharmonic" )
  {
    // Get initial setup
    MatrixXr initial_V;
    std::vector<MatrixXr> initial_gradients;
    std::vector<std::vector<MatrixXr>> layout;
    //build_quadratic_biharmonic_square(initial_V, initial_gradients, layout);
    build_quadratic_torus_square(initial_V, initial_gradients, layout);
    
    // Build independent variables
    DynamicMatrixXv V;
    std::vector<DynamicVectorXv> positions;
    std::vector<DynamicMatrixXv> gradients;
    initialize_six_split_variables(
      initial_V,
      initial_gradients,
      positions,
      gradients
    );

    // Perturb gradient and optimize
    gradients[4](0, 0) += DynamicVar(1);
    gradients[4](0, 1) += DynamicVar(2);
    gradients[4](0, 2) += DynamicVar(3);
    gradients[4](1, 0) += DynamicVar(4);
    gradients[4](1, 1) += DynamicVar(5);
    gradients[4](1, 2) += DynamicVar(6);

    // Sum quadratic surface mapping energies for the three triangles
    DynamicVar energy(0);
    TriangleEnergy triangle_energy;
    energy += triangle_energy.triangle_quadratic_surface_mapping_energy(
      V.row(0),
      V.row(1),
      V.row(4),
      gradients[0],
      gradients[1],
      gradients[4],
      layout[0][0],
      layout[0][1],
      layout[0][2]
    );
    energy += triangle_energy.triangle_quadratic_surface_mapping_energy(
      V.row(1),
      V.row(2),
      V.row(4),
      gradients[1],
      gradients[2],
      gradients[4],
      layout[1][0],
      layout[1][1],
      layout[1][2]
    );
    energy += triangle_energy.triangle_quadratic_surface_mapping_energy(
      V.row(2),
      V.row(3),
      V.row(4),
      gradients[2],
      gradients[3],
      gradients[4],
      layout[2][0],
      layout[2][1],
      layout[2][2]
    );
    energy += triangle_energy.triangle_quadratic_surface_mapping_energy(
      V.row(3),
      V.row(0),
      V.row(4),
      gradients[3],
      gradients[0],
      gradients[4],
      layout[3][0],
      layout[3][1],
      layout[3][2]
    );

    // Compute energy with gradient and hessian
    double energy_value = compute_variable_value(energy);
    compute_variable_gradient(energy, gradient);
    compute_variable_hessian(energy, hessian);
    spdlog::trace("Energy: {}", energy_value);
    spdlog::trace("Gradient:\n{}", gradient);

    // Compute newton direction for the center gradient
    VectorXr center_gradient = gradient.segment(39, 6);
    MatrixXr center_hessian = hessian.block(39, 39, 6, 6);
    VectorXr center_newton_direction = center_hessian.ldlt().solve(center_gradient);
    spdlog::trace("Newton direction:\n{}", center_newton_direction);

    // Newton direction is the perturbation
    REQUIRE ( float_equal(center_newton_direction[0], 1.0) );
    REQUIRE ( float_equal(center_newton_direction[1], 2.0) );
    REQUIRE ( float_equal(center_newton_direction[2], 3.0) );
    REQUIRE ( float_equal(center_newton_direction[3], 4.0) );
    REQUIRE ( float_equal(center_newton_direction[4], 5.0) );
    REQUIRE ( float_equal(center_newton_direction[5], 6.0) );
  }
}
*/

// FIXME Redundant triangle energy test
//TEST_CASE ( "Corner least squares energy can be computed", "[optimize_spline_surface]")
//{
//  MatrixXr P(3, 3);
//  MatrixXr G(2, 3);
//  DynamicVectorXv pi, pj, pk;
//  VectorXr uj(2);
//  VectorXr uk(2);
//  DynamicMatrixXv Gi;
//
//  SECTION ( "Zero" )
//  {
//    P << 
//      0, 0, 0,
//      0, 0, 0,
//      0, 0, 0;
//    G <<
//      0, 0, 0,
//      0, 0, 0;
//    uj << 1, 0;
//
//    build_corner_variables(P, G, pi, pj, pk, Gi);
//    TriangleEnergy triangle_energy;
//    DynamicVar energy = triangle_energy.triangle_corner_to_corner_least_squares_energy(pi, pj, Gi, uj);
//    spdlog::trace("Zero edge energy: {}", compute_variable_value(energy));
//    //spdlog::trace("Zero edge gradient: {}", energy.derivatives());
//    REQUIRE ( float_equal(compute_variable_value(energy), 0.0) );
//  }
//
//  SECTION ( "Constant" )
//  {
//    P << 
//      0, 0, 0,
//      1, 0, 0,
//      0, 1, 0;
//    G <<
//      1, 0, 0,
//      0, 1, 0;
//    uj << 1, 0;
//    uk << 0, 1; 
//
//    build_corner_variables(P, G, pi, pj, pk, Gi);
//    TriangleEnergy triangle_energy;
//    DynamicVar energy = triangle_energy.triangle_corner_to_corner_least_squares_energy(pi, pj, Gi, uj);
//    spdlog::trace("Constant edge energy: {}", compute_variable_value(energy));
//    //spdlog::trace("Constant edge gradient: {}", energy.derivatives());
//    REQUIRE ( float_equal(compute_variable_value(energy), 0.0) );
//  }
//}