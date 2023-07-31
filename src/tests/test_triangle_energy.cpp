#include <catch2/catch_test_macros.hpp>
#include "common.h"

#include "triangle_energy.h"
#include "generate_shapes.h"

typedef DScalar2<double, DynamicGradient, DynamicHessian> DynamicVar;
typedef DScalar2<double, SixSplitGradient, SixSplitHessian> SixSplitVar;
typedef DScalar2<double, TwelveSplitGradient, TwelveSplitHessian>
  TwelveSplitVar;

typedef Eigen::Matrix<DynamicVar, 1, 3> DynamicVector3v;
typedef Eigen::Matrix<DynamicVar, 2, 3> DynamicMatrix2x3v;
typedef Eigen::Matrix<DynamicVar, Eigen::Dynamic, 1> DynamicVectorXv;
typedef Eigen::Matrix<DynamicVar, 1, Eigen::Dynamic> DynamicOneFormXv;
typedef Eigen::Matrix<DynamicVar, Eigen::Dynamic, Eigen::Dynamic>
  DynamicMatrixXv;

TEST_CASE ( "The position difference triangle energy can be computed" )
{
  TriangleEnergy<DynamicVar>  triangle_energy;

  SECTION ( "Vertex zero energy" )
  {
    DynamicVector3v p(3);
    p <<
      DynamicVar(1), DynamicVar(2), DynamicVar(3);
    SpatialVector p_0(3);
    p_0 << 1, 2, 3;
    DynamicVar energy = triangle_energy.vertex_position_difference_energy(p, p_0);
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Vertex energy" )
  {
    DynamicVector3v p(3);
    p <<
      DynamicVar(1), DynamicVar(2), DynamicVar(3);
    SpatialVector p_0(3);
    p_0 << 0, 3, 5;
    DynamicVar energy = triangle_energy.vertex_position_difference_energy(p, p_0);
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 6.0) );
  }

  SECTION ( "Triangle zero energy" )
  {
    DynamicVector3v pi(3);
    DynamicVector3v pj(3);
    DynamicVector3v pk(3);
    pi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(0);
    pj <<
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    pk <<
      DynamicVar(0), DynamicVar(0), DynamicVar(1);
    SpatialVector target_pi(3);
    SpatialVector target_pj(3);
    SpatialVector target_pk(3);
    target_pi <<
      1, 0, 0;
    target_pj <<
      0, 1, 0;
    target_pk <<
      0, 0, 1;
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<SpatialVector, 3> target_face_vertex_positions = { target_pi, target_pj, target_pk };
    DynamicVar energy = triangle_energy.triangle_position_difference_energy(
      face_vertex_positions,
      target_face_vertex_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Triangle energy" )
  {
    DynamicVector3v pi(3);
    DynamicVector3v pj(3);
    DynamicVector3v pk(3);
    pi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(0);
    pj <<
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    pk <<
      DynamicVar(0), DynamicVar(0), DynamicVar(1);
    SpatialVector target_pi(3);
    SpatialVector target_pj(3);
    SpatialVector target_pk(3);
    target_pi <<
      1, 0, 0;
    target_pj <<
      0, 2, 0;
    target_pk <<
      2, 3, 1;
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<SpatialVector, 3> target_face_vertex_positions = { target_pi, target_pj, target_pk };
    DynamicVar energy = triangle_energy.triangle_position_difference_energy(
      face_vertex_positions,
      target_face_vertex_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 14.0) );
  }
}

TEST_CASE ( "The least squares triangle energy can be computed" )
{
  TriangleEnergy<DynamicVar> triangle_energy;

  // Positions
  DynamicVector3v pi(3);
  DynamicVector3v pj(3);
  DynamicVector3v pk(3);
  pi <<
    DynamicVar(0), DynamicVar(0), DynamicVar(0);
  pj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0);
  pk <<
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };

  // Gradients
  DynamicMatrix2x3v Gi(2, 3);
  Gi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gj(2, 3);
  Gj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gk(2, 3);
  Gk <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };

  // UV
  Matrix2x2r UVi(2, 2);
  UVi <<
    1, 0,
    0, 1;
  Matrix2x2r UVj(2, 2);
  UVj <<
    -1, 1,
    -1, 0;
  Matrix2x2r UVk(2, 2);
  UVk <<
    0, -1,
    1, -1;
  std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };

  SECTION ( "Edge zero energy" )
  {
    PlanarPoint uvj(2);
    uvj <<
      1, 0;
    DynamicVar energy = triangle_energy.triangle_corner_to_corner_least_squares_energy(pi, pj, Gi, uvj);
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Edge energy" )
  {
    PlanarPoint uvj(2);
    uvj <<
      3, 1;
    DynamicVar energy = triangle_energy.triangle_corner_to_corner_least_squares_energy(pi, pj, Gi, uvj);
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 5.0) );
  }

  SECTION ( "Triangle zero energy" )
  {
    DynamicVar energy = triangle_energy.triangle_vertex_least_squares_energy(
       face_vertex_positions,
       face_vertex_gradients,
       corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Triangle energy" )
  {
    // Perturb first gradients for new edge derivatives [2, 1, 0] and [0, 1, 0]
    // Change in energy: +2
    face_vertex_gradients[0](0, 0) = 2;
    face_vertex_gradients[0](0, 1) = 1;

    // Perturb second UV coordinates for new edge derivatives [0, 1, 0] and [-1, 1, 0]
    // Target edge derivatives are [-1, 1, 0] and [-1, 0, 0] (before pk perturbation)
    // Change in energy: +2
    corner_to_corner_uv_positions[1] <<
       0, 1,
      -1, 1;
    
    // Perturb third positions for new position [0, 1, 1]
    // Change in energy: +4
    face_vertex_positions[2][2] = DynamicVar(1);

    DynamicVar energy = triangle_energy.triangle_vertex_least_squares_energy(
      face_vertex_positions,
      face_vertex_gradients,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 8.0) );
  }
}


TEST_CASE ( "The quadratic surface mapping triangle energy can be computed" )
{
  spdlog::set_level(spdlog::level::debug);
  TriangleEnergy<DynamicVar>  triangle_energy;

  // Positions
  // Note: Order is important here. Swapping the order changes which two coordinates
  // of u, v, and w are used
  DynamicVector3v pi(3);
  DynamicVector3v pj(3);
  DynamicVector3v pk(3);
  pi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0);
  pj <<
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  pk <<
    DynamicVar(0), DynamicVar(0), DynamicVar(0);
  std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };

  // Gradients
  DynamicMatrix2x3v Gi(2, 3);
  Gi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gj(2, 3);
  Gj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gk(2, 3);
  Gk <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };

  // UV
  Matrix2x2r UVi(2, 2);
  UVi <<
    -1, 1,
    -1, 0;
  Matrix2x2r UVj(2, 2);
  UVj <<
    0, -1,
    1, -1;
  Matrix2x2r UVk(2, 2);
  UVk <<
    1, 0,
    0, 1;
  std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };

  SECTION ( "Triangle zero energy" )
  {
    DynamicVar energy = triangle_energy.triangle_six_split_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2 + 2y^2 + 3xy
    // Derivatives are dz/dx = 2x + 3y and dz/dy = 3x + 4y with ||H||^2 = 38
    face_vertex_positions[0] <<
      DynamicVar(1), DynamicVar(0), DynamicVar(1);

    face_vertex_positions[1] <<
      DynamicVar(0), DynamicVar(1), DynamicVar(2);
    
    face_vertex_gradients[0] <<
      DynamicVar(1), DynamicVar(0), DynamicVar(2),
      DynamicVar(0), DynamicVar(1), DynamicVar(3);
    
    face_vertex_gradients[1] <<
      DynamicVar(1), DynamicVar(0), DynamicVar(3),
      DynamicVar(0), DynamicVar(1), DynamicVar(4);
    
    DynamicVar energy = triangle_energy.triangle_six_split_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 114.0) );
  }

}

TEST_CASE ( "The parametrized quadratic surface mapping triangle energy can be computed" )
{
  spdlog::set_level(spdlog::level::debug);
  TriangleEnergy<DynamicVar>  triangle_energy;

  // Positions
  // Note: Order is important here. Swapping the order changes which two coordinates
  // of u, v, and w are used
  DynamicVector3v pi(3);
  DynamicVector3v pj(3);
  DynamicVector3v pk(3);
  pi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0);
  pj <<
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  pk <<
    DynamicVar(0), DynamicVar(0), DynamicVar(0);

  // Gradients
  DynamicMatrix2x3v Gi(2, 3);
  Gi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gj(2, 3);
  Gj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gk(2, 3);
  Gk <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);

  // uv
  PlanarPoint uvi(2);
  uvi <<
    1, 0;
  PlanarPoint uvj(2);
  uvj <<
    0, 1;
  PlanarPoint uvk(2);
  uvk <<
    0, 0;

  // UV
  Matrix2x2r UVi(2, 2);
  UVi <<
    -1, 1,
    -1, 0;
  Matrix2x2r UVj(2, 2);
  UVj <<
    0, -1,
    1, -1;
  Matrix2x2r UVk(2, 2);
  UVk <<
    1, 0,
    0, 1;

  SECTION ( "Triangle zero energy" )
  {
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

  SECTION ( "Triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2 + 2y^2 + 3xy
    // Derivatives are dz/dx = 2x + 3y and dz/dy = 3x + 4y with ||H||^2 = 38
    pi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(1);
    pj <<
      DynamicVar(0), DynamicVar(1), DynamicVar(2);
    Gi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(2),
      DynamicVar(0), DynamicVar(1), DynamicVar(3);
    Gj <<
      DynamicVar(1), DynamicVar(0), DynamicVar(3),
      DynamicVar(0), DynamicVar(1), DynamicVar(4);

    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 19.0) );
  }

  SECTION ( "Scaled triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2 + 2y^2 + 3xy
    // Derivatives are dz/dx = 2x + 3y and dz/dy = 3x + 4y with ||H||^2 = 38
    pi <<
      DynamicVar(2), DynamicVar(0), DynamicVar(4);
    pj <<
      DynamicVar(0), DynamicVar(1), DynamicVar(2);
    Gi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(4),
      DynamicVar(0), DynamicVar(1), DynamicVar(6);
    Gj <<
      DynamicVar(1), DynamicVar(0), DynamicVar(3),
      DynamicVar(0), DynamicVar(1), DynamicVar(4);
    uvi <<
      2, 0;
    UVi <<
      -2, 1,
      -2, 0;
    UVj <<
      0, -1,
      2, -1;
    UVk <<
      2, 0,
      0, 1;

    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 38.0) );
  }

  SECTION ( "Nonstandard triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2 + 2y^2 + 3xy
    // Derivatives are dz/dx = 2x + 3y and dz/dy = 3x + 4y with ||H||^2 = 38
    pi <<
      DynamicVar(2), DynamicVar(1), DynamicVar(12);
    pj <<
      DynamicVar(0), DynamicVar(1), DynamicVar(2);
    Gi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(7),
      DynamicVar(0), DynamicVar(1), DynamicVar(10);
    Gj <<
      DynamicVar(1), DynamicVar(0), DynamicVar(3),
      DynamicVar(0), DynamicVar(1), DynamicVar(4);
    uvi <<
      2, 1;
    UVi <<
      -2,  0,
      -2, -1;
    UVj <<
      0, -1,
      2,  0;
    UVk <<
      2, 1,
      0, 1;

    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 38.0) );
  }

  SECTION ( "Perturbed triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2
    // Derivatives are dz/dx = 2x and dz/dy = 0 over an oblique triangle
    pi <<
      DynamicVar(0.1), DynamicVar(0), DynamicVar(0.01);
    pj <<
      DynamicVar(1), DynamicVar(0), DynamicVar(1);
    pk <<
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    Gi <<
      DynamicVar(1), DynamicVar(0), DynamicVar(0.2),
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    Gj <<
      DynamicVar(1), DynamicVar(0), DynamicVar(2),
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    Gk <<
      DynamicVar(1), DynamicVar(0), DynamicVar(0),
      DynamicVar(0), DynamicVar(1), DynamicVar(0);
    uvi <<
      0.1, 0;
    uvj <<
      1, 0;
    uvk <<
      0, 1;
    
    UVi.row(0) = uvj - uvi;
    UVi.row(1) = uvk - uvi;
    
    UVj.row(0) = uvk - uvj;
    UVj.row(1) = uvi - uvj;

    UVk.row(0) = uvi - uvk;
    UVk.row(1) = uvj - uvk;
    
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.45 * 4) );
  }

  SECTION ( "Second perturbed triangle energy" )
  {
    // Sample positions and gradients from quadratic surface z = x^2
    // Derivatives are dz/dx = 2x and dz/dy = 0 over an oblique triangle
    pi <<
      DynamicVar(-0.4), DynamicVar(-0.5), DynamicVar(0.16);
    pj <<
      DynamicVar(-0.5), DynamicVar(0.5), DynamicVar(0.25);
    pk <<
      DynamicVar(-1.5), DynamicVar(0.5), DynamicVar(2.25);
    Gi <<
      DynamicVar(-0.1), DynamicVar(-1), DynamicVar(0.0809476),
      DynamicVar(1), DynamicVar(-0.1), DynamicVar(-0.805201);
    Gj <<
      DynamicVar(-1), DynamicVar(0), DynamicVar(1),
      DynamicVar(0), DynamicVar(-1), DynamicVar(0);
    Gk <<
      DynamicVar(0), DynamicVar(-1), DynamicVar(0),
      DynamicVar(1), DynamicVar(0), DynamicVar(-3);
    uvi <<
      -0.4, -0.5;
    uvj <<
      -0.5, 0.5;
    uvk <<
      -1.5, 0.5;
    
    UVi <<
    -0.980198, -0.19802,
    -0.881188, -1.18812;
    UVj <<
      1, 0,
      -0.1, 1;
    UVk <<
      1, 1.1,
      0, 1;

    // Expected: 1.9945806273025122 = 0.49999999999999994 * 3.989161254605025
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    //REQUIRE( float_equal(energy_value, 0.45 * 4) );
  }
}

TEST_CASE ( "The parametrized quadratic surface mapping triangle energy is permutation invariant" )
{
  spdlog::set_level(spdlog::level::debug);
  TriangleEnergy<DynamicVar>  triangle_energy;

  // Positions
  // Note: Order is important here. Swapping the order changes which two coordinates
  // of u, v, and w are used
  DynamicVector3v pi(3);
  DynamicVector3v pj(3);
  DynamicVector3v pk(3);
  pi <<
    DynamicVar(0), DynamicVar(0), DynamicVar(0);
  pj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0);
  pk <<
    DynamicVar(0), DynamicVar(1), DynamicVar(0);

  // Gradients
  DynamicMatrix2x3v Gi(2, 3);
  Gi <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gj(2, 3);
  Gj <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);
  DynamicMatrix2x3v Gk(2, 3);
  Gk <<
    DynamicVar(1), DynamicVar(0), DynamicVar(0),
    DynamicVar(0), DynamicVar(1), DynamicVar(0);

  // uv
  PlanarPoint uvi(2);
  uvi <<
    0, 0;
  PlanarPoint uvj(2);
  uvj <<
    1, 0;
  PlanarPoint uvk(2);
  uvk <<
    0, 1;

  // UV
  Matrix2x2r UVi(2, 2);
  UVi <<
    1, 0,
    0, 1;
  Matrix2x2r UVj(2, 2);
  UVj <<
    -1, 1,
    -1, 0;
  Matrix2x2r UVk(2, 2);
  UVk <<
    0, -1,
    1, -1;

  SECTION ( "Triangle zero energy" )
  {
    std::array<DynamicVector3v, 3> face_vertex_positions = { pi, pj, pk };
    std::array<DynamicMatrix2x3v, 3> face_vertex_gradients = { Gi, Gj, Gk };
    std::array<PlanarPoint, 3> face_vertex_uv_positions = { uvi, uvj, uvk };
    std::array<Matrix2x2r, 3> corner_to_corner_uv_positions = { UVi, UVj, UVk };
    DynamicVar energy = triangle_energy.triangle_six_split_parametrized_quadratic_surface_mapping_energy(
      face_vertex_positions,
      face_vertex_gradients,
      face_vertex_uv_positions,
      corner_to_corner_uv_positions
    );
    double energy_value = compute_variable_value(energy);
    REQUIRE( float_equal(energy_value, 0.0) );
  }

//  SECTION ( "Triangle energy" )
//  {
//    // Sample positions and gradients from quadratic surface z = x^2 + 2y^2 + 3xy
//    // Derivatives are dz/dx = 2x + 3y and dz/dy = 3x + 4y with ||H||^2 = 38
//    pj <<
//      DynamicVar(1), DynamicVar(0), DynamicVar(1);
//    pk <<
//      DynamicVar(0), DynamicVar(1), DynamicVar(2);
//    Gj <<
//      DynamicVar(1), DynamicVar(0), DynamicVar(2),
//      DynamicVar(0), DynamicVar(1), DynamicVar(3);
//    Gk <<
//      DynamicVar(1), DynamicVar(0), DynamicVar(3),
//      DynamicVar(0), DynamicVar(1), DynamicVar(4);
//    DynamicVar energy = triangle_energy.triangle_parametrized_quadratic_surface_mapping_energy(
//       pi,  pj,  pk,
//       Gi,  Gj,  Gk,
//      uvi, uvj, uvk,
//      UVi, UVj, UVk
//      );
//    double energy_value = compute_variable_value(energy);
//    REQUIRE( float_equal(energy_value, 19.0) );
//  }

}