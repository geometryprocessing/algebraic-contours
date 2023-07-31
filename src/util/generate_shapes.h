// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <Eigen/Core>

#include "common.h"
#include "conic.h"
#include <vector>
#include <random>
#include "position_data.h"
#include "affine_manifold.h"



int flatten(
  int i,
  int j,
  int n
);

// Generate circle of given radius
//
// param[in] radius: radius of the circle
// return: parametrized circle (missing point at bottom)
Conic generate_circle(
  double radius
);

SpatialVector generate_torus_point(
  double major_radius,
  double minor_radius,
  int i,
  int j,
  int resolution,
  double angle_offset = 0.0
);

double generate_angle(
  double i,
  int resolution,
  double angle_offset = 0.0
);

double generate_angle_derivative(
  int resolution
);

/// Generate a quadratic surface with an ellipse as the parametric contour.
///
/// @param[out] surface_mapping_coeffs: Coefficients for the quadratic surface
/// @param[out] normal_mapping_coeffs: Coefficients for the quadratic surface normal
void generate_elliptic_contour_quadratic_surface(
  Matrix6x3r &surface_mapping_coeffs,
  Matrix6x3r &normal_mapping_coeffs
);


// ***************
// VF construction
// ***************

void generate_equilateral_triangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double length=1
);


void generate_right_triangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double width=1,
  double height=1
);


void generate_rectangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double width=1,
  double height=1
);


void generate_square_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double length=1
);

// Generate simple tetrahedron mesh
//
// param[out] V: Tetrahedron vertices
// param[out] F: Tetrahedron facets
void generate_tetrahedron_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F
);


// Generate simple torus mesh
//
// param[out] V: Torus vertices
// param[out] F: Torus facets
void generate_minimal_torus_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double major_radius=3,
  double minor_radius=1
);


// ********************
// Polygon construction
// ********************

// Generate rectangle with lower left corner (x0, y0) and upper right corner (x0, y0)
//
// param[in] x0, y0: lower left coordinates
// param[in] x1, y1: upper right coordinates
// param[out] rectangle_boundary_coeffs: half-plane coefficients for the rectangle
void generate_rectangle(
  double x0,
  double y0,
  double x1,
  double y1,
  std::vector<Eigen::Matrix<double, 3, 1>> &rectangle_boundary_coeffs
);


// Generate square centered at the origin with given side length
//
// param[in] length: square side length
// param[out] rectangle_boundary_coeffs: half-plane coefficients for the square
void generate_square(
  double length,
  std::vector<Eigen::Matrix<double, 3, 1>> &square_boundary_coeffs
);


// ***********************
// Point grid construction
// ***********************


// Generate a grid of uniformly spaced points in the xy plane.
//
// param[in] resolution: number of points in the x and y directions
// param[in] delta: spacing between points
// param[in] x0: x coordinate of lower left corner of grid
// param[in] y0: y coordinate of lower left corner of grid
// param[out] point_grid: output point grid
void generate_plane_grid(
  int resolution,
  double delta,
  double x0,
  double y0,
  std::vector< std::vector<SpatialVector> > &point_grid
);

// Generate a grid of points on a torus. The spacing is uniform in the usual
// angular parametrization of the torus.
//
// param[in] resolution: number of points in the x and y directions
// param[in] major_radius: radius of the major circle of the torus
// param[in] minor_radius: radius of the minor circle of the torus
// param[out] point_grid: output point grid
void generate_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  std::vector< std::vector<SpatialVector> > &point_grid
);

// Generate a grid of points on a plane. 
//
// param[in] resolution: number of points in the x and y directions
// param[in] u_slope: slope in the u direction
// param[in] v_slope: slope in the v direction
// param[out] point_grid: output point grid
void generate_plane_grid(
  int resolution,
  double u_slope,
  double v_slope,
  std::vector< std::vector<SpatialVector> > &point_grid
);

// Generate a grid of points on a quadratic surface. 
//
// param[in] resolution: number of points in the x and y directions
// param[in] u_curvature: curvature in the u direction
// param[in] v_curvature: curvature in the v direction
// param[out] point_grid: output point grid
void generate_quadratic_grid(
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  double u_curvature,
  double v_curvature,
  double uv_curvature,
  std::vector< std::vector<SpatialVector> > &point_grid
);


/// Generate a grid of points on a torus with sinusoidal variation in the major radius.
///
/// The spacing is uniform in the usual angular parametrization of the torus.
///
/// @param[in] resolution: number of points in the x and y directions
/// @param[in] major_radius: radius of the major circle of the torus
/// @param[in] minor_radius: radius of the minor circle of the torus
/// @param[in] amplitude: amplitude of the sinusoidal variations
/// @param[in] frequency: number of waves
/// @param[out] point_grid: output point grid
void generate_sinusoidal_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  double amplitude,
  double frequency,
  std::vector< std::vector<SpatialVector> > &point_grid
);

// Generate a grid of points on a bumpy torus. The spacing is uniform in the usual
// angular parametrization of the torus.
//
// param[in] resolution: number of points in the x and y directions
// param[in] major_radius: radius of the major circle of the torus
// param[in] minor_radius: radius of the minor circle of the torus
// param[in] stddev: standard deviation for random perturbations
// param[out] point_grid: output point grid
void generate_bumpy_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  double stddev,
  std::vector< std::vector<SpatialVector> > &point_grid
);

// Generate a grid of points on a quadratic surface with a unit circle as
// its contour.
//
// param[in] length: number of points in the x and y directions
// param[in] x0: x translation of the center of the ellipse
// param[in] y0: x translation of the center of the ellipse
// param[out] point_grid: output point grid
void generate_ellipse_contour_quadratic_grid(
  int length,
  double x0,
  double y0,
  std::vector< std::vector<SpatialVector> > &point_grid
);


// Generate a grid of points on a quadratic surface with a hyperbola
// its contour.
//
// param[in] length: number of points in the x and y directions
// param[in] x0: x translation of the center of the hyperbola
// param[in] y0: x translation of the center of the hyperbola
// param[out] point_grid: output point grid
void generate_hyperbola_contour_quadratic_grid(
  int length,
  double x0,
  double y0,
  std::vector< std::vector<SpatialVector> > &point_grid
);


void generate_perturbed_quadratic_grid(
  int resolution,
  double u_curvature,
  double v_curvature,
  std::vector< std::vector<SpatialVector> > &point_grid
);


// ***********************
// General utility methods
// ***********************


// Generate a VF mesh for a point grid for visualization.
//
// param[in] point_grid: input point grid
// param[out] V: vertices of the output mesh (same as the point grid)
// param[out] F: triangulation of the point grid
// param[out] l: uv parametrizatioin lengths for the grid
template <typename Point>
void generate_mesh_from_grid(
  const std::vector< std::vector<Point> > &point_grid,
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  std::vector<std::vector<double>> &l,
  bool closed_surface
) {
  int n = point_grid.size();
  assert(n > 0);
  assert(n == point_grid[0].size());
  V.resize(n * n, point_grid[0][0].size());

  // Flatten vertices in the point grid to a standard vector V of vertices
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      V.row(flatten(i,j,n)) = point_grid[i][j];
    }
  }

  // Use periodic boundary triangulation if closed surface
  int N;
  if (closed_surface) {
    N = n;
  } else {
    N = n - 1;
  }

  // Create triangulation F for the grid
  F.resize(2 * N * N, 3);
  l.resize(2 * N * N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      // Create first face
      F(2*flatten(i, j, N), 0) = flatten(i, j, n);
      F(2*flatten(i, j, N), 1) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N), 2) = flatten(i, (j + 1) % n, n);
      l[2*flatten(i, j, N)].resize(3);
      l[2*flatten(i, j, N)][0] = std::sqrt(2);
      l[2*flatten(i, j, N)][1] = 1.0;
      l[2*flatten(i, j, N)][2] = 1.0;

      // Create second face
      F(2*flatten(i, j, N) + 1, 0) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N) + 1, 1) = flatten((i + 1) % n, (j + 1) % n, n);
      F(2*flatten(i, j, N) + 1, 2) = flatten(i, (j + 1) % n, n);
      l[2*flatten(i, j, N) + 1].resize(3);
      l[2*flatten(i, j, N) + 1][0] = 1.0;
      l[2*flatten(i, j, N) + 1][1] = std::sqrt(2);
      l[2*flatten(i, j, N) + 1][2] = 1.0;
    }
  }

  spdlog::debug("Faces:\n{}", F);
}

void generate_mesh_from_grid(
  const std::vector< std::vector<SpatialVector> > &point_grid,
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  std::vector<std::vector<double>> &l,
  bool closed_surface
);


void generate_interior_faces(
  int resolution,
  std::vector<int> &faces
);


void generate_gradients_from_grid(
  const std::vector< std::vector<Matrix2x3r> > &gradient_grid,
  std::vector<Matrix2x3r> &gradients,
  std::vector<int> &boundary_vertices
);

void generate_quadratic_gradients_grid(
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  double u_curvature,
  double v_curvature,
  double uv_curvature,
  std::vector< std::vector<Matrix2x3r> > &gradient_grid
);

void generate_global_layout_grid(
  int resolution,
  std::vector< std::vector<PlanarPoint> > &layout_grid
);

