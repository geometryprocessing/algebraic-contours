// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "generate_shapes.h"


// VF meshes

void generate_equilateral_triangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double length
) {
  V.resize(3, 3);
  F.resize(1, 3);

  V <<
    length,      0,      0,
         0, length,      0,
         0,      0, length;

  F << 
    0, 1, 2;
}


void generate_right_triangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double width,
  double height
) {
  V.resize(3, 3);
  F.resize(1, 3);

  V <<
        0,      0, 0,
    width,      0, 0,
        0, height, 0;

  F << 
    0, 1, 2;
}


void generate_rectangle_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double width,
  double height
) {
  V.resize(4, 3);
  F.resize(2, 3);

  V <<
        0,      0, 0,
    width,      0, 0,
    width, height, 0,
        0, height, 1;

  F << 
    0, 1, 2,
    1, 3, 2;
}

void generate_square_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double length
) {
  generate_rectangle_VF(V, F, length, length);
}


void generate_tetrahedron_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F
) {
  V.resize(4, 3);
  F.resize(4, 3);

  V <<
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1;

  F << 
    0, 2, 1,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3;
}

void generate_minimal_torus_VF(
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  double major_radius,
  double minor_radius
) {
  V.resize(9, 3);
  F.resize(18, 3);
  int resolution = 3;

  V.row(0) = generate_torus_point(major_radius, minor_radius, 0, 0, resolution);
  V.row(1) = generate_torus_point(major_radius, minor_radius, 0, 1, resolution);
  V.row(2) = generate_torus_point(major_radius, minor_radius, 0, 2, resolution);
  V.row(3) = generate_torus_point(major_radius, minor_radius, 1, 0, resolution);
  V.row(4) = generate_torus_point(major_radius, minor_radius, 1, 1, resolution);
  V.row(5) = generate_torus_point(major_radius, minor_radius, 1, 2, resolution);
  V.row(6) = generate_torus_point(major_radius, minor_radius, 2, 0, resolution);
  V.row(7) = generate_torus_point(major_radius, minor_radius, 2, 1, resolution);
  V.row(8) = generate_torus_point(major_radius, minor_radius, 2, 2, resolution);

  F << 
    0, 1, 3,
    1, 4, 3,
    1, 2, 4,
    2, 5, 4,
    2, 0, 5,
    0, 3, 5,
    3, 4, 6,
    4, 7, 6,
    4, 5, 7,
    5, 8, 7,
    5, 3, 8,
    3, 6, 8,
    6, 7, 0,
    7, 1, 0,
    7, 8, 1,
    8, 2, 1,
    8, 6, 2,
    6, 0 , 2;
}



int flatten(
  int i,
  int j,
  int n
) {
  return j * n + i;
}

void generate_plane_grid(
  int resolution,
  double delta,
  double x0,
  double y0,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      point_grid[i][j] = SpatialVector(3);
      point_grid[i][j](0) = x0 + (delta * i);
      point_grid[i][j](1) = y0 + (delta * j);
      point_grid[i][j](2) = 0.0;
    }
  }
}

double generate_angle(
  double i,
  int resolution,
  double angle_offset
) {
  return angle_offset + 2 * M_PI * i / float(resolution);
}

double generate_angle_derivative(
  int resolution
) {
  return 2 * M_PI / float(resolution);
}

SpatialVector generate_torus_point(
  double major_radius,
  double minor_radius,
  int i,
  int j,
  int resolution,
  double angle_offset
) {
  SpatialVector point;
  double theta = generate_angle(i, resolution, angle_offset);
  double phi = generate_angle(j, resolution, angle_offset);
  point(0) = (major_radius + minor_radius * std::cos(theta)) * std::cos(phi);
  point(1) = (major_radius + minor_radius * std::cos(theta)) * std::sin(phi);
  point(2) = minor_radius * std::sin(theta);

  return point;
}



void generate_interior_faces(
  int resolution,
  std::vector<int> &faces
) {
  faces.reserve(2 * resolution * resolution);
  for (int i = 2; i < resolution-2; ++i)
  {
    for (int j = 2; j < resolution-2; ++j)
    {
        faces.push_back(2*flatten(i, j, resolution) + 0);
        faces.push_back(2*flatten(i, j, resolution) + 1);
    }
  }
}

void generate_plane_grid(
  int resolution,
  double u_slope,
  double v_slope,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  double center = double(resolution) / 2.0;
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      double u = float(i) - center;
      double v = float(j) - center;
      point_grid[i][j] = SpatialVector(u, v, u_slope * u + v_slope * v);
    }
  }
}

void generate_perturbed_quadratic_grid(
  int resolution,
  double u_curvature,
  double v_curvature,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  double center = double(resolution) / 2.0;
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      double u = float(i) - center;
      double v = float(j) - center;

      point_grid[i][j] = SpatialVector(
        u + 0.1 * u * v,
        v + 0.1 * u * u + 0.1 * v * v,
        0.5 * u_curvature * u * u + 0.5 * v_curvature * v * v);
    }
  }
}

void generate_quadratic_grid(
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  double u_curvature,
  double v_curvature,
  double uv_curvature,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  point_grid.resize(layout_grid.size());
  for (int i = 0; i < layout_grid.size(); ++i) {
    point_grid[i].resize(layout_grid[i].size());
    for (int j = 0; j < layout_grid[i].size(); ++j) {
      double u = layout_grid[i][j][0];
      double v = layout_grid[i][j][1];
      point_grid[i][j] = SpatialVector(
        u,
        v,
        0.5 * u_curvature * u * u + 0.5 * v_curvature * v * v + uv_curvature * u * v);
    }
  }
}

void generate_quadratic_gradients_grid(
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  double u_curvature,
  double v_curvature,
  double uv_curvature,
  std::vector< std::vector<Matrix2x3r> > &gradient_grid
) {
  gradient_grid.resize(layout_grid.size());
  for (int i = 0; i < layout_grid.size(); ++i) {
    gradient_grid[i].resize(layout_grid[i].size());
    for (int j = 0; j < layout_grid[i].size(); ++j) {
      double u = layout_grid[i][j][0];
      double v = layout_grid[i][j][1];

      gradient_grid[i][j].resize(2, 3);
      gradient_grid[i][j] <<
        1, 0, u_curvature * u + uv_curvature * v,
        0, 1, v_curvature * v + uv_curvature * u;
    }
  }
}


void generate_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  double angle_offset = 0.1;
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      point_grid[i][j] = generate_torus_point(
        major_radius,
        minor_radius,
        i,
        j,
        resolution,
        angle_offset
      );
    }
  }
}


void generate_sinusoidal_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  double amplitude,
  double frequency,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  double angle_offset = 0.1;
  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      // Get sinusoidal major radius
      double phi = generate_angle(j, resolution, angle_offset);
      double sinusoidal_major_radius = major_radius + amplitude * std::cos(frequency * phi);
      point_grid[i][j] = generate_torus_point(
        sinusoidal_major_radius,
        minor_radius,
        i,
        j,
        resolution,
        angle_offset
      );
    }
  }
}


void generate_bumpy_torus_grid(
  int resolution,
  double major_radius,
  double minor_radius,
  double stddev,
  std::vector< std::vector<SpatialVector> > &point_grid
) {
  double angle_offset = 0.1;

  // Create random number generator
  std::random_device rd;
  std::mt19937 e2(rd());
  std::normal_distribution<double> gaussian(0, stddev);

  point_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i) {
    point_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j) {
      double pert = gaussian(e2);
      point_grid[i][j] = generate_torus_point(
        major_radius,
        minor_radius + pert,
        i,
        j,
        resolution,
        angle_offset
      );
    }
  }
}

void generate_global_layout_grid(
  int resolution,
  std::vector< std::vector<PlanarPoint> > &layout_grid
) {
  double center = double(resolution) / 2.0;
  layout_grid.resize(resolution);
  for (int i = 0; i < resolution; ++i)
  {
    layout_grid[i].resize(resolution);
    for (int j = 0; j < resolution; ++j)
    {
      layout_grid[i][j] = PlanarPoint(i - center, j - center);
    }
  }

}




void generate_gradients_from_grid(
  const std::vector< std::vector<Matrix2x3r> > &gradient_grid,
  std::vector<Matrix2x3r> &gradients,
  std::vector<int> &interior_vertices
) {
  int n = gradient_grid.size();
  assert(n > 0);
  assert(n == gradient_grid[0].size());
  gradients.resize(n * n);
  interior_vertices.reserve(n * n);

  // Flatten vertices in the point grid to a standard vector V of vertices
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      gradients[flatten(i,j,n)] = gradient_grid[i][j];

      // Record boundary vertices
      if (
           (i > 0)
        && (i < n - 1)
        && (j > 0)
        && (j < n - 1)
      ) {
        interior_vertices.push_back(flatten(i, j, n));
      }
    }
  }
  std::sort(interior_vertices.begin(), interior_vertices.end());
}


void generate_mesh_from_grid(
  const std::vector< std::vector<SpatialVector> > &point_grid,
  const std::vector< std::vector<PlanarPoint> > &layout_grid,
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  std::vector<std::vector<double>> &l,
  bool closed_surface
) {
  int n = point_grid.size();
  assert(n > 0);
  assert(n == point_grid[0].size());
  V.resize(n * n, 3);

  // Flatten vertices in the point grid to a standard vector V of vertices
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      V(flatten(i,j,n),0) = point_grid[i][j][0];
      V(flatten(i,j,n),1) = point_grid[i][j][1];
      V(flatten(i,j,n),2) = point_grid[i][j][2];
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
      VectorXr layout_00 = layout_grid[i][j];
      VectorXr layout_10 = layout_grid[(i + 1) % n][j];
      VectorXr layout_01 = layout_grid[i][(j + 1) % n];
      VectorXr layout_11 = layout_grid[(i + 1) % n][(j + 1) % n];

      // Create first face
      F(2*flatten(i, j, N), 0) = flatten(i, j, n);
      F(2*flatten(i, j, N), 1) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N), 2) = flatten(i, (j + 1) % n, n);
      l[2*flatten(i, j, N)].resize(3);
      l[2*flatten(i, j, N)][0] = (layout_01 - layout_10).norm();
      l[2*flatten(i, j, N)][1] = (layout_00 - layout_01).norm();
      l[2*flatten(i, j, N)][2] = (layout_10 - layout_00).norm();

      // Create second face
      F(2*flatten(i, j, N) + 1, 0) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N) + 1, 1) = flatten((i + 1) % n, (j + 1) % n, n);
      F(2*flatten(i, j, N) + 1, 2) = flatten(i, (j + 1) % n, n);
      l[2*flatten(i, j, N) + 1].resize(3);
      l[2*flatten(i, j, N) + 1][0] = (layout_01 - layout_11).norm();
      l[2*flatten(i, j, N) + 1][1] = (layout_10 - layout_01).norm();
      l[2*flatten(i, j, N) + 1][2] = (layout_11 - layout_10).norm();
    }
  }

  spdlog::debug("Faces:\n{}", F);
}



void generate_boundary_mesh_from_grid(
  const std::vector< std::vector<SpatialVector> > &point_grid,
  Eigen::MatrixXd &V,
  Eigen::MatrixXi &F,
  std::vector<std::vector<double>> &l,
  bool closed_surface
) {
  int n = point_grid.size();
  assert(n > 0);
  assert(n == point_grid[0].size());
  V.resize(n * n, 3);

  // Flatten vertices in the point grid to a standard vector V of vertices
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      V(flatten(i,j,n),0) = point_grid[i][j][0];
      V(flatten(i,j,n),1) = point_grid[i][j][1];
      V(flatten(i,j,n),2) = point_grid[i][j][2];
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

void generate_ellipse_contour_quadratic_grid(
  int length,
  double x0,
  double y0,
  std::vector< std::vector<VectorXr> > &point_grid
) {
  point_grid.resize(length);
  for (int i = 0; i < length; ++i) {
    point_grid[i].resize(length);
    for (int j = 0; j < length; ++j) {
      // Build point sample for torus
      double u = (i - x0);
      double v = (j - y0);
      point_grid[i][j] = VectorXr(3);
      point_grid[i][j](0) = 0.5 * u * u - 0.5 * v * v + v;
      point_grid[i][j](1) = u * v + u;
      point_grid[i][j](2) = v;
    }
  }

}

void generate_hyperbola_contour_quadratic_grid(
  int length,
  double x0,
  double y0,
  std::vector< std::vector<VectorXr> > &point_grid
) {
  point_grid.resize(length);
  for (int i = 0; i < length; ++i) {
    point_grid[i].resize(length);
    for (int j = 0; j < length; ++j) {
      // Build point sample for torus
      double u = (i - x0);
      double v = (j - y0);
      point_grid[i][j] = VectorXr(3);
      point_grid[i][j](0) = v + u * u;
      point_grid[i][j](1) = u + 0.5 * v * v;
      point_grid[i][j](2) = v;
    }
  }
}

// Halfplane convex polygons

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
) {
  rectangle_boundary_coeffs.resize(4);
  rectangle_boundary_coeffs[0] = Eigen::Matrix<double, 3, 1>(-x0,  1.0,  0.0);
  rectangle_boundary_coeffs[1] = Eigen::Matrix<double, 3, 1>(-y0,  0.0,  1.0);
  rectangle_boundary_coeffs[2] = Eigen::Matrix<double, 3, 1>( x1, -1.0,  0.0);
  rectangle_boundary_coeffs[3] = Eigen::Matrix<double, 3, 1>( y1,  0.0, -1.0);
}

void generate_square(
  double length,
  std::vector<Eigen::Matrix<double, 3, 1>> &square_boundary_coeffs
) {
  generate_rectangle(
    -length/2.0,
    -length/2.0,
    length/2.0,
    length/2.0,
    square_boundary_coeffs
  );
}

// Conics 

Conic generate_circle(
  double radius
) {
  MatrixXr P_coeffs(3, 2);
  VectorXr Q_coeffs(3);
  P_coeffs <<          0.0,  radius,
              2.0 * radius,     0.0,
                       0.0, -radius;
  Q_coeffs << 1.0,
              0.0,
              1.0;

  return Conic(P_coeffs, Q_coeffs);
}

// Quadratic surface mappings

void generate_elliptic_contour_quadratic_surface(
  Matrix6x3r &surface_mapping_coeffs,
  Matrix6x3r &normal_mapping_coeffs
) {
  // Build surface mapping coefficients
  surface_mapping_coeffs.setZero(6, 3);

  // x coordinate
  v_coeff<3>(surface_mapping_coeffs, 0) =  1.0;
  uv_coeff<3>(surface_mapping_coeffs, 0) =  0.5;
  vv_coeff<3>(surface_mapping_coeffs, 0) = -0.5;

  // y coordinate
  u_coeff<3>(surface_mapping_coeffs, 1) =  1.0;
  uv_coeff<3>(surface_mapping_coeffs, 1) =  1.0;

  // z coordinate
  v_coeff<3>(surface_mapping_coeffs, 2) =  1.0;

  // Build normal mapping coefficients
  normal_mapping_coeffs.setZero(6, 3);

  // x coordinate
  const_coeff<3>(normal_mapping_coeffs, 0) =  1.0;
  v_coeff<3>(normal_mapping_coeffs, 0) =  1.0;

  // y coordinate
  u_coeff<3>(normal_mapping_coeffs, 1) = -1.0;

  // z coordinate
  const_coeff<3>(normal_mapping_coeffs, 2) = -1.0;
  uu_coeff<3>(normal_mapping_coeffs, 2) =  1.0;
  vv_coeff<3>(normal_mapping_coeffs, 2) =  1.0;
}

void generate_uv_square(
  int resolution,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &F
) {
  int n = resolution;
  assert(n > 0);
  UV.resize(n * n, 2);

  // Build uv points
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double u = i * 1.0;
      double v = j * 1.0;
      UV.row(flatten(i,j,n)) << u, v;
    }
  }

  // Create triangulation F for the grid
  int N = n - 1;
  F.resize(2 * N * N, 3);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      // Create first face
      F(2*flatten(i, j, N), 0) = flatten(i, j, n);
      F(2*flatten(i, j, N), 1) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N), 2) = flatten(i, (j + 1) % n, n);

      // Create second face
      F(2*flatten(i, j, N) + 1, 0) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, N) + 1, 1) = flatten((i + 1) % n, (j + 1) % n, n);
      F(2*flatten(i, j, N) + 1, 2) = flatten(i, (j + 1) % n, n);
    }
  }
}


void generate_uv_torus(
  int resolution,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &F
) {
  int n = resolution;
  assert(n > 0);
  UV.resize(n * n, 2);

  // Build uv points
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double u = i * 1.0;
      double v = j * 1.0;
      UV.row(flatten(i,j,n)) << u, v;
    }
  }

  // Create triangulation F for the grid
  F.resize(2 * n * n, 3);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // Create first face
      F(2*flatten(i, j, n), 0) = flatten(i, j, n);
      F(2*flatten(i, j, n), 1) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, n), 2) = flatten(i, (j + 1) % n, n);

      // Create second face
      F(2*flatten(i, j, n) + 1, 0) = flatten((i + 1) % n, j, n);
      F(2*flatten(i, j, n) + 1, 1) = flatten((i + 1) % n, (j + 1) % n, n);
      F(2*flatten(i, j, n) + 1, 2) = flatten(i, (j + 1) % n, n);
    }
  }
}