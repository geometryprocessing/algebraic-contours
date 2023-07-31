// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "sample_general_contour_frame.h"

#include "discretize.h"
#include "evaluate_general_contour_frame.h"
#include "polynomial_function.h"

// void
// sample_spline_surface_general_contour_frame(
//   const QuadraticSplineSurface& spline_surface,
//   const MatrixXr& frame,
//   const SurfaceDiscretizationParameters& surface_disc_params,
//   MatrixXr& tangents,
//   MatrixXr& normals,
//   MatrixXr& tangent_normals)
//{
//   return;
/*
// Discretize parametric domain
Eigen::MatrixXd V;
Eigen::MatrixXi F;
discretize_parametric_domain(surface_disc_params, V, F);

// Evaluate frame at all parametric points
tangents.resize(V.rows(), V.cols());
normals.resize(V.rows(), V.cols());
tangent_normals.resize(V.rows(), V.cols());
for (int vi = 0; vi < V.rows(); ++vi) {
  double x = V(vi, 0);
  double y = V(vi, 1);
  VectorXr point(2);
  point << x, y;
  VectorXr tangent;
  VectorXr normal;
  VectorXr tangent_normal;
  evaluate_spline_surface_general_contour_frame(
    spline_surface,
    frame,
    point,
    tangent,
    normal,
    tangent_normal
  );
  tangents.row(vi) = tangent;
  normals.row(vi) = normal;
  tangent_normals.row(vi) = tangent_normal;

  // Validity frame check
  assert( float_equal(tangent.dot(normal), 0.0) );
  assert( float_equal(tangent.dot(tangent_normal), 0.0) );
  assert( float_equal(tangent_normal.dot(normal), 0.0) );
}
*/
//}
