// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "generate_mesh.h"

void generate_zwart_powell_mesh(int dimension)
{
  Eigen::VectorXd V;
  Eigen::VectorXi F;

  // Use trivial vertices
  V.resize(6 * dimension * dimension);
  V.setZero();


  // Build face matrix
  F.resize(12 * dimension * dimension);
  for (int i = 0; i < dimension; ++i)
  {
    for (int j = 0; j < dimension; ++j)
    {
      // Build vertex indices
      int v00 = flatten(i, j, dimension);
      int v10 = flatten((i + 1) % dimension, j, dimension);
      int v01 = flatten(i, (j + 1) % dimension, dimension);
      int v11 = flatten((i + 1) % dimension, (j + 1) % dimension, dimension);
      int v_center = flatten(i, j, dimension) + dimension * dimension;

      // Face 0
      F(3*(flatten(i, j, dimension) + 0) + 0) = v_center;
      F(3*(flatten(i, j, dimension) + 0) + 1) = v01;
      F(3*(flatten(i, j, dimension) + 0) + 2) = v00;

      // Face 1
      F(3*(flatten(i, j, dimension) + 1) + 0) = v_center;
      F(3*(flatten(i, j, dimension) + 1) + 1) = v10;
      F(3*(flatten(i, j, dimension) + 1) + 2) = v11;

      // Face 2
      F(3*(flatten(i, j, dimension) + 2) + 0) = v_center;
      F(3*(flatten(i, j, dimension) + 2) + 1) = v00;
      F(3*(flatten(i, j, dimension) + 2) + 2) = v10;
      
      // Face 3
      F(3*(flatten(i, j, dimension) + 3) + 0) = v_center;
      F(3*(flatten(i, j, dimension) + 3) + 1) = v11;
      F(3*(flatten(i, j, dimension) + 3) + 2) = v01;
    }
  }

}
