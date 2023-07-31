// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once
#include "common.h"

void generate_uv_square(
  int resolution,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &F
);


void generate_uv_torus(
  int resolution,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &F
);
