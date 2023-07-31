// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include <igl/is_vertex_manifold.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>
// #include <gmm/gmm.h>
// #include <igl/copyleft/comiso/nrosy.h>
#include "affine_manifold.h"
// #include "conformal_ideal_delaunay/ConformalInterface.hh"

/// Given a VF mesh with a global uv layout, remove the cones where the layout
/// does not result in a flat metric.
///
/// @param[in] V: input mesh vertices
/// @param[in] F: input mesh faces
/// @param[in] uv: input mesh uv coordinates
/// @param[in] F_uv: input mesh uv layout faces
/// @param[in] V_flat: mesh vertices with cones removed
/// @param[in] F_flat: mesh faces with cones removed
/// @param[in] uv_flat: mesh uv coordinates with cones removed
/// @param[in] F_uv_flat: mesh uv layout faces with cones removed
void
remove_cones_from_uv(const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const Eigen::MatrixXd& uv,
                     const Eigen::MatrixXi& F_uv,
                     Eigen::MatrixXd& V_flat,
                     Eigen::MatrixXi& F_flat,
                     Eigen::MatrixXd& uv_flat,
                     Eigen::MatrixXi& F_uv_flat);

void
generate_metric_from_uv(const Eigen::MatrixXi& F,
                        const MatrixXr& uv,
                        std::vector<std::vector<double>>& l);
