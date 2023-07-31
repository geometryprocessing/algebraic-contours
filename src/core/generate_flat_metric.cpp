// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "generate_flat_metric.h"

#include <igl/facet_components.h>

void
generate_metric_from_uv(const Eigen::MatrixXi& F,
                        const MatrixXr& uv,
                        std::vector<std::vector<double>>& l)
{
  Eigen::Index num_faces = F.rows();
  Eigen::Index face_size = F.cols();
  assert(face_size == 3);
  l.resize(num_faces);

  // Iterate over faces
  for (Eigen::Index i = 0; i < num_faces; ++i) {
    l[i].resize(face_size);

    // Iterate over vertices in face i
    for (Eigen::Index j = 0; j < face_size; ++j) {
      // Get the length of the edge opposite face vertex j
      PlanarPoint prev_uv = uv.row(F(i, (j + 2) % face_size));
      PlanarPoint next_uv = uv.row(F(i, (j + 1) % face_size));
      PlanarPoint edge_vector = prev_uv - next_uv;
      l[i][j] = edge_vector.norm();
    }
  }
}

void
remove_cones_from_uv(const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const Eigen::MatrixXd& uv,
                     const Eigen::MatrixXi& F_uv,
                     Eigen::MatrixXd& V_flat,
                     Eigen::MatrixXi& F_flat,
                     Eigen::MatrixXd& uv_flat,
                     Eigen::MatrixXi& F_uv_flat)
{
  // Get the cones of the metric
  std::vector<AffineManifold::Index> cones;
  AffineManifold cone_manifold(F, uv, F_uv);
  cone_manifold.compute_cones(cones);
  SPDLOG_TRACE("Removing cones at {}", formatted_vector(cones));

  // Restrict the manifold to the flat vertices
  std::vector<Eigen::Index> removed_faces;
  remove_mesh_vertices(V, F, cones, V_flat, F_flat, removed_faces);

  // Remove the faces around cones from the uv layout as well
  remove_mesh_faces(uv, F_uv, removed_faces, uv_flat, F_uv_flat);
}
