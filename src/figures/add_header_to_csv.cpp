// Copyright 2023 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include <fstream>
#include <cassert>
#include "common.h"

// Helper executable to add a header to the timing csv file

int main(int argc, char *argv[]) {
  if (argc < 1) return 0;
  std::string output_dir = argv[1];

  std::ofstream out_view_independent(join_path(output_dir, "view_independent.csv"),
                                     std::ios::app);
  std::ofstream out_per_view(join_path(output_dir, "per_view.csv"), std::ios::app);

  out_view_independent << "mesh name, num triangles, time spline surface, time "
                          "patch boundary edges,\n";
  out_per_view
      << "mesh name, rotation matrix, total time per view, surface update, "
         "compute contour, compute cusps, compute intersections, compute "
         "visibility, graph building, num segments, num int cusps, num bound "
         "cusps, num intersection call, num ray inter call, num patches\n";

  out_view_independent.close();
  out_per_view.close();

  return 0;
}