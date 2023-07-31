add_library(ContourNetworkLib
  compute_contours.cpp
  compute_closed_contours.cpp
  compute_curve_frame.cpp
  compute_cusps.cpp
  compute_intersections.cpp
  compute_rational_bezier_curve_intersection.cpp
  compute_ray_intersections.cpp
  compute_ray_intersections_pencil_method.cpp
  contour_network.cpp
  discretize.cpp
  evaluate_general_contour_frame.cpp
  intersection_heuristics.cpp
  project_curves.cpp
  projected_curve_network.cpp
  ray_intersection_heuristics.cpp
  sample_general_contour_frame.cpp
  split_contours.cpp
  validity.cpp
  write_output.cpp
)
target_include_directories(ContourNetworkLib PUBLIC .)
target_link_libraries(ContourNetworkLib PUBLIC
  QuadraticContoursCoreLib
  QuadraticSplineSurfaceLib
)
target_compile_definitions(ContourNetworkLib PUBLIC
  SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
)

target_compile_options(ContourNetworkLib PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)