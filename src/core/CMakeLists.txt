# Build core library
add_library(QuadraticContoursCoreLib
  abstract_curve_network.cpp
  apply_transformation.cpp
  affine_manifold.cpp
  bivariate_quadratic_function.cpp
  compute_boundaries.cpp
  conic.cpp
  convert_conic.cpp
  convex_polygon.cpp
  evaluate_surface_normal.cpp
  generate_flat_metric.cpp
  halfedge.cpp
  generate_colormap.cpp
  generate_transformation.cpp
  globals.cpp
  intersect_conic.cpp
  interval.cpp
  line_segment.cpp
  parametrize_conic.cpp
  polynomial_function.cpp
  vertex_circulator.cpp
)
target_include_directories(QuadraticContoursCoreLib PUBLIC .)
target_link_libraries(QuadraticContoursCoreLib PUBLIC
  Eigen3::Eigen
  igl::core
  polyscope
  spdlog::spdlog
  SuiteSparse::CHOLMOD
  SVGLib
)

# Optionall include GoTools dependencies
if (USE_GOTOOLS)
  target_include_directories(QuadraticContoursCoreLib PUBLIC
    ${gotools_SOURCE_DIR}/gotools-core/include
    ${gotools_SOURCE_DIR}/intersections/include
  )
  target_link_libraries(QuadraticContoursCoreLib PUBLIC
    GoToolsCore
    GoIntersections
  )
endif()

# Set compile time debugger level
target_compile_definitions(QuadraticContoursCoreLib PUBLIC
  SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
  # SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_OFF
)
target_compile_options(QuadraticContoursCoreLib PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

