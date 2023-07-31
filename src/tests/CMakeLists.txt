# Create executable for tests
Include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1
)
FetchContent_MakeAvailable(Catch2)

# Make testing library
add_library(QuadraticContoursTestLib
  test_affine_manifold.cpp
  test_apply_transformation.cpp
  test_assemble_matrix.cpp
  test_convert_conic.cpp
  test_compute_ray_intersections.cpp
  test_conic.cpp
  test_convert_conic.cpp
  test_evaluate_surface.cpp
  test_halfedge.cpp
  test_optimize_spline_surface.cpp
  test_parametrize_conic.cpp
  test_polynomial_function.cpp
  test_projected_curve_network.cpp
  test_quadratic_spline_surface.cpp
  test_rational_function.cpp
  test_twelve_split_spline.cpp
  test_vertex_circulator.cpp
#  test_position_data.cpp
#  test_intersect_conic.cpp
#  test_six_split_spline.cpp
#  test_triangle_energy.cpp
#  test_compute_intersections.cpp
)
target_link_libraries(QuadraticContoursTestLib PRIVATE
  Catch2::Catch2WithMain
  QuadraticContoursLib
  QuadraticContoursUtilLib
)

# Make testing executable
add_executable(QuadraticContoursTests test.cpp)
target_link_libraries(QuadraticContoursTests PRIVATE
  Catch2::Catch2WithMain
  QuadraticContoursLib
  QuadraticContoursUtilLib
  QuadraticContoursTestLib
)

# Append all testing extras
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
