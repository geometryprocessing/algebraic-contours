add_executable(generate_algebraic_contours generate_algebraic_contours.cpp)
target_link_libraries(generate_algebraic_contours PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_algebraic_contours PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(view_quadratic_spline view_quadratic_spline.cpp)
target_link_libraries(view_quadratic_spline PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(view_quadratic_spline PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_hidden_cones generate_hidden_cones.cpp)
target_link_libraries(generate_hidden_cones PRIVATE
  QuadraticContoursLib
  ConeGenerationLib
  CLI11::CLI11
)
target_compile_options(generate_hidden_cones PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(parameterize_mesh parameterize_mesh.cpp)
target_link_libraries(parameterize_mesh PRIVATE
  QuadraticContoursLib
  conformal_cpp
  Eigen3::Eigen
  spdlog::spdlog
  igl::core
  CLI11::CLI11
)
target_compile_options(parameterize_mesh PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)
