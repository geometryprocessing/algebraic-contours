add_executable(generate_algebraic_contours generate_algebraic_contours.cpp)
target_link_libraries(generate_algebraic_contours PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_algebraic_contours PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)