add_executable(add_header_to_csv add_header_to_csv.cpp)
target_link_libraries(add_header_to_csv PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(add_header_to_csv PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(animate_approximation_error animate_approximation_error.cpp)
target_link_libraries(animate_approximation_error PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(animate_approximation_error PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(animate_rotation animate_rotation.cpp)
target_link_libraries(animate_rotation PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(animate_rotation PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_example_figure generate_example_figure.cpp)
target_link_libraries(generate_example_figure PRIVATE
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_example_figure PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_perspective_figure generate_perspective_figure.cpp)
target_link_libraries(generate_perspective_figure PRIVATE 
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_perspective_figure PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_perspective_distortion_figure generate_perspective_distortion_figure.cpp)
target_link_libraries(generate_perspective_distortion_figure PRIVATE 
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_perspective_distortion_figure PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_smoothness_figure generate_smoothness_figure.cpp)
target_link_libraries(generate_smoothness_figure PRIVATE 
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_smoothness_figure PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)

add_executable(generate_timing_data generate_timing_data.cpp)
target_link_libraries(generate_timing_data PRIVATE 
  QuadraticContoursLib
  CLI11::CLI11
)
target_compile_options(generate_timing_data PRIVATE
  -Wall -Wpedantic -Wextra -Werror
)
