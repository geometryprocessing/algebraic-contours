# Build utility library
add_library(QuadraticContoursUtilLib
	generate_position_data.cpp
	generate_shapes.cpp
	generate_mesh.cpp
	generate_parametric_domains.cpp
)
target_include_directories(QuadraticContoursUtilLib PUBLIC .)
target_link_libraries(QuadraticContoursUtilLib PUBLIC
  QuadraticContoursLib
)
