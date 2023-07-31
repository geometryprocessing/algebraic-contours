cmake_minimum_required(VERSION 3.20)
project(SVG VERSION 0.0.1)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# Get required packages
list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(eigen)

# Build SVG library
add_library(SVGLib
  src/svg.cpp
)
target_include_directories(SVGLib PUBLIC src)
target_link_libraries(SVGLib PUBLIC
  Eigen3::Eigen
)
