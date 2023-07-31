# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/build-rel-with-deb-info/tbb-build"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/tmp"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src"
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp${cfgdir}") # cfgdir has leading slash
endif()
