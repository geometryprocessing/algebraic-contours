# Distributed under the OSI-approved MIT License.  See accompanying
# file LICENSE or https://github.com/Crascit/DownloadProject for details.

cmake_minimum_required(VERSION 2.8.2)

project(tbb-download NONE)

include(ExternalProject)
ExternalProject_Add(tbb-download
                    GIT_CONFIG advice.detachedHead=false;GIT_REPOSITORY;https://github.com/wjakob/tbb.git;GIT_TAG;20357d83871e4cb93b2c724fe0c337cd999fd14f
                    SOURCE_DIR          "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb"
                    BINARY_DIR          "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/build-rel-with-deb-info/tbb-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
)
