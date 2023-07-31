# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitclone-lastrun.txt" AND EXISTS "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitinfo.txt" AND
  "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitclone-lastrun.txt" IS_NEWER_THAN "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" 
            clone --no-checkout --config "advice.detachedHead=false" --config "advice.detachedHead=false" "https://github.com/wjakob/tbb.git" "tbb"
    WORKING_DIRECTORY "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/wjakob/tbb.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git" 
          checkout "20357d83871e4cb93b2c724fe0c337cd999fd14f" --
  WORKING_DIRECTORY "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '20357d83871e4cb93b2c724fe0c337cd999fd14f'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/tbb'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitinfo.txt" "/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/rcapouel/work/research/2022_quadratic_contours/quadratic_contours/ext/.cache/tbb/tbb-download-prefix/src/tbb-download-stamp/tbb-download-gitclone-lastrun.txt'")
endif()
