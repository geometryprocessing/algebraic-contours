if(TARGET SmoothSplineLib)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    smooth-spline
    SYSTEM
    GIT_REPOSITORY https://github.com/rjc8237/smooth-spline.git
    GIT_TAG main
)
FetchContent_MakeAvailable(smooth-spline)