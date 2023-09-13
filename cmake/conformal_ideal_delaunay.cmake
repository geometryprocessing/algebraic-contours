include(FetchContent)
FetchContent_Declare(
    conformal_ideal_delaunay
    SYSTEM
    GIT_REPOSITORY https://github.com/rjc8237/ConformalIdealDelaunay.git
    GIT_TAG clean_layout
)
FetchContent_MakeAvailable(conformal_ideal_delaunay)

