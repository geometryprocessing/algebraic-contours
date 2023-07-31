if(TARGET polyscope)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    polyscope
    SYSTEM
    GIT_REPOSITORY https://github.com/rjc8237/polyscope.git
)
#    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
FetchContent_MakeAvailable(polyscope)
