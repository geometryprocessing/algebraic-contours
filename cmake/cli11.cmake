if(TARGET CLI11::CLI11)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    cli11
    SYSTEM
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.3.2
)
FetchContent_MakeAvailable(cli11)