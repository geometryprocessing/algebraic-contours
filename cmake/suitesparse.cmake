if(TARGET SuiteSparse::SuiteSparseConfig)
  return()
endif()

include(FetchContent)
FetchContent_Declare(
    suitesparse
	SYSTEM
    GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse.git
	GIT_TAG stable
)
FetchContent_MakeAvailable(suitesparse)
