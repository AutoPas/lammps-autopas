if (PKG_USER-AUTOPAS)

    if (CMAKE_VERSION VERSION_LESS "3.16")
        message(FATAL_ERROR "For the AutoPas package you need at least cmake-3.16")
    endif ()

    enable_language(C)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    ######### 1. Clone and build AutoPas library #########

    # Enable FetchContent CMake module
    include(FetchContent)

    # Select https (default) or ssh path.
    set(autopasRepoPath https://github.com/ssauermann/AutoPas.git)
    if (GIT_SUBMODULES_SSH)
        set(autopasRepoPath git@github.com:ssauermann/AutoPas.git)
    endif ()

    # Download and install autopas
    FetchContent_Declare(
            autopas
            GIT_REPOSITORY ${autopasRepoPath}
            # GIT_TAG f639d8b77eb62b84ffb3717ca4a3e25f1caaea86
            GIT_TAG origin/cmake-fetchcontent
            #GIT_TAG origin/feature/regionParticleIteratorIncrease
            #BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/autopas/build
            #BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/autopas/build/src/autopas/libautopas.a
            #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/autopas
    )

    option(AUTOPAS_BUILD_TESTS "" OFF)
    option(AUTOPAS_BUILD_EXAMPLES "" OFF)
    option(AUTOPAS_BUILD_TARGET_DOC "" OFF)
    option(AUTOPAS_ENABLE_ADDRESS_SANITIZER "" OFF)
    option(AUTOPAS_OPENMP "" ${BUILD_OMP})
    option(spdlog_ForceBundled "" ON)
    option(Eigen3_ForceBundled "" ON)
    option(yaml-cpp_ForceBundled "" ON)

    mark_as_advanced(AUTOPAS_BUILD_TESTS)
    mark_as_advanced(AUTOPAS_BUILD_EXAMPLES)
    mark_as_advanced(AUTOPAS_ENABLE_ADDRESS_SANITIZER)
    mark_as_advanced(AUTOPAS_OPENMP)

    FetchContent_MakeAvailable(autopas)

    list(APPEND LAMMPS_LINK_LIBS "autopas")


    ######### 2. AutoPas package settings #########

    set(USER-AUTOPAS_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-AUTOPAS)

    # Sources that are not a style
    set(USER-AUTOPAS_SOURCES #)${USER-AUTOPAS_SOURCES_DIR}/fix_autopas.cpp
            # ${USER-AUTOPAS_SOURCES_DIR}/thr_data.cpp
            # ${USER-AUTOPAS_SOURCES_DIR}/thr_omp.cpp
            )

    set_property(GLOBAL PROPERTY "AUTOPAS_SOURCES" "${USER-AUTOPAS_SOURCES}")

    # detects styles which have USER-AUTOPAS version
    RegisterStylesExt(${USER-AUTOPAS_SOURCES_DIR} autopas AUTOPAS_SOURCES)

    get_property(USER-AUTOPAS_SOURCES GLOBAL PROPERTY AUTOPAS_SOURCES)

    list(APPEND LIB_SOURCES ${USER-AUTOPAS_SOURCES})
    include_directories(${USER-AUTOPAS_SOURCES_DIR})
endif ()
