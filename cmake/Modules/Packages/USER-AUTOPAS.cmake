option(PKG_USER-AUTOPAS "Enable AutoPas support" ON)

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
    set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
    if (GIT_SUBMODULES_SSH)
        set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
    endif ()

    # Final version of 2 Body AutoPas
    set(AUTOPAS_TAG v2.0.0 CACHE STRING "AutoPas Git tag or commit id to use")

    # Download and install autopas
    FetchContent_Declare(
            autopas
            GIT_REPOSITORY ${autopasRepoPath}
            GIT_TAG ${AUTOPAS_TAG}
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
    # Link detached library which contains the particle properties library
    list(APPEND LAMMPS_LINK_LIBS "molecularDynamicsLibrary")


    ######### 2. AutoPas package settings #########

    set(USER-AUTOPAS_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-AUTOPAS)

    add_definitions(-DLMP_AUTOPAS)

    # Sources that are not a style
    set(USER-AUTOPAS_SOURCES
            ${USER-AUTOPAS_SOURCES_DIR}/autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/atom_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/atom_vec_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/comm_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/domain_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/modify_autopas.cpp
            #${USER-AUTOPAS_SOURCES_DIR}/nbin_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/neighbor_autopas.cpp
            ${USER-AUTOPAS_SOURCES_DIR}/output_autopas.cpp
            )

    set_property(GLOBAL PROPERTY "AUTOPAS_SOURCES" "${USER-AUTOPAS_SOURCES}")

    # detects styles which have USER-AUTOPAS version
    RegisterStylesExt(${USER-AUTOPAS_SOURCES_DIR} autopas AUTOPAS_SOURCES)

    # register autopas-only styles
    # RegisterNBinStyle(${USER-AUTOPAS_SOURCES_DIR}/nbin_autopas.h)

    get_property(USER-AUTOPAS_SOURCES GLOBAL PROPERTY AUTOPAS_SOURCES)

    list(APPEND LIB_SOURCES ${USER-AUTOPAS_SOURCES})
    include_directories(${USER-AUTOPAS_SOURCES_DIR})
endif ()
