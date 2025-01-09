cmake_minimum_required(VERSION 3.18)
project(BloodFlow LANGUAGES CXX)

###
# BASIC SETTINGS
###
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "CMake Generator:            ${CMAKE_GENERATOR}")
message(STATUS "CMake Build Types:          ${CMAKE_CONFIGURATION_TYPES}")
message(STATUS "Current Build Type:         ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###
# EXTERNAL PROJECT SETUP
###
include(ExternalProject)
include(FetchContent)

option(FETCH_PALABOS "Download and manage Palabos via FetchContent" ON)
option(FETCH_LAMMPS  "Download and manage LAMMPS via ExternalProject" ON)
option(BUILD_EXAMPLE_SINGLECELL  "Build SingleCell Example" ON)
option(ENABLE_ASCENT  "Ascent ENABLED" ON)

# ----------------------------------------------------------------------------
# Palabos: fetch only, do not build
# ----------------------------------------------------------------------------
if(FETCH_PALABOS)
    FetchContent_Declare(
        palabos
        GIT_REPOSITORY https://gitlab.com/unigespc/palabos.git
        GIT_TAG        e498e8ad7f24fd7ff87313670db7873703c1fd3f
    )
    FetchContent_Populate(palabos)

    if(NOT palabos_POPULATED)
        message(FATAL_ERROR "Failed to fetch Palabos")
    endif()

    set(PALABOS_ROOT ${palabos_SOURCE_DIR} CACHE PATH
        "Path to automatically fetched Palabos source directory" FORCE)
else()
    set(PALABOS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../palabos" CACHE PATH
        "Path to manually installed Palabos directory")
endif()

# ----------------------------------------------------------------------------
# LAMMPS: fetch and build with 'make yes-MOLECULE' then 'make mpi mode=lib'
# ----------------------------------------------------------------------------
if(FETCH_LAMMPS)
    ExternalProject_Add(
        LAMMPSExt
        GIT_REPOSITORY  https://github.com/lammps/lammps.git
        GIT_TAG         cd6ab17d36523b08c6a8c9a1a224a4f2b70d8a46
        UPDATE_COMMAND  ""
        PATCH_COMMAND
            ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/rbc"       # Adjust path if needed
                "${CMAKE_BINARY_DIR}/lammps-src/src"
        CONFIGURE_COMMAND ""  # Not using CMake for LAMMPS config
        BUILD_COMMAND
            ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_BINARY_DIR}/lammps-src/src yes-MOLECULE
            && ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_BINARY_DIR}/lammps-src/src mpi mode=lib
        INSTALL_COMMAND ""
        SOURCE_DIR       "${CMAKE_BINARY_DIR}/lammps-src"
        BINARY_DIR       "${CMAKE_BINARY_DIR}/lammps-build"
    )
    set(LAMMPS_DIR "${CMAKE_BINARY_DIR}/lammps-src" CACHE PATH
        "Path to automatically fetched LAMMPS source directory" FORCE)
else()
    set(LAMMPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../lammps" CACHE PATH
        "Path to manually installed LAMMPS directory")
endif()

###
# FIND ASCENT (REQUIRED)
###
if(ENABLE_ASCENT)
  find_package(Ascent REQUIRED)
  if(ASCENT_FOUND)
      message(STATUS "FOUND ASCENT")
  endif()
endif()

###
# OTHER ENVIRONMENT PATHS
###
# For example, if BloodFlow root is the top-level dir, you can keep this:
set(BLOODFLOW_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH
    "Path to BloodFlow root directory")

# Attempt to locate the coupling directory if needed
find_path(COUPLING_DIR palabos-lammps)

# LAMMPS library path
set(LAMMPS_LIB "${LAMMPS_DIR}/src")

###
# COMPILER FLAGS
###
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "Detected compiler: ${CMAKE_CXX_COMPILER_ID}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wnon-virtual-dtor -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -DPLB_DEBUG -O0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Detected compiler: MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG   "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message(FATAL_ERROR "CXX compiler not recognized: ${CMAKE_CXX_COMPILER_ID}")
endif()

###
# MPI SUPPORT
###
option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message(STATUS "Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

###
# POSIX SUPPORT
###
if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message(STATUS "Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif()

###
# INCLUDE DIRECTORIES
###
# Palabos (headers)
include_directories(palabos PUBLIC "${PALABOS_ROOT}/src")
include_directories(palabos PUBLIC "${PALABOS_ROOT}/externalLibraries")
include_directories(palabos PUBLIC "${PALABOS_ROOT}/externalLibraries/Eigen3")

# IBM / RBC code
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/ibm")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/rbc")

# LAMMPS
include_directories("${LAMMPS_LIB}")
include_directories("${COUPLING_DIR}/src")

# Ascent
if(ENABLE_ASCENT)
  include_directories(
      ${ASCENT_INSTALL_PREFIX}/include/ascent
      ${CONDUIT_DIR}/include/conduit
  )
  add_definitions(-DASCENT_ENABLED)
endif()

###
# CREATE PALABOS STATIC LIBRARY
# (If you need to compile Palabos code from .cpp files)
###
file(GLOB_RECURSE PALABOS_SRC
    "${PALABOS_ROOT}/src/*.cpp"
)
file(GLOB_RECURSE EXT_SRC
    "${PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp"
)
add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###
# CREATE IMPORTED LAMMPS LIBRARY
###
add_library(lammps STATIC IMPORTED)
add_dependencies(lammps LAMMPSExt)
set_property(TARGET lammps PROPERTY
    IMPORTED_LOCATION "${LAMMPS_LIB}/liblammps_mpi.a"
)

###
# ADD SUBDIRECTORIES FOR EXAMPLES
###
if(BUILD_EXAMPLE_SINGLECELL)
  add_subdirectory(examples/singleCell)
endif()

message(STATUS "Top-level CMake configuration complete.")
