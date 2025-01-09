cmake_minimum_required(VERSION 3.18)
project(palabosExample LANGUAGES CXX)

###
# BASIC SETTINGS
###
# Default to Release if no build type is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "CMake Generator:            ${CMAKE_GENERATOR}")
message(STATUS "CMake Build Types:          ${CMAKE_CONFIGURATION_TYPES}")
message(STATUS "Current Build Type:         ${CMAKE_BUILD_TYPE}")

# Use at least C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###
# EXTERNAL PROJECT SETUP
###
include(ExternalProject)

option(FETCH_PALABOS "Download and manage Palabos via ExternalProject" ON)
option(FETCH_LAMMPS  "Download and manage LAMMPS via ExternalProject" ON)

# ----------------------------------------------------------------------------
# Palabos: fetch only, do not build
# ----------------------------------------------------------------------------
if(FETCH_PALABOS)
    ExternalProject_Add(
        PalabosExt
        GIT_REPOSITORY  https://gitlab.com/unigespc/palabos.git
        GIT_TAG         e498e8ad7f24fd7ff87313670db7873703c1fd3f
        UPDATE_COMMAND  ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND  ""
        # Source/Build directories
        SOURCE_DIR       "${CMAKE_BINARY_DIR}/palabos-src"
        BINARY_DIR       "${CMAKE_BINARY_DIR}/palabos-build"
    )
    set(PALABOS_ROOT "${CMAKE_BINARY_DIR}/palabos-src" CACHE PATH
        "Path to automatically fetched Palabos source directory" FORCE)
endif()

# ----------------------------------------------------------------------------
# LAMMPS: fetch and build with 'make yes-MOLECULE' then 'make mpi mode=lib -j16'
# ----------------------------------------------------------------------------
if(FETCH_LAMMPS)
    ExternalProject_Add(
        LAMMPSExt
        GIT_REPOSITORY  https://github.com/lammps/lammps.git
        GIT_TAG         cd6ab17d36523b08c6a8c9a1a224a4f2b70d8a46
        UPDATE_COMMAND  ""
        PATCH_COMMAND 
        ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/../../rbc"
            "${CMAKE_BINARY_DIR}/lammps-src/src"
        CONFIGURE_COMMAND ""  # Not using CMake for LAMMPS config in this example
        # Build in the cloned src directory:
        BUILD_COMMAND
            ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_BINARY_DIR}/lammps-src/src yes-MOLECULE
            && ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_BINARY_DIR}/lammps-src/src mpi mode=lib -j16
        INSTALL_COMMAND  ""
        SOURCE_DIR       "${CMAKE_BINARY_DIR}/lammps-src"
        BINARY_DIR       "${CMAKE_BINARY_DIR}/lammps-build"
    )

    # We expect the library to appear under lammps-src/src after build
    set(LAMMPS_DIR "${CMAKE_BINARY_DIR}/lammps-src" CACHE PATH
        "Path to automatically fetched LAMMPS source directory" FORCE)
endif()

###
# FIND ASCENT (Required)
###
find_package(Ascent REQUIRED)
if(ASCENT_FOUND)
    message(STATUS "FOUND ASCENT")
endif()

###
# OTHER ENVIRONMENT PATHS
###
set(BLOODFLOW_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH
    "Path to BloodFlow root directory")

# Locate the coupling directory if it exists
find_path(COUPLING_DIR palabos-lammps)

# LAMMPS library path
set(LAMMPS_LIB "${LAMMPS_DIR}/src")

# Executable name
set(EXECUTABLE_NAME "cellFlow")

###
# COMPILER FLAGS
###
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR 
    CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR 
    CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "Detected compiler: ${CMAKE_CXX_COMPILER_ID}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wnon-virtual-dtor -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Detected compiler: MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message(FATAL_ERROR "CXX compiler not recognized: ${CMAKE_CXX_COMPILER_ID}. CMake will quit.")
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
# Palabos
include_directories("${PALABOS_ROOT}/src")
include_directories("${PALABOS_ROOT}/externalLibraries")
include_directories("${PALABOS_ROOT}/externalLibraries/Eigen3")

# BloodFlow
include_directories("${BLOODFLOW_ROOT}/ibm")

# LAMMPS
include_directories("${LAMMPS_LIB}")
include_directories("${COUPLING_DIR}/src")

# Ascent
include_directories(
    ${ASCENT_INSTALL_PREFIX}/include/ascent
    ${CONDUIT_DIR}/include/conduit
)

###
# PALABOS LIBRARY (not built, just a static lib from .cpp files)
###
file(GLOB_RECURSE PALABOS_SRC
    "${PALABOS_ROOT}/src/*.cpp"
)
file(GLOB_RECURSE EXT_SRC
    "${PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp"
)

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###
# LAMMPS LIBRARY (IMPORTED)
###
# After 'make mpi mode=lib', we expect liblammps_mpi.a in ${LAMMPS_LIB}.
add_library(lammps STATIC IMPORTED)
set_property(TARGET lammps PROPERTY
    IMPORTED_LOCATION "${LAMMPS_LIB}/liblammps_mpi.a"
)

###
# MAIN EXECUTABLE
###
set(SRC
    ${EXECUTABLE_NAME}.cpp
    Bridge.cpp
)

add_executable(${EXECUTABLE_NAME} ${SRC})

target_link_libraries(${EXECUTABLE_NAME}
    palabos
    lammps
    ascent::ascent_mpi
)

if(ENABLE_MPI)
    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()

###
# COPY RESOURCE FILES
###
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/in.lmp4cell"
    "${CMAKE_CURRENT_BINARY_DIR}/in.lmp4cell"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/in.tencells"
    "${CMAKE_CURRENT_BINARY_DIR}/in.tencells"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ascent_actions.yaml"
    "${CMAKE_CURRENT_BINARY_DIR}/ascent_actions.yaml"
    COPYONLY
)

message(STATUS "CMake configuration complete.")
