# Optionally, you can set a project name:
project(BidirMultiCellsExample LANGUAGES CXX)

# All code for singlecell example
set(SINGLECELL_SOURCES
    cellFlow.cpp
)

# Create the executable
add_executable(bidirMultiCells ${SINGLECELL_SOURCES})

# Prepare a list of libraries to link against
set(LIBS
    palabos
    lammps
)

# Conditionally add Ascent to the list
if(ENABLE_ASCENT)
    list(APPEND LIBS insitu)
endif()

# If MPI is enabled, link it here as well
if(ENABLE_MPI)
    list(APPEND LIBS ${MPI_CXX_LIBRARIES})
endif()

# Link libraries from top-level
target_link_libraries(bidirMultiCells PRIVATE ${LIBS})

# Resource files (copied locally to the build folder)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/in.lmp4cell"
    "${CMAKE_CURRENT_BINARY_DIR}/in.lmp4cell"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/in.cells"
    "${CMAKE_CURRENT_BINARY_DIR}/in.cells"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/rbc.molecule"
    "${CMAKE_CURRENT_BINARY_DIR}/rbc.molecule"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ascent_actions.yaml"
    "${CMAKE_CURRENT_BINARY_DIR}/ascent_actions.yaml"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/trame.yaml"
    "${CMAKE_CURRENT_BINARY_DIR}/trame.yaml"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cinema.yaml"
    "${CMAKE_CURRENT_BINARY_DIR}/cinema.yaml"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ascent/ascent_trame_bridge.py"
    "${CMAKE_CURRENT_BINARY_DIR}/ascent/ascent_trame_bridge.py"
    COPYONLY
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/trame/trame_app.py"
    "${CMAKE_CURRENT_BINARY_DIR}/trame/trame_app.py"
    COPYONLY
)
