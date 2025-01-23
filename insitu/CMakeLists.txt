# insitu/CMakeLists.txt
cmake_minimum_required(VERSION 3.18)

project(insitu LANGUAGES CXX)

# Collect all your .cpp or .cxx files
file(GLOB INSITU_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    # Add more subdirectories if needed
)

# Create a static (or shared) library
add_library(insitu STATIC ${INSITU_SOURCES})

# Include paths so that code linking to this library can see its headers
target_include_directories(insitu
    PUBLIC
        # Tells dependent targets to look for headers in insitu/ directly
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        ${ASCENT_INSTALL_PREFIX}/include/ascent
)

# If your library depends on Ascent or other libraries, link them here

target_compile_definitions(insitu PUBLIC -DENABLE_ASCENT)
target_link_libraries(insitu PUBLIC ascent::ascent_mpi) 

# Or link any other dependencies as needed
# target_link_libraries(insitu PUBLIC someOtherLib)
