cmake_minimum_required(VERSION 3.27)
project(dna_repair_plus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

# Find required packages
find_package(Protobuf REQUIRED)
find_package(OpenMP REQUIRED)

# Include directories for protobuf generated files
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROTOBUF_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


###############
### INDEXER ###
###############

# Source files for the make_index executable
set(MAKE_INDEX_SOURCES
        kmerization.cpp
        genome.cpp
        core/component/matrix.cpp
        core/component/index.cpp
        core/component/proto/matrix.pb.cc
        core/component/proto/index.pb.cc
)

# Create the make_index executable
add_executable(make_index ${MAKE_INDEX_SOURCES})
target_link_libraries(make_index PUBLIC ${PROTOBUF_LIBRARY} OpenMP::OpenMP_CXX)


###############
#### TEST #####
###############

# Google Test configuration
include(FetchContent)
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_subdirectory(tests)

