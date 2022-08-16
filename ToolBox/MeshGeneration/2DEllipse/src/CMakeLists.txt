cmake_minimum_required(VERSION 3.10)

# define the project we are building
project(MeshGeneration VERSION 1.0)

# we are using c++11, set standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file
# this will pass the location of triangle to the files
configure_file(EllipseFromOutline.hpp.in EllipseFromOutline.hpp)

# the executable we want to compile
add_executable(MeshGeneration EllipseFromOutline.cpp)

target_include_directories(MeshGeneration PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )