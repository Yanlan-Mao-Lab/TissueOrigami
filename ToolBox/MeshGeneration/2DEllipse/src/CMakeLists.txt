cmake_minimum_required(VERSION 3.10)

# define the project we are building
project(MeshGeneration VERSION 1.0)

# we are using c++11, set standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set triangle path - change this later to auto-detect triangle's location, or make it a subroutine
# but for now, we can set this variable here, could even try to get auto-detection working?
#set(TRIANGLE_PATH "../triangle/triangle")
find_program(TRIANGLE_PATH
            triangle triangle.o
            "../../triangle") # this causes a link to the local version of triangle - should make triangle a prequiste (callable from the command line) and drop this extra argument

# configure a header file
# this will pass the location of triangle to the files
configure_file(EllipseFromOutline.hpp.in EllipseFromOutline.hpp)

# the executable we want to compile
add_executable(MeshGeneration EllipseFromOutline.cpp)

target_include_directories(MeshGeneration PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )