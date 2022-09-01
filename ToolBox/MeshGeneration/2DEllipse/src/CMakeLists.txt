cmake_minimum_required(VERSION 3.10)

# define the project we are building
project(EllipseFromOutline VERSION 1.0)

# we are using c++11, set standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set triangle path - change this later to auto-detect triangle's location, or make it a subroutine
# but for now, we can set this variable here, could even try to get auto-detection working?
find_program(TRIANGLE_PATH
            triangle triangle.o
            "../../triangle") # add optional search path for triangle here

# configure a header file
# this will pass the location of triangle to the files
configure_file(EllipseFromOutline.hpp.in EllipseFromOutline.hpp)

# the executable we want to compile
add_executable(EllipseFromOutline 
    EllipseFromOutline.cpp
    ArgumentParser.hpp
    ArgumentParser.cpp
    )

target_include_directories(EllipseFromOutline PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )