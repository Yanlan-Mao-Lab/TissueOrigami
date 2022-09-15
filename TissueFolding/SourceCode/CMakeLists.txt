cmake_minimum_required(VERSION 3.20)

# define the project we are building
project(TissueFolding VERSION 1.0)

# for this build, we are using C++14, set standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# populate the source files
set(TF_SOURCES 
    main.cpp
    Prism.cpp
    ReferenceShapeBase.cpp
    ShapeBase.cpp
    Simulation.cpp
    Node.cpp
    ModelInputObject.cpp
    RandomGenerator.cpp
    NewtonRaphsonSolver.cpp
    Analysis.cpp
    Lumen.cpp
    TimeSeriesPhysicalProperties.cpp
    YoungsModulusModifier.cpp
)

# if building with PARDISO, we will need to append the PARDISO library name and pthread too
set(TF_LINK_LIBS
    gsl
    gslcblas
    gomp
    lapack
)

# additional options to pass to the compiler
# openmp is needed regardless, we will append the PARDISO macro if we are building with PARDISO
set(TF_COMPILE_OPTS
    -fopenmp
)

# the executable we want to compile
add_executable(TissueFolding ${TF_SOURCES})

# we need to know if we are trying to build with Pardiso or not
# this command means the -DPARDISO option can be passed from the command line
# the default behaviour is to build _without_ PARDISO however
option(PARDISO "If ON, TissueFolding will be built with PARDISO. Default is OFF." OFF)
message(STATUS "Pardiso is: ${PARDISO}")

# if we want to build with PARDISO, we also should allow the user to pass "hints" to cmake on where to find PARDISO
# as such, we will allow the passing of the variables
# PLOC : The /path/to/the/pardiso.so file
# PNAME : The name of the .so file containing the PARDISO library
# The command-line PARDISO=ON option can still be passed, however without these hints we expect it to fail to locate the library

# if building with pardiso, locate the library using the information we've been given
if (PARDISO)
    message(STATUS "Pardiso location provided: ${PLOC}")
    message(STATUS "Pardiso library provided : ${PNAME}")

    # attempt to locate, place location into variable PARDISO_LIBRARY
    find_library(PARDISO_LIBRARY 
                NAMES "${PNAME}" 
                PATHS "${PLOC}"
                REQUIRED)

    # check that the PARDISO library was found, error if not
    if(NOT PARDISO_LIBRARY)
        message(FATAL_ERROR "Pardiso not found under the location/name specified")
    else()
        message(STATUS "Found PARDISO at ${PARDISO_LIBRARY}")
        # NOTE: GCC compilers require linking in the order in which libraries depend on each other
        # In the PARDISO docs, pardiso must be linked to before opnblas, gfortran, openmp, pthread, and m (math library)
        # Ergo, we must prepend rather than append here. [See https://pardiso-project.org/manual/manual.pdf, PG35]
        list(PREPEND TF_LINK_LIBS ${PARDISO_LIBRARY} gfortran pthread m)
    endif()

    # ensure the preprocessor macro for BUILD_WITH_PARDISO is defined,
    # by adding this to the compile options
    list(APPEND TF_COMPILE_OPTS "-DBUILD_WITH_PARDISO")

else()
    message(STATUS "Building without PARDISO")
endif()

# include directories are the same for both PARDISO and non-PARDISO builds
# since we are planning to build in an empty directory, no need to inlude PROJECT_BINARY_DIR here
target_include_directories(TissueFolding PUBLIC
                ${PROJECT_SOURCE_DIR}
)

# pass in additional compile options
message(STATUS "Passing in compiler options ${TF_COMPILE_OPTS}")
target_compile_options(TissueFolding PUBLIC
                ${TF_COMPILE_OPTS}
)

# tell the compiler to link libaries
message(STATUS "Linking with the following libraries: ${TF_LINK_LIBS}")
target_link_libraries(TissueFolding PUBLIC
                ${TF_LINK_LIBS}
)

#g++ -std=c++14 -lgsl -lgslcblas -lgomp -llapack -fopenmp -l/home/ccaegra/Documents/libpardiso600-GNU720-X86-64.so -lpthread -fopenmp -DBUILD_WITH_PARDISO -o TissueFolding -I/home/caegra/Documents/TissueOrigami/TissueFolding/SourceCode/ main.cpp Prism.cpp ReferenceShapeBase.cpp ShapeBase.cpp Simulation.cpp Node.cpp ModelInputObject.cpp RandomGenerator.cpp NewtonRaphsonSolver.cpp Analysis.cpp Lumen.cpp TimeSeriesPhysicalProperties.cpp YoungsModulusModifier.cpp 