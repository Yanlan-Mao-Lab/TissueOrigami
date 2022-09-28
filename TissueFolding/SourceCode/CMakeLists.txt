cmake_minimum_required(VERSION 3.20)

# define the project we are building
project(TissueFolding VERSION 1.0)

# for this build, we are using C++14, set standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# we need to know if we are trying to build with Pardiso or not
# this command means the -DPARDISO option can be passed from the command line
# the default behaviour is to build _without_ PARDISO
option(PARDISO "If ON, TissueFolding will be built with PARDISO. Default is OFF." OFF)
# if we want to build with PARDISO, we also should allow the user to pass "hints" to cmake on where to find PARDISO
# as such, we will allow the passing of the variables
# PLOC : The /path/to/the/pardiso.so file
# PNAME : The name of the .so file containing the PARDISO library
# The command-line PARDISO=ON option can still be passed, however without these hints we expect it to fail to locate the library
message(STATUS "Pardiso is: ${PARDISO}")

# find libraries we always need
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# populate the source files
set(SIMULATION_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Prism.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceShapeBase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ShapeBase.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Simulation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Node.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ModelInputObject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RandomGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NewtonRaphsonSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Analysis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Lumen.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeSeriesPhysicalProperties.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/YoungsModulusModifier.cpp
)

# the executable we want to compile
add_executable(TissueFolding ${SIMULATION_SOURCES})

# include directories are the same for both PARDISO and non-PARDISO builds
target_include_directories(TissueFolding PUBLIC
                ${PROJECT_SOURCE_DIR}
                ${Boost_INCLUDE_DIRS}
)

# building with PARDISO requires additional linked libraries
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
        # [See https://pardiso-project.org/manual/manual.pdf, PG35]

        # find the additional libraries that we need
        # (also need -pthread option)
        set(THREADS_PREFER_PTHREAD_FLAG ON)
        find_package(Threads REQUIRED)
    endif()

    # ensure the preprocessor macro for BUILD_WITH_PARDISO is defined,
    # by adding this to the compile options
    target_compile_options(TissueFolding PUBLIC "-DBUILD_WITH_PARDISO")

    # add the required libraries (including PARDISO)
    target_link_libraries(TissueFolding PUBLIC
                ${PARDISO_LIBRARY}          # PARDISO (must be 1st))
                ${Boost_LIBRARIES}          # Boost
                ${GSL_LIBRARIES}            # GSL
                ${BLAS_LIBRARIES}           # BLAS
                ${LAPACK_LIBRARIES}         # LAPACK (before gfortran)
                ${OpenMP_CXX_LIBRARIES}     # OpenMP (after PARDISO)
                gfortran                    # gfortran (after PARDISO)
                Threads::Threads            # pthread (after PARDISO)
                m                           # math.h (after PARDISO)
                )
else()
    message(STATUS "Building without PARDISO")

    # add the required libraries (including PARDISO)
    target_link_libraries(TissueFolding PUBLIC
                ${Boost_LIBRARIES}          # Boost
                ${GSL_LIBRARIES}            # GSL
                ${BLAS_LIBRARIES}           # BLAS
                ${LAPACK_LIBRARIES}         # LAPACK
                ${OpenMP_CXX_LIBRARIES}     # OpenMP
                )
endif()
