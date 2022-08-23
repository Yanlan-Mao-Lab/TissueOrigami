cmake_minimum_required(VERSION 3.10)

# define the project we are building
project(TissueFolding VERSION 1.0)

# for this build, we are using C++14, set standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# need to find pardiso (if compiling with Pardiso): https://pardiso-project.org/manual/manual.pdf

# FOR RUNNING THE CODE, the following environment variables must be set!
# LD_LIBRARY_PATH must point to the Pardiso library (.so file)
# OMP_NUM_THREADS (if using multi-threaded version)

# FOR LINKING, THE FOLLOWING COMMANDS ARE USED ON LINUX ARCHITECTURES!
# gcc <source/objects files> -o <executable>
# -L <Path to directory of PARDISO> -lpardiso700-GNU461-X86-64
# -L <Path to directory of LAPACK/BLAS>
# -l <Fast LAPACK and BLAS libraries> -lgfortran -fopenmp -lpthread -lm

# WINDOWS LINKING IS DIFFERENT STILL!

# the executable we want to compile
add_executable(TissueFolding
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

target_include_directories(TissueFolding PUBLIC
                "${PROJECT_SOURCE_DIR}"
)

target_compile_options(TissueFolding PUBLIC
                -fopenmp
)

target_link_libraries(TissueFolding PUBLIC
                gsl
                gslcblas
                gomp
                lapack
)

target_include_directories(TissueFolding PUBLIC
                            "${PROJECT_BINARY_DIR}"
)