cmake_minimum_required(VERSION 3.20)

# define the project we are building
project(VisualiseTissueFolding VERSION 1.0)

# for this build, we are using C++14, set standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(OpenGL REQUIRED) # find OpenGL library...
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED) # which is different to the Qt5 OpenGL component here!

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

set(SIMULATION_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../TissueFolding/SourceCode/")

set(UI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GLWidget.cpp
)

set(SIMULATION_SOURCES
${SIMULATION_SOURCE_DIR}/Prism.cpp
${SIMULATION_SOURCE_DIR}/ReferenceShapeBase.cpp
${SIMULATION_SOURCE_DIR}/ShapeBase.cpp
${SIMULATION_SOURCE_DIR}/Simulation.cpp
${SIMULATION_SOURCE_DIR}/Node.cpp
${SIMULATION_SOURCE_DIR}/ModelInputObject.cpp
${SIMULATION_SOURCE_DIR}/RandomGenerator.cpp
${SIMULATION_SOURCE_DIR}/NewtonRaphsonSolver.cpp
${SIMULATION_SOURCE_DIR}/Analysis.cpp
${SIMULATION_SOURCE_DIR}/Lumen.cpp
${SIMULATION_SOURCE_DIR}/TimeSeriesPhysicalProperties.cpp
${SIMULATION_SOURCE_DIR}/YoungsModulusModifier.cpp
)


add_executable(VisualiseTissueFolding ${UI_SOURCES} ${SIMULATION_SOURCES})

target_include_directories(VisualiseTissueFolding PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../TissueFolding/SourceCode/" ${Boost_INCLUDE_DIRS})
target_link_libraries(VisualiseTissueFolding PRIVATE Qt5::OpenGL Qt5::Widgets ${OPENGL_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

