#curr path for ubuntu
CurrPath = /TissueFolding/

TARGET = $$CurrPath/Debug/TissueFolding

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++14
QMAKE_LFLAGS +=  -fopenmp
 
# Input
OBJECTS_DIR += $$CurrPath/Debug/

#HEADERS += $$CurrPath/SourceCode/*.h
HEADERS += $$CurrPath/SourceCode/Prism.h \
	$$CurrPath/SourceCode/ReferenceShapeBase.h \
 	$$CurrPath/SourceCode/ShapeBase.h \
    $$CurrPath/SourceCode/Simulation.h \
    $$CurrPath/SourceCode/Node.h \
    $$CurrPath/SourceCode/ModelInputObject.h \
	$$CurrPath/SourceCode/RandomGenerator.h \
	$$CurrPath/SourceCode/NewtonRaphsonSolver.h \
	$$CurrPath/SourceCode/Analysis.h \
	$$CurrPath/SourceCode/Lumen.h \
	$$CurrPath/SourceCode/YoungsModulusModifier.h \
	$$CurrPath/SourceCode/TimeSeriesPhysicalProperties.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
	$$CurrPath/SourceCode/Prism.cpp \
	$$CurrPath/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/SourceCode/ShapeBase.cpp \
    $$CurrPath/SourceCode/Simulation.cpp \
    $$CurrPath/SourceCode/Node.cpp \
    $$CurrPath/SourceCode/ModelInputObject.cpp \
	$$CurrPath/SourceCode/RandomGenerator.cpp \
	$$CurrPath/SourceCode/NewtonRaphsonSolver.cpp \
	$$CurrPath/SourceCode/Analysis.cpp \
	$$CurrPath/SourceCode/Lumen.cpp \
	$$CurrPath/SourceCode/YoungsModulusModifier.cpp \
	$$CurrPath/SourceCode/TimeSeriesPhysicalProperties.cpp

LIBS += -lgsl -lgslcblas -lgomp -fopenmp -llapack 

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFolding-Docker.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
