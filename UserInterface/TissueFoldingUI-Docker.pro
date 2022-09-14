CurrPath = /UserInterface

TARGET = $$CurrPath/Debug/TissueFoldingUI

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++17
QMAKE_LFLAGS +=  -fopenmp

# #for linux - include if you want to run without pardiso
# QMAKE_CXXFLAGS += -fopenmp -std=c++17 -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += opengl
 
# Input
HEADERS += 	$$CurrPath/SourceCode/MainWindow.h \
			$$CurrPath/SourceCode/GLWidget.h \
			$$CurrPath/../TissueFolding/SourceCode/*.h

SOURCES +=	$$CurrPath/SourceCode/main.cpp \
			$$CurrPath/SourceCode/MainWindow.cpp \
			$$CurrPath/SourceCode/GLWidget.cpp \
			$$CurrPath/../TissueFolding/SourceCode/Prism.cpp \
			$$CurrPath/../TissueFolding/SourceCode/ReferenceShapeBase.cpp \
 			$$CurrPath/../TissueFolding/SourceCode/ShapeBase.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/Simulation.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/Node.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/ModelInputObject.cpp \
			$$CurrPath/../TissueFolding/SourceCode/RandomGenerator.cpp \
			$$CurrPath/../TissueFolding/SourceCode/NewtonRaphsonSolver.cpp \
			$$CurrPath/../TissueFolding/SourceCode/Analysis.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/Lumen.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/TimeSeriesPhysicalProperties.cpp \
        	$$CurrPath/../TissueFolding/SourceCode/YoungsModulusModifier.cpp

LIBS += -L/usr/local/include -lgsl -lgslcblas -fopenmp -llapack -lgomp -lpthread -lgfortran -lm

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFoldingUI.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
