#Use the command line:
#qmake -o Makefile TissueFoldingUI.pro
#To generate the Makefile. Then you can build the project either with eclipse, or with
#the commands on the header of the Makefile


#the line necessary for the compilation of non-visual interface model is below. This should be in the makefile inside TissueFolding/Debug/:
#g++ -L/usr/include/  -L/opt/Qt5.2.1/5.2.1/gcc_64/lib  -o "TissueFolding" $(OBJS) $(USER_OBJS) $(LIBS) $(INCLUDEPATH) -lgsl  -lgslcblas -lpardiso500-GNU461-X86-64 -fopenmp -llapack -lpthread
	
#FOR LEGION, when you add a new source/header file to the sourcecode, you do need to update the *.mk files. 
#These include the files under "/home/melda/Documents/TissueFolding/TissueFolding/Debug/" and 
#                              "/home/melda/Documents/TissueFolding/TissueFolding/Debug/SourceCode"
#The procedure is, you need to build the latest version (updates the *.mk files). Then go to the makefile of the non-visual version
# under  "/home/melda/Documents/TissueFolding/TissueFolding/Debug/". Change the compiler line with the line given below. Make the latest non-visual version.
# These steps will ensure the the subdirectory repositories are updated. Then copy all (as listed above) to Legion. The make file of Legion 
# is not the same as the makfile of Tethys. The compiler line is:
#   	g++  -o "TissueFolding" $(OBJS) $(USER_OBJS) $(LIBS) -L${OPENBLASROOT}/lib -lopenblas -lpardiso500-GNU481-X86-64  -fopenmp  -I/shared/ucl/apps/gsl/1.16/gcc/include  -lgsl -lgslcblas
# Then you can make on legion.


#curr path for ubuntu
#CurrPath = /home/melda/Documents/TissueFolding/UserInterface/

#curr path for MacOS:
CurrPath = ./

TARGET = $$CurrPath/Debug/TissueFoldingUI

#qt creator on mac can read config not QMAKE_CXXFLAGS, ubuntu can read QMAKE_CXXFLAGS, not CONFIG

#for ubuntu
QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++17
QMAKE_LFLAGS +=  -fopenmp

QMAKE_CXXFLAGS += -fopenmp -std=c++17 -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS

#for mac
#CONFIG += c++17 -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS
#QMAKE_CXXFLAGS += -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += opengl
 
# Input


HEADERS       +=  	$$CurrPath/SourceCode/MainWindow.h \
			$$CurrPath/SourceCode/GLWidget.h \
			$$CurrPath/../TissueFolding/SourceCode/*.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
	$$CurrPath/SourceCode/MainWindow.cpp \
	$$CurrPath/SourceCode/GLWidget.cpp \
	$$CurrPath/../TissueFolding/SourceCode/Prism.cpp \
	$$CurrPath/../TissueFolding/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/../TissueFolding/SourceCode/ShapeBase.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Simulation.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Node.cpp \
        $$CurrPath/../TissueFolding/SourceCode/ModelInputObject.cpp \
	$$CurrPath/../TissueFolding/SourceCode/RandomGenerator.cpp \
#	$$CurrPath/../TissueFolding/SourceCode/NewtonRaphsonSolver.cpp \
	$$CurrPath/../TissueFolding/SourceCode/Analysis.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Lumen.cpp \
        $$CurrPath/../TissueFolding/SourceCode/TimeSeriesPhysicalProperties.cpp \
        $$CurrPath/../TissueFolding/SourceCode/YoungsModulusModifier.cpp

#libs and includes for linux for independent license pardiso:
#LIBS += -L/usr/include -lgsl -lgslcblas -lpardiso600-GNU720-X86-64  -fopenmp -llapack -lgomp -lpthread -lgfortran -lm
LIBS += -L/usr/local/include -lgsl -lgslcblas -fopenmp -llapack -lgomp -lpthread -lgfortran -lm

# libs and includes for MacOS Intel processor
#LIBS += -L/usr/include -L/usr/local/lib/ -lgsl -lgslcblas -L/usr/local/Cellar/boost/1.58.0/include -lpardiso600-MACOS-X86-64
#LIBS += -L/usr/include -L/usr/local/lib/ -lgsl -lgslcblas -L/usr/local/Cellar/boost/1.58.0/include
#INCLUDEPATH += /usr/local/Cellar/boost/1.71.0/include /usr/local/include/

# libs and includes for MacOS M1 Chip processor
#LIBS += -L/opt/include -L/opt/homebrew/lib/ -lgsl -lgslcblas -L/opt/homebrew/Cellar/boost/1.78.0_1/include
#INCLUDEPATH += /opt/homebrew/Cellar/boost/1.78.0_1/include /opt/homebrew/include/

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFoldingUI.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
