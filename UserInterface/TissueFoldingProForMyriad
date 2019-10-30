#Use the command line:
# For the server, oyu need the correct modules, 
# ignore the module lines on ubuntu. 
# On mac, qt creator should help you out anyway and you do not use terminal for this stage.

# -- Install the qt module to use correct qmake first (below s for QT 4, QT 5 is different if you ever need it, :
# module load qt/4.8.6/gnu-4.9.2
# -- To generate the Makefile, do a qmake:
# qmake TissueFolding.pro 
# -- Then you will need to load all the modules that are needed to compile and run the simulation. These include all
#    the libraries (such as boost or gsl). You can find them in Job submitter scripts. Because you would need to load
#    these at each job initiation as well, therefore they are automated in Melda's scripts. Look for 
#    ../Run00001/Job00001 first, if you dont have any ofg these, check inside /Scriptwriters/Writer_JobSubmission.cpp
# module unload compilers mpi mkl
# module load compilers/gnu/4.9.2
# module load gsl/1.16/gnu-4.9.2 
# module load openblas/0.2.14/gnu-4.9.2
# module load python/2.7.9
# module load boost/1_54_0/gnu-4.9.2
# -- You also need the export library paths for pardiso. Find them exactly where the modules are, as mentioned above.
# export PATH=/home/ucbpnkh/LocalLibs/:$PATH
# export LD_LIBRARY_PATH=/home/ucbpnkh/LocalLibs/:$LD_LIBRARY_PATH
# export LIBRARY_PATH=/home/ucbpnkh/LocalLibs/:$LIBRARY_PATH
# export OMP_NUM_THREADS=4
# export OPENBLAS_NUM_THREADS=1
# -- Now make:
# make
# -- Copy the executable under Debug (/Debug/TissueFolding) to the bin.
# scp -r ../TissueFolding/Debug/TissueFolding ./   
#(this will work is you are in the bin folder - change the paths if you are not).
	
#FOR LEGION, when you add a new source/header file to the sourcecode, you do need to update the *.mk files. 
#These include the files under "/home/melda/Documents/TissueFolding/TissueFolding/Debug/" and check legion page
# (https://wiki.rc.ucl.ac.uk/wiki/Development_Tools)
#                              "/home/melda/Documents/TissueFolding/TissueFolding/Debug/SourceCode"
#Install the qt module to use correct qmake first (below s for QT 4, QT 5 is different if you ever need it, :
#module load qt/4.8.6/gnu-4.9.2
#then qmake: 
#qmake -o Makefile TissueFoldingUI.pro
#The procedure is, you need to build the latest version (updates the *.mk files). Then go to the makefile of the non-visual version
# under  "/home/melda/Documents/TissueFolding/TissueFolding/Debug/". Change the compiler line with the line given below. Make the latest non-visual version.
# These steps will ensure the the subdirectory repositories are updated. Then copy all (as listed above) to Legion. The make file of Legion 
# is not the same as the makfile of Tethys. The compiler line is:
#   	g++  -o "TissueFolding" $(OBJS) $(USER_OBJS) $(LIBS) -L${OPENBLASROOT}/lib -lopenblas -lpardiso500-GNU481-X86-64  -fopenmp  -I/shared/ucl/apps/gsl/1.16/gcc/include  -lgsl -lgslcblas
# Then you can make on legion.

#curr path for ubuntu
#CurrPath = /home/melda/Documents/TissueFolding/TissueFolding/

#curr path for MacOS:
#CurrPath = ./

#CurrPath for Legion/Myriad:
CurrPath = /home/ucbpnkh/Scratch/TissueFolding/

TARGET = $$CurrPath/Debug/TissueFolding

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++14
QMAKE_LFLAGS +=  -fopenmp
#QMAKE_CXXFLAGS += -fopenmp -std=c++14 -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS
 
# Input
OBJECTS_DIR += $$CurrPath/Debug/

HEADERS       +=  	$$CurrPath/SourceCode/*.h

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
	$$CurrPath/SourceCode/TimeSeriesPhysicalProperties.cpp \
	$$CurrPath/SourceCode/YoungsModulusModifier.cpp

LIBS += -L/usr/include -lgsl -lgslcblas -lgomp -fopenmp -lpardiso600-GNU720-X86-64   -llapack 

# libs and includes for MacOS
# LIBS += -L/usr/include -L/usr/local/lib/ -lgsl -lgslcblas -L/usr/local/Cellar/boost/1.58.0/include -lpardiso500-MACOS-X86-64
# INCLUDEPATH += /usr/local/Cellar/boost/1.58.0/include /usr/local/include/
#CONFIG += c++17 -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS
#QMAKE_CXXFLAGS += -D DO_NOT_USE_OMP -D DO_NOT_SOLVE_SYSTEM_OF_EQUATIONS

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFolding.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
