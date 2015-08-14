#Use the command line:
#qmake -o Makefile TissueFoldingUI.pro
#To generate the Makefile. Then you can build the project either with eclipse, or with
#the commands on the header of the Makefile

#curr path for ubuntu
CurrPath = /home/melda/Documents/TissueFolding/UserInterface/

#curr path for MacOS:
#CurrPath = ./

TARGET = $$CurrPath/Debug/TissueFoldingUI

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
	$$CurrPath/../TissueFolding/SourceCode/Triangle.cpp \
	$$CurrPath/../TissueFolding/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/../TissueFolding/SourceCode/ShapeBase.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Simulation.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Node.cpp \
        $$CurrPath/../TissueFolding/SourceCode/ModelInputObject.cpp \

#libs and includes for linux:
LIBS += -L/usr/include -lgsl -lgslcblas -lpardiso500-GNU461-X86-64  -fopenmp  -llapack

# libs and includes for MacOS
# LIBS += -L/usr/include -L/usr/local/lib/ -lgsl -lgslcblas -L/usr/local/Cellar/boost/1.58.0/include -lpardiso500-MACOS-X86-64
# INCLUDEPATH += /usr/local/Cellar/boost/1.58.0/include /usr/local/include/


# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFoldingUI.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
