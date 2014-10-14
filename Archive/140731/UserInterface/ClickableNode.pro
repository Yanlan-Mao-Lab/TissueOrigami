CurrPath = /home/melda/Documents/TissueFolding/UserInterface/

TARGET = $$CurrPath/Debug/UserInterface

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += opengl
 
# Input


HEADERS       +=  	$$CurrPath/SourceCode/MainWindow.h \
			$$CurrPath/SourceCode/GLWidget.h \
			$$CurrPath/SourceCode/TriangleFromModel.h
			$$CurrPath/../TissueFolding/SourceCode/*.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
	$$CurrPath/SourceCode/MainWindow.cpp \
	$$CurrPath/SourceCode/GLWidget.cpp \
	$$CurrPath/SourceCode/TriangleFromModel.cpp \
	$$CurrPath/../TissueFolding/SourceCode/Prism.cpp \
	$$CurrPath/../TissueFolding/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/../TissueFolding/SourceCode/ShapeBase.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Simulation.cpp \
        $$CurrPath/../TissueFolding/SourceCode/Node.cpp \

# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES ClickableNode.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
