/*
 * man.cc
 *
 *  Created on: 18 Mar 2014
 *      Author: melda
 */

#include <QGraphicsWidget>
#include <QtGui>
#include <QtWidgets>
#include "MainWindow.h"
#include "GLWidget.h"

#include "../TissueFolding/SourceCode/Simulation.h"
#include <vector>

class MainWindow;
class SurfaceBase;
class GLWidget;

Simulation* Sim01;
int main(int argc, char **argv)
{
	bool Success = false;
	Sim01 = new Simulation();
	Sim01->displayIsOn = true;
	if (argc<2){
		Sim01->DisplaySave = false;
		cerr<<"Using default settings"<<endl;
		Success = true;
	}
	else{
		Success = Sim01->readExecutableInputs(argc, argv);
	}
	if (Success == 0 ){
		cout<<"Error in input to executable"<<endl;
		return true;
	}

	QApplication app(argc, argv);

	if (Sim01->DisplaySave){
		cout<<"Initiating simulation display"<<endl;
		Success = Sim01->initiateSavedSystem();
	}
	else{
		Success = Sim01->initiateSystem();
		cout<<"system initiated"<<endl;
		for (int i=0; i<Sim01->Elements.size(); ++i){
			//This is the initial setup, the elements should take the actual positions of the nodes, this corresponds to RK step 4, RKId= 3
			Sim01->Elements[i]->updatePositions(3,Sim01->Nodes);
		}
	}

	if (Success == 0 ){
		cout<<"System is not initiated successfully, terminating"<<endl;
		return true;
	}

	MainWindow mw(Sim01);
	mw.show();
	mw.MainGLWidget->show();
	mw.raise();
	mw.setGeometry(50, 50, 1800, 1100);
	mw.MainScene->update();

	return app.exec();
}




