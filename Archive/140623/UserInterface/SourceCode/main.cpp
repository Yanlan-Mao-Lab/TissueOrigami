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
	QApplication app(argc, argv);

	bool DisplayOnTheGo = true;
	Sim01 = new Simulation();
	if (DisplayOnTheGo){
		Sim01->initiateNodes();
		Sim01->initiateElements();
		Sim01->initiateSystemForces();
		Sim01->calculateSystemCentre();
		cout<<"finalised initiation"<<endl;
	}

	MainWindow mw(Sim01);

	mw.show();
	mw.MainGLWidget->show();
	mw.raise();
	mw.setGeometry(100, 100, 1200, 900);
	mw.MainScene->update();


	return app.exec();
}




