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
	bool Success = false;
	Sim01 = new Simulation();
	if (DisplayOnTheGo){
		Success = Sim01->initiateSystem();
	}
	if (Success == 0 ){
		cout<<"System is not initiated successfully, terminating"<<endl;
		return true;
	}
	MainWindow mw(Sim01);

	mw.show();
	mw.MainGLWidget->show();
	mw.raise();
	mw.setGeometry(100, 100, 1200, 900);
	mw.MainScene->update();

	return app.exec();
}




