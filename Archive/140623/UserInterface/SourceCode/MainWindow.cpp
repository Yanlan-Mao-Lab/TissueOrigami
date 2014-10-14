/*
 * MainWindow.cpp
 *
 *  Created on: 18 Mar 2014
 *      Author: melda
 */
#include <iostream>
#include <QtGui>
#include <QtWidgets>
#include <QGraphicsWidget>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
//#include <QLabel>
#include "MainWindow.h"
#include "GLWidget.h"
#include <sstream>

using namespace std;

class MainWindow;

MainWindow::MainWindow(Simulation* Sim01)
 {
	MainScene = new QGraphicsScene;
	MainScene->setSceneRect(0, 0, 800, 480);
    MainGrid = new QGridLayout;
    CentralWidget = new QWidget;
    this->Sim01 = Sim01;
    nCoordBox = 6;

    setWindowTitle(tr("Hello!!"));
    GenerateControlPanel();
    SetUpGLWidget();
    SetUpCentralWidget();
    SetViewBackgroundColour();


    MainScene->update();

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(timerSimulationStep()));
	timer->start(0);
 };

MainWindow::~MainWindow(){
	cerr<<"called the destructor"<<endl;
	delete MainGLWidget;
	delete MainGrid;
	delete MainScene;
};

void MainWindow::SetViewBackgroundColour(){
}

void MainWindow::GenerateControlPanel(){
	ControlPanelMainHBox = new QVBoxLayout();
	ControlPanelMainHBox->setSpacing(2);

	//Generating Selection Display Panel:
	QGridLayout *SelectionDisplayGrid = new QGridLayout;
	SetUpSelectionDisplayGrid(SelectionDisplayGrid);
	ControlPanelMainHBox->addLayout(SelectionDisplayGrid,Qt::AlignTop);

	//Generating the quit button:
	QPushButton	*QuitButton = new QPushButton("Quit",this);
	QuitButton->setFixedWidth(100);
	//Connecting the quit button to close slot of main window, this will call the destructor.
	connect(QuitButton, SIGNAL(clicked()), this,SLOT(close()));
	ControlPanelMainHBox->addWidget(QuitButton,0,Qt::AlignBottom| Qt::AlignRight);

    //Adding the control panel vertical box to the main grind of the main window.
	MainGrid->addLayout(ControlPanelMainHBox,0,1,Qt::AlignLeft);
	MainGrid->setColumnStretch(1,-2);
}

void MainWindow::SetUpGLWidget(){
	MainGLWidget = new GLWidget();
	MainGLWidget->Sim01 = Sim01;
	cerr<<"From MainGLWIDGET: Element list size: "<<MainGLWidget->Sim01->Elements.size()<<endl;
	MainGrid->addWidget(MainGLWidget,0,0,Qt::AlignLeft);
	MainGrid->setColumnStretch(0,10);
}

void MainWindow::SetUpCentralWidget(){
    setCentralWidget(CentralWidget);
    CentralWidget->setParent(this);
    //CentralWidget->setStyleSheet("QWidget { background-color: LightGrey; }");
    CentralWidget->setStyleSheet("QWidget { background-color: rgb(233,229,243) }");
    CentralWidget->setLayout(MainGrid);
    connect(MainGLWidget, SIGNAL(SelectedItemChanged()), this, SLOT(SelectedItemChange()));
}

void MainWindow::SetUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid){
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);

	setItemSelectionTitles(font, boldFont, SelectionDisplayGrid);
	setCoordBoxes(font, boldFont, SelectionDisplayGrid);

	//Adding a last row with high stretch to push the upper columns to the top of the window:
	SelectionDisplayGrid->setRowStretch(4+nCoordBox,10);
}

void MainWindow::setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid){
	QLabel *PanelTitle = new QLabel("SelectedItemProperties");
	PanelTitle->setFont(boldFont);
	QLabel *NameTitle = new QLabel("Name:");
	NameTitle->setFont(boldFont);
	NameBox = new QLineEdit();
	NameBox->setPlaceholderText ( "No input" );
	NameBox->setReadOnly(true);
	NameBox->setFont(font);

	QLabel *CoordTitle = new QLabel("Coordinates");
	CoordTitle ->setFont(boldFont);
	QLabel *CoordTitlex = new QLabel("x");
	CoordTitlex ->setFont(boldFont);
	QLabel *CoordTitley = new QLabel("y");
	CoordTitley ->setFont(boldFont);
	QLabel *CoordTitlez = new QLabel("z");
	CoordTitlez ->setFont(boldFont);

	SelectionDisplayGrid->addWidget(PanelTitle,0,0,1,3,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameTitle,1,0,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameBox,1,1,1,2,Qt::AlignLeft);

	SelectionDisplayGrid->addWidget(CoordTitle,2,0,1,3,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlex,3,1,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitley,3,2,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlez,3,3,1,1,Qt::AlignHCenter);
}

void MainWindow::setCoordBoxes(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid){
	CoordLabel_n[0] = new QLabel("Node 1");
	CoordLabel_n[1] = new QLabel("Node 2");
	CoordLabel_n[2] = new QLabel("Node 3");
	CoordLabel_n[3] = new QLabel("Node 4");
	CoordLabel_n[4] = new QLabel("Node 5");
	CoordLabel_n[5] = new QLabel("Node 6");
	for (int i = 0 ;i<nCoordBox; ++i){
		CoordLabel_n[i] ->setFont(boldFont);
		CoordBox_x[i] = new QLineEdit();
		CoordBox_y[i] = new QLineEdit();
		CoordBox_z[i] = new QLineEdit();
		CoordBox_x[i]->setPlaceholderText( "No input" );
		CoordBox_y[i]->setPlaceholderText( "No input" );
		CoordBox_z[i]->setPlaceholderText( "No input" );
		CoordBox_x[i]->setReadOnly(true);
		CoordBox_y[i]->setReadOnly(true);
		CoordBox_z[i]->setReadOnly(true);
		CoordBox_x[i]->setFont(font);
		CoordBox_y[i]->setFont(font);
		CoordBox_z[i]->setFont(font);
		CoordBox_x[i]->setFixedWidth(70);
		CoordBox_y[i]->setFixedWidth(70);
		CoordBox_z[i]->setFixedWidth(70);
		SelectionDisplayGrid->addWidget(CoordLabel_n[i],i+4,0,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_x[i],i+4,1,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_y[i],i+4,2,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_z[i],i+4,3,1,1,Qt::AlignLeft);
	}
}

void MainWindow::SelectedItemChange(){
	//cerr<<"Main window saw the selection change"<<endl;
    QString tmpstring = QString::fromStdString(MainGLWidget->SelectedItemName);
    NameBox->setText(tmpstring);
    for (int i = 0 ;i<nCoordBox; ++i){
    	 CoordBox_x[i]->setText ( "" );
    	 CoordBox_y[i]->setText ( "" );
    	 CoordBox_z[i]->setText ( "" );
    	 CoordBox_x[i]->setEnabled(false);
    	 CoordBox_y[i]->setEnabled(false);
    	 CoordBox_z[i]->setEnabled(false);
    	 if (MainGLWidget->SelectedPos.size()>i*3){
			 CoordBox_x[i]->setText ( MainGLWidget->SelectedPos[i*3] );
			 CoordBox_y[i]->setText ( MainGLWidget->SelectedPos[i*3+1] );
			 CoordBox_z[i]->setText ( MainGLWidget->SelectedPos[i*3+2] );
			 CoordBox_x[i]->setEnabled(true);
			 CoordBox_y[i]->setEnabled(true);
			 CoordBox_z[i]->setEnabled(true);

		}
    }
 };

void MainWindow::timerSimulationStep(){
	//cout<<"Called the function via timer"<<endl;
	//while (Sim01->time<0.01){
		Sim01->RunOneStep();
		MainGLWidget->update();
	//}
}
