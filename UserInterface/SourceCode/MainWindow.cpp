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

#include "ElementBasicDisplay.h"
#include "ElementPropertySelection.h"

using namespace std;

class MainWindow;

MainWindow::MainWindow(Simulation* Sim01)
 {

	interatorForPressure = 0;

	MainScene = new QGraphicsScene;
    MainScene->setSceneRect(0, 0, 1000, 480);
    MainGrid = new QGridLayout;
    CentralWidget = new QWidget;
    this->Sim01 = Sim01;
    Sim01->calculateBoundingBox();
	double boundingBoxWidth  = Sim01->boundingBox[1][1] - Sim01->boundingBox[0][1];
    this->analyser01 = new Analysis(3, Sim01->saveDirectoryToDisplayString, Sim01->Nodes, boundingBoxWidth);

    nCoordBox = n_nodes_per_element;
    setWindowTitle(tr("Tissue Origami"));
    generateControlPanel();
    setUpGLWidget();
    setUpCentralWidget();
    setViewBackgroundColour();


    //QLabel *SimTime = new QLabel("SimTime");
    //SimTime->setText( "aaa" );
    simulationStartClock = std::clock();
    simulationStartTime = time(0);
    displayedSimulationLength = false;

    MainScene->update();

    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(timerSimulationStep()));
    cout<<"starting timer"<<endl;
	timer->start(0);
	cout<<"finalised mainwindow initiation"<<endl;
 };

MainWindow::~MainWindow(){
	//cerr<<"called the destructor"<<endl;
	delete Sim01;
	delete analyser01;
	delete MainGLWidget;
	delete MainGrid;
	delete MainScene;
};

void MainWindow::setViewBackgroundColour(){
}

void MainWindow::generateControlPanel(){
	ControlPanelMainHBox = new QVBoxLayout();
	ControlPanelMainHBox->setSpacing(2);

	// perpare the basic element info display using the default constructor
	ElementProps = new ElementBasicDisplay;
	// initialise the validators for the node and element selection boxes
	ElementProps->setNodeSelectionValidator(Sim01->Nodes.size() - 1, this);
	ElementProps->setElementSelectionValidator(Sim01->Elements.size() - 1, this);
	// create connections for the node and element selection boxes
	connect(&(ElementProps->node_selection_box), SIGNAL(textChanged(const QString &)), this, SLOT(manualNodeSelection(const QString &)));
	connect(&(ElementProps->element_selection_box), SIGNAL(textChanged(const QString &)), this, SLOT(manualElementSelection(const QString &)));
	// connect basic element display to update when lookingAtNewElement signal is sent out
	connect(this, SIGNAL(lookingAtNewElement(std::unique_ptr<ShapeBase> *)), ElementProps, SLOT(updateDisplayValues(std::unique_ptr<ShapeBase> *)));
	// connect to the main display
	ControlPanelMainHBox->addLayout(ElementProps,Qt::AlignTop);

	// prepare the user-specified element property information display
	PropertySelection = new ElementPropertySelection;
	// connect element property selection updates to lookingAtNewElement signal
	connect(this, SIGNAL(lookingAtNewElement(std::unique_ptr<ShapeBase> *)), PropertySelection, SLOT(updatePropertyValues(std::unique_ptr<ShapeBase> *)));
	// connect to the main display
	ControlPanelMainHBox->addLayout(PropertySelection,Qt::AlignCenter);

	// connect the export selected element properties button to the node information display
    connect(PropertySelection, SIGNAL(writeNodePositionsToFile(QString, std::unique_ptr<ShapeBase> *)), ElementProps, SLOT(writeNodePositions(QString, std::unique_ptr<ShapeBase> *)));

	// Generating project display options panel:
	QGridLayout *ProjectDisplayOptionsGrid = new QGridLayout;
    setUpProjectDisplayOptionGrid(ProjectDisplayOptionsGrid);
	ControlPanelMainHBox->addLayout(ProjectDisplayOptionsGrid,Qt::AlignTop);

	//Generating view options Panel:
	QGridLayout *ViewOptionsGrid = new QGridLayout;
	setUpViewOptionsGrid(ViewOptionsGrid);
	ControlPanelMainHBox->addLayout(ViewOptionsGrid,Qt::AlignBottom);

	//Generating the quit button:
	QHBoxLayout *BottomLineBox = new QHBoxLayout; // the bottom line will include quit button only for now
	QPushButton	*QuitButton = new QPushButton("Quit",this);
	QuitButton->setFixedWidth(100);
	//Connecting the quit button to close slot of main window, this will call the destructor.
	connect(QuitButton, SIGNAL(clicked()), this,SLOT(close()));
	BottomLineBox->addWidget(QuitButton,Qt::AlignBottom| Qt::AlignRight);
	ControlPanelMainHBox->addLayout(BottomLineBox,Qt::AlignBottom);

	//Adding the control panel vertical box to the main grid of the main window.
	MainGrid->addLayout(ControlPanelMainHBox,0,1,Qt::AlignLeft);
	MainGrid->setColumnStretch(1,-2);
}

void MainWindow::setUpGLWidget(){
	MainGLWidget = new GLWidget();
	MainGLWidget->Sim01 = Sim01;
	MainGLWidget->analyser01 = analyser01;
	MainGLWidget->currNodeNumber = Sim01->nNodes;
	MainGLWidget->DisplayStrainRange[0] = StrainSpinBoxes[0]->value();
	MainGLWidget->DisplayStrainRange[1] = StrainSpinBoxes[1]->value();
    for (int i =0; i<4; ++i){
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[0]->value();
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[1]->value();
    }
	//cerr<<"From MainGLWIDGET: Element list size: "<<MainGLWidget->Sim01->Elements.size()<<endl;
	//For background analysis, do not add the gl widget to the window, and make it tiny. It will be closed at the end of analuysis
    MainGrid->addWidget(MainGLWidget,0,0,Qt::AlignLeft);
    //MainGLWidget->resize(50,5);
	//setting up time display:
	QFont boldFont("SansSerif", 9, QFont::Bold,true);
	TimeTitle = new QLabel("Simulation Time: 0 hr");
	TimeTitle->setFont(boldFont);
	MainGrid->addWidget(TimeTitle,0,0,Qt::AlignBottom| Qt::AlignRight);
	MainGrid->setColumnStretch(0,10);
}

void MainWindow::setUpCentralWidget(){
    setCentralWidget(CentralWidget);
    CentralWidget->setParent(this);
    //CentralWidget->setStyleSheet("QWidget { background-color: LightGrey; }");
    CentralWidget->setStyleSheet("QWidget { background-color: rgb(233,229,243) }");
    CentralWidget->setLayout(MainGrid);
    connect(MainGLWidget, SIGNAL(SelectedItemChanged(bool)), this, SLOT(SelectedItemChange(bool)));
    connect(MainGLWidget, SIGNAL(NeedToClearManualElementSelection()), this, SLOT(ManualElementSelectionReset()));
    connect(MainGLWidget, SIGNAL(NeedToClearManualNodeSelection()), this, SLOT(ManualNodeSelectionReset()));
}

void MainWindow::setUpProjectDisplayOptionGrid(QGridLayout *ProjectDisplayOptionsGrid){
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	setStrainDisplayMenu(ProjectDisplayOptionsGrid);
	setPysPropDisplayMenu(ProjectDisplayOptionsGrid);
	setDisplayPreferences(ProjectDisplayOptionsGrid);
	int widthPhysProp = PysPropComboBox->minimumSizeHint().width();
	int widthStrain = StrainComboBox->minimumSizeHint().width();
	if (widthPhysProp>widthStrain){
		StrainComboBox->setMinimumWidth(widthPhysProp);
	}
	else{
		PysPropComboBox->setMinimumWidth(widthStrain);
	}
	//Adding a last row with high stretch to push the upper columns to the top of the window:
	ProjectDisplayOptionsGrid->setRowStretch(10,10);
}

void MainWindow::setUpViewOptionsGrid(QGridLayout *ViewOptionsGrid){
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	//Adding first row with high stretch to push the upper columns to the bottom of the window:
	ViewOptionsGrid->setRowStretch(0,10);
	//Adding last column with high stretch to push the previous columns to the left of the window:
	ViewOptionsGrid->setColumnStretch(10,10);

	ClippingSliders[0] = new QSlider(Qt::Horizontal);
	ClippingSliders[0]->setValue(99);
	connect(ClippingSliders[0], SIGNAL(valueChanged(int)), this, SLOT(xClipChange(int)));
	ClippingSliders[1] = new QSlider(Qt::Horizontal);
	ClippingSliders[1]->setValue(99);
	connect(ClippingSliders[1], SIGNAL(valueChanged(int)), this, SLOT(yClipChange(int)));
	ClippingSliders[2] = new QSlider(Qt::Horizontal);
	ClippingSliders[2]->setValue(99);
	connect(ClippingSliders[2], SIGNAL(valueChanged(int)), this, SLOT(zClipChange(int)));
	ViewOptionsGrid->addWidget(ClippingSliders[0],1,0,1,1);
	ViewOptionsGrid->addWidget(ClippingSliders[1],2,0,1,1);
	ViewOptionsGrid->addWidget(ClippingSliders[2],3,0,1,1);
	PerspectiveButton = new QPushButton("Switch To \n Orthagonal View",this);
	PerspectiveButton->setFixedWidth(150);
	//Connecting the button to toggle function
	connect(PerspectiveButton, SIGNAL(clicked()), this,SLOT(updateOrthagonalPerspectiveViewToggle()));
	ViewOptionsGrid->addWidget(PerspectiveButton,4,0,1,1);
	//Add button to switch to top view:
	TopViewButton = new QPushButton("Top \n View",this);
	connect(TopViewButton, SIGNAL(clicked()), this,SLOT(updateToTopView()));
	TopViewButton->setFixedWidth(60);
	ViewOptionsGrid->addWidget(TopViewButton,4,1,1,1);
	//Add button to switch to front view:
	FrontViewButton = new QPushButton("Front \n View",this);
	connect(FrontViewButton, SIGNAL(clicked()), this,SLOT(updateToFrontView()));
	FrontViewButton->setFixedWidth(60);
	ViewOptionsGrid->addWidget(FrontViewButton,4,2,1,1);
	//Add button to switch to side view:
	SideViewButton = new QPushButton("Side \n View",this);
	connect(SideViewButton, SIGNAL(clicked()), this,SLOT(updateToSideView()));
	SideViewButton->setFixedWidth(60);
	ViewOptionsGrid->addWidget(SideViewButton,4,3,1,1);
	//Add button to switch symmetricity view:
	SymmetricityDisplayButton = new QPushButton("Hide \n Symmetric",this);
	connect(SymmetricityDisplayButton, SIGNAL(clicked()), this,SLOT(updateDrawSymmetricityViewToggle()));
	SideViewButton->setFixedWidth(60);
	ViewOptionsGrid->addWidget(SymmetricityDisplayButton,4,4,1,1);
	//Add button to perspective view:
	PerspectiveViewButton = new QPushButton("Perspective \n View",this);
	connect(PerspectiveViewButton, SIGNAL(clicked()), this,SLOT(updateToPerspectiveView()));
	PerspectiveViewButton->setFixedWidth(60);
	ViewOptionsGrid->addWidget(PerspectiveViewButton,3,4,1,1);
}

void MainWindow::setStrainDisplayMenu(QGridLayout *ProjectDisplayOptionsGrid){
	DisplayCheckBoxes[0] = new QCheckBox("Strain");
	DisplayCheckBoxes[0]->setChecked(false);
	connect(DisplayCheckBoxes[0] , SIGNAL(stateChanged(int)),this,SLOT(updateStrainCheckBox(int)));

	StrainComboBox = new QComboBox();
	StrainComboBox->addItem("Volumetric Strain (via Fe)");
	StrainComboBox->addItem("Strain in DV");
	StrainComboBox->addItem("Strain in AP");
	StrainComboBox->addItem("Strain in AB");
	StrainComboBox->addItem("Shear in xy");
	StrainComboBox->addItem("Shear in yz");
	StrainComboBox->addItem("Shear in xz");
	StrainComboBox->setEnabled(false);
	connect(StrainComboBox , SIGNAL(currentIndexChanged(int)),this,SLOT(updateStrain(int)));

    StrainSpinBoxes[0] = new  QDoubleSpinBox();
    StrainSpinBoxes[1] = new  QDoubleSpinBox();
    StrainSpinBoxes[0]->setRange( -1.0, -0.00 );
    StrainSpinBoxes[0]->setSingleStep( 0.005 );
    StrainSpinBoxes[0]->setValue ( -0.35 );
    StrainSpinBoxes[0]->setDecimals(5);
    StrainSpinBoxes[0]->setEnabled(false);
    StrainSpinBoxes[1]->setRange( 0.0, 1.0 );
    StrainSpinBoxes[1]->setSingleStep( 0.005 );
    StrainSpinBoxes[1]->setValue( 0.35 );
    StrainSpinBoxes[1]->setDecimals(3);
    StrainSpinBoxes[1]->setEnabled(false);
	connect(StrainSpinBoxes[0] , SIGNAL(valueChanged(double)),this,SLOT(updateStrainSpinBoxes()));
    connect(StrainSpinBoxes[1], SIGNAL(valueChanged(double)), this, SLOT(updateStrainSpinBoxes()));

    ProjectDisplayOptionsGrid->addWidget(DisplayCheckBoxes[0],0,0,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainComboBox,1,0,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainSpinBoxes[0],2,0,1,1,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(StrainSpinBoxes[1],2,1,1,1,Qt::AlignLeft);
}

void MainWindow::setPysPropDisplayMenu(QGridLayout *ProjectDisplayOptionsGrid){

	DisplayCheckBoxes[1] = new QCheckBox("Physical Properties");
	DisplayCheckBoxes[1]->setChecked(false);
	connect(DisplayCheckBoxes[1] , SIGNAL(stateChanged(int)),this,SLOT(updatePysCheckBox(int)));

	PysPropComboBox = new QComboBox();
	PysPropComboBox->addItem("External Viscosity");
	PysPropComboBox->addItem("Internal Viscosity");
	PysPropComboBox->addItem("Young Modulus");
	PysPropComboBox->addItem("Poisson Ratio");
	PysPropComboBox->addItem("Volume (xyz) Growth Rate (fold per 24hr)");
	PysPropComboBox->addItem("Volume Growth (fold total)");
	PysPropComboBox->addItem("Emergent Size & Shape");
	PysPropComboBox->addItem("ShapeChangeRate_z");
	PysPropComboBox->setEnabled(false);
	connect(PysPropComboBox , SIGNAL(currentIndexChanged(int)),this,SLOT(updatePysProp(int)));

	PysPropSpinBoxes[0] = new  QDoubleSpinBox();
	PysPropSpinBoxes[1] = new  QDoubleSpinBox();
	PysPropSpinBoxes[0]->setRange ( 0, 10 );
	PysPropSpinBoxes[0]->setSingleStep( 0.1 );
	PysPropSpinBoxes[0]->setValue ( 3 );
	PysPropSpinBoxes[0]->setEnabled(false);
	PysPropSpinBoxes[1]->setRange( 0, 10.0 );
	PysPropSpinBoxes[1]->setSingleStep( 0.1 );
	PysPropSpinBoxes[1]->setValue( 5.0 );
	PysPropSpinBoxes[1]->setEnabled(false);
    connect(PysPropSpinBoxes[0], SIGNAL(valueChanged (double)), this, SLOT(updatePysPropSpinBoxes()));
    connect(PysPropSpinBoxes[1], SIGNAL(valueChanged (double)), this, SLOT(updatePysPropSpinBoxes()));


    //SelectionDisplayGrid->addWidget(DisplayCheckBoxes[1],4+nCoordBox,2,1,2,Qt::AlignLeft);
	//SelectionDisplayGrid->addWidget(PysPropComboBox,5+nCoordBox,2,1,2,Qt::AlignLeft);
	//SelectionDisplayGrid->addWidget(PysPropSpinBoxes[0],6+nCoordBox,2,1,1,Qt::AlignLeft);
	//SelectionDisplayGrid->addWidget(PysPropSpinBoxes[1],6+nCoordBox,3,1,1,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(DisplayCheckBoxes[1],0,2,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropComboBox,1,2,1,2,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropSpinBoxes[0],2,2,1,1,Qt::AlignLeft);
    ProjectDisplayOptionsGrid->addWidget(PysPropSpinBoxes[1],2,3,1,1,Qt::AlignLeft);
}

void MainWindow::setDisplayPreferences(QGridLayout *ProjectDisplayOptionsGrid){
    //draw pipette aspiration CheckBox
    DisplayPreferencesCheckBoxes[0] = new QCheckBox("Display Pipette");
	DisplayPreferencesCheckBoxes[0]->setChecked(false);
	//draw net forces checkbox
	DisplayPreferencesCheckBoxes[1] = new QCheckBox("Net Forces");
	DisplayPreferencesCheckBoxes[1]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[1] , SIGNAL(stateChanged(int)),this,SLOT(updateNetForceCheckBox(int)));
	//draw fixed Nodes checkbox
	DisplayPreferencesCheckBoxes[2] = new QCheckBox("Fixed Nodes");
	DisplayPreferencesCheckBoxes[2]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[2] , SIGNAL(stateChanged(int)),this,SLOT(updateFixedNodesCheckBox(int)));
	DisplayPreferencesCheckBoxes[3] = new QCheckBox("ScaleCube");
	DisplayPreferencesCheckBoxes[3]->setChecked(true);
	connect(DisplayPreferencesCheckBoxes[3] , SIGNAL(stateChanged(int)),this,SLOT(updateScaleBarCheckBox(int)));
	DisplayPreferencesCheckBoxes[4] = new QCheckBox("Display Peripodial Membrane");
	DisplayPreferencesCheckBoxes[4]->setChecked(true);
	connect(DisplayPreferencesCheckBoxes[4] , SIGNAL(stateChanged(int)),this,SLOT(updatePeripodialDisplayCheckBox(int)));
	DisplayPreferencesCheckBoxes[5] = new QCheckBox("Display Columnar Layer");
	DisplayPreferencesCheckBoxes[5]->setChecked(true);
	connect(DisplayPreferencesCheckBoxes[5] , SIGNAL(stateChanged(int)),this,SLOT(updateColumnarLayerDisplayCheckBox(int)));
	//draw packing forces checkbox
	DisplayPreferencesCheckBoxes[6] = new QCheckBox("Packing Forces");
	DisplayPreferencesCheckBoxes[6]->setChecked(false);
    connect(DisplayPreferencesCheckBoxes[6] , SIGNAL(stateChanged(int)),this,SLOT(updatePackingForceCheckBox(int)));
	//draw Bounding Box checkbox
	DisplayPreferencesCheckBoxes[7] = new QCheckBox("Bounding Box");
    DisplayPreferencesCheckBoxes[7]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[7] , SIGNAL(stateChanged(int)),this,SLOT(updateBoundingBoxCheckBox(int)));
	//draw net forces checkbox
	//draw marking ellipses
	DisplayPreferencesCheckBoxes[8] = new QCheckBox("Marking Ellipses");
	DisplayPreferencesCheckBoxes[8]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[8] , SIGNAL(stateChanged(int)),this,SLOT(updateMarkingEllipseCheckBox(int)));

	//Draw volume redistribution checkbox:
	DisplayPreferencesCheckBoxes[9] = new QCheckBox("Growth redistribution");
	DisplayPreferencesCheckBoxes[9]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[9] , SIGNAL(stateChanged(int)),this,SLOT(updateGrowthRedistributionCheckBox(int)));

	//Draw node binding
	DisplayPreferencesCheckBoxes[10] = new QCheckBox("Node binding");
	DisplayPreferencesCheckBoxes[10]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[10] , SIGNAL(stateChanged(int)),this,SLOT(updateDrawNodeBindingCheckBox(int)));

	//Draw Lumen:
	DisplayPreferencesCheckBoxes[11] = new QCheckBox("Mark Lumen");
	DisplayPreferencesCheckBoxes[11]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[11] , SIGNAL(stateChanged(int)),this,SLOT(updateLumenDisplayCheckBox(int)));


    ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[0],3,0,1,2,Qt::AlignLeft);  // display pipette
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[7],3,2,1,2,Qt::AlignLeft);  // display bounding box
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[1],4,0,1,2,Qt::AlignLeft);  // Net Forces
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[6],4,2,1,2,Qt::AlignLeft);  // Packing Forces
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[2],5,0,1,1,Qt::AlignLeft);  // Fixed Nodes
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[3],6,0,1,2,Qt::AlignLeft); // Scale Bar
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[4],7,0,1,2,Qt::AlignLeft); // Display Peripodial Membrane
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[5],8,0,1,2,Qt::AlignLeft); // Display Columnar Layer
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[8],7,2,1,2,Qt::AlignLeft); // Display Marked Ellipses
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[9],8,2,1,2,Qt::AlignLeft); // Display Volume redistribution
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[10],9,2,1,2,Qt::AlignLeft); // Display node binding
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[11],9,0,1,2,Qt::AlignLeft); // Display Lumen

}



void  MainWindow::xClipChange(int k){
	MainGLWidget->xClip = Sim01->boundingBox[0][0] +( ( Sim01->boundingBox[1][0] - Sim01->boundingBox[0][0] ) * (double) (k+10)/100.0 );
	MainGLWidget->updateClipping();
	//cout<<"x:" <<k<<" "<<MainGLWidget->xClip<<endl;
}

void  MainWindow::yClipChange(int k){
	MainGLWidget->yClip = Sim01->boundingBox[0][1] +( ( Sim01->boundingBox[1][1] - Sim01->boundingBox[0][1] ) * (double) (90-k)/100.0 );
	MainGLWidget->updateClipping();
	//cout<<"y:" <<k<<" "<<MainGLWidget->yClip<<endl;
}

void  MainWindow::zClipChange(int k){
	MainGLWidget->zClip = Sim01->boundingBox[0][2] +( ( Sim01->boundingBox[1][2] - Sim01->boundingBox[0][2] ) * (double) (k+10)/100.0 );
	//MainGLWidget->zClip = Sim01->boundingBox[0][0] +( ( Sim01->boundingBox[1][0] - Sim01->boundingBox[0][0] ) * (double) (k+10)/100.0 );
	MainGLWidget->updateClipping();
	//cout<<"z:" <<k<<" "<<MainGLWidget->zClip<<endl;
}

void  MainWindow::updateOrthagonalPerspectiveViewToggle(){
	//cout<<"button clicked, perspevctiveView: "<<MainGLWidget->PerspectiveView<<endl;
	if (MainGLWidget->PerspectiveView){
		//the view was perspective, toggling to orthagonal, and changing the text in button for future toggle choice:
		MainGLWidget->PerspectiveView = false;
		PerspectiveButton->setText("Switch To \n Perspective View");
	}
	else{
		MainGLWidget->PerspectiveView = true;
		PerspectiveButton->setText("Switch To \n Orthagonal View");
	}
}

void  MainWindow::updateDrawSymmetricityViewToggle(){
	//cout<<"button clicked, perspevctiveView: "<<MainGLWidget->PerspectiveView<<endl;
	if (MainGLWidget->drawSymmetricity){
		//the view was perspective, toggling to orthagonal, and changing the text in button for future toggle choice:
		MainGLWidget->drawSymmetricity = false;
		SymmetricityDisplayButton->setText("Show \n Symmetric");
	}
	else{
		MainGLWidget->drawSymmetricity = true;
		SymmetricityDisplayButton->setText("Hide \n Symmetric");
	}
}

void  MainWindow::updateToTopView(){
	MainGLWidget->updateToTopView();
}

void  MainWindow::updateToFrontView(){
	MainGLWidget->updateToFrontView();
}

void  MainWindow::updateToSideView(){
	MainGLWidget->updateToSideView();
}

void  MainWindow::updateToPerspectiveView(){
	MainGLWidget->updateToPerspectiveView();
}

void  MainWindow::updateNetForceCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawNetForces = true;
	else
		MainGLWidget->drawNetForces = false;
}

void  MainWindow::updateMarkingEllipseCheckBox(int s){
	if (s ==2 ){
		MainGLWidget->drawMarkingEllipses = true;
	}
	else{
		MainGLWidget->drawMarkingEllipses = false;
	}
}

void  MainWindow::updateGrowthRedistributionCheckBox(int s){
	if (s ==2 ){
		MainGLWidget->drawGrowthRedistribution = true;
	}
	else{
		MainGLWidget->drawGrowthRedistribution = false;
	}
}

void  MainWindow::updateDrawNodeBindingCheckBox(int s){
	if (s ==2 ){
		MainGLWidget->drawNodeBinding = true;
	}
	else{
		MainGLWidget->drawNodeBinding = false;
	}
}


void  MainWindow::updateDisplayPipette(int s){
    if ( s == 2 )
        MainGLWidget->displayPipette = true;
    else
        MainGLWidget->displayPipette = false;
}

void  MainWindow::updatePackingForceCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawPackingForces = true;
	else
		MainGLWidget->drawPackingForces = false;
}


void  MainWindow::updateFixedNodesCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->DisplayFixedNodes = true;
	else
		MainGLWidget->DisplayFixedNodes = false;
}

void  MainWindow::updateScaleBarCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawTissueScaleBar = true;
	else
		MainGLWidget->drawTissueScaleBar = false;
}

void  MainWindow::updatePeripodialDisplayCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawPeripodialMembrane = true;
	else
		MainGLWidget->drawPeripodialMembrane = false;
}

void  MainWindow::updateColumnarLayerDisplayCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawColumnar = true;
	else
		MainGLWidget->drawColumnar = false;
}

void  MainWindow::updateLumenDisplayCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->drawLumen = true;
	else
		MainGLWidget->drawLumen = false;
}

void  MainWindow::updateBoundingBoxCheckBox(int s){
	if ( s == 2 )
		MainGLWidget->displayBoundingBox  = true;
	else
		MainGLWidget->displayBoundingBox  = false;
}

void MainWindow::updateStrain(int s){
	MainGLWidget->StrainToDisplay = s;
	MainGLWidget->update();
	cout<<"Strain to display: "<<MainGLWidget->StrainToDisplay <<endl;
}


void MainWindow::updatePysProp(int s){
	MainGLWidget->PysPropToDisplay = s;
	float low = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0];
	float high = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1];
	float min[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][0],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][2]};
	float max[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][1],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][3]};
	int   decimals = MainGLWidget->DisplayPysPropDecimals[MainGLWidget->PysPropToDisplay];
	float step = MainGLWidget->DisplayPysPropSteps[MainGLWidget->PysPropToDisplay];
	PysPropSpinBoxes[0]->setRange( min[0], max[0] );
	PysPropSpinBoxes[0]->setValue (low);
	PysPropSpinBoxes[0]->setDecimals(decimals);
	PysPropSpinBoxes[0]->setSingleStep( step );
	PysPropSpinBoxes[1]->setRange( min[1], max[1] );
	PysPropSpinBoxes[1]->setValue (high);
	PysPropSpinBoxes[1]->setDecimals(decimals);
	PysPropSpinBoxes[1]->setSingleStep( step );
	//MainGLWidget->update();
	//cout<<"Pys prop to display: "<<MainGLWidget->StrainToDisplay <<endl;

}
void MainWindow::updateStrainCheckBox(int s){
	if (s == 0){
		MainGLWidget->StrainToDisplay = -1;
		StrainComboBox->setEnabled(false);
		StrainSpinBoxes[0]->setEnabled(false);
		StrainSpinBoxes[1]->setEnabled(false);
		MainGLWidget->DisplayStrains = false;
	}
	else{
		DisplayCheckBoxes[1]->setChecked(false);
		MainGLWidget->StrainToDisplay = StrainComboBox->currentIndex ();
		StrainComboBox->setEnabled(true);
		StrainSpinBoxes[0]->setEnabled(true);
		StrainSpinBoxes[1]->setEnabled(true);
		MainGLWidget->DisplayStrains = true;
	}
	//MainGLWidget->update();
	cout<<"Strain Check Box updated, active/nonactive: "<<s <<endl;

}

void MainWindow::updatePysCheckBox(int s){
	if (s == 0){
		//cout<<"updating phys prop checkbox, signal is false"<<endl;
		MainGLWidget->PysPropToDisplay = -1;
		PysPropComboBox->setEnabled(false);
		PysPropSpinBoxes[0]->setEnabled(false);
		PysPropSpinBoxes[1]->setEnabled(false);
		MainGLWidget->DisplayPysProp = false;
	}
	else{
		//cout<<"updating phy prop checkbox, signal is true"<<endl;
		DisplayCheckBoxes[0]->setChecked(false);
		MainGLWidget->PysPropToDisplay = PysPropComboBox->currentIndex();
		float low = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0];
		float high = MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1];
		float min[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][0],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][2]};
		float max[2] = {MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][1],MainGLWidget->DisplayPysPropBounds[MainGLWidget->PysPropToDisplay][3]};
		int decimals = MainGLWidget->DisplayPysPropDecimals[MainGLWidget->PysPropToDisplay];
		float step =  MainGLWidget->DisplayPysPropSteps[MainGLWidget->PysPropToDisplay];

		PysPropSpinBoxes[0]->setRange( min[0], max[0] );
		PysPropSpinBoxes[0]->setValue (low);
		PysPropSpinBoxes[0]->setDecimals(decimals);
		PysPropSpinBoxes[0]->setSingleStep( step );
		PysPropSpinBoxes[1]->setRange( min[1], max[1] );
		PysPropSpinBoxes[1]->setValue (high);
		PysPropSpinBoxes[1]->setDecimals(decimals);
		PysPropSpinBoxes[1]->setSingleStep( step );
		PysPropComboBox->setEnabled(true);
		PysPropSpinBoxes[0]->setEnabled(true);
		PysPropSpinBoxes[1]->setEnabled(true);
		MainGLWidget->DisplayPysProp = true;
	}
	//MainGLWidget->update();
}

void MainWindow::updateStrainSpinBoxes(){
	MainGLWidget->DisplayStrainRange[0] = StrainSpinBoxes[0]->value();
	MainGLWidget->DisplayStrainRange[1] = StrainSpinBoxes[1]->value();
	//MainGLWidget->update();
	//cout<<"Strain Spin Box updated: "<<StrainSpinBoxes[0]->value()<<" "<<StrainSpinBoxes[1]->value()<<endl;
}


void MainWindow::updatePysPropSpinBoxes(){
	MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][0] = PysPropSpinBoxes[0]->value();
	MainGLWidget->DisplayPysPropRange[MainGLWidget->PysPropToDisplay][1] = PysPropSpinBoxes[1]->value();
	//MainGLWidget->update();
	//cout<<"Physical property Spin Box updated: "<<StrainSpinBoxes[0]->value()<<" "<<StrainSpinBoxes[1]->value()<<endl;

}

void MainWindow::updateTimeText(){
	//cout<<"inside updateTimeText"<<endl;
	//round to two digits:
	QString timeStringInSec = QString::number(Sim01->currSimTimeSec);
	QString timeStringInHr = QString::number(Sim01->currSimTimeSec/3600.0, 'f', 2);
	QString timeStringInHrAEL = QString::number((Sim01->currSimTimeSec/3600.0 + 48), 'f', 2);
	timeStringInSec = "Simulation Time: " + timeStringInHr + " hr ("+timeStringInHrAEL+" hr AEL) - ( "+ timeStringInSec + " sec, " + QString::number(Sim01->timestep) + " steps )"; ;
	TimeTitle->setText(timeStringInSec);
	//cout<<"finalised updateTimeText"<<endl;
}

void MainWindow::SelectedItemChange(bool element_found){
	std::unique_ptr<ShapeBase> *new_element = nullptr;
	// fetch the name of the selected element, if we selected an element
	// otherwise, nullptr indicates that deselection has occurred
    if (element_found) {
		int new_element_index = MainGLWidget->SelectedItemIndex;	// the new element index
		new_element = &Sim01->Elements[new_element_index]; 			// pointer to the new element
	}
	emit lookingAtNewElement(new_element);
 };

void MainWindow::manualNodeSelection(const QString &newValue){
	MainGLWidget->manualNodeSelection(newValue.toInt());
	//cerr<<"Manual Node Selection Update"<<newValue.toInt()<<endl;
}

void MainWindow::manualElementSelection(const QString &newValue){
	MainGLWidget->manualElementSelection(newValue.toInt());
	//cerr<<"Manual Element Selection Update: "<<newValue.toInt()<<endl;
}

void MainWindow::ManualElementSelectionReset(){
	MainGLWidget->ManualNodeSelection = false;
	MainGLWidget->ManualSelectedNodeId = -100;
	// block signals whilst resetting
	ElementProps->element_selection_box.blockSignals(true);
	// reset the text in the selection box
	ElementProps->setElementSelectionValidator(Sim01->Elements.size()-1, this);
	ElementProps->element_selection_box.setText("");
	// reopen to user input
	ElementProps->element_selection_box.blockSignals(false);
}

void MainWindow::ManualNodeSelectionReset(){
	// block signals whilst resetting
	ElementProps->node_selection_box.blockSignals(true);
	// reset the text in the selection box
	ElementProps->setNodeSelectionValidator(Sim01->Nodes.size()-1, this);
	ElementProps->node_selection_box.setText("");
	// reopen to user input
	ElementProps->node_selection_box.blockSignals(false);
}

void MainWindow::testAdhesionsAndCurveConstruction(){
	//testing adhesion!!
	Sim01->thereIsEmergentEllipseMarking = true;
	Sim01->detectPacingNodes();
	int n = Sim01->pacingNodeCouples0.size();
	for(int nodeCoupleIterator = 0 ; nodeCoupleIterator<n; ++nodeCoupleIterator){
		int slaveNodeId = Sim01->pacingNodeCouples0[nodeCoupleIterator];
		int masterNodeId = Sim01->pacingNodeCouples1[nodeCoupleIterator];
		cout<<"checking node pair: "<<slaveNodeId<<" "<<masterNodeId<<endl;
		for (int i=0;i<3; ++i){
			Sim01->Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
			Sim01->Nodes[masterNodeId]->isMaster[i] = true;
		}
	}
	for (const auto& itNode : Sim01->Nodes){
		if (itNode->onFoldInitiation){
			int n = itNode->connectedElementIds.size();
			for (int i=0; i<n; ++i){
				int currElementId = itNode->connectedElementIds[i];
				//check if at least two nodes of the element are at curves:
				bool changeEllipseId = Sim01->Elements[currElementId]->hasEnoughNodesOnCurve(Sim01->Nodes);
				if (changeEllipseId){
					//cout<<"changing Id for element "<<currElementId<<endl;
					Sim01->Elements[currElementId]->insideEllipseBand = true;
					Sim01->Elements[currElementId]->coveringEllipseBandId = 100;
					//Sim01->Elements[currElementId]->assignEllipseBandIdToWholeTissueColumn(Sim01->TissueHeightDiscretisationLayers,Sim01->Nodes,Sim01->Elements);
				}
			}
		}
	}
}
void MainWindow::timerSimulationStep(){
    //cout<<"Called the function via timer"<<endl;

    bool 	automatedSave = false;
    int		viewSelection = -1; //0: top, 1: cross, 2: perspective 3: side, 4: clone cross-section
    int 	displayAutomatedStrain = -1; //-1 no strain display, 1 display DV strains, 2 AP strains;
    bool 	analyseResults = false;
    bool 	slowstepsOnDisplay = false;
	bool 	slowstepsOnRun = false;
	int 	slowWaitTime = 0;
	//display DV strains
	if (displayAutomatedStrain>-1){
		MainGLWidget->DisplayStrains = true;
		MainGLWidget->StrainToDisplay =displayAutomatedStrain;
	}
	float zClipValueForCloneCross =99;
	if (Sim01->DisplaySave){
		if (Sim01->timestep == 0){
			if( automatedSave ){
				Sim01->assignTips();
				if (viewSelection == 0){
					MainGLWidget->updateToTopView(); //display tissue from the top view
					MainGLWidget->drawSymmetricity = true; //show symmetric
				}
				else if (viewSelection == 1){
					MainGLWidget->updateToFrontView(); //display tissue from the front view
					MainGLWidget->drawSymmetricity = false; //hide symmetric
				}
				else if (viewSelection == 2){
					MainGLWidget->updateToPerspectiveView(); //display tissue from the tilted view
					MainGLWidget->drawSymmetricity = true; //show symmetric
					MainGLWidget->drawTissueScaleBar= false;
				}
				else if (viewSelection == 3){
					MainGLWidget->updateToSideView(); //display tissue from the tilted view
					MainGLWidget->drawSymmetricity = true; //show symmetric
				}
				else if (viewSelection == 4){
					//MainGLWidget->updateToCloneCrossView();
					MainGLWidget->updateToTopView();
					MainGLWidget->drawSymmetricity = false;

					double *c = Sim01->Elements[14415]->getCentre();
					zClipValueForCloneCross = c[0]+c[1]+12.5;
					MainGLWidget->zClip=zClipValueForCloneCross;
					xClipChange(99);
					yClipChange(99);
					MainGLWidget->updateClipping();
				}
				//MainGLWidget->drawPeripodialMembrane= false;
				MainGLWidget->PerspectiveView = false; //switch to orthogoanal view type.
				Sim01->saveImages = true;
				Sim01->saveDirectory = Sim01->saveDirectoryToDisplayString;
			}
		}
		if( automatedSave && viewSelection == 4){
			//

			double *c = Sim01->Elements[14415]->getCentre();
			zClipValueForCloneCross = c[0]+c[1]+12.5;
			MainGLWidget->zClip=zClipValueForCloneCross;
			delete c;

			//
			xClipChange(99);
			yClipChange(99);
			MainGLWidget->updateClipping();
		}
		if(!Sim01->reachedEndOfSaveFile){
			cout<<" updating step"<<endl;
			Sim01->updateOneStepFromSave();
			//testAdhesionsAndCurveConstruction();

			//end of testing adhesion
			if (Sim01->timestep >= 0){
				if (analyseResults){
					double boundingBoxLength = Sim01->boundingBox[1][0] - Sim01->boundingBox[0][0];
					double boundingBoxWidth  = Sim01->boundingBox[1][1] - Sim01->boundingBox[0][1];
					cout<<" calculateBoundingBoxSizeAndAspectRatio"<<endl;
					analyser01->calculateBoundingBoxSizeAndAspectRatio(Sim01->currSimTimeSec,boundingBoxLength,boundingBoxWidth);
					cout<<" calculateContourLineLengthsDV"<<endl;
					analyser01->calculateContourLineLengthsDV(Sim01->Nodes);
					analyser01->saveApicalCircumferencePosition(Sim01->currSimTimeSec,Sim01->Nodes);

					cout<<" finished analysis"<<endl;
					Sim01->thereIsEmergentEllipseMarking = true;
					Sim01->detectPacingNodes();
					analyser01->saveNodesOnFold(Sim01->currSimTimeSec,Sim01->Nodes);
				}
			}
			Sim01->calculateDVDistance();
			for (int a = 0; a<0; a++){ //11 for 6 hours with 1800 sec time step
				Sim01->updateOneStepFromSave();
				Sim01->calculateDVDistance();
			}
			cout<<" updating time text"<<endl;
			//Sim01->fixNode0InPosition(36,0,0);
			updateTimeText();
			QTime dieTime= QTime::currentTime().addSecs(0.01);
            while( QTime::currentTime() < dieTime ){
			    QCoreApplication::processEvents(QEventLoop::AllEvents, 1);
            }
			if (Sim01->saveImages){
				takeScreenshot();
			}
			if (slowstepsOnDisplay){
				QTime dieTime= QTime::currentTime().addSecs(slowWaitTime);
				while( QTime::currentTime() < dieTime ){
					QCoreApplication::processEvents(QEventLoop::AllEvents, slowWaitTime);
				}
			}
			cout<<"end of loop"<<endl;
		}
		else{
			//close();
			if (automatedSave){
				MainGLWidget->close();
				close();
			}
		}
	}
	else{
		//cout<<" step: "<<Sim01->timestep<<" currSimTimeSec: "<<Sim01->currSimTimeSec<<" simLength: "<< Sim01->SimLength<<endl;
		if (Sim01->currSimTimeSec <= Sim01->SimLength){
			cout<<"started step: "<<Sim01->currSimTimeSec<<" length: "<<Sim01->SimLength<<endl;
			bool Success = Sim01->runOneStep();
			updateTimeText();
			if (slowstepsOnRun){
				QTime dieTime= QTime::currentTime().addSecs(slowWaitTime);
				while( QTime::currentTime() < dieTime ){
					QCoreApplication::processEvents(QEventLoop::AllEvents, slowWaitTime);
				}
			}
			if (Sim01->saveImages && Sim01->timestep%Sim01->imageSaveInterval == 0){
				takeScreenshot();
			}
			if (Success == false ){
				//there is a flipped element, I am not continuing simulation
				cout<<"there is a flipped element, I am not continuing simulation"<<endl;
				Sim01->timestep = 2*Sim01->SimLength;
			}
		}
        else if(!displayedSimulationLength){
        	displayedSimulationLength = true;
            Sim01->wrapUpAtTheEndOfSimulation();
            Sim01->writeRelaxedMeshFromCurrentState();
            //Sim01->writeMeshRemovingAblatedRegions();
            //double durationClock = ( std::clock() - simulationStartClock ) / (double) CLOCKS_PER_SEC;
            //double durationTime = std::difftime(std::time(0), simulationStartTime);
            //cout<<"Simulation time: "<<durationTime<<" sec, Simulation clock: "<<durationClock<<" sec"<<endl;
            //cout<<" Final lumen volume: "<<Sim01->tissueLumen->currentIdealVolume<<endl;
            //close();
        }
    }
    MainGLWidget->update();
}

void MainWindow::takeScreenshot(){
	cout<<"taking screenshot"<<endl;
	QPixmap originalPixmap;
	QScreen *screen = qApp->primaryScreen();
	int x = MainGLWidget->geometry().x();
	int y = MainGLWidget->geometry().y();
	int w = MainGLWidget->width();
	int h = MainGLWidget->height();
	originalPixmap = screen->grabWindow(this->winId(),x,y,w,h);
	//QString timepoint = QString::number(timestep);
	QString timepoint = QString("%1").arg(Sim01->timestep,6,10,QChar('0'));
	QString timeinsec = QString("%1").arg(Sim01->timestep*Sim01->dt);
	QString directory = QString(Sim01->saveDirectory.c_str());
	QString fileName = directory+"/ScreenShots/frame"+ timepoint +"-"+timeinsec+"sec.png";
	cout<<fileName.toStdString()<<endl;
	originalPixmap.save(fileName, "png");
}
