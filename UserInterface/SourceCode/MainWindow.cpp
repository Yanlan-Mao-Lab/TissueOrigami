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

	interatorForPressure = 0;

	MainScene = new QGraphicsScene;
    MainScene->setSceneRect(0, 0, 1000, 480);
    MainGrid = new QGridLayout;
    CentralWidget = new QWidget;
    this->Sim01 = Sim01;
    nCoordBox = 6;

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
	delete MainGLWidget;
	delete MainGrid;
	delete MainScene;
};

void MainWindow::setViewBackgroundColour(){
}

void MainWindow::generateControlPanel(){
	ControlPanelMainHBox = new QVBoxLayout();
	ControlPanelMainHBox->setSpacing(2);

	//Generating Selection Display Panel:
	QGridLayout *SelectionDisplayGrid = new QGridLayout;
	setUpSelectionDisplayGrid(SelectionDisplayGrid);
	ControlPanelMainHBox->addLayout(SelectionDisplayGrid,Qt::AlignTop);

	//Generating project display options panel:
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
	MainGLWidget->currNodeNumber = Sim01->nNodes;
	MainGLWidget->DisplayStrainRange[0] = StrainSpinBoxes[0]->value();
	MainGLWidget->DisplayStrainRange[1] = StrainSpinBoxes[1]->value();
    for (int i =0; i<4; ++i){
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[0]->value();
    	MainGLWidget->DisplayPysPropRange[i][0] = PysPropSpinBoxes[1]->value();
    }
	//cerr<<"From MainGLWIDGET: Element list size: "<<MainGLWidget->Sim01->Elements.size()<<endl;
	//MainGrid->setStyleSheet("border: 1px solid red");
	MainGrid->addWidget(MainGLWidget,0,0,Qt::AlignLeft);
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
    connect(MainGLWidget, SIGNAL(SelectedItemChanged()), this, SLOT(SelectedItemChange()));
    connect(MainGLWidget, SIGNAL(NeedToClearManualElementSelection()), this, SLOT(ManualElementSelectionReset()));
    connect(MainGLWidget, SIGNAL(NeedToClearManualNodeSelection()), this, SLOT(ManualNodeSelectionReset()));
}

void MainWindow::setSelectionByIdSection(QFont font1, QGridLayout *SelectionDisplayGrid){
	//QFont boldFont("SansSerif", 9, QFont::Bold,true);
	//QFont font("SansSerif", 9);
	QLabel *NodeSelectTitle = new QLabel("Select <br> Node:");
	NodeSelectTitle->setFont(font1);
	//setWordWrap(true);
	NodeSelectBox = new QLineEdit();
	NodeSelectBox->setPlaceholderText ( QString("# 0-%1").arg(Sim01->Nodes.size()-1) );
	NodeSelectBox->setFont(font1);
	NodeSelectBox->setStyleSheet("background-color: white");
	NodeSelectBox->setFixedWidth(70);
	NodeSelectBox->setValidator( new QIntValidator(0, Sim01->Nodes.size()-1, this) );

	QLabel *ElementSelectTitle = new QLabel("Select <br> Element:");
	ElementSelectTitle->setFont(font1);
	ElementSelectBox = new QLineEdit();
	ElementSelectBox->setPlaceholderText ( QString("# 0-%1").arg(Sim01->Elements.size()-1) );
	ElementSelectBox->setFont(font1);
	ElementSelectBox->setStyleSheet("background-color: white");
	ElementSelectBox->setFixedWidth(70);
	ElementSelectBox->setValidator( new QIntValidator(0, Sim01->Elements.size()-1, this) );
	connect(NodeSelectBox, SIGNAL(textChanged(const QString &)), this, SLOT(manualNodeSelection(const QString &)));

	SelectionDisplayGrid->addWidget(NodeSelectTitle,0,3,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(ElementSelectTitle,0,4,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NodeSelectBox,1,3,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(ElementSelectBox,1,4,1,1,Qt::AlignLeft);
	connect(ElementSelectBox, SIGNAL(textChanged(const QString &)), this, SLOT(manualElementSelection(const QString &)));
}

void MainWindow::setUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid){
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	setItemSelectionTitles(font, boldFont, SelectionDisplayGrid);
	setCoordBoxes(font, boldFont, SelectionDisplayGrid);
	setSelectionByIdSection(font, SelectionDisplayGrid);
}

void MainWindow::setUpProjectDisplayOptionGrid(QGridLayout *ProjectDisplayOptionsGrid){
	QFont boldFont("SansSerif", 10, QFont::Bold,true);
	QFont font("SansSerif", 10);
	setStrainDisplayMenu(ProjectDisplayOptionsGrid);
	setPysPropDisplayMenu(ProjectDisplayOptionsGrid);
	setDisplayPreferences(ProjectDisplayOptionsGrid);
	setMyosinComboBox(ProjectDisplayOptionsGrid);
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

void MainWindow::setMyosinComboBox(QGridLayout *ProjectDisplayOptionsGrid){
	MyosinComboBox = new QComboBox();
	MyosinComboBox->addItem("Uniform Myosin");
	MyosinComboBox->addItem("Uni-polar myosin");
	MyosinComboBox->setEnabled(true);
	connect(MyosinComboBox , SIGNAL(currentIndexChanged(int)),this,SLOT(updateMyosinComboBox(int)));
	ProjectDisplayOptionsGrid->addWidget(MyosinComboBox,6,2,1,2,Qt::AlignLeft);
}

void MainWindow::setStrainDisplayMenu(QGridLayout *ProjectDisplayOptionsGrid){
	DisplayCheckBoxes[0] = new QCheckBox("Strain");
	DisplayCheckBoxes[0]->setChecked(false);
	connect(DisplayCheckBoxes[0] , SIGNAL(stateChanged(int)),this,SLOT(updateStrainCheckBox(int)));

	StrainComboBox = new QComboBox();
	StrainComboBox->addItem("Average(DV-AP-AB) Strain");
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

void MainWindow::setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid){
	QLabel *PanelTitle = new QLabel("Selected Item Properties");
	PanelTitle->setFont(boldFont);
	QLabel *NameTitle = new QLabel("Name:");
	NameTitle->setFont(boldFont);
	NameBox = new QLineEdit();
	NameBox->setPlaceholderText ( "No input" );
	NameBox->setReadOnly(true);
	NameBox->setFont(font);

	//QLabel *CoordTitle = new QLabel("Coordinates");
	//CoordTitle ->setFont(boldFont);
	QLabel *NodeIdTitle = new QLabel("id");
	NodeIdTitle ->setFont(boldFont);
	QLabel *CoordTitlex = new QLabel("x");
	CoordTitlex ->setFont(boldFont);
	QLabel *CoordTitley = new QLabel("y");
	CoordTitley ->setFont(boldFont);
	QLabel *CoordTitlez = new QLabel("z");
	CoordTitlez ->setFont(boldFont);

	SelectionDisplayGrid->addWidget(PanelTitle,0,0,1,2,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameTitle,1,0,1,1,Qt::AlignLeft);
	SelectionDisplayGrid->addWidget(NameBox,1,1,1,2,Qt::AlignLeft);

	//SelectionDisplayGrid->addWidget(CoordTitle,2,0,1,3,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(NodeIdTitle,2,1,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlex,2,2,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitley,2,3,1,1,Qt::AlignHCenter);
	SelectionDisplayGrid->addWidget(CoordTitlez,2,4,1,1,Qt::AlignHCenter);
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
		CoordBox_id[i] = new QLineEdit();
		CoordBox_x[i] = new QLineEdit();
		CoordBox_y[i] = new QLineEdit();
		CoordBox_z[i] = new QLineEdit();
		CoordBox_id[i]->setPlaceholderText( "No input" );
		CoordBox_x[i]->setPlaceholderText( "No input" );
		CoordBox_y[i]->setPlaceholderText( "No input" );
		CoordBox_z[i]->setPlaceholderText( "No input" );
		CoordBox_id[i]->setReadOnly(true);
		CoordBox_x[i]->setReadOnly(true);
		CoordBox_y[i]->setReadOnly(true);
		CoordBox_z[i]->setReadOnly(true);
		CoordBox_id[i]->setFont(font);
		CoordBox_x[i]->setFont(font);
		CoordBox_y[i]->setFont(font);
		CoordBox_z[i]->setFont(font);
		CoordBox_id[i]->setFixedWidth(70);
		CoordBox_x[i]->setFixedWidth(70);
		CoordBox_y[i]->setFixedWidth(70);
		CoordBox_z[i]->setFixedWidth(70);
		SelectionDisplayGrid->addWidget(CoordLabel_n[i],i+3,0,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_id[i],i+3,1,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_x[i],i+3,2,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_y[i],i+3,3,1,1,Qt::AlignLeft);
		SelectionDisplayGrid->addWidget(CoordBox_z[i],i+3,4,1,1,Qt::AlignLeft);
	}
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
	DisplayPreferencesCheckBoxes[3]->setChecked(false);
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
	DisplayPreferencesCheckBoxes[8] = new QCheckBox("Myosin");
	DisplayPreferencesCheckBoxes[8]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[8] , SIGNAL(stateChanged(int)),this,SLOT(updateMyosinCheckBox(int)));

	//draw marking ellipses
	DisplayPreferencesCheckBoxes[9] = new QCheckBox("MarkingEllipses");
	DisplayPreferencesCheckBoxes[9]->setChecked(false);
	connect(DisplayPreferencesCheckBoxes[9] , SIGNAL(stateChanged(int)),this,SLOT(updateMarkingEllipseCheckBox(int)));

    ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[0],3,0,1,2,Qt::AlignLeft);  // display pipette
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[7],3,2,1,2,Qt::AlignLeft);  // display bounding box
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[1],4,0,1,2,Qt::AlignLeft);  // Net Forces
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[6],4,2,1,2,Qt::AlignLeft);  // Packing Forces
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[2],5,0,1,1,Qt::AlignLeft);  // Fixed Nodes
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[8],5,2,1,2,Qt::AlignLeft);  // display myosin box
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[3],6,0,1,2,Qt::AlignLeft); // Scale Bar
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[4],7,0,1,2,Qt::AlignLeft); // Display Peripodial Membrane
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[5],8,0,1,2,Qt::AlignLeft); // Display Columnar Layer
	ProjectDisplayOptionsGrid->addWidget(DisplayPreferencesCheckBoxes[9],7,2,1,2,Qt::AlignLeft); // Display Marked Ellipses
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

void  MainWindow::updateMyosinCheckBox(int s){
	cout<<"updating myosin checkbox, s: "<<s<<endl;
	if (s == 0){
		MainGLWidget->MyosinToDisplay = -1;
		MyosinComboBox->setEnabled(false);
		MainGLWidget->drawMyosinForces = false;
	}
	else{
		MainGLWidget->MyosinToDisplay = MyosinComboBox->currentIndex ();
		MyosinComboBox->setEnabled(true);
		MainGLWidget->drawMyosinForces = true;
		cout<<" drawMyosinForces: "<<MainGLWidget->drawMyosinForces<<endl;
		cout<<" MyosinToDisplay: "<<MainGLWidget->MyosinToDisplay <<endl;
	}
}

void  MainWindow::updateMarkingEllipseCheckBox(int s){
	if (s ==2 ){
		MainGLWidget->drawMarkingEllipses = true;
	}
	else{
		MainGLWidget->drawMarkingEllipses = false;
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

void MainWindow::updateMyosinComboBox(int s){
	MainGLWidget->MyosinToDisplay = s;
	MainGLWidget->update();
	//cout<<"Strain to display: "<<MainGLWidget->StrainToDisplay <<endl;
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
	cout<<"inside updateTimeText"<<endl;
	//round to two digits:
	QString timeStringInSec = QString::number(Sim01->currSimTimeSec);
	QString timeStringInHr = QString::number(Sim01->currSimTimeSec/3600.0, 'f', 2);
	QString timeStringInHrAEL = QString::number((Sim01->currSimTimeSec/3600.0 + 48), 'f', 2);
	timeStringInSec = "Simulation Time: " + timeStringInHr + " hr ("+timeStringInHrAEL+" hr AEL) - ( "+ timeStringInSec + " sec, " + QString::number(Sim01->timestep) + " steps )"; ;
	TimeTitle->setText(timeStringInSec);
	cout<<"finalised updateTimeText"<<endl;
}

void MainWindow::SelectedItemChange(){
	//cerr<<"Main window saw the selection change"<<endl;
    QString tmpstring = QString::fromStdString(MainGLWidget->SelectedItemName);
    NameBox->setText(tmpstring);
    for (int i = 0 ;i<nCoordBox; ++i){
    	 CoordBox_id[i]->setText ( "" );
    	 CoordBox_x[i]->setText ( "" );
    	 CoordBox_y[i]->setText ( "" );
    	 CoordBox_z[i]->setText ( "" );
    	 CoordBox_id[i]->setEnabled(false);
    	 CoordBox_x[i]->setEnabled(false);
    	 CoordBox_y[i]->setEnabled(false);
    	 CoordBox_z[i]->setEnabled(false);
    	 if ((signed int)MainGLWidget->SelectedPos.size()>i*3){
    		 CoordBox_id[i]->setText ( MainGLWidget->SelectedId[i] );
			 CoordBox_x[i]->setText ( MainGLWidget->SelectedPos[i*3] );
			 CoordBox_y[i]->setText ( MainGLWidget->SelectedPos[i*3+1] );
			 CoordBox_z[i]->setText ( MainGLWidget->SelectedPos[i*3+2] );
			 CoordBox_id[i]->setEnabled(true);
			 CoordBox_x[i]->setEnabled(true);
			 CoordBox_y[i]->setEnabled(true);
			 CoordBox_z[i]->setEnabled(true);
		}
    }
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
	ElementSelectBox->blockSignals(true);
	ElementSelectBox->setValidator( new QIntValidator(0, Sim01->Elements.size()-1, this) );
	ElementSelectBox->setPlaceholderText ( QString("# 0 to %1").arg(Sim01->Elements.size()-1) );
	ElementSelectBox->setText("");
	ElementSelectBox->blockSignals(false);
}

void MainWindow::ManualNodeSelectionReset(){
	NodeSelectBox->blockSignals(true);
	NodeSelectBox->setValidator( new QIntValidator(0, Sim01->Nodes.size()-1, this) );
	NodeSelectBox->setPlaceholderText ( QString("# 0 to %1").arg(Sim01->Nodes.size()-1) );
	NodeSelectBox->setText("");
	NodeSelectBox->blockSignals(false);
}

void MainWindow::timerSimulationStep(){
    //cout<<"Called the function via timer"<<endl;
    //for (int j=0; j<3; ++j){
    //    cout<<"Node 2025 at circmference: "<<Sim01->Nodes[2025]->atCircumference<<" fixed pos "<<j<<" "<<Sim01->Nodes[2025]->FixedPos[j]<<endl;
    //}
	bool 	automatedSave = false;
	bool 	analyseResults = true;
	bool 	slowstepsOnDisplay = false;
	bool 	slowstepsOnRun = false;
	int 	slowWaitTime = 10;
	if (Sim01->DisplaySave){
		if (Sim01->timestep == 0){
			if (analyseResults){
				cout<<"Analysin results"<<endl;
				analyser01 = new Analysis(3, Sim01->nElements, Sim01->saveDirectoryToDisplayString, Sim01->Nodes);
			}
			if( automatedSave){
				Sim01->assignTips();
				//MainGLWidget->updateToTopView(); //display tissue from the top view
				//MainGLWidget->updateToFrontView(); //display tissue from the front view
				MainGLWidget->updateToPerspectiveView(); //display tissue from the tilted view
				MainGLWidget->drawSymmetricity = true; //hide symmetric
				//MainGLWidget->drawPeripodialMembrane= false;
				MainGLWidget->PerspectiveView = false; //switch to orthogoanal view type.
				Sim01->saveImages = true;
				Sim01->saveDirectory = Sim01->saveDirectoryToDisplayString;
			}
		}


		//Sim01->PipetteSuction = true;
		//MainGLWidget->displayPipette = true;
        //Sim01->pipetteInnerRadius = 12.4;
        //Sim01->pipetteCentre[0] = 0;
        //Sim01->pipetteCentre[1] = 0;
        //Sim01->pipetteCentre[2] = 25;
        //Sim01->pipetteThickness = 2.8;
		//double pressureArray[117]= {2.5,5,5,7.5,7.5,10,10,12.5,12.5,15,15,17.5,17.5,20,20,22.5,22.5,24.525,24.525,27.5,27.5,30,30,32.5,32.5,35,35,37.5,37.5,40,40,42.5,42.5,45,45,47.5,47.5,49.05,49.05,52.5,52.5,55,55,57.5,57.5,60,60,62.5,62.5,65,65,67.5,67.5,70,70,73.575,73.575,75,75,77.5,77.5,80,80,82.5,82.5,85,85,87.5,87.5,90,90,92.5,92.5,95,95,98.01,98.01,100,100,102.5,102.5,105,105,107.5,107.5,110,110,112.5,112.5,115,115,117.5,117.5,120,120,122.625,122.625,125,125,127.5,127.5,130,130,132.5,132.5,135,135,137.5,137.5,140,140,142.5,142.5,145,145,147.15,147.15};
		if(!Sim01->reachedEndOfSaveFile){
			Sim01->updateOneStepFromSave();
			if (Sim01->timestep >= 0){
				double boundingBoxLength = Sim01->boundingBox[1][0] - Sim01->boundingBox[0][0];
				double boundingBoxWidth  = Sim01->boundingBox[1][1] - Sim01->boundingBox[0][1];


				/*Sim01->SuctionPressure[2] = pressureArray[interatorForPressure];
				double zMax = -10000;
				int idMax = -10;
				for (vector<Node*>::iterator itNode = Sim01->Nodes.begin(); itNode < Sim01->Nodes.end(); ++itNode){
					if ((*itNode)->Position[2] > zMax){
						zMax = (*itNode)->Position[2];
						idMax = (*itNode)->Id;
					}
				}
				cout<<"Pipette suction: "<<Sim01->SuctionPressure[2]<<" max suction: "<<zMax<<" from node "<<idMax<<endl;
				interatorForPressure++;
				 */



				if (analyseResults){
					analyser01->calculateBoundingBoxSizeAndAspectRatio(Sim01->currSimTimeSec,boundingBoxLength,boundingBoxWidth);
					analyser01->calculateContourLineLengthsDV(Sim01->Nodes);
					analyser01->findApicalKinkPointsDV(Sim01->currSimTimeSec,Sim01->boundingBox[0][0], boundingBoxLength, Sim01->Nodes);
					analyser01->calculateTissueVolumeMap(Sim01->Elements, Sim01->currSimTimeSec,Sim01->boundingBox[0][0],Sim01->boundingBox[0][1],boundingBoxLength,boundingBoxWidth);
				}
			}
			Sim01->calculateDVDistance();
			for (int a = 0; a<0; a++){ //11 for 6 hours with 1800 sec time step
				Sim01->updateOneStepFromSave();
				Sim01->calculateDVDistance();
			}
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

			/*for (int i =0 ; i<Sim01->nElements; i++){
				if (i == 345 || i == 174 || i == 302 || i == 1041){
					cout.precision(6);
					cout<<"After display frame: Time: "<<Sim01->currSimTimeSec<<" Element : "<<Sim01->Elements[i]->Id<<" grownVolume: "<<Sim01->Elements[i]->GrownVolume<<" reference volume: "<<Sim01->Elements[i]->ReferenceShape->Volume<<endl;
				}
			}*/
			//sleep(60);
			//spitting coordinates:
			//Sim01->CoordinateDisplay();
			//Sim01->TissueAxisPositionDisplay();
		}
		else{
			//close();
			if (automatedSave){
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
			/*for(vector<ShapeBase*>::iterator itElement=Sim01->Elements.begin(); itElement<Sim01->Elements.end(); ++itElement){
				if((*itElement)->tissuePlacement == 1){ //apical element
						double* growthRate =(*itElement)->getGrowthRate();
						gsl_matrix* elementFG = (*itElement)->getFg();
						double detFg = (*itElement)->determinant3by3Matrix(elementFG);
						double* ReletivePos = new double[2];
						(*itElement)->getRelativePosInBoundingBox(ReletivePos);
						cout<<" TheGrowthRateOfElement: "<<(*itElement)->Id<<" "<<growthRate[0]<<" "<<growthRate[1]<<" "<<growthRate[2]<<" curent detFg: "<<detFg<<" grownVolume: "<<(*itElement)->GrownVolume<<" ReferenceVolume "<<(*itElement)->ReferenceShape->Volume<<" Relative pos: "<<ReletivePos[0]<<" "<<ReletivePos[1]<<endl;
						gsl_matrix_free(elementFG);
						delete[] ReletivePos;
				}
			}*/
            //Sim01->writeMeshRemovingAblatedRegions();
            double durationClock = ( std::clock() - simulationStartClock ) / (double) CLOCKS_PER_SEC;
            double durationTime = std::difftime(std::time(0), simulationStartTime);
            cout<<"Simulation time: "<<durationTime<<" sec, Simulation clock: "<<durationClock<<" sec"<<endl;
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
