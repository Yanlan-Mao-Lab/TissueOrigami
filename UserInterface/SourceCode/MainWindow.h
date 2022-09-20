/*
 * MainWindow.h
 *
 *  Created on: 18 Mar 2014
 *      Author: melda
 */

#ifndef MAINWINDOW_H_
#define MAINWINDOW_H_

#include <QtGui>
#include <QtWidgets>
#include <QGraphicsWidget>
#include <QMainWindow>
#include <QLineEdit>
#include <QTimer>
#include <QSlider>
//Ubuntu version:
//#include <time.h>
#include <ctime>
class GLWidget;

#include "../TissueFolding/SourceCode/Simulation.h"
#include "../TissueFolding/SourceCode/Analysis.h"

#include "ElementPropertiesUI.h"


using namespace std;

class MainWindow : public QMainWindow
 {
     Q_OBJECT

 public:
    MainWindow(Simulation* Sim01);
    ~MainWindow();
    QGraphicsScene	*MainScene;
    QVBoxLayout		*ControlPanelMainHBox;
    QGridLayout		*MainGrid;
    GLWidget 		*MainGLWidget;
    Simulation* Sim01;
    Analysis* analyser01;
	int interatorForPressure;


public slots:
    void 	SelectedItemChange(bool element_found);
    void	manualNodeSelection(const QString &);
    void 	manualElementSelection(const QString &);
    void 	ManualElementSelectionReset();
    void 	ManualNodeSelectionReset();
    void	testAdhesionsAndCurveConstruction();
    void 	timerSimulationStep();
    void 	updateStrain(int);
    void 	updateStrainCheckBox(int);
    void 	updateStrainSpinBoxes();
    void 	updatePysProp(int s);
    void 	updatePysCheckBox(int);
    void 	updatePysPropSpinBoxes();
    void    updateDisplayPipette(int);
    void 	updateNetForceCheckBox(int);
    void	updateMarkingEllipseCheckBox(int);
    void   	updateGrowthRedistributionCheckBox(int s);
    void	updateDrawNodeBindingCheckBox(int s);
    void 	updatePackingForceCheckBox(int);
    void 	updateFixedNodesCheckBox(int);
    void  	updateScaleBarCheckBox(int);
    void  	updatePeripodialDisplayCheckBox(int s);
    void  	updateColumnarLayerDisplayCheckBox(int s);
    void    updateLumenDisplayCheckBox(int s);
    void  	updateBoundingBoxCheckBox(int s);
    void  	updateOrthagonalPerspectiveViewToggle();
    void  	updateToTopView();
    void  	updateToFrontView();
    void  	updateToSideView();
    void  	updateToPerspectiveView();
    void  	updateDrawSymmetricityViewToggle();
    void    updateSelectElementPropertyDisplay(int option);
    void 	xClipChange(int);
    void 	yClipChange(int);
    void 	zClipChange(int);

//signals:
 //   void StrainComboBoxCanged();
 private:
    void setViewBackgroundColour();
    void generateControlPanel();
    void setUpView();
    void setUpGLWidget();
    void setUpCentralWidget();
    void setUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid);
    void setUpProjectDisplayOptionGrid(QGridLayout *ProjectDisplayOptionsGrid);
    void setUpViewOptionsGrid(QGridLayout *ViewOptionsGrid);
    void setSelectionByIdSection(QFont font, QGridLayout *SelectionDisplayGrid);
    void setCoordBoxes(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid);
    void setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid);
    void setStrainDisplayMenu(QGridLayout *DisplayOptionsGrid);
    void setPysPropDisplayMenu(QGridLayout *DisplayOptionsGrid);
    void setDisplayPreferences(QGridLayout *SelectionDisplayGrid);
    void updateTimeText();
    void takeScreenshot();
    QWidget		*CentralWidget;

    /**
     * @brief ElementPropertiesUI handles the display of basic element information, and selection of nodes and elements by passing their ID
     * 
     * This panel renders:
     * - The "selected item properties" header
     * - The node selection and element selection boxes
     * - The printout of node information associated to a given element
     */
    ElementPropertiesUI *ElementProps;

    QTimer *timer;
    int nCoordBox;
    QCheckBox		*DisplayCheckBoxes[2];
    QComboBox   	*StrainComboBox;
    QDoubleSpinBox	*StrainSpinBoxes[2];
    QComboBox   	*PysPropComboBox;
    QDoubleSpinBox 	*PysPropSpinBoxes[2];
    QGroupBox   	*ColourCodingBox;
    QCheckBox		*DisplayPreferencesCheckBoxes[12];
    QLabel			*SimTime;
    QPushButton		*PerspectiveButton;
    QPushButton		*TopViewButton;
    QPushButton		*FrontViewButton;
    QPushButton		*SideViewButton;
    QPushButton		*PerspectiveViewButton;
    QPushButton		*SymmetricityDisplayButton;
    QSlider			*ClippingSliders[3];
    QLineEdit 		*TimeBox;
    QLabel 			*TimeTitle;
    //MacOS version:
    std::clock_t    simulationStartClock;	//simulation clock as in cpu clock time (will be the same regardless of using single or multi-processors. sleep during process etc.
    std::time_t 	simulationStartTime;	//simulation time in real time, given in seconds

    bool            displayedSimulationLength;
 };


#endif /* MAINWINDOW_H_ */
