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

#include "Simulation.h"
#include "Analysis.h"


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
    void 	xClipChange(int);
    void 	yClipChange(int);
    void 	zClipChange(int);

signals:
   void deselectedElement();                                            // emitted when an element is deselected
   void lookingAtNewElement(std::unique_ptr<ShapeBase> *element);       // emitted when the selected element changes

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
     * @brief This layout contains the widgets that handle the display of basic element information (name, associated nodes, their positions, etc) and the selection of nodes/elements by providing a manual index.
     */
    ElementBasicDisplay *ElementProps;

    /**
     * @brief Allows the user to select a particular element property (growth, force, etc) and displays the value of this property for the currently selected element.
     * 
     * Auto-updates the display based on the property selected, and disables when an element is deselected.
     */
    ElementPropertySelection *PropertySelection;

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
