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
class GLWidget;

#include "../TissueFolding/SourceCode/Simulation.h"

class MainWindow : public QMainWindow
 {
     Q_OBJECT

 public:
    MainWindow(Simulation* SIm01);
    ~MainWindow();
    QGraphicsScene	*MainScene;
    QVBoxLayout		*ControlPanelMainHBox;
    QGridLayout		*MainGrid;
    GLWidget 		*MainGLWidget;
    Simulation* Sim01;

public slots:
    void SelectedItemChange();
    void timerSimulationStep();
    void updateStrain(int);
    void updateStrainCheckBox(int);
    void updateStrainSpinBoxes(double);
    void updatePysProp(int s);
    void updatePysCheckBox(int);
    void updatePysPropSpinBoxes(double d);
    void updateNormalCheckBox(int);
    void updateNetForceCheckBox(int);
    void updateVelocityCheckBox(int);
//signals:
 //   void StrainComboBoxCanged();
 private:
    void setViewBackgroundColour();
    void generateControlPanel();
    void setUpView();
    void setUpGLWidget();
    void setUpCentralWidget();
    void setUpSelectionDisplayGrid(QGridLayout *SelectionDisplayGrid);
    void setCoordBoxes(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid);
    void setItemSelectionTitles(QFont font, QFont boldFont, QGridLayout *SelectionDisplayGrid);
    void setStrainDisplayMenu(QGridLayout *DisplayOptionsGrid);
    void setPysPropDisplayMenu(QGridLayout *DisplayOptionsGrid);
    void setDisplayPreferences(QGridLayout *SelectionDisplayGrid);
    void takeScreenshot();
    QWidget		*CentralWidget;
    QLineEdit 	*NameBox;
    QTimer *timer;
    int nCoordBox;
    QLineEdit 		*CoordBox_x[6];
    QLineEdit 		*CoordBox_y[6];
    QLineEdit 		*CoordBox_z[6];
    QLabel  		*CoordLabel_n[6];
    QCheckBox		*DisplayCheckBoxes[2];
    QComboBox   	*StrainComboBox;
    QDoubleSpinBox	*StrainSpinBoxes[2];
    QComboBox   	*PysPropComboBox;
    QDoubleSpinBox 	*PysPropSpinBoxes[2];
    QGroupBox   	*ColourCodingBox;
    QCheckBox		*DisplayPreferencesCheckBoxes[3];
    QLabel			*SimTime;
 };


#endif /* MAINWINDOW_H_ */
