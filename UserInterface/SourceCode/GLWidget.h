/*
 * GLWidget.h
 *
 *  Created on: 20 Mar 2014
 *      Author: melda
 */

#ifndef GLWIDGET_H_
#define GLWIDGET_H_

#include <QGLWidget>
#include <QtOpenGL>
//#include <iostream>
#include <vector>
#include <string.h>
#include <math.h>
using namespace std;

#include "../TissueFolding/SourceCode/Simulation.h"

 class GLWidget : public QGLWidget
 {
     Q_OBJECT

 public:
     GLWidget(QWidget *parent = 0);
     ~GLWidget();

     Simulation* Sim01;
     QSize minimumSizeHint() const;
     QSize sizeHint() const;
     bool ItemSelected;
     string SelectedItemName;
     int SelectedItemIndex;
     bool DisplayStrains;
     float DisplayStrainRange[2];
     int StrainToDisplay;
     bool DisplayPysProp;
     int PysPropToDisplay;
     float DisplayPysPropRange[5][2];	//current range
     float DisplayPysPropBounds[5][4];  //the minimum and maximum they can get
     vector <QString> SelectedPos;
     bool drawNormals;
     bool drawNetForces;
     bool drawVelocities;


 signals:
 	 void SelectedItemChanged();

 protected:
     void initializeGL();
     void paintGL();
     void resizeGL(int width, int height);
     void mousePressEvent(QMouseEvent *event);
     void mouseReleaseEvent(QMouseEvent *event);
     void mouseMoveEvent(QMouseEvent *event);
     void wheelEvent(QWheelEvent *event);
     void ObjectSelection(QPoint LastPos);
     void resetItemSelectionInfo();
     void findElement();
     void GetColourOfPoint(QPoint LastPos);
     void DrawForPicking ();
     void generate3DObject();

     //Element drawing functions
     void drawElement(int i, bool picking);
     void drawReferenceElement(int i);
     void drawPrism(int i);
     void drawPrismLateral(int i);
     void drawPrismForPicking(int i);
     void drawReferencePrism(int i);
     void drawReferencePrismLateral(int i);
     void fillItemSelectionInfo(int i);


 private:
     QPoint lastPos;
     QPoint InitialClickPos;
     int 	MouseButton;
     float 	obj_pos[3];
     float 	aspectratio;
     int 	PickedColour[4];
     float 	ReferenceLineThickness, MainShapeLineThickness;
     bool 	DisplayFixedNodes;

    // QColor qtGreen;
     QColor 	qtPurple;
     double 	Qcurr[4], Qlast[4];
     float 		MatRot[16];
     bool 		checkPickedColour(int* ElementColour);
     void 		rotateByQuaternians(double* Qrot);
     void 		normaliseCurrentRotationAngle (double* Qrot);
     void 		rotateCurrentRotationQuaternian(double* Qrot);
     void 		rotateMatrix();
     void 		getDisplayColour(float* StrainColour, float Data);
     float** 	getElementColourList(int i);
     void 		drawColourbar();
     void 		drawFixedNodes();
     void		drawAxesArrows();
     void 		drawNormalToPrism(int i);
     void 		drawNormalToReferencePrism(int i);
     void 		drawTissueCoordSystemPrism(int i);
     void 		drawNormalToPrismLateral(int i);
     void 		drawNormalToReferencePrismLateral(int i);
     void 		drawForces();
     void 		drawNodeVelocities();
     void 		drawArrow3D(double* pos, double* endPoint, double r, double g, double b);
 };


#endif /* GLWIDGET_H_ */
