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

/*! GLWidget class */
 class GLWidget : public QGLWidget
 {
     Q_OBJECT

 public:
     GLWidget(QWidget *parent = 0);
     ~GLWidget();

     Simulation* Sim01;
     QSize		minimumSizeHint() const;
     QSize 		sizeHint() const;
     void 		manualElementSelection(int i);
     void 		manualNodeSelection(int i);
     void		updateClipping();
     bool 		ItemSelected;
     string 	SelectedItemName;
     int 		SelectedItemIndex;
     bool 		ManualNodeSelection;
     int	 	ManualSelectedNodeId;
     bool 		DisplayStrains;
     float 		DisplayStrainRange[2];
     int 		StrainToDisplay;
     bool 		DisplayPysProp;
     int 		PysPropToDisplay;
     float 		DisplayPysPropRange[5][2];	//current range
     float 		DisplayPysPropBounds[5][4];  //the minimum and maximum they can get
     int		DisplayPysPropDecimals[5];	//the decimal points for spin boxes
     float 		DisplayPysPropSteps[5];	//current step of the spinbox
     vector <QString> SelectedPos;
     vector <QString> SelectedId;
     bool 		drawTissueCoordinates;
     bool 		drawNetForces;
     bool 		drawPackingForces;
     bool 		drawVelocities;
     bool 		drawTissueScaleBar;
     bool 		drawPeripodialMembrane;
     bool 		drawColumnar;
     bool 		PerspectiveView;
     bool		displayBoundingBox;
     double  	xClip, yClip, zClip;

 signals:
 	 void SelectedItemChanged();
 	 void NeedToClearManualElementSelection();
 	 void NeedToClearManualNodeSelection();

 protected:
     void initializeGL();
     void paintGL();
     void resizeGL(int width, int height);
     void mousePressEvent(QMouseEvent *event);
     void mouseReleaseEvent(QMouseEvent *event);
     void mouseMoveEvent(QMouseEvent *event);
     void wheelEvent(QWheelEvent *event);
     void ObjectSelection(QPoint LastPos);
     void resetItemSelectionInfo(int source);
     void findElement();
     bool findElement(int i);
     bool findNode(int i);
     void getColourOfPoint(QPoint LastPos);
     void drawForPicking ();
     void generate3DObject();
     void initialiseNodeColourList();

     //Element drawing functions
     bool checkIfDrawingElement(int i);
     bool checkIfDrawingNode(int i);
     void drawElement(int i, bool picking);
     void highlightElement(int i);
     void highlightNode(int i);
     void drawReferenceElement(int i);
     void drawPrism(int i);
     void drawTriangle(int i);
     void drawPrismForPicking(int i);
     void drawTriangleForPicking(int i);
     void drawReferencePrism(int i);
     void drawReferenceTriangle(int i);
     void highlightPrism(int i);
     void highlightTriangle(int i);
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
     float** NodeColourList;
     float orthoViewLimits[6];
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
     void 		getVelocityColour(float* OutputColour, float Data);
     void 		getForceColour(float* OutputColour, float Data);
     void 		constructNodeColourList();
     float** 	getElementColourList(int i);
     void 		drawColourbar();
     void 		drawFixedNodes();
     void		drawAxesArrows();
     void 		drawScaleBar();
     void 		drawTissueCoordSystemPrism(int i);
     void 		drawTissueCoordSystemTetrahedron(int i);
     void 		drawTissueCoordSystemTriangle(int i);
     void 		drawForces();
     void 		drawPackForces();
     void 		drawNodeVelocities();
     void 		drawBoundingBox();
     void 		drawArrow3D(double* pos, double* endPoint, double r, double g, double b);
 };


#endif /* GLWIDGET_H_ */
