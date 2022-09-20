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
#include <array>
using namespace std;

#include "../TissueFolding/SourceCode/Simulation.h"
#include "../TissueFolding/SourceCode/Analysis.h"

/*! GLWidget class */
 class GLWidget : public QGLWidget
 {
     Q_OBJECT

 public:
     GLWidget(QWidget *parent = 0);
     ~GLWidget();

     Simulation* Sim01;
     Analysis* 	analyser01;
     QSize		minimumSizeHint() const;
     QSize 		sizeHint() const;
     void 		manualElementSelection(int i);
     void 		manualNodeSelection(int i);
     void		updateClipping();
     void		updateToTopView();
     void		updateToFrontView();
     void		updateToSideView();
     void		updateToPerspectiveView();
     void		updateToCloneCrossView();
     void		drawPointsForDisplay();
     void		drawEnclosingShell();
     bool 		ItemSelected;
     string 	SelectedItemName;
     int 		SelectedItemIndex;
     bool 		ManualNodeSelection;
     int	 	ManualSelectedNodeId;
     bool 		DisplayFixedNodes;
     bool 		DisplayStrains;
     float 		DisplayStrainRange[2];
     int 		StrainToDisplay;
     bool 		DisplayPysProp;
     int 		PysPropToDisplay;
     float 		DisplayPysPropRange[8][2];	//current range
     float 		DisplayPysPropBounds[8][4];  //the minimum and maximum they can get
     int		DisplayPysPropDecimals[8];	//the decimal points for spin boxes
     float 		DisplayPysPropSteps[8];	//current step of the spinbox
     vector <QString> SelectedPos;
     vector <QString> SelectedId;
     bool 		drawNetForces;
     bool		drawMarkingEllipses;
     bool		drawGrowthRedistribution;
     bool 		drawNodeBinding;
     bool 		drawPackingForces;
     bool 		drawTissueScaleBar;
     bool 		drawPeripodialMembrane;
     bool 		drawColumnar;
     bool		drawLumen;
     bool 		PerspectiveView;
     bool		displayBoundingBox;
     bool       displayPipette;
     double  	xClip, yClip, zClip;
     bool 		drawSymmetricity;
     size_t		currNodeNumber;
     float 		obj_pos[3];

 signals:
 	 void SelectedItemChanged(bool found_element=true);
 	 void NeedToClearManualElementSelection();
 	 void NeedToClearManualNodeSelection();

 protected:
     void initializeGL();
     void 	paintGL();
     void 	resizeGL(int width, int height);
     void 	mousePressEvent(QMouseEvent *event);
     void 	mouseReleaseEvent(QMouseEvent *event);
     void 	mouseMoveEvent(QMouseEvent *event);
     void 	wheelEvent(QWheelEvent *event);
     void 	ObjectSelection(QPoint LastPos);
     void 	resetItemSelectionInfo(int source);
     bool 	findElement();
     bool 	findElement(int i);
     bool 	findNode(int i);
     void 	getColourOfPoint(QPoint LastPos);
     void 	drawForPicking ();
     void 	generate3DObject();
     void 	initialiseNodeColourList();
     void	reInitialiseNodeColourList(size_t oldNodeNumber);

     //Element drawing functions
     bool 	checkIfDrawingElement(int i);
     std::array<bool,7> 	checkIfDrawingElementSymmetric(int i);
     bool 	checkIfDrawingNode(int i);
     void 	drawElement(size_t i, bool picking);
     void 	highlightElement(int i);
     void 	highlightNode(int i);
     void 	drawReferenceElement(int i);
     void 	drawPrism(int i, bool symmetricX, bool symmetricY, bool symmetricZ);
     void 	drawTriangle(int i);
     void 	drawPrismForPicking(int i);
     void 	drawTriangleForPicking(int i);
     void 	drawReferencePrism(int i);
     void 	drawReferenceTriangle(int i);
     void 	highlightPrism(int i);
     void 	highlightTriangle(int i);
     void 	fillItemSelectionInfo(int i);



 private:
     QPoint lastPos;
     QPoint InitialClickPos;
     int 	MouseButton;
     
     float 	aspectratio;
     int 	PickedColour[4];
     float 	ReferenceLineThickness, MainShapeLineThickness;
     float** NodeColourList;
     float orthoViewLimits[6];
    // QColor qtGreen;
     QColor 	qtPurple;
     double 	Qcurr[4], Qlast[4];
     float 		MatRot[16];
     bool 		checkPickedColour(std::array<int,3> ElementColour);
     void 		rotateByQuaternians(double* Qrot);
     void 		normaliseCurrentRotationAngle (double* Qrot);
     void 		rotateCurrentRotationQuaternian(double* Qrot);
     void 		rotateMatrix();
     void 		getDisplayColour(float* StrainColour, float Data);
     void 		getVelocityColour(float* OutputColour, float Data);
     void 		getForceColour(float* OutputColour, float Data);
     void 		getConcentrationColour(float* OutputColour, float concentration);
     void 		constructNodeColourList();
     float** 	getElementColourList(int i);
     void 		drawColourbar();
     void 		drawFixedNodes();
     void		drawBoundNodes();
     void		drawAxesArrows();
     void 		drawScaleBar();
     void 		drawForces();
     void 		drawPackForces();
     void 		drawBoundingBox();
     void 		drawPipette();
     void		drawAFMBead();
     void 		drawArrow3D(const std::array<double,3> &pos, const std::array<double,3> &endPoint, double r, double g, double b);
 };


#endif /* GLWIDGET_H_ */
