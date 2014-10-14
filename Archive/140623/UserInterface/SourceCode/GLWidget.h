/*
 * GLWidget.h
 *
 *  Created on: 20 Mar 2014
 *      Author: melda
 */

#ifndef GLWIDGET_H_
#define GLWIDGET_H_

#include <QGLWidget>
//#include <iostream>
#include <vector>
#include <string.h>
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
     vector <QString> SelectedPos;



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
     void drawPrismForPicking(int i);
     void drawReferencePrism(int i);
     void fillItemSelectionInfo(int i);


 private:
     QPoint lastPos;
     QPoint InitialClickPos;
     int MouseButton;
     float obj_pos[3];
     float Rotate[3];
     float aspectratio;
     int PickedColour[4];
     float ReferenceLineThickness, MainShapeLineThickness;
    // QColor qtGreen;
     QColor qtPurple;

     bool checkPickedColour(int* ElementColour);
 };


#endif /* GLWIDGET_H_ */
