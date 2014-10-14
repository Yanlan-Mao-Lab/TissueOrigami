#include <iostream>
#include <math.h>
#include <vector>
#include <string.h>

#include <QtGui>
#include <QtOpenGL>
#include <QWheelEvent>

#include "GLWidget.h"
#include "../TissueFolding/SourceCode/ShapeBase.h"

using namespace std;




 //GLWidget::GLWidget(QWidget *parent)
  //   : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
 GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
 {
	 obj_pos[0] = -4.0f;
	 obj_pos[1] =  2.0f;
	 obj_pos[2] =  20.0f;
	 Rotate[0] = -50.0;
	 Rotate[1] = -2.0;
	 Rotate[2] = 80.0;
	 PickedColour[0] = 255;
	 PickedColour[1] = 255;
	 PickedColour[2] = 255;
	 PickedColour[3] = 1;
	 ItemSelected = false;
	 SelectedItemIndex = -1;
	 SelectedItemName = "N/A";
     qtPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
     aspectratio =1.0;
     ReferenceLineThickness = 1.0;
     MainShapeLineThickness = 1.0;

  	 setSizePolicy(QSizePolicy ::Expanding , QSizePolicy ::Expanding );
 }

 GLWidget::~GLWidget()
 {
 }

 QSize GLWidget::minimumSizeHint() const
 {
     return QSize(50, 50);
 }

 QSize GLWidget::sizeHint() const
 {
     return QSize(6000, 6000);
 }

 void GLWidget::initializeGL()
 {
	 qglClearColor(QColor::fromRgbF(1, 1, 1, 1));
     glEnable(GL_DEPTH_TEST);
     setAutoBufferSwap(true);
     GLfloat LineRange[2];
     glGetFloatv(GL_LINE_WIDTH_RANGE,LineRange);
     ReferenceLineThickness = (LineRange[1]-LineRange[0])/2.0;
     MainShapeLineThickness = (LineRange[1]-LineRange[0])/10.0;
 }

 void GLWidget::paintGL()
 {
	 QSize viewport_size = size();
	 glViewport(0, 0, viewport_size.width(), viewport_size.height());

	 glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();
	 glFrustum(-1*aspectratio, 1*aspectratio, -1, 1, 1, 500); // near and far match your triangle Z distance

	 glMatrixMode(GL_MODELVIEW);

	 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	 glLoadIdentity();

	 glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );

	 glTranslatef( Sim01->SystemCentre[0], Sim01->SystemCentre[1], -Sim01->SystemCentre[2]);
	 glRotatef(Rotate[0], 1.0, 0.0, 0.0);
	 glRotatef(Rotate[1], 0.0, 1.0, 0.0);
	 glRotatef(Rotate[2], 0.0, 0.0, 1.0);
	 glTranslatef( -Sim01->SystemCentre[0], -Sim01->SystemCentre[1], Sim01->SystemCentre[2]);

	 //cout<<"obj_pos: "<<obj_pos[0]<<" "<<obj_pos[1]<<" "<<obj_pos[2]<<endl;
	 //cout<<"Rotate: "<<Rotate[0]<<" "<<Rotate[1]<<" "<<Rotate[2]<<endl;
	 for (int i =0; i<Sim01->Elements.size();i++){
		 drawElement(i,false);
	 }
	 if (ItemSelected){
		 //cout<<"item is selected calling reference shape drawing!"<<endl;
		 drawReferenceElement(SelectedItemIndex);
	 }
	 swapBuffers();
 }

 void GLWidget::drawElement(int i, bool picking){
	 int ShapeType = Sim01->Elements[i]->getShapeType();
	 if (ShapeType == 1 ){
		 if (picking){
			drawPrismForPicking(i);
		 }
		 else{
			 drawPrism(i);
		 }
	 }
 }

 void  GLWidget::drawReferenceElement(int i){
	 int ShapeType = Sim01->Elements[i]->getShapeType();
	 if (ShapeType == 1 ){
		 drawReferencePrism(i);
	 }
 }

 void GLWidget::drawPrism(int i){
	//Drawing the surfaces
	const int nTriangle = 8; //top + bottom + 2 for each side.
	int TriangleConnectivity[nTriangle][3] = {{0,1,2},{3,4,5},{0,2,3},{2,3,5},{0,1,3},{1,3,4},{1,2,5},{1,5,4}};
	const int nLineStrip = 12;
	int BorderConnectivity[nLineStrip] = {0,2,5,3,0,1,4,3,5,4,1,2};

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
	glBegin(GL_TRIANGLES);
		glColor3f(0.75,1.0,1.0);
		for (int j =0; j<nTriangle;++j){
			for (int k =0; k<3; ++k){
				int pointId = TriangleConnectivity[j][k];
				float x = Sim01->Elements[i]->Positions[pointId][0];
				float y = Sim01->Elements[i]->Positions[pointId][1];
				float z = Sim01->Elements[i]->Positions[pointId][2];
				glVertex3f( x, y, z);
			}
		}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);

 	glLineWidth(MainShapeLineThickness);

	//Drawing the borders
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3f(0,0,0);
	glBegin(GL_LINE_STRIP);
		for (int j =0; j<nLineStrip;++j){
			int pointId = BorderConnectivity[j];
			//cout<<j<<" point id: "<<pointId<<endl;
			float x = Sim01->Elements[i]->Positions[pointId][0];
			float y = Sim01->Elements[i]->Positions[pointId][1];
			float z = Sim01->Elements[i]->Positions[pointId][2];
			//cout<<"x,y,z: "<<x<<" "<<y<<" "<<z<<endl;
			glVertex3f( x, y, z);
		}
	glEnd();
 }

 void GLWidget::drawReferencePrism(int i){
	const int nLineStrip = 12;
 	int BorderConnectivity[nLineStrip] = {0,2,5,3,0,1,4,3,5,4,1,2};

 	glLineWidth(ReferenceLineThickness);
 	//Drawing the borders
 	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 	glColor3f(1,0,0);
 	glBegin(GL_LINE_STRIP);
 		for (int j =0; j<nLineStrip;++j){
 			int pointId = BorderConnectivity[j];
 			double** pos = Sim01->Elements[i]->getReferencePos();
 			//cout<<j<<" point id: "<<pointId<<endl;
 			float x = pos[pointId][0];
 			float y = pos[pointId][1];
 			float z = pos[pointId][2];
 			//cout<<"x,y,z: "<<x<<" "<<y<<" "<<z<<endl;
 			glVertex3f( x, y, z);
 		}
 	glEnd();
  }

 void GLWidget::resizeGL(int width, int height)
 {
     int  side = qMin(width, height);
     glViewport((width - side) / 2, (height - side) / 2, width, height);
     aspectratio = float(width) /float (height);
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     glMatrixMode(GL_MODELVIEW);
 }

 void GLWidget::mousePressEvent(QMouseEvent *event)
 {
     lastPos = event->pos();
     if(event->button()==Qt::LeftButton){
         MouseButton = 0;
         InitialClickPos = event->pos();
         //ObjectSelection(lastPos);
     }
     else if(event->button()==Qt::MidButton){
    	 MouseButton = 1;
     }
     else if(event->button()==Qt::RightButton){
         MouseButton = 2;
     }
 }

 void GLWidget::mouseReleaseEvent(QMouseEvent *event)
  {
      lastPos = event->pos();
      if(event->button()==Qt::LeftButton){
          MouseButton = 0;
          int dx = InitialClickPos.x() - lastPos.x();
          if (dx<0) {
        	  dx = dx*(-1);
          }
          if (dx<5){
        	  ObjectSelection(lastPos);
          }
      }
  }

 void GLWidget::mouseMoveEvent(QMouseEvent *event)
 {
	 if (MouseButton==0){
		 int dx = event->x() - lastPos.x();
		 float speed = 0.3;
		 Rotate[2] += dx*speed; //movement in x rotates the scene around x axis (move mouse up to push the object up)
		 while (Rotate[2] >180) {Rotate[2] -=360;}
		 while (Rotate[2] <-180) {Rotate[2] +=360;}
		 lastPos = event->pos();
		 updateGL();
	 }
	 else if(MouseButton==1){
    	 int dx = event->x() - lastPos.x();
    	 int dy = event->y() - lastPos.y();
    	 float speed = 0.01;
    	 obj_pos[0] += dx*speed;
    	 obj_pos[1] += -dy*speed;
    	 //cerr<<"right button "<<dx<<" "<<dy<<" "<<obj_pos[0]<<" "<<obj_pos[1]<<endl;
    	 lastPos = event->pos();
    	 updateGL();
     }
	 else if(MouseButton==2){
    	 int dx = event->x() - lastPos.x();
    	 int dy = event->y() - lastPos.y();
    	 float speed = 0.3;
    	 if(!GL_DEPTH_TEST){
    		 speed *=-1.0;
    	 }
    	 Rotate[0] += dy*speed; //movement in y rotates the scene around x axis (move mouse up to push the object up)
    	 Rotate[1] += dx*speed; //same as above: movement in x causing rotation in y.
    	 while (Rotate[0] >180) {Rotate[0] -=360;}
    	 while (Rotate[1] >180) {Rotate[1] -=360;}
    	 while (Rotate[0] <-180) {Rotate[0] +=360;}
    	 while (Rotate[1] <-180) {Rotate[1] +=360;}
    	 lastPos = event->pos();
    	 updateGL();
     }
 }

 void GLWidget::wheelEvent(QWheelEvent *event)
  {
	 float numDegrees = event->delta() / 8;
	 float numSteps = numDegrees / 30;
	 obj_pos[2] += -numSteps;
	 updateGL();
  }

 void GLWidget::ObjectSelection(QPoint LastPos){
	 DrawForPicking();
	 GetColourOfPoint(LastPos);
	 resetItemSelectionInfo();
	 findElement();
	 emit SelectedItemChanged();
 }

 void GLWidget::resetItemSelectionInfo(){
	 ItemSelected = false;
	 SelectedItemName = "";
	 SelectedItemIndex = -1;
	 while ( SelectedPos.size()>0){
		 SelectedPos.pop_back();
 	 }
 }

 void GLWidget::GetColourOfPoint(QPoint LastPos){
	 QSize viewport_size = size();
	 unsigned char pixels[4];
	 glReadPixels(LastPos.x(), viewport_size.height() - LastPos.y(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &pixels);
	 //cerr << "pos: "<<LastPos.x()<<" "<<LastPos.x()<<" rgba: " << (int)pixels[0] <<" "<< (int)pixels[1] <<" "<< (int)pixels[2] << " " << (int)pixels[3] << endl;
	 PickedColour[0] = (int)pixels[0];
	 PickedColour[1] = (int)pixels[1];
	 PickedColour[2] = (int)pixels[2];
	 PickedColour[3] = (int)pixels[3];
 }

 void GLWidget::findElement(){
	 for (int i =0; i<Sim01->Elements.size();i++){
		 int* ElementColour = Sim01->Elements[i]->getIdentifierColour();
		 ItemSelected = checkPickedColour(ElementColour);
		 if (ItemSelected){
			fillItemSelectionInfo(i);
			SelectedItemIndex = i;
		    update();
			break;
		}
	}
 }

 void GLWidget::fillItemSelectionInfo(int i){
	SelectedItemName = Sim01->Elements[i]->getName();
	int nNodes = Sim01->Elements[i]->getNodeNumber();
	int nDim = Sim01->Elements[i]->getDim();
	for (int j=0;j<nNodes;j++){
		for (int k =0 ;k<nDim; k++){
			QString tmpstring = QString::number(Sim01->Elements[i]->Positions[j][k], 'f', 2);
			SelectedPos.push_back(tmpstring);
			//cout<<"j: "<<j<<"k: "<<k<<" string: "<<tmpstring.toStdString()<<endl;
		}
	}
 }

 bool GLWidget::checkPickedColour(int* ElementColour){
	 for (int i=0; i<3; ++i){
		 if (ElementColour[i] != PickedColour[i]){
			  return false;
		 }
	 }
	 return true;

 }

 void GLWidget::DrawForPicking(){
	 for (int i =0; i<Sim01->Elements.size();i++){
		drawElement(i,true);
	 }
	  //To debug, you can actually draw the colour buffer to the screen, and see the change in behaviour
	 //swapBuffers();
	 glEnable(GL_DITHER);
 }

 void GLWidget::drawPrismForPicking(int i){
 	//Drawing the surfaces
 	const int nTriangle = 8; //top + bottom + 2 for each side.
 	int TriangleConnectivity[nTriangle][3] = {{0,1,2},{3,4,5},{0,2,3},{2,3,5},{0,1,3},{1,3,4},{1,2,5},{1,5,4}};

 	int* ElementColour;
 	ElementColour = Sim01->Elements[i]->getIdentifierColour();
 	//cout<<"Element "<<i<<" Color: "<<ElementColour[0]<<" "<<ElementColour[1]<<" "<<ElementColour[2]<<endl;
	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);glBegin(GL_TRIANGLES);

	float r = ElementColour[0];
	float g = ElementColour[1];
	float b = ElementColour[2];
	glColor3f(r/255,g/255,b/255);

	for (int j =0; j<nTriangle;++j){
		for (int k =0; k<3; ++k){
			int pointId = TriangleConnectivity[j][k];
			float x = Sim01->Elements[i]->Positions[pointId][0];
			float y = Sim01->Elements[i]->Positions[pointId][1];
			float z = Sim01->Elements[i]->Positions[pointId][2];
			glVertex3f( x, y, z);
		}
	}
	glEnd();
  }

