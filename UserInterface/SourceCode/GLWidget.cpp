#include <iostream>
#include <math.h>
#include <vector>
#include <string.h>

#include <QtGui>
#include <QtOpenGL>
#include <QWheelEvent>

#include "GLWidget.h"
#include "../TissueFolding/SourceCode/ShapeBase.h"

//needed for MacOS version:
//#include <OpenGL/OpenGL.h>
//#include <OpenGL/gl3.h>

using namespace std;

 GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
 {
	 //cout<<"initiating gl widget"<<endl;
	 obj_pos[0] = -10.0f;
	 obj_pos[1] =  4.0f;
     obj_pos[2] =  250.0f;//250.0f;
	 MatRot[0]  = 1.0; MatRot[1]  = 0.0; MatRot[2]  = 0.0; MatRot[3]  = 0.0;
	 MatRot[4]  = 0.0; MatRot[5]  = 1.0; MatRot[6]  = 0.0; MatRot[7]  = 0.0;
	 MatRot[8]  = 0.0; MatRot[9]  = 0.0; MatRot[10] = 1.0; MatRot[11] = 0.0;
	 MatRot[12] = 0.0; MatRot[13] = 0.0; MatRot[14] = 0.0; MatRot[15] = 1.0;
	 Qcurr[0] =  1; Qcurr[1] =  0; Qcurr[2] =  0; Qcurr[3] =  0;
	 Qlast[0] =  1; Qlast[1] =  0; Qlast[2] =  0; Qlast[3] =  0;
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
     DisplayStrains = false;
     DisplayPysProp = false;
     //current ranges:
     DisplayPysPropRange[0][0] = 1.0; DisplayPysPropRange[0][1] = 100.0;
     DisplayPysPropRange[1][0] = 1.0; DisplayPysPropRange[1][1] = 5000.0;
     DisplayPysPropRange[2][0] = 0.0; DisplayPysPropRange[2][1] = 0.5;
     DisplayPysPropRange[3][0] = -1.0; DisplayPysPropRange[3][1] = 10.0;
     DisplayPysPropRange[4][0] = -1.0; DisplayPysPropRange[3][1] = 10.0;
     //the minimum and maximum they can get:
     DisplayPysPropBounds[0][0] = 0.0; DisplayPysPropBounds[0][1] = 50.0;
     DisplayPysPropBounds[0][2] = 51.0; DisplayPysPropBounds[0][3] = 250.0;
     DisplayPysPropBounds[1][0] = 1.0; DisplayPysPropBounds[1][1] = 50.0;
     DisplayPysPropBounds[1][2] = 51.0; DisplayPysPropBounds[1][3] = 10000.0;
     DisplayPysPropBounds[2][0] = 0.0; DisplayPysPropBounds[2][1] = 0.1;
     DisplayPysPropBounds[2][2] = 0.11; DisplayPysPropBounds[2][3] = 0.5;
     DisplayPysPropBounds[3][0] = -4.0; DisplayPysPropBounds[3][1] = 0.0;
     DisplayPysPropBounds[3][2] = 0.0; DisplayPysPropBounds[3][3] = 4.0;
     DisplayPysPropBounds[4][0] = -6.0; DisplayPysPropBounds[4][1] = 0.0;
     DisplayPysPropBounds[4][2] = 0.0; DisplayPysPropBounds[4][3] = 6.0;
     //the decimals to display:
     DisplayPysPropDecimals[0] = 0;
     DisplayPysPropDecimals[1] = 0;
	 DisplayPysPropDecimals[2] = 2;
	 DisplayPysPropDecimals[3] = 2;
	 DisplayPysPropDecimals[4] = 2;
	 DisplayPysPropSteps[0] = 1;
	 DisplayPysPropSteps[1] = 10;
	 DisplayPysPropSteps[2] = 0.05;
	 DisplayPysPropSteps[3] = 0.05;
	 DisplayPysPropSteps[4] = 0.05;

  	 setSizePolicy(QSizePolicy ::Expanding , QSizePolicy ::Expanding );
  	 drawNetForces = false;
  	 drawPackingForces = false;
     drawVelocities = false;
     drawPeripodialMembrane = true;
     drawColumnar = true;
     ManualNodeSelection = false;
     ManualSelectedNodeId = -100;
     PerspectiveView = true;
     orthoViewLimits[0] = -250;
     orthoViewLimits[1] =  250;
     orthoViewLimits[2] = -250;
     orthoViewLimits[3] =  250;
     orthoViewLimits[4] = -1000;
     orthoViewLimits[5] =  1000;
     displayBoundingBox = false;
     xClip = 1000.0;
     yClip = 1000.0;
     zClip = 1000.0;
     cout<<"gl initiated"<<endl;
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
     MainShapeLineThickness = (LineRange[1]-LineRange[0])/5.0;
     initialiseNodeColourList();
     DisplayStrains = false;
     glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );
     glMultMatrixf(MatRot);

    //On MacOS version:
	//int a = 1;
    //    int *swap_interval = &a;
    //    CGLContextObj cgl_context = CGLGetCurrentContext();
    //    CGLSetParameter(cgl_context, kCGLCPSwapInterval, swap_interval);
 }

 void GLWidget::initialiseNodeColourList(){
	 const int n = Sim01->Nodes.size();
     NodeColourList = new float*[n];
     for (int i=0; i<n; ++i){
    	 NodeColourList[i]=new float[3];
    	 NodeColourList[i][0]=0.0;
    	 NodeColourList[i][1]=0.0;
    	 NodeColourList[i][2]=0.0;
     }
 }

 void GLWidget::paintGL()
 {
     //glScalef(this->devicePixelRatio(),this->devicePixelRatio(),this->devicePixelRatio());
     QSize viewport_size = size();
     glViewport(0, 0, viewport_size.width()*this->devicePixelRatio(), viewport_size.height()*this->devicePixelRatio());

	 glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();

	 if (PerspectiveView) {
         glFrustum(-1*aspectratio, 1*aspectratio, -1, 1, 1, 10000); // near and far match your triangle Z distance
	 }
	 else{
		 glOrtho(orthoViewLimits[2]*aspectratio, orthoViewLimits[3]*aspectratio, orthoViewLimits[2], orthoViewLimits[3], orthoViewLimits[4],orthoViewLimits[5]);
	 }

	 glMatrixMode(GL_MODELVIEW);
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	 glLoadIdentity();
	 constructNodeColourList();
	 drawColourbar();
	 drawAxesArrows();
	 if (drawTissueScaleBar){
		 drawScaleBar();
	 }
	 glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] );
	 glTranslatef( Sim01->SystemCentre[0], Sim01->SystemCentre[1], -Sim01->SystemCentre[2]);
	 glMultMatrixf(MatRot);
	 glTranslatef( -Sim01->SystemCentre[0], -Sim01->SystemCentre[1], Sim01->SystemCentre[2]);


	 if (ItemSelected){
		 drawReferenceElement(SelectedItemIndex);
		 highlightElement(SelectedItemIndex);
	 }
	 if (ManualNodeSelection){
		 highlightNode(ManualSelectedNodeId);
	 }
	 int n = Sim01->Elements.size();
	 for (int i =0; i<n;i++){
		 drawElement(i,false);
	 }
	 drawForces();
	 drawPackForces();
	 drawNodeVelocities();
	 drawBoundingBox();
     drawPipette();

/*
	glColor3f(0.8,0,0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces

	cout<<"Drawing triangles now: "<<endl;
	for (int i =0; i<Sim01->TrianglesToDraw.size();i++){
		glBegin(GL_TRIANGLES);
		for (int k=0;k<3;++k){
			float x = Sim01->NodesToDraw [Sim01->TrianglesToDraw[i][k]][0];
			float y = Sim01->NodesToDraw [Sim01->TrianglesToDraw[i][k]][1];
			float z = Sim01->NodesToDraw [Sim01->TrianglesToDraw[i][k]][2];
			cout<<Sim01->TrianglesToDraw[i][k]<<" ";
			glVertex3f( x, y, z);
		}
		cout<<endl;
		glEnd();
	}
*/

	 DisplayFixedNodes= true;
	 if (DisplayFixedNodes){
		 drawFixedNodes();
	 }

     //swapBuffers();
 }


 void GLWidget::updateClipping(){
	 //cout<<"updating the clipping"<<endl;
	 int n = Sim01->Elements.size();
	 for (int i=0; i<n; ++i){
		 Sim01->Elements[i]->checkDisplayClipping(xClip, yClip, zClip);
	 }
 }

 bool GLWidget::checkIfDrawingElement(int i){
	 bool drawthisElement = true;
	 if (Sim01->Elements[i]->IsClippedInDisplay){
		 drawthisElement = false;
	 }
	 if (Sim01->Elements[i]->IsAblated){
		 drawthisElement = false;
	 }
	 if (!drawPeripodialMembrane && Sim01->Elements[i]->tissueType == 1){	//I am NOT drawing peripodial membrane and this element is a peripodial element
	 	 drawthisElement = false;
	 }
	 if (!drawColumnar && Sim01->Elements[i]->tissueType == 0){		//I am NOT drawing columnar layer and this element is a columnar element
		 drawthisElement = false;
	 }
	 return drawthisElement;
 }

 bool GLWidget::checkIfDrawingNode(int i){
 	 bool drawthisNode = true;
 	 if (!drawPeripodialMembrane && Sim01->Nodes[i]->tissueType == 1){	//I am NOT drawing peripodial membrane and this element is a peripodial node
 	 	 drawthisNode = false;
 	 }
	 if (!drawColumnar && Sim01->Nodes[i]->tissueType == 0){	//I am NOT drawing columnar layer and this element is a columnar node
		 drawthisNode = false;
	 }
 	 return drawthisNode;
 }

 void GLWidget::drawElement(int i, bool picking){
	 bool drawCurrentElement = checkIfDrawingElement(i);
	 if (drawCurrentElement){
		 int ShapeType = Sim01->Elements[i]->getShapeType();
		 if (ShapeType == 1 ){
			 if (picking){
				drawPrismForPicking(i);
			 }
			 else{
				 drawPrism(i);
			 }
		 }
		 else if (ShapeType == 4 ){
			 if (picking){
				drawTriangleForPicking(i);
			 }
			 else{
				drawTriangle(i);
			 }
		 }
	 }
 }


void GLWidget::highlightNode(int i){
	float markerSize = 1.0;
	double x = Sim01->Nodes[i]->Position[0];
	double y = Sim01->Nodes[i]->Position[1];
	double z[3] = {Sim01->Nodes[i]->Position[2]-markerSize,Sim01->Nodes[i]->Position[2]+markerSize,Sim01->Nodes[i]->Position[2]};
	glColor3f(0.8,0.8,0.0);
	glBegin(GL_TRIANGLES);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x+markerSize, y, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x-markerSize, y, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x, y+markerSize, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x, y-markerSize, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x+markerSize, y+markerSize, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x+markerSize, y-markerSize, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x-markerSize, y+markerSize, z[2]);
		for (int k =0; k<2; ++k){glVertex3f( x, y, z[k]);}
		glVertex3f( x-markerSize, y-markerSize, z[2]);
	glEnd();
}

 void GLWidget::highlightElement(int i){
	int ShapeType = Sim01->Elements[i]->getShapeType();
	if (ShapeType == 1){
		highlightPrism(i);
	}
	else if (ShapeType == 4){
		highlightTriangle(i);
	}
 }

 void GLWidget::drawReferenceElement(int i){
	 //glPushMatrix();
	 	// glTranslatef( 0, -20, 0 );
		 int ShapeType = Sim01->Elements[i]->getShapeType();
		 if (ShapeType == 1 ){
			 drawReferencePrism(i);
		 }
		 else if (ShapeType == 4 ){
			 drawReferenceTriangle(i);
		 }
	// glPopMatrix();
 }

 void GLWidget::constructNodeColourList(){
	float threshold = 1E-10;
	vector<Node*>::iterator itNode;
	for (itNode=Sim01->Nodes.begin(); itNode<Sim01->Nodes.end(); ++itNode){
		float* currColour;
		currColour = new float[3];
		if(!DisplayStrains && !DisplayPysProp && !drawNetForces && !drawPackingForces && !drawVelocities){
			//I am not displaying ant data on the colours, therefore I do not need any calculaitons on element basis, the colour is constant
			if ((*itNode)->tissueType == 0){ // columnar layer
				NodeColourList[(*itNode)->Id][0]=0.75;
				NodeColourList[(*itNode)->Id][1]=1.0;
				NodeColourList[(*itNode)->Id][2]=1.0;
			}
			else if ((*itNode)->tissueType == 1){ // Peripodial membrane
				NodeColourList[(*itNode)->Id][0]=1.0;
				NodeColourList[(*itNode)->Id][1]=1.0;
				NodeColourList[(*itNode)->Id][2]=0.75;
			}
		}
		else{
			if(DisplayStrains){
				float StrainMag = 0.0;
				int nConnectedElements = (*itNode)->connectedElementIds.size();
				for (int i=0;i<nConnectedElements; ++i){
					float TmpStrainMag =0.0, PlasticStrainMag=0;
					Sim01->Elements[(*itNode)->connectedElementIds[i]]->getStrain(StrainToDisplay, TmpStrainMag);
                    StrainMag += (TmpStrainMag)*(*itNode)->connectedElementWeights[i];
				}
				getDisplayColour(currColour, StrainMag);
			}
			else if (DisplayPysProp){
				float PysPropMag = 0.0;
				//float* PysPropColour;
				//PysPropColour = new float[3];
				//If the physical property is viscosity, then get the colour directly
				if (PysPropToDisplay == 0){
					PysPropMag = (*itNode)->Viscosity;
				}
				else{
					int nConnectedElements = (*itNode)->connectedElementIds.size();
					for (int i=0;i<nConnectedElements; ++i){
						float TmpPysPropMag =0.0;
                        Sim01->Elements[(*itNode)->connectedElementIds[i]]->getPysProp(PysPropToDisplay, TmpPysPropMag, Sim01->dt);
						PysPropMag += TmpPysPropMag*(*itNode)->connectedElementWeights[i];
					}
				}
				getDisplayColour(currColour,PysPropMag);
			}
			else if(drawNetForces){
				float ForceMag = 0.0;
				double F[3];
				F[0] = Sim01->SystemForces[0][(*itNode)->Id][0];
				F[1] = Sim01->SystemForces[0][(*itNode)->Id][1];
				F[2] = Sim01->SystemForces[0][(*itNode)->Id][2];
				ForceMag = F[0]* F[0] + F[1]*F[1] + F[2]* F[2];
				ForceMag = pow(ForceMag,(float)0.5);
				if (ForceMag>threshold){
					getForceColour(currColour,ForceMag);
				}
				else{
					currColour[0] = 1.0;
					currColour[1] = 1.0;
					currColour[2] = 0.8;
				}
			}
			else if(drawPackingForces){
                //cout<<"Drawing Packing Forces"<<endl;
				float ForceMag = 0.0;
				double F[3];
				F[0] = Sim01->PackingForces[0][(*itNode)->Id][0];
				F[1] = Sim01->PackingForces[0][(*itNode)->Id][1];
				F[2] = Sim01->PackingForces[0][(*itNode)->Id][2];
				ForceMag = F[0]* F[0] + F[1]*F[1] + F[2]* F[2];
				ForceMag = pow(ForceMag,(float)0.5);
                //cout<<"PAcking Force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
                if (ForceMag>threshold){
					getForceColour(currColour,ForceMag);
				}
				else{
					currColour[0] = 1.0;
					currColour[1] = 1.0;
					currColour[2] = 0.8;
				}
			}
			else if(drawVelocities){
				float VelMag = 0.0;
				double v[3] = {(*itNode)->Velocity[0][0],(*itNode)->Velocity[0][1],(*itNode)->Velocity[0][2]};
				VelMag = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
				VelMag = pow(VelMag,(float)0.5);
				if (VelMag>threshold){
					getVelocityColour(currColour, VelMag);
				}
				else{
					currColour[0] = 1.0;
					currColour[1] = 1.0;
					currColour[2] = 0.8;
				}
			}
			NodeColourList[(*itNode)->Id][0]=currColour[0];
			NodeColourList[(*itNode)->Id][1]=currColour[1];
			NodeColourList[(*itNode)->Id][2]=currColour[2];
			delete[] currColour;
		}
	}
}

 float** GLWidget::getElementColourList(int i){
	const int n = Sim01->Elements[i]->getNodeNumber();
	int* NodeIds = Sim01->Elements[i]->getNodeIds();
	float** NodeColours;
	NodeColours = new float*[n];
	for (int j = 0; j<n; j++){		NodeColours[j] = new float[3];
		NodeColours[j][0]=NodeColourList[NodeIds[j]][0];
		NodeColours[j][1]=NodeColourList[NodeIds[j]][1];
		NodeColours[j][2]=NodeColourList[NodeIds[j]][2];
	}
	return NodeColours;
 }


 void GLWidget::drawPrism(int i){
	 //Drawing the surfaces
	const int nTriangle = 8; //top + bottom + 2 for each side.
	int TriangleConnectivity[nTriangle][3] = {{0,1,2},{3,4,5},{0,2,3},{2,3,5},{0,1,3},{1,3,4},{1,2,5},{1,5,4}};
	const int nLineStrip = 12;
	int BorderConnectivity[nLineStrip] = {0,2,5,3,0,1,4,3,5,4,1,2};
	float** NodeColours;
	NodeColours = getElementColourList(i);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
	glBegin(GL_TRIANGLES);
		for (int j =0; j<nTriangle;++j){
			for (int k =0; k<3; ++k){
				int pointId = TriangleConnectivity[j][k];
				glColor3f(NodeColours[pointId][0],NodeColours[pointId][1],NodeColours[pointId][2]);
				float x = Sim01->Elements[i]->Positions[pointId][0];
				float y = Sim01->Elements[i]->Positions[pointId][1];
				float z = Sim01->Elements[i]->Positions[pointId][2];
				//cout<<"triangle : "<<j<<" point: "<<k<<" position : "<<x<<" y "<<y<<" z "<<z<<" colour: "<<NodeColours[pointId][0]<<" "<<NodeColours[pointId][1]<<" "<<NodeColours[pointId][2]<<endl;
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
			float x = Sim01->Elements[i]->Positions[pointId][0];
			float y = Sim01->Elements[i]->Positions[pointId][1];
			float z = Sim01->Elements[i]->Positions[pointId][2];
			glVertex3f( x, y, z);
		}
	glEnd();
	delete[] NodeColours;
 }

 void GLWidget::drawTriangle(int i){
 	//Drawing the surfaces
 	const int nTriangle = 1; //a triangle with 3 points needs 1 actual triangle to draw
 	int TriangleConnectivity[nTriangle][3] = {{0,1,2}};
 	const int nLineStrip = 4;
 	int BorderConnectivity[nLineStrip] = {0,1,2,0};
 	float** NodeColours;
 	NodeColours = getElementColourList(i);

 	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
 	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
 	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
 	glBegin(GL_TRIANGLES);
 		for (int j =0; j<nTriangle;++j){
 			for (int k =0; k<3; ++k){
 				int pointId = TriangleConnectivity[j][k];
 				glColor3f(NodeColours[pointId][0],NodeColours[pointId][1],NodeColours[pointId][2]);
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
 			float x = Sim01->Elements[i]->Positions[pointId][0];
 			float y = Sim01->Elements[i]->Positions[pointId][1];
 			float z = Sim01->Elements[i]->Positions[pointId][2];
 			glVertex3f( x, y, z);
 		}
 	glEnd();
 	delete[] NodeColours;
  }

 void GLWidget::getForceColour(float* OutputColour, float Data){
	 double scale2[2] = {0,1000.0};
	 double r = (Data- scale2[0])/(scale2[1]-scale2[0]);
	 //OutputColour[0] = r;
	 //OutputColour[1] = 0.0;
	 //OutputColour[2] = 0.0;
	 OutputColour[0] = 1.0;
	 OutputColour[1] = 1.0-r;
	 OutputColour[2] = 0.8*(1.0-r);
 }

 void GLWidget::getVelocityColour(float* OutputColour, float Data){
	 double scale2[2] = {0,10.0};
	 double b = (Data- scale2[0])/(scale2[1]-scale2[0]);
	 //OutputColour[0] = 0.0;
	 //OutputColour[1] = b;
	 //OutputColour[2] = 0.0;
	 OutputColour[0] = 1.0-b;
	 OutputColour[1] = 1.0;
	 OutputColour[2] = 0.8*(1.0-b);
 }

 void GLWidget::getDisplayColour(float* OutputColour, float Data){
	 float DataMin=0, DataMax=0;
	 if(DisplayStrains){
		 DataMin = DisplayStrainRange[0];
		 DataMax = DisplayStrainRange[1];
	 }
	 else if (DisplayPysProp){
		 DataMin = DisplayPysPropRange[PysPropToDisplay][0];
		 DataMax = DisplayPysPropRange[PysPropToDisplay][1];
	 }
	 float segment = (DataMax - DataMin)/5.0;
	 float r,g,b;
	 float minR = 0.4;
	 float minB = 0.4;
	 OutputColour[0] = 0.0;
	 OutputColour[1] = 0.0;
	 OutputColour[2] = 0.0;
	 if ((Data - DataMin) < segment/2){
		 float d = (Data - DataMin);
		 r = (1.0 - minR)/ (segment/2.0) * d + minR;
		 g = 0;
		 b = 0;
	 }
	 else if ((Data - DataMin) <segment*1.5){
		 float d = (Data - DataMin) - segment/2.0;
		 r = 1.0;
		 g = d/segment;
		 b = 0;
	 }
	 else if ((Data - DataMin)< segment*2.5){
		 float d = (Data - DataMin) - 1.5*segment;
		 r = 1.0 - d/segment;
		 g = 1.0;
		 b = 0;
	 }
	 else if ((Data - DataMin)< segment*3.5){
		 float d = (Data - DataMin) - 2.5* segment;
		 r = 0;
		 g = 1.0;
		 b = d/segment;
	 }
	 else if ((Data - DataMin)< segment*4.5){
		 float d = (Data - DataMin) - 3.5* segment;
		 r = 0;
		 g = 1.0 - d/segment;
		 b = 1;
	 }
	 else{
		 float d = (Data - DataMin) - 4.5* segment;
		 r = 0;
		 g = 0;
		 b = 1.0 - d*(1-minB)/(segment/2.0);
	 }
	 //invert the display from blue to red for growth!:
	 if (DisplayPysProp &&  PysPropToDisplay==3){
		 double tmpred = r;
		 r=b;
		 b=tmpred;
	 }
	 //cout<<"data: "<<Data<<" colour: "<<r<<" "<<g<<" "<<b<<"Datamin: "<<DataMin<<"DataMax: "<<DataMin<<endl;
	 OutputColour[0] = r;
	 OutputColour[1] = g;
	 OutputColour[2] = b;
 }

 void GLWidget::highlightPrism(int i){
	const int nTriangle = 8; //top + bottom + 2 for each side.
	int TriangleConnectivity[nTriangle][3] = {{0,1,2},{3,4,5},{0,2,3},{2,3,5},{0,1,3},{1,3,4},{1,2,5},{1,5,4}};
	const int nLineStrip = 12;
	int BorderConnectivity[nLineStrip] = {0,2,5,3,0,1,4,3,5,4,1,2};
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
	glBegin(GL_TRIANGLES);
		for (int j =0; j<nTriangle;++j){
			for (int k =0; k<3; ++k){
				int pointId = TriangleConnectivity[j][k];
				glColor3f(0.30,0.30,0.30);
				float x = Sim01->Elements[i]->Positions[pointId][0];
				float y = Sim01->Elements[i]->Positions[pointId][1];
				float z = Sim01->Elements[i]->Positions[pointId][2];
				glVertex3f( x, y, z);
			}
		}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);
	glLineWidth(2*MainShapeLineThickness);
	//Drawing the borders
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3f(0,1.0,0);
	glBegin(GL_LINE_STRIP);
		for (int j =0; j<nLineStrip;++j){
			int pointId = BorderConnectivity[j];
			float x = Sim01->Elements[i]->Positions[pointId][0];
			float y = Sim01->Elements[i]->Positions[pointId][1];
			float z = Sim01->Elements[i]->Positions[pointId][2];
			glVertex3f( x, y, z);
		}
	glEnd();
 }


 void GLWidget::highlightTriangle(int i){
	const int nTriangle = 1; //a triangle with 3 points needs 1 actual triangle to draw
	int TriangleConnectivity[nTriangle][3] = {{0,1,2}};
	const int nLineStrip = 4;
	int BorderConnectivity[nLineStrip] = {0,1,2,0};

 	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
 	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
 	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
 	glBegin(GL_TRIANGLES);
 		for (int j =0; j<nTriangle;++j){
 			for (int k =0; k<3; ++k){
 				int pointId = TriangleConnectivity[j][k];
 				glColor3f(0.30,0.30,0.30);
 				float x = Sim01->Elements[i]->Positions[pointId][0];
 				float y = Sim01->Elements[i]->Positions[pointId][1];
 				float z = Sim01->Elements[i]->Positions[pointId][2];
 				glVertex3f( x, y, z);
 			}
 		}
 	glEnd();
 	glDisable(GL_POLYGON_OFFSET_FILL);
 	glLineWidth(2*MainShapeLineThickness);
 	//Drawing the borders
 	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 	glColor3f(0,1.0,0);
 	glBegin(GL_LINE_STRIP);
 		for (int j =0; j<nLineStrip;++j){
 			int pointId = BorderConnectivity[j];
 			float x = Sim01->Elements[i]->Positions[pointId][0];
 			float y = Sim01->Elements[i]->Positions[pointId][1];
 			float z = Sim01->Elements[i]->Positions[pointId][2];
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
 	//cout<<"getting the reference pos:"<<endl;
 	double** pos = Sim01->Elements[i]->getReferencePos();
 	//cout<<"drawing the reference pos:"<<endl;
 	//drawing the posiions of reference elementin red
 	glColor3f(1,0,0);
 	glBegin(GL_LINE_STRIP);
 		for (int j =0; j<nLineStrip;++j){
 			int pointId = BorderConnectivity[j];
 			//cout<<j<<" point id: "<<pointId<<endl;
 			if(j==0){glColor3f(0,0,1);}else{glColor3f(1,0,0);}
 			float x = pos[pointId][0];
 			float y = pos[pointId][1];
 			float z = pos[pointId][2];
 			//cout<<"x,y,z: "<<x<<" "<<y<<" "<<z<<endl;
 			glVertex3f( x, y, z);
 		}
 	glEnd();
}

 void GLWidget::drawReferenceTriangle(int i){
	const int nLineStrip = 4;
	int BorderConnectivity[nLineStrip] = {0,1,2,0};

 	glLineWidth(ReferenceLineThickness);
 	//Drawing the borders
 	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 	double** pos = Sim01->Elements[i]->getReferencePos();
	//drawing the posiions of reference elementin red
 	glColor3f(1,0,0);
 	glBegin(GL_LINE_STRIP);
 		for (int j =0; j<nLineStrip;++j){
 			int pointId = BorderConnectivity[j];
 			//cout<<j<<" point id: "<<pointId<<endl;
 			if(j==0){glColor3f(0,0,1);}else{glColor3f(1,0,0);}
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
     width *= this->devicePixelRatio();
     height *= this->devicePixelRatio();
     int  side = qMin(width, height);
     glViewport((width - side) / 2, (height - side) / 2, width, height);
     aspectratio = float(width) /float (height);
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
     //cout<<(width - side) / 2<<"  "<<(height - side) / 2<<"  "<<width<<"  "<<height<<"  "<<side<<endl;
 }

 void GLWidget::mousePressEvent(QMouseEvent *event)
 {
     lastPos = event->pos();
     if(event->button()==Qt::LeftButton){
         MouseButton = 0;
         InitialClickPos = event->pos();
         //cout<<"Read an initial click at: "<<InitialClickPos.x()<<" "<<InitialClickPos.y()<<endl;
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
          //cout<<"read a mouse relaese event at: "<<lastPos.x()<<" "<<lastPos.y()<<endl;
          int dx = InitialClickPos.x() - lastPos.x();
          int dy = InitialClickPos.x() - lastPos.x();
          dx *=this->devicePixelRatio();
          dy *=this->devicePixelRatio();
          if (dx<0) {
              dx *= -1;
          }
          if (dy<0) {
              dy *= -1;
          }
          if (dx<5 && dy<5){
        	  ObjectSelection(lastPos);
          }
          else{
              //cout<<"dx is above 5, the mouse has been moved I will move the scene"<<endl;
              //this is an alternative view movement option added as mac OS X does not have a default middle button click,
              //replica of MouseButton == 1 event:
              int dx = lastPos.x() - InitialClickPos.x();
              int dy = lastPos.y() - InitialClickPos.y();
              dx *=this->devicePixelRatio();
              dy *=this->devicePixelRatio();
              cout<<"dx: "<<dx<<" dy: "<<dy<<endl;
              float speed = 0.01;
              obj_pos[0] +=  dx*speed;
              obj_pos[1] += -dy*speed;
              lastPos = event->pos();
              updateGL();
          }
      }
  }

 void GLWidget::mouseMoveEvent(QMouseEvent *event)
 {
	 if (MouseButton==0){
		 lastPos = event->pos();
	 }
	 else if(MouseButton==1){
         //eliminating this for mac:
         int dx = event->x() - lastPos.x();
    	 int dy = event->y() - lastPos.y();
         dx *=this->devicePixelRatio();
         dy *=this->devicePixelRatio();
    	 float speed = 0.01;
    	 obj_pos[0] +=  dx*speed;
    	 obj_pos[1] += -dy*speed;
         //cerr<<"right button "<<dx<<" "<<dy<<" "<<obj_pos[0]<<" "<<obj_pos[1]<<endl;
    	 lastPos = event->pos();
         updateGL();
     }
	 else if(MouseButton==2){

		 float speed = 50;
		 QSize viewport_size = size();
         float width  = viewport_size.width()/this->devicePixelRatio();
         float height = viewport_size.height()/this->devicePixelRatio();
		 double initialPos[3] = {lastPos.x()/(width*2.0) - 1.0, lastPos.y()/(height*2.0) - 1.0,  0.0};
		 double finalPos[3]   = {event->x()/(width*2.0) - 1.0, event->y()/(height*2.0)  - 1.0,  0.0};
         initialPos[1] = (-1.0) * initialPos[1];
         finalPos[1]   = (-1.0) * finalPos[1];
         float r = 20;
		 float r2 = r*r;
		 //Projecting event location:
		 double lengthSq = finalPos[0]*finalPos[0]+finalPos[1]*finalPos[1];
		 if (lengthSq  <= r2){
			 finalPos[2] = pow((r2 - lengthSq),0.5);
		 }
		 else
		 {
			 double length = pow(lengthSq,0.5);
			 finalPos[0] = finalPos[0]/length*r;
			 finalPos[1] = finalPos[1]/length*r;
			 finalPos[2] = finalPos[2]/length*r;
		 }
		 double mag = pow((finalPos[0]*finalPos[0] + finalPos[1]*finalPos[1] + finalPos[2]*finalPos[2]),0.5);
		 if (fabs(mag) > 0.001 && fabs(mag - 1.0f) > 0.001) {
			 finalPos[0] = finalPos[0]/mag;
			 finalPos[1] = finalPos[1]/mag;
			 finalPos[2] = finalPos[2]/mag;
		 }
		 //projecting last pos:
		 lengthSq = initialPos[0]*initialPos[0]+initialPos[1]*initialPos[1];
		 if (lengthSq  <= r2)
			 initialPos[2] = pow((r2 - lengthSq),0.5);
		 else
		 {
			 double  length = pow(lengthSq ,0.5);
			 initialPos[0] = initialPos[0]/length*r;
			 initialPos[1] = initialPos[1]/length*r;
			 initialPos[2] = initialPos[2]/length*r;
		 }
		 mag = pow((initialPos[0]*initialPos[0] + initialPos[1]*initialPos[1] + initialPos[2]*initialPos[2]),0.5);
		 if (fabs(mag) > 0.001 && fabs(mag - 1.0f) > 0.001) {
			 initialPos[0] = initialPos[0]/mag;
			 initialPos[1] = initialPos[1]/mag;
			 initialPos[2] = initialPos[2]/mag;
		 }

		 double cross[3];
		 cross[0] = initialPos[1]*finalPos[2] - initialPos[2]*finalPos[1];
		 cross[1] = initialPos[2]*finalPos[0] - initialPos[0]*finalPos[2];
		 cross[2] = initialPos[0]*finalPos[1] - initialPos[1]*finalPos[0];
		 mag = pow((cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]),0.5);
		 cross[0] /= mag;
		 cross[1] /= mag;
		 cross[2] /= mag;

		 double dot;
		 dot = initialPos[0]*finalPos[0] + initialPos[1]*finalPos[1] + initialPos[2]*finalPos[2];
		 double angle = acosf(min(1.0,dot));
		 angle *= speed;

         if (angle>1E-4 || angle < -1E-4){
             double* Qrot;
             Qrot = new double[4];
             Qrot[0]= cosf(angle); Qrot[1] = sinf(angle)* cross[0]; Qrot[2] = sinf(angle)*cross[1]; Qrot[3] = sinf(angle)*cross[2];
             //cout<<"QRot:		"<<Qrot[0]<<" "<<Qrot[1]<<" "<<Qrot[2]<<" "<<Qrot[3]<<endl;
             rotateByQuaternians(Qrot);
         }
         //cout<<"initialpos: "<<initialPos[0]<<" "<<initialPos[1]<<" "<<initialPos[2]<<endl;
         //cout<<"fianlpos:   "<<finalPos[0]<<" "<<finalPos[1]<<" "<<finalPos[2]<<endl;
         //cout<<"angle:      "<<angle<<endl;

		 lastPos = event->pos();
		 updateGL();
     }
 }

 void GLWidget::rotateByQuaternians (double* Qrot){
	//ArcBall rotation
	normaliseCurrentRotationAngle(Qrot);
	rotateCurrentRotationQuaternian(Qrot);
	rotateMatrix();
 }

 void GLWidget::normaliseCurrentRotationAngle (double* Qrot){
	//Normalising the input:
	double mag2 = Qrot[0] * Qrot[0] + Qrot[1] * Qrot[1] + Qrot[2] * Qrot[2] + Qrot[3] * Qrot[3];
	if (fabs(mag2) > 0.00000001 && fabs(mag2 - 1.0f) > 0.00000001) {
		double mag = pow(mag2,0.5);
		Qrot[0] /= mag;
		Qrot[1] /= mag;
		Qrot[2] /= mag;
		Qrot[3] /= mag;
	}
 }

 void GLWidget::rotateCurrentRotationQuaternian(double* Qrot){
	//Multiplying current rotation quaternian with the rotation I want:
	Qcurr[0] = 	Qlast[0] * Qrot[0] - Qlast[1] * Qrot[1] - Qlast[2] * Qrot[2] - Qlast[3] * Qrot[3];
	Qcurr[1] =  Qlast[0] * Qrot[1] + Qlast[1] * Qrot[0] - Qlast[2] * Qrot[3] + Qlast[3] * Qrot[2];
	Qcurr[2] =  Qlast[0] * Qrot[2] + Qlast[1] * Qrot[3] + Qlast[2] * Qrot[0] - Qlast[3] * Qrot[1];
	Qcurr[3] =  Qlast[0] * Qrot[3] - Qlast[1] * Qrot[2] + Qlast[2] * Qrot[1] + Qlast[3] * Qrot[0];

	Qlast[0] = Qcurr[0];
	Qlast[1] = Qcurr[1];
	Qlast[2] = Qcurr[2];
	Qlast[3] = Qcurr[3];
 }

 void GLWidget::rotateMatrix(){
	double x2 = Qcurr[1] * Qcurr[1];  double y2 = Qcurr[2] * Qcurr[2];  double z2 = Qcurr[3] * Qcurr[3];
	double xy = Qcurr[1] * Qcurr[2];  double xz = Qcurr[1] * Qcurr[3];  double yz = Qcurr[2] * Qcurr[3];
	double wx = Qcurr[0] * Qcurr[1];  double wy = Qcurr[0] * Qcurr[2];  double wz = Qcurr[0] * Qcurr[3];

	MatRot[0]  = 1.0f - 2.0f * (y2 + z2);
	MatRot[1]  = 2.0f * (xy + wz);
	MatRot[2]  = 2.0f * (xz - wy);
	MatRot[3]  = 0.0f;

	MatRot[4]  = 2.0f * (xy - wz);
	MatRot[5]  = 1.0f - 2.0f * (x2 + z2);
	MatRot[6]  = 2.0f * (yz + wx);
	MatRot[7]  = 0.0f;

	MatRot[8]  = 2.0f * (xz + wy);
	MatRot[9]  = 2.0f * (yz - wx);
	MatRot[10] = 1.0f - 2.0f * (x2 + y2);
	MatRot[11] = 0.0f;

	MatRot[12] = 0.0f;
	MatRot[13] = 0.0f;
	MatRot[14] = 0.0f;
	MatRot[15] = 1.0f;
 }

 void GLWidget::wheelEvent(QWheelEvent *event)
  {
	 float numDegrees = event->delta() / 8;
	 float numSteps = numDegrees / 30;
	 float speed = 5.0;
	 //if (PerspectiveView){
		 obj_pos[2] += -speed*numSteps;
	 //}
	 //else{
		 //ortagonal view, the zoomiing should change the x & y clipping
		// speed *= 2.0;
		 orthoViewLimits[2] += speed*numSteps;
		 orthoViewLimits[3] -= speed*numSteps;
		 orthoViewLimits[0] = orthoViewLimits[2]*aspectratio;
		 orthoViewLimits[1] = orthoViewLimits[3]*aspectratio;
	 //}
	 updateGL();
  }

 void GLWidget::ObjectSelection(QPoint LastPos){
	 drawForPicking();
	 getColourOfPoint(LastPos);
	 resetItemSelectionInfo(1);
	 findElement();
	 emit SelectedItemChanged();
 }


 void GLWidget::manualElementSelection(int i){
	 resetItemSelectionInfo(2);
	 bool validElementId = findElement(i);
	 if (validElementId){
		emit SelectedItemChanged();
	 }
 }

 void GLWidget::manualNodeSelection(int i){
	 resetItemSelectionInfo(3);
	 bool validNodeId = findNode(i);
	 if (validNodeId){
		emit SelectedItemChanged();
	 }
 }

 void GLWidget::resetItemSelectionInfo(int source){
	 if (source == 1){
		 //The source is item selection via screen click, I will clear up the changes of both manual inputs
		 emit NeedToClearManualElementSelection();
		 emit NeedToClearManualNodeSelection();
	 }
	 else if (source == 2){
		 //source is a manual element selection, I will clear up the manual node selection input
		 emit NeedToClearManualNodeSelection();
	 }
	 else if (source == 3){
		 //source is manual node selection, I will clear the manual element selection
		 emit NeedToClearManualElementSelection();
	 }
	 ItemSelected = false;
	 SelectedItemName = "";
	 SelectedItemIndex = -1;
	 while ( SelectedPos.size()>0){
		 SelectedPos.pop_back();
 	 }
	 while ( SelectedId.size()>0){
		 SelectedId.pop_back();
 	 }
 }

 void GLWidget::getColourOfPoint(QPoint LastPos){
	 QSize viewport_size = size();
	 unsigned char pixels[4];
     glReadPixels(LastPos.x()*this->devicePixelRatio(), viewport_size.height()*this->devicePixelRatio() - LastPos.y()*this->devicePixelRatio(), 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &pixels);
	 //cerr << "pos: "<<LastPos.x()<<" "<<LastPos.x()<<" rgba: " << (int)pixels[0] <<" "<< (int)pixels[1] <<" "<< (int)pixels[2] << " " << (int)pixels[3] << endl;
	 PickedColour[0] = (int)pixels[0];
	 PickedColour[1] = (int)pixels[1];
	 PickedColour[2] = (int)pixels[2];
	 PickedColour[3] = (int)pixels[3];
 }

 void GLWidget::findElement(){
	 int n = Sim01->Elements.size();
	 for (int i =0; i<n;i++){
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

bool GLWidget::findElement(int i){
	int n = Sim01->Elements.size();
	if (i<n){
		if (checkIfDrawingElement(i)){
			ItemSelected = true;
			fillItemSelectionInfo(i);
			SelectedItemIndex = i;
			//update();
			return true;
		}
	}
	return false;
 }

bool GLWidget::findNode(int i){
	int n = Sim01->Nodes.size();
	if (i<n){
		//Node id is small enough
		n = Sim01->Elements.size();
		for (int j =0; j<n; ++j){
			if(Sim01->Elements[j]->DoesPointBelogToMe(i)){
				//The node belongs to the element:
				if (checkIfDrawingElement(j)){
					//We are drawing the element
					ItemSelected = true;
					fillItemSelectionInfo(j);
					SelectedItemIndex = j;
					ManualNodeSelection = true;
					ManualSelectedNodeId = i;
					//update();
					return true;
				}
			}
		}
	}
	return false;
 }

 void GLWidget::fillItemSelectionInfo(int i){
	SelectedItemName = Sim01->Elements[i]->getName();
	int nNodes = Sim01->Elements[i]->getNodeNumber();
	int nDim = Sim01->Elements[i]->getDim();
	int* NodeIds =  Sim01->Elements[i]->getNodeIds();
	for (int j=0;j<nNodes;j++){
		for (int k =0 ;k<nDim; k++){
			QString tmpstring = QString::number(Sim01->Elements[i]->Positions[j][k], 'f', 2);
			SelectedPos.push_back(tmpstring);
			//cout<<"j: "<<j<<"k: "<<k<<" string: "<<tmpstring.toStdString()<<endl;
		}
		QString tmpstring = QString::number(Sim01->Elements[i]->NodeIds[j], 'f', 0);
		SelectedId.push_back(tmpstring);
	}
	//cout<<"SelectedItemName: "<<SelectedItemName<<endl;
 }

 void GLWidget::drawColourbar(){
	 glPushMatrix();
	 	 if (PerspectiveView) {
	 		 glTranslatef( 0, 0, -20.0f);
	 	 }
	 	 else{
	 		 glScalef(13.0,13.0,1.0);
	 	 }
		 bool draw = false;
		 float DataMin, DataMax;
		 if(DisplayStrains){
			 DataMin = DisplayStrainRange[0];
			 DataMax = DisplayStrainRange[1];
			 draw = true;
		 }
		 else if (DisplayPysProp){
			 DataMin = DisplayPysPropRange[PysPropToDisplay][0];
			 DataMax = DisplayPysPropRange[PysPropToDisplay][1];
			 draw = true;
		 }
		 if (draw){
			 float barlength = 30.0;
			 float xinit = -barlength/2.0;
			 float y = 17;
			 float dy = 1.0;
			 int slices = 400;
			 float dx = barlength/slices;
			 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			 float dData = (DataMax - DataMin)/slices;
			 float Data = DataMin;
			 float currX = xinit;
			 //cout<<"starins visible: datarange: "<<DataMin<<" "<<DataMax<<" data: "<<Data<<endl;
			 while (Data < DataMax){
				 float* DisplayColour;
				 DisplayColour = new float[3];
				 getDisplayColour(DisplayColour, Data);
				 glColor3f(DisplayColour[0],DisplayColour[1],DisplayColour[2]);
				 glRectf(currX,y,currX+dx,y+dy);
				 Data += dData;
				 currX = currX + dx;
			 }
		 }
	 glPopMatrix();
 }

 void GLWidget::drawFixedNodes(){
 	int n = Sim01->Nodes.size();
 	for(int i = 0; i<n; ++i){
 		bool PositionFixed=false;
 		float colour[3]={0.0,0.0,0.0};
 		if(Sim01->Nodes[i]->FixedPos[0]){	//x-dim fixed on node
 			colour[0]=0.5;
 			PositionFixed=true;
 		}
 		if(Sim01->Nodes[i]->FixedPos[1]){	//y-dim fixed on node
 			colour[1]=0.5;
 			PositionFixed = true;
 		}
 		if(Sim01->Nodes[i]->FixedPos[2]){	//z-dim fixed on node
 			colour[2]=0.5;
 			PositionFixed = true;
 		}

 		if (PositionFixed){
 			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
 			glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
 			glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
 			glDisable(GL_POLYGON_OFFSET_FILL);
 			glColor3f(colour[0],colour[1],colour[2]);
 			if (Sim01->Nodes[i]->nDim == 2){
 				float trianglepoints[3][2];
 				trianglepoints[0][0]=Sim01->Nodes[i]->Position[0]- 0.1;
 				trianglepoints[0][1]=Sim01->Nodes[i]->Position[1]- 0.1;
 				trianglepoints[1][0]=Sim01->Nodes[i]->Position[0]+ 0.1;
 				trianglepoints[1][1]=Sim01->Nodes[i]->Position[1]- 0.1;
 				trianglepoints[2][0]=Sim01->Nodes[i]->Position[0];
 				trianglepoints[2][1]=Sim01->Nodes[i]->Position[1]+ 0.1;
 				glBegin(GL_TRIANGLES);
 					for (int j = 0; j < 3; ++j){
 						glVertex3f( trianglepoints[j][0],  trianglepoints[j][1], 0);
 					}
 				glEnd();
 			}
 			else{
 				//this is the 3D case:
 				float trianglepoints[4][3];
 				trianglepoints[0][0]=Sim01->Nodes[i]->Position[0]- 0.1;
 				trianglepoints[0][1]=Sim01->Nodes[i]->Position[1]- 0.1;
 				trianglepoints[0][2]=Sim01->Nodes[i]->Position[2]- 0.1;
 				trianglepoints[1][0]=Sim01->Nodes[i]->Position[0]+ 0.1;
 				trianglepoints[1][1]=Sim01->Nodes[i]->Position[1]- 0.1;
 				trianglepoints[1][2]=Sim01->Nodes[i]->Position[2];
 				trianglepoints[2][0]=Sim01->Nodes[i]->Position[0];
 				trianglepoints[2][1]=Sim01->Nodes[i]->Position[1]+ 0.1;
 				trianglepoints[2][2]=Sim01->Nodes[i]->Position[2];
 				trianglepoints[3][0]=Sim01->Nodes[i]->Position[0]- 0.1;
 				trianglepoints[3][1]=Sim01->Nodes[i]->Position[1]- 0.1;
 				trianglepoints[3][2]=Sim01->Nodes[i]->Position[2]+ 0.1;
 				int order[4][3]={{0,1,2},{0,2,3},{0,1,3},{1,2,3}};
 				glBegin(GL_TRIANGLES);
 					for (int j = 0; j < 4; ++j){
 						for (int k = 0; k < 3; ++k){
 							glVertex3f( trianglepoints[order[j][k]][0],  trianglepoints[order[j][k]][1], trianglepoints[order[j][k]][2]);
 						}
 					}
 				glEnd();
 			}
 		}
 	}
  }

 void GLWidget::drawScaleBar(){
	glPushMatrix();
		//glTranslatef( -15.0f, +15.0f, -obj_pos[2]);
		glTranslatef( 0.0f, 0.0f, -obj_pos[2]);
		glMultMatrixf(MatRot);
		float size = 5.0; //one side of the cube is 5 microns
		float Points[8][3]={{0,0,0},{size,0,0},{size,size,0},{0,size,0},{0,0,size},{size,0,size},{size,size,size},{0,size,size}};
		int FaceConnectivity[12][3] = {{0,1,2},{0,2,3},{1,2,6},{1,6,5},{4,5,6},{4,6,7},{4,7,3},{0,3,4},{0,1,4},{1,4,5},{2,3,7},{2,7,6}};
		int BorderConnectivity[16] = {0,1,2,3,0,4,5,1,2,6,5,4,7,3,7,6};
		glColor3f(0,0,0);
		glDisable(GL_DITHER);
		glEnable(GL_DEPTH_TEST);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_TRIANGLES);
			for (int j =0; j<12;++j){
				for (int k =0; k<3; ++k){
					int pointId = FaceConnectivity[j][k];
					float x = Points[pointId][0];
					float y = Points[pointId][1];
					float z = Points[pointId][2];
					glVertex3f( x, y, z);
				}
			}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);

		glColor3f(0.4,0.4,0.4);
		glLineWidth(ReferenceLineThickness);
		//Drawing the borders
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBegin(GL_LINE_STRIP);
			for (int j =0; j<16;++j){
				int pointId = BorderConnectivity[j];
				float x = Points[pointId][0];
				float y = Points[pointId][1];
				float z = Points[pointId][2];
				glVertex3f( x, y, z);
			}
		glEnd();
	glPopMatrix();
 }

 void GLWidget::drawAxesArrows(){
	glPushMatrix();
		if (PerspectiveView) {
			glTranslatef( -15.0f, -15.0f, -20.0f);	 }
		else{
			glTranslatef( -200.0f, -200.0f, 0.0f);
			glScalef(15.0,15.0,1.0);
		}
		glMultMatrixf(MatRot);
		glLineWidth(ReferenceLineThickness);
		glBegin(GL_LINES);
			glColor3f(1,0,0);
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(2.0f, 0.0f, 0.0f);
			glColor3f(0,1,0);
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(0.0f, 2.0f, 0.0f);
			glColor3f(0,0,1);
			glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(0.0f, 0.0f, 2.0f);
		glEnd();
	glPopMatrix();
 }

 bool GLWidget::checkPickedColour(int* ElementColour){
	 for (int i=0; i<3; ++i){
		 if (ElementColour[i] != PickedColour[i]){
			  return false;
		 }
	 }
	 return true;

 }

 void GLWidget::drawForPicking(){
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

 void GLWidget::drawTriangleForPicking(int i){
   	//Drawing the surfaces
	const int nTriangle = 1; //a triangle with 3 points needs 1 actual triangle to draw
	int TriangleConnectivity[nTriangle][3] = {{0,1,2}};

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


 void GLWidget::drawForces(){
	 if (drawNetForces){
		 double threshold2 = 1E-16;
		 double minlength = 0.3, maxlength = 2;
		 double minlength2 = minlength*minlength, maxlength2 = maxlength*maxlength;
		 double scale2[2] = {0,10.0}, scale = 10.0;
		 double scalesq = scale*scale;
		 int n = Sim01->Nodes.size();
		 for (int i =0; i<n; ++i){
			 bool drawCurrentNode = checkIfDrawingNode(i);
			 if (drawCurrentNode){
				 double* F;
				 F = new double[3];
				 F[0] = Sim01->SystemForces[0][i][0];
				 F[1] = Sim01->SystemForces[0][i][1];
				 F[2] = Sim01->SystemForces[0][i][2];
				 //cout<<"Force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
				 //check if the force is large enough to display:
				 double mag2 = F[0]* F[0] + F[1]*F[1] + F[2]* F[2];
				 if (mag2 > threshold2){
					 double mag = pow(mag2,0.5);
					 double r = (mag- scale2[0])/(scale2[1]-scale2[0]);
					 double a = mag2/scalesq;
					 double currscale = scale;
					 if (a < minlength2 ){
						 currscale = mag/minlength;
					 }
					 else if ( a > maxlength2){
						 currscale = mag/maxlength;
					 }
					 //cout<<"Element: "<<i<<" F: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<" Fmag: "<<mag<<" scale: "<<currscale<<" r: "<<r <<endl;
					 F[0] =  F[0]/currscale + Sim01->Nodes[i]->Position[0];
					 F[1] =  F[1]/currscale + Sim01->Nodes[i]->Position[1];
					 F[2] =  F[2]/currscale + Sim01->Nodes[i]->Position[2];
					 drawArrow3D(Sim01->Nodes[i]->Position, F, r, 0.0, 0.0);
				 }
				 delete[] F;
			 }
		 }
	 }
 }

 void GLWidget::drawPackForces(){
	 if (drawPackingForces){
		 double threshold2 = 1E-16;
		 double minlength = 0.3, maxlength = 2;
		 double minlength2 = minlength*minlength, maxlength2 = maxlength*maxlength;
		 double scale2[2] = {0,0.1}, scale = 10.0;
		 double scalesq = scale*scale;
		 int n = Sim01->Nodes.size();
		 for (int i =0; i<n; ++i){
			 bool drawCurrentNode = checkIfDrawingNode(i);
			 if (drawCurrentNode){
				 double* F;
				 F = new double[3];
				 F[0] = Sim01->PackingForces[0][i][0];
				 F[1] = Sim01->PackingForces[0][i][1];
				 F[2] = Sim01->PackingForces[0][i][2];
                 //cout<<"Force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
				 //check if the force is large enough to display:
				 double mag2 = F[0]* F[0] + F[1]*F[1] + F[2]* F[2];
				 if (mag2 > threshold2){
					 double mag = pow(mag2,0.5);
					 double r = (mag- scale2[0])/(scale2[1]-scale2[0]);
					 double a = mag2/scalesq;
					 double currscale = scale;
					 if (a < minlength2 ){
						 currscale = mag/minlength;
					 }
					 else if ( a > maxlength2){
						 currscale = mag/maxlength;
					 }
					 F[0] =  F[0]/currscale + Sim01->Nodes[i]->Position[0];
					 F[1] =  F[1]/currscale + Sim01->Nodes[i]->Position[1];
					 F[2] =  F[2]/currscale + Sim01->Nodes[i]->Position[2];
					 drawArrow3D(Sim01->Nodes[i]->Position, F, r, 0.0, 0.0);
				 }
				 delete[] F;
			 }
		 }
	 }

 }
 void GLWidget::drawNodeVelocities(){
	 if (drawVelocities){
		 double threshold2 = 1E-16;
		 double scale2[2] = {0,0.01}, scale = 0.5;
		 double minlength = 0.3, maxlength = 2;
		 double minlength2 = minlength*minlength, maxlength2 = maxlength*maxlength;
		 double scalesq = scale*scale;
		 int n = Sim01->Nodes.size();
		 for (int i =0; i<n; ++i){
			 bool drawCurrentNode = checkIfDrawingNode(i);
			 if (drawCurrentNode){
				 double* v;
				 v = new double[3];
				 v[0] = Sim01->Nodes[i]->Velocity[0][0];
				 v[1] = Sim01->Nodes[i]->Velocity[0][1];
				 v[2] = Sim01->Nodes[i]->Velocity[0][2];
				 //check if the force is large enough to display:
				 double mag2 = v[0]* v[0] + v[1]*v[1] + v[2]* v[2];
				 if (mag2 > threshold2){
					 double mag = pow(mag2,0.5);
					 double a = mag2/scalesq;
					 if (a < minlength2 ){
						 scale = mag /minlength;
					 }
					 else if ( a > maxlength2){
						 scale = mag /maxlength;
					 }
					 double b = (mag- scale2[0])/(scale2[1]-scale2[0]);
					 //cout<<"Colour coding: "<<b<<endl;
					 v[0] = v[0]/scale + Sim01->Nodes[i]->Position[0];
					 v[1] = v[1]/scale + Sim01->Nodes[i]->Position[1];
					 v[2] = v[2]/scale + Sim01->Nodes[i]->Position[2];
					 drawArrow3D(Sim01->Nodes[i]->Position, v, 0.0, 0.0, b);
				 }
				 delete[] v;
			 }
		 }
	 }
 }

 void GLWidget::drawPipette(){
     if (Sim01->PipetteSuction && displayPipette){
         for (int i=0; i<2; i++){
             //drawing inner and outer borders:
             double distance  = Sim01->pipetteRadius;
             glColor3f(0.5,0.5,0.6);
             if (i == 1){
                 distance  += Sim01->pipetteThickness;
                 glColor3f(0.6,0.6,0.7);
             }
             float DEG2RAD = 3.14159/180;
             glLineWidth(ReferenceLineThickness);
             //glLineStipple(1, 0x3F07);
             //glEnable(GL_LINE_STIPPLE);
             glBegin(GL_LINE_LOOP);
                for (int i=0; i < 180; i++){
                   float degInRad = i*DEG2RAD;
                   float x = cos(degInRad)*distance + Sim01->pipetteCentre[0];
                   float y = sin(degInRad)*distance + Sim01->pipetteCentre[1];
                   float z = Sim01->pipetteCentre[2];
                   glVertex3f(x,y,z);
                   glBegin(GL_LINES);
                          glVertex3f(x,y,z);
                          glVertex3f(x,y,z+100);
                   glEnd();
                }
             glEnd();
         }
     }
 }

 void GLWidget::drawBoundingBox(){
	 if (displayBoundingBox){
		 glLineWidth(ReferenceLineThickness);
		 glColor3f(1.0,0.0,1.0);
		 //glLineStipple(1, 0x3F07);
		 //glEnable(GL_LINE_STIPPLE);
		 glBegin(GL_LINE_LOOP);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[0][1],Sim01->boundingBox[0][2]);
		 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[0][1],Sim01->boundingBox[0][2]);
		 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[1][1],Sim01->boundingBox[0][2]);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[1][1],Sim01->boundingBox[0][2]);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[0][1],Sim01->boundingBox[0][2]);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[0][1],Sim01->boundingBox[1][2]);
		 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[0][1],Sim01->boundingBox[1][2]);
		 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[1][1],Sim01->boundingBox[1][2]);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[1][1],Sim01->boundingBox[1][2]);
		 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[0][1],Sim01->boundingBox[1][2]);
		 glEnd();
		 glBegin(GL_LINES);
	 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[0][1],Sim01->boundingBox[0][2]);
	 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[0][1],Sim01->boundingBox[1][2]);
	 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[1][1],Sim01->boundingBox[0][2]);
	 		glVertex3f(Sim01->boundingBox[1][0],Sim01->boundingBox[1][1],Sim01->boundingBox[1][2]);
	 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[1][1],Sim01->boundingBox[0][2]);
	 		glVertex3f(Sim01->boundingBox[0][0],Sim01->boundingBox[1][1],Sim01->boundingBox[1][2]);
		 glEnd();
		 //glDisable(GL_LINE_STIPPLE);
	 }
 }

 void GLWidget::drawArrow3D(double* pos, double* endPoint, double r, double g, double b){
	glColor3f(r,g,b);
	double v[3],u[3],cross[3];
	v[0]= 0.2 * (endPoint[0]-pos[0]);
	v[1]= 0.2 * (endPoint[1]-pos[1]);
	v[2]= 0.2 * (endPoint[2]-pos[2]);
	//this is the 3D case:
	glLineWidth(ReferenceLineThickness);
	glBegin(GL_LINES);
		glVertex3f(pos[0],pos[1],pos[2]);
		glVertex3f(endPoint[0] - v[0], endPoint[1] - v[1], endPoint[2] - v[2]);
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL); // Avoid Stitching!
	glPolygonOffset(1.0, 1.0); 	//These are necessary so the depth test can keep the lines above surfaces
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (v[0]>1E-8){
		u[0] = v[0];
		u[1] = -v[2];
		u[2] = v[1];
	}
	else if (v[1]>1E-8){
		u[0] = v[0];
		u[1] = -v[2];
		u[2] = v[1];
	}
	else {
		u[0] = -v[1];
		u[1] = v[0];
		u[2] = v[2];
	}
	//get the normal to the two vectors you found:
	cross[0] = u[1]*v[2] - u[2]*v[1];
	cross[1] = u[2]*v[0] - u[0]*v[2];
	cross[2] = u[0]*v[1] - u[1]*v[0];
	//now I will define the pyramid with these three vectors.
	float trianglepoints[5][3];
	trianglepoints[0][0]=endPoint[0];
	trianglepoints[0][1]=endPoint[1];
	trianglepoints[0][2]=endPoint[2];
	trianglepoints[1][0]=endPoint[0] - v[0] - u[0];
	trianglepoints[1][1]=endPoint[1] - v[1] - u[1];
	trianglepoints[1][2]=endPoint[2] - v[2] - u[2];
	trianglepoints[2][0]=endPoint[0] - v[0] + u[0];
	trianglepoints[2][1]=endPoint[1] - v[1] + u[1];
	trianglepoints[2][2]=endPoint[2] - v[2] + u[2];
	trianglepoints[3][0]=endPoint[0] - v[0] - cross[0];
	trianglepoints[3][1]=endPoint[1] - v[1] - cross[1];
	trianglepoints[3][2]=endPoint[2] - v[2] - cross[2];
	trianglepoints[4][0]=endPoint[0] - v[0] + cross[0];
	trianglepoints[4][1]=endPoint[1] - v[1] + cross[1];
	trianglepoints[4][2]=endPoint[2] - v[2] + cross[2];
	int order[6][3]={{0,1,3},{0,1,4},{0,2,3},{0,2,4},{1,3,4},{2,3,4}};
	glBegin(GL_TRIANGLES);
		for (int j = 0; j < 6; ++j){
			for (int k = 0; k < 3; ++k){
				glVertex3f( trianglepoints[order[j][k]][0],  trianglepoints[order[j][k]][1], trianglepoints[order[j][k]][2]);
			}
		}
	glEnd();
	glColor3f(1.0,1.0,1.0);
	glLineWidth(MainShapeLineThickness);
	glBegin(GL_LINES);
		glVertex3f(trianglepoints[0][0], trianglepoints[0][1], trianglepoints[0][2]);
		glVertex3f(trianglepoints[1][0], trianglepoints[1][1], trianglepoints[1][2]);
		glVertex3f(trianglepoints[0][0], trianglepoints[0][1], trianglepoints[0][2]);
		glVertex3f(trianglepoints[2][0], trianglepoints[2][1], trianglepoints[2][2]);
		glVertex3f(trianglepoints[0][0], trianglepoints[0][1], trianglepoints[0][2]);
		glVertex3f(trianglepoints[3][0], trianglepoints[3][1], trianglepoints[3][2]);
		glVertex3f(trianglepoints[0][0], trianglepoints[0][1], trianglepoints[0][2]);
		glVertex3f(trianglepoints[4][0], trianglepoints[4][1], trianglepoints[4][2]);
	glEnd();
 }



