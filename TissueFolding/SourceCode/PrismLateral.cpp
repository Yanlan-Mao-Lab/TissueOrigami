/*
 * PrismLateral.cpp
 *
 *  Created on: 23 Sep 2014
 *      Author: melda
 */

#include "PrismLateral.h"
#include "ReferenceShapeBase.h"


using namespace std;

PrismLateral::PrismLateral(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId): Prism(tmpNodeIds, Nodes, CurrId){
	//cout<<"constructed lateral prism"<<endl;
}

PrismLateral::~PrismLateral(){
	//cout<<"called the destructor for prism class"<<endl;
	//Only delete what you have constructed in the child element's constructor
}


void  PrismLateral::setElasticProperties(double E, double v){
	this -> E = E;
	this -> v = v; //poisson ratio
	if (v>0.5){v = 0.5;}
	else if (v<0.0){v = 0.0;}
	using namespace boost::numeric::ublas;
	D = zero_matrix<double>(6,6);
	double multiplier = E/((1+v)*(1-2*v));
	D(0,0)= multiplier*(1-v);	D(0,1)=	multiplier*v;		D(0,2)=	multiplier*v;
	D(1,0)= multiplier*v;		D(1,1)= multiplier*(1-v);	D(1,2)= multiplier*v;
	D(2,0)= multiplier*v;		D(2,1)= multiplier*v;		D(2,2)= multiplier*(1-v);
	D(3,3)= multiplier*(1-2*v)/2;
	D(4,4)= multiplier*(1-2*v)/2;
	D(5,5)= multiplier*(1-2*v)/2;

	//making it more resistive in Ex,y
	D(1,1) *= 3.0;
	D(2,2) *= 3.0;
}

void PrismLateral::setViscosity(double ApicalVisc,double BasalVisc, vector <Node*>& Nodes){
	Nodes[NodeIds[0]]->Viscosity = BasalVisc;
	Nodes[NodeIds[3]]->Viscosity = BasalVisc;
	Nodes[NodeIds[1]]->Viscosity = ApicalVisc;
	Nodes[NodeIds[4]]->Viscosity = ApicalVisc;
	double avrVisc = (ApicalVisc + BasalVisc)/2.0;
	Nodes[NodeIds[2]]->Viscosity = avrVisc;
	Nodes[NodeIds[5]]->Viscosity = avrVisc;
}

void PrismLateral::calculateNormals(){
	if (alingmentTurn == 0){
		calculateNormalToBasalSide();
		calculateReferenceNormalToBasalSide();
	}
	else{
		calculateNormalToApicalSide();
		calculateReferenceNormalToApicalSide();
	}
}

void PrismLateral::calculateNormalToBasalSide(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = PositionsAlignedToReference[3][0] - PositionsAlignedToReference[0][0];
	vec1[1] = PositionsAlignedToReference[3][1] - PositionsAlignedToReference[0][1];
	vec1[2] = PositionsAlignedToReference[3][2] - PositionsAlignedToReference[0][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = PositionsAlignedToReference[5][0] - PositionsAlignedToReference[0][0];
	vec2[1] = PositionsAlignedToReference[5][1] - PositionsAlignedToReference[0][1];
	vec2[2] = PositionsAlignedToReference[5][2] - PositionsAlignedToReference[0][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	//now I have the normal, but the normal must look into the prism
	//For the bottom case, the vector should have a less than 90 degree angle
	//with the vector towards the top. lets say from node 0 to 3:
	double* vec3;
	vec3 = new double[3];
	vec3[0] = PositionsAlignedToReference[1][0] - PositionsAlignedToReference[0][0];
	vec3[1] = PositionsAlignedToReference[1][1] - PositionsAlignedToReference[0][1];
	vec3[2] = PositionsAlignedToReference[1][2] - PositionsAlignedToReference[0][2];
	double dotp = dotProduct3D(vec3, CurrentNormal);
	if(dotp < 0){
		CurrentNormal[0] *= (-1.0);
		CurrentNormal[1] *= (-1.0);
		CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(CurrentNormal);
	//cout<<"Normal calculated: "<<CurrentNormal[0]<<" "<<CurrentNormal[1]<<" "<<CurrentNormal[2]<<endl;
}

void PrismLateral::calculateReferenceNormalToBasalSide(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
	vec1[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
	vec1[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[0][0];
	vec2[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[0][1];
	vec2[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[0][2];

	crossProduct3D(vec1,vec2,ReferenceShape->CurrentNormal);
	//now I have the normal, but the normal must look into the prism
	//For the bottom case, the vector should have a less than 90 degree angle
	//with the vector towards the top. lets say from node 0 to 3:
	double* vec3;
	vec3 = new double[3];
	vec3[0] = ReferenceShape->Positions[1][0] - ReferenceShape->Positions[0][0];
	vec3[1] = ReferenceShape->Positions[1][1] - ReferenceShape->Positions[0][1];
	vec3[2] = ReferenceShape->Positions[1][2] - ReferenceShape->Positions[0][2];
	double dotp = dotProduct3D(vec3, ReferenceShape->CurrentNormal);
	if(dotp < 0 ){
		ReferenceShape->CurrentNormal[0] *= (-1.0);
		ReferenceShape->CurrentNormal[1] *= (-1.0);
		ReferenceShape->CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(ReferenceShape->CurrentNormal);
}

void PrismLateral::calculateNormalToApicalSide(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = PositionsAlignedToReference[4][0] - PositionsAlignedToReference[1][0];
	vec1[1] = PositionsAlignedToReference[4][1] - PositionsAlignedToReference[1][1];
	vec1[2] = PositionsAlignedToReference[4][2] - PositionsAlignedToReference[1][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = PositionsAlignedToReference[5][0] - PositionsAlignedToReference[1][0];
	vec2[1] = PositionsAlignedToReference[5][1] - PositionsAlignedToReference[1][1];
	vec2[2] = PositionsAlignedToReference[5][2] - PositionsAlignedToReference[1][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	double* vec3;
	vec3 = new double[3];
	vec3[0] = PositionsAlignedToReference[0][0] - PositionsAlignedToReference[1][0];
	vec3[1] = PositionsAlignedToReference[0][1] - PositionsAlignedToReference[1][1];
	vec3[2] = PositionsAlignedToReference[0][2] - PositionsAlignedToReference[1][2];
	double dotp = dotProduct3D(vec3, CurrentNormal);
	if(dotp < 0 ){
		CurrentNormal[0] *= (-1.0);
		CurrentNormal[1] *= (-1.0);
		CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(CurrentNormal);
	//cout<<"Normal calculated: "<<CurrentNormal[0]<<" "<<CurrentNormal[1]<<" "<<CurrentNormal[2]<<endl;
}

void PrismLateral::calculateReferenceNormalToApicalSide(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = ReferenceShape->Positions[4][0] - ReferenceShape->Positions[1][0];
	vec1[1] = ReferenceShape->Positions[4][1] - ReferenceShape->Positions[1][1];
	vec1[2] = ReferenceShape->Positions[4][2] - ReferenceShape->Positions[1][2];
	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[1][0];
	vec2[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[1][1];
	vec2[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[1][2];

	crossProduct3D(vec1,vec2,ReferenceShape->CurrentNormal);
	double* vec3;
	vec3 = new double[3];
	vec3[0] = ReferenceShape->Positions[0][0] - ReferenceShape->Positions[1][0];
	vec3[1] = ReferenceShape->Positions[0][1] - ReferenceShape->Positions[1][1];
	vec3[2] = ReferenceShape->Positions[0][2] - ReferenceShape->Positions[1][2];
	double dotp = dotProduct3D(vec3, ReferenceShape->CurrentNormal);
	if(dotp < 0 ){
		ReferenceShape->CurrentNormal[0] *= (-1.0);
		ReferenceShape->CurrentNormal[1] *= (-1.0);
		ReferenceShape->CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(ReferenceShape->CurrentNormal);
}

void PrismLateral::getCurrentAlignmentSides(double* RefSide, double* ShapeSide){
	//side vector for alignment:
	if (alingmentTurn == 0){
		ShapeSide[0] = PositionsAlignedToReference[3][0] - PositionsAlignedToReference[0][0];
		ShapeSide[1] = PositionsAlignedToReference[3][1] - PositionsAlignedToReference[0][1];
		ShapeSide[2] = PositionsAlignedToReference[3][2] - PositionsAlignedToReference[0][2];
		RefSide[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
		RefSide[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
		RefSide[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	}
	else{
		ShapeSide[0] = PositionsAlignedToReference[4][0] - PositionsAlignedToReference[1][0];
		ShapeSide[1] = PositionsAlignedToReference[4][1] - PositionsAlignedToReference[1][1];
		ShapeSide[2] = PositionsAlignedToReference[4][2] - PositionsAlignedToReference[1][2];
		RefSide[0] = ReferenceShape->Positions[4][0] - ReferenceShape->Positions[1][0];
		RefSide[1] = ReferenceShape->Positions[4][1] - ReferenceShape->Positions[1][1];
		RefSide[2] = ReferenceShape->Positions[4][2] - ReferenceShape->Positions[1][2];
	}
	normaliseVector3D(RefSide);
	normaliseVector3D(ShapeSide);
}

void PrismLateral::getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace){
	double* Ref2ndSide;
	Ref2ndSide = new double[3];
	double* Shape2ndSide;
	Shape2ndSide = new double[3];
	//side vector for alignment:
	if (alingmentTurn == 0){
		Shape2ndSide[0] = PositionsAlignedToReference[5][0] - PositionsAlignedToReference[0][0];
		Shape2ndSide[1] = PositionsAlignedToReference[5][1] - PositionsAlignedToReference[0][1];
		Shape2ndSide[2] = PositionsAlignedToReference[5][2] - PositionsAlignedToReference[0][2];
		Ref2ndSide[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[0][0];
		Ref2ndSide[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[0][1];
		Ref2ndSide[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[0][2];
	}
	else{
		Shape2ndSide[0] = PositionsAlignedToReference[5][0] - PositionsAlignedToReference[1][0];
		Shape2ndSide[1] = PositionsAlignedToReference[5][1] - PositionsAlignedToReference[1][1];
		Shape2ndSide[2] = PositionsAlignedToReference[5][2] - PositionsAlignedToReference[1][2];
		Ref2ndSide[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[1][0];
		Ref2ndSide[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[1][1];
		Ref2ndSide[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[1][2];
	}
	crossProduct3D(ShapeSide,Shape2ndSide,ShapeFace);
	crossProduct3D(RefSide,Ref2ndSide,RefFace);
	cout<<"Shape2ndSide: "<<Shape2ndSide[0]<<" "<<Shape2ndSide[1]<<" "<<Shape2ndSide[2]<<endl;
	cout<<"Ref2ndSide: "<<Ref2ndSide[0]<<" "<<Ref2ndSide[1]<<" "<<Ref2ndSide[2]<<endl;
	cout<<"ShapeFace: "<<ShapeFace[0]<<" "<<ShapeFace[1]<<" "<<ShapeFace[2]<<endl;
	cout<<"RefFace: "<<RefFace[0]<<" "<<RefFace[1]<<" "<<RefFace[2]<<endl;
}

void PrismLateral::calculateXVecForTissueCoordAlignment(double* u ){
	u[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
	u[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
	u[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	normaliseVector3D(u);
}
