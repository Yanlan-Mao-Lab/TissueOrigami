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
	cout<<"constructed lateral prism"<<endl;
}
/*	nNodes = 6;
	nDim = 3;
	Id = CurrId;
	NodeIds = new int[6];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[3];
	CurrGrowthStrainAddition = new double[3];
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
	}
	updatedReference = false;
	TissueCoordinateSystemUpToDate = false;
	CurrShapeChangeStrainsUpToDate = false;
	CurrGrowthStrainsUpToDate = false;
	IsGrowing = false;
	IsChangingShape = false;
	GrewInThePast = false;
	ChangedShapeInThePast = false;

	setIdentificationColour();
	setShapeType("PrismLateral");
	ReferenceShape = new ReferenceShapeBase("PrismLateral");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	alingmentTurn=0;
	setNormals();
	setRefShapePosBuffers();


	Strain = boost::numeric::ublas::zero_vector<double>(6);
	StrainTissueMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	PlasticStrain = boost::numeric::ublas::zero_vector<double>(6);
	LocalGrowthStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalShapeChangeStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalPlasticStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrGrowthStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrShapeChangeStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrPlasticStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);

	CurrShapeChangeToAdd[0] = 0;
	CurrShapeChangeToAdd[1] = 0;
	CurrShapeChangeToAdd[2] = 0;
	TissueCoordinateSystem = new double[9];
	TissueCoordinateSystem[0]=1.0;
	TissueCoordinateSystem[1]=0.0;
	TissueCoordinateSystem[2]=0.0;
	TissueCoordinateSystem[3]=0.0;
	TissueCoordinateSystem[4]=1.0;
	TissueCoordinateSystem[5]=0.0;
	TissueCoordinateSystem[6]=0.0;
	TissueCoordinateSystem[7]=0.0;
	TissueCoordinateSystem[8]=1.0;
	RefToTissueRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatT = boost::numeric::ublas::zero_matrix<double>(3,3);
	TissueToWorldRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	setTissueCoordsRotationsBuffers();
}
*/

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
	vec1[0] = Positions[3][0] - Positions[0][0];
	vec1[1] = Positions[3][1] - Positions[0][1];
	vec1[2] = Positions[3][2] - Positions[0][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = Positions[5][0] - Positions[0][0];
	vec2[1] = Positions[5][1] - Positions[0][1];
	vec2[2] = Positions[5][2] - Positions[0][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	//now I have the normal, but the normal must look into the prism
	//For the bottom case, the vector should have a less than 90 degree angle
	//with the vector towards the top. lets say from node 0 to 3:
	double* vec3;
	vec3 = new double[3];
	vec3[0] = Positions[1][0] - Positions[0][0];
	vec3[1] = Positions[1][1] - Positions[0][1];
	vec3[2] = Positions[1][2] - Positions[0][2];
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
	vec1[0] = Positions[4][0] - Positions[1][0];
	vec1[1] = Positions[4][1] - Positions[1][1];
	vec1[2] = Positions[4][2] - Positions[1][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = Positions[5][0] - Positions[1][0];
	vec2[1] = Positions[5][1] - Positions[1][1];
	vec2[2] = Positions[5][2] - Positions[1][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	double* vec3;
	vec3 = new double[3];
	vec3[0] = Positions[0][0] - Positions[1][0];
	vec3[1] = Positions[0][1] - Positions[1][1];
	vec3[2] = Positions[0][2] - Positions[1][2];
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
		ShapeSide[0] = Positions[3][0] - Positions[0][0];
		ShapeSide[1] = Positions[3][1] - Positions[0][1];
		ShapeSide[2] = Positions[3][2] - Positions[0][2];
		RefSide[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
		RefSide[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
		RefSide[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	}
	else{
		ShapeSide[0] = Positions[4][0] - Positions[1][0];
		ShapeSide[1] = Positions[4][1] - Positions[1][1];
		ShapeSide[2] = Positions[4][2] - Positions[1][2];
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
		Shape2ndSide[0] = Positions[5][0] - Positions[0][0];
		Shape2ndSide[1] = Positions[5][1] - Positions[0][1];
		Shape2ndSide[2] = Positions[5][2] - Positions[0][2];
		Ref2ndSide[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[0][0];
		Ref2ndSide[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[0][1];
		Ref2ndSide[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[0][2];
	}
	else{
		Shape2ndSide[0] = Positions[5][0] - Positions[1][0];
		Shape2ndSide[1] = Positions[5][1] - Positions[1][1];
		Shape2ndSide[2] = Positions[5][2] - Positions[1][2];
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

void PrismLateral::calculateZVecForTissueCoordAlignment(double* u){
	u[0] = ReferenceShape->Positions[1][0] - ReferenceShape->Positions[0][0];
	u[1] = ReferenceShape->Positions[1][1] - ReferenceShape->Positions[0][1];
	u[2] = ReferenceShape->Positions[1][2] - ReferenceShape->Positions[0][2];
	normaliseVector3D(u);
}

void PrismLateral::calculateXVecForTissueCoordAlignment(double* u ){
	u[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
	u[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
	u[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	normaliseVector3D(u);
}
