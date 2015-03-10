#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include </usr/include/gsl/gsl_matrix.h>
#include </usr/include/gsl/gsl_linalg.h>

#include "Node.h"
#include "ReferenceShapeBase.h"

using namespace std;

class ShapeBase{
private:
	void ParentErrorMessage();
protected:
	int 	tissuePlacement; //1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
	int 	ShapeType;

	int 	nNodes;

	int 	nDim;
	int* 	IdentifierColour;
	double* GrowthRate;
	double* ShapeChangeRate;
	double* CurrGrowthStrainAddition;
	bool 	CurrShapeChangeStrainsUpToDate;
	bool 	CurrGrowthStrainsUpToDate;
	bool 	GrewInThePast;
	bool 	ChangedShapeInThePast;
	void 	setShapeType(string TypeName);
	void 	readNodeIds(int* tmpNodeIds);
	void 	setPositionMatrix(vector<Node*>& Nodes);
	void 	setTissuePlacement(vector<Node*>& Nodes);
	void 	setTissueType(vector<Node*>& Nodes);
	void 	setReferencePositionMatrix();
	void 	setIdentificationColour();
	void 	rotateReferenceElementByRotationMatrix(double* rotMat);
	bool 	InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);
	int 	determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm);
	void	crossProduct3D(double* u, double* v, double* cross);
	double  dotProduct3D(double* u, double* v);
	void 	updateNodeIdsFromSave(ifstream& file);
	void 	updateReferencePositionMatrixFromSave(ifstream& file);
	virtual void calculateReferenceVolume(){ParentErrorMessage();};
	//bool	areSidesFacingSameDirection(double* RefSide, double* ShapeSide);
	void 	updateTissueCoordStrain();
	void 	updateTissueCoordPlasticStrain();
	void 	updateTissueCoordinateSystem();
	void 	updatePositionsAlignedToReferenceWithBuffers();
	void 	calculateGrowthInLocalCoordinates(double* strainsToAdd);
	void 	calculateReferenceCoordSysAlignedToTissue(double* RefCoords);
	bool 	calculateWorldToTissueRotMat(double* v);
	bool 	calculateGrowthStrainsRotMat(double* v);
	void 	updateTissueCoordinateSystem(double* TissueCoords);
	void	calculateForces2D(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile);
	void	calculateForces3D(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile);
	void 	calculatedudEdXde(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde);
	void 	dudEdXde3D(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde);
	void 	dudEdXde2D(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde);

	boost::numeric::ublas::matrix<double> D;
	boost::numeric::ublas::matrix<int> CoeffMat;
	boost::numeric::ublas::matrix<double> k;
	boost::numeric::ublas::matrix<double> B;
	boost::numeric::ublas::matrix<double> BE;
	boost::numeric::ublas::matrix<double> Bo;
	boost::numeric::ublas::vector<double> Forces;

	double E, v;




public:
	int 	Id;
	boost::numeric::ublas::matrix<double> LocalGrowthStrainsMat;
	int* 	NodeIds;
	virtual ~ShapeBase(){
			//while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
			};
	double** Positions;
	ReferenceShapeBase* ReferenceShape;
	boost::numeric::ublas::vector<double> Strain;
	boost::numeric::ublas::vector<double> RK1Strain;
	boost::numeric::ublas::matrix<double> StrainTissueMat;
	boost::numeric::ublas::vector<double> PlasticStrain;
	boost::numeric::ublas::matrix<double> CurrPlasticStrainsInTissueCoordsMat;
	bool 	IsGrowing;
	bool 	IsChangingShape;
	bool 	WorldToTissueRotMatUpToDate;
	bool 	GrowthStrainsRotMatUpToDate;
	bool	NormalForPackingUpToDate;
	int 	tissueType;	//Columnar layer = 0, peripodium = 1
	bool	IsAblated;
	double 	CurrShapeChangeToAdd[3];
	//int alingmentTurn;
	//double* CurrentNormal;
	//bool updatedReference;
	double* TissueCoordinateSystem;
	double* normalForPacking;
	int 	getId();
	string 	getName();
	int 	getShapeType();
	int 	getNodeNumber();
	int* 	getNodeIds();
	int 	getDim();
	int* 	getIdentifierColour();
	double* getCentre();
	void 	getStrain(int type, float &StrainMag);
	void 	getPlasticStrain(int type, float &StrainMag);
	void	getTissueCoordinaSystem(double* TissueCoords);
	void 	getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag);
	void 	getPysProp(int type, float &PysPropMag);
	double 	getYoungModulus();
	double 	getPoissonRatio();
	double* getGrowthRate();
	double* getShapeChangeRate();
	double** getReferencePos();
	double*	 getReferenceNormal();
	void 	displayName();
	void	displayNodeIds();
	void 	displayPositions();
	void 	displayReferencePositions();
	void 	displayIdentifierColour();
	virtual void setElasticProperties(double EApical,double EBasal, double EMid,double v){ParentErrorMessage();};
	virtual void calculateBasalNormal(double * normal){ParentErrorMessage();};
	virtual void AlignReferenceBaseNormalToZ(){ParentErrorMessage();};
	void 	updateGrowthRate(double scalex, double scaley, double scalez);
	virtual void calculateReferenceStiffnessMatrix(){ParentErrorMessage();};
	void 	calculateForces(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile);
	void 	updatePositions(int RKId, vector<Node*>& Nodes);
	void 	setGrowthRate(double x, double y, double z);
	void 	setShapeChangeRate(double x, double y, double z);
	void 	updateGrowthToAdd(double* growthscale);
	void 	updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes);

	virtual void  checkHealth(){ParentErrorMessage();};
	void 	resetCurrStepGrowthData();
	void 	resetCurrStepShapeChangeData();
	//void 	calculatePlasticStrain();
	//void 	changeShape(double shapechangescale, int axis);
	void 	updateShapeFromSave(ifstream& file);
	void 	displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
	void	normaliseVector3D(double* v);
	double 	determinant3by3Matrix(double* rotMat);
	double 	determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat);
	double 	determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat);
	void	calculateRotationAngleSinCos(double* u, double* v, double& c, double& s);
	void	calculateRotationAxis(double* u, double* v,double* rotAx);
	void	constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,double* rotMat);


	void alignElementOnReference();
	virtual void correctFor2DAlignment(){ParentErrorMessage();};
	virtual double getApicalSideLengthAverage(){ParentErrorMessage();};
	virtual void getApicalTriangles(vector <int> &ApicalTriangles){ParentErrorMessage();};
	virtual int getCorrecpondingApical(int currNodeId){ParentErrorMessage();};
	virtual bool IsThisNodeMyBasal(int currNodeId){ParentErrorMessage();};
	virtual double getElementHeight(){ParentErrorMessage();};
	virtual bool IsPointCloseEnoughForPacking(double* Pos, float threshold){ParentErrorMessage();};
	virtual void calculateNormalForPacking(){ParentErrorMessage();};
	virtual void AddPackingToApicalSurface(double Fx, double Fy,double Fz, int RKId,  double ***SystemForces, double ***PackingForces, vector<Node*> &Nodes){ParentErrorMessage();};
	virtual void getApicalNodePos(double* posCorner){ParentErrorMessage();};
	virtual bool IspointInsideApicalTriangle(double x, double y,double z){ParentErrorMessage();};
	bool DoesPointBelogToMe(int IdNode);
	void updatePositionsAlignedToReferenceForRK();
	void growShape();
	void assignVolumesToNodes(vector <Node*>& Nodes);
	void assignElementToConnectedNodes(vector <Node*>& Nodes);
	void removeMassFromNodes(vector <Node*>& Nodes);
	bool RotatedElement;
	boost::numeric::ublas::matrix<double> WorldToTissueRotMat;
	double **PositionsInTissueCoord;
	double **PositionsAlignedToReference;
	boost::numeric::ublas::matrix<double> WorldToReferenceRotMat;
	boost::numeric::ublas::matrix<double> GrowthStrainsRotMat;

	//bool 	calculateAlignmentRotationMatrix(double** RefNormalised, double* rotMat);
	bool 	calculateAlignmentScore(double** RefNormalised);
	void 	bringShapePositionsToOrigin(double** RefNormalised, double* refCentre);
	void 	bringPositionAlignedToReferenceToOrigin(double* refCentre);
	bool	calculateDisplacementGradientRotationMatrix(double** RefNormalised, double* rotMat);
	void 	updateElementsNodePositions(int RKId, double ***SystemForces, vector <Node*>& Nodes, double dt);
	void 	updateReferencePositionMatrixFromMeshInput(ifstream& file);
	void	fillNodeNeighbourhood(vector<Node*>& Nodes);
};

#endif
