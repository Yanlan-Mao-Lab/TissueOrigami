#ifndef Prism_H
#define Prism_H


#include "ShapeBase.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class Prism : public ShapeBase{

protected:
    void setTissueCoordsRotationsBuffers();
    void getCurrRelaxedShape(gsl_matrix * CurrRelaxedShape);
    void setShapeFunctionDerivatives(gsl_matrix * ShapeFuncDer,double eta, double zeta, double nu);
    void setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack);
	void setCoeffMat();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double zeta, double nu);
    void calculateCurrNodalForces(gsl_matrix *gslcurrg, gsl_matrix *gslcurrF, int pointNo);

    void calculateNormalToBottom();
	void calculateReferenceNormalToBottom();
	void calculateNormalToTop();
	void calculateReferenceNormalToTop();
	void getCurrentAlignmentSides(double*, double*);
	void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	void updateAlignmentTurn();
	//void updateReferenceShapeBaseFromBuffer();
	//void resetBuffersAfterGrowth();
	//void calculateZVecForTissueCoordAlignment(double* u);
	//void calculateXVecForTissueCoordAlignment(double* u);
	void calculateReferenceVolume();
	void calculatePlaneNormals(double** normals);
	void assignNodalVector(double* vec, int id0, int id1);
	bool checkNodePlaneConsistency(double** normals);
	double getApicalSideLengthAverage();

public:
	Prism(int* NodeIds,vector<Node*>& Nodes, int CurrId);
	~Prism();
	void  setElasticProperties(double EApical, double EBasal, double EMid,double v);
	void  calculateBasalNormal(double * normal);
	void  AlignReferenceBaseNormalToZ();
    void  calculateElementShapeFunctionDerivatives();

	//void  calculateForces(int RKid, double** SystemForces);
	void  checkHealth();
	void getApicalTriangles(vector <int> &ApicalTriangles);
	int getCorrecpondingApical(int currNodeId);
	bool IsThisNodeMyBasal(int currNodeId);
	double getElementHeight();
	void AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, int RKId,  double ***SystemForces,  double ***PackingForces,vector<Node*> &Nodes);
	bool IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode);
	void calculateNormalForPacking(int tissuePlacement);
	void calculateApicalArea();
	void calculateBasalArea();
	void calculateMyosinForces();
	void distributeMyosinForce(bool isIsotropic, bool apical);

	void getApicalNodePos(double* posCorner);
	void getBasalNodePos(double* posCorner);
	bool IspointInsideTriangle(int tissueplacement,double x, double y,double z);
	void checkRotationConsistency3D();

};

#endif
