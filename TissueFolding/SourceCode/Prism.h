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
	boost::numeric::ublas::matrix<double> kTopAlignedBuffer;
	boost::numeric::ublas::matrix<double> BTopAlignedBuffer;
	boost::numeric::ublas::matrix<double> BETopAlignedBuffer;
	boost::numeric::ublas::matrix<double> kBottomAlignedBuffer;
	boost::numeric::ublas::matrix<double> BBottomAlignedBuffer;
	boost::numeric::ublas::matrix<double> BEBottomAlignedBuffer;
	double** RefShapePosBottomAlignedBuffer;
	double** RefShapePosTopAlignedBuffer;



	boost::numeric::ublas::matrix<double> TissueToWorldRotMatTopAlignedBuffer;
	boost::numeric::ublas::matrix<double> RefToTissueRotMatTopAlignedBuffer;
	boost::numeric::ublas::matrix<double> RefToTissueRotMatTTopAlignedBuffer;
	double* TissueCoordinateSystemTopAlignedBuffer;

	boost::numeric::ublas::matrix<double> TissueToWorldRotMatBottomAlignedBuffer;
	boost::numeric::ublas::matrix<double> RefToTissueRotMatBottomAlignedBuffer;
	boost::numeric::ublas::matrix<double> RefToTissueRotMatTBottomAlignedBuffer;
	double* TissueCoordinateSystemBottomAlignedBuffer;


	void setNormals();
	void setRefShapePosBuffers();
	void setTissueCoordsRotationsBuffers();
	void getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape);
	void setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer,double eta, double zeta, double nu);
	void setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack);
	void setCoeffMat();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, double eta, double zeta, double nu);
	void calculateNormals();
	void calculateNormalToBottom();
	void calculateReferenceNormalToBottom();
	void calculateNormalToTop();
	void calculateReferenceNormalToTop();
	void getCurrentAlignmentSides(double*, double*);
	void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	void updateAlignmentTurn();
	void updateReferenceShapeBaseFromBuffer();
	void setStiffnessMatrixBuffers();
	void resetBuffersAfterGrowth();
	void calculateZVecForTissueCoordAlignment(double* u);
	void calculateXVecForTissueCoordAlignment(double* u);
public:
	Prism(int* NodeIds,vector<Node*>& Nodes, int CurrId);
	~Prism();
	void  setElasticProperties(double E,double v);
	void  setViscosity(double ApicalVisc,double BasalVisc, vector <Node*>& Nodes);
	void  calculateReferenceStiffnessMatrix();
	void  calculateForces(int RKid, double** SystemForces);

};

#endif
