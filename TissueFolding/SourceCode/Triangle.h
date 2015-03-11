/*
 * Triangle.h
 *
 *  Created on: 6 Jan 2015
 *      Author: melda
 */

#ifndef TRIANGLE_H
#define TRIANGLE_H


#include "ShapeBase.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class Triangle : public ShapeBase{

protected:
	int normalCrossOrder[2];
	double slabHeight;
	//double 	apicalZDir;
	void setTissueCoordsRotationsBuffers();
	void getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape);
	void setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer,double eta, double nu);
	void setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack);
	void setCoeffMat();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double nu);
	//void calculateNormalToBottom();
	//void calculateReferenceNormalToBottom();
	//void calculateNormalToTop();
	//void calculateReferenceNormalToTop();
	//void getCurrentAlignmentSides(double*, double*);
	//void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	//void updateAlignmentTurn();
	void calculateReferenceVolume();
	double getApicalSideLengthAverage();
	double getElementHeight();

public:
	Triangle(int* NodeIds,vector<Node*>& Nodes, int CurrId, double h);
	~Triangle();
	void  setElasticProperties(double EApical, double EBasal, double EMid,double v);
	void  calculateApicalNormalCrossOrder(double* SystemCentre);
	void  AlignReferenceApicalNormalToZ(double* SystemCentre);
	void  correctFor2DAlignment();
	void  calculateReferenceStiffnessMatrix();
	void  checkHealth();
	void AddPackingToApicalSurface(double Fx, double Fy,double Fz, int RKId,  double ***SystemForces, double ***PackingForces, vector<Node*> &Nodes);
	bool IsPointCloseEnoughForPacking(double* Pos, float threshold);
	void getApicalNodePos(double* posCorner);
	void calculateNormalForPacking();
	bool IspointInsideApicalTriangle(double x, double y,double z);
};

#endif /* TRIANGLE_H_ */
