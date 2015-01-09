/*
 * Tetrahedron.h
 *
 *  Created on: 19 Dec 2014
 *      Author: melda
 */

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "ShapeBase.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class Tetrahedron : public ShapeBase{

protected:
	void setTissueCoordsRotationsBuffers();
	void getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape);
	void setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer,double eta, double zeta, double nu);
	void setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack);
	void setCoeffMat();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double zeta, double nu);
	void calculateNormalToBottom();
	void calculateReferenceNormalToBottom();
	void calculateNormalToTop();
	void calculateReferenceNormalToTop();
	void getCurrentAlignmentSides(double*, double*);
	void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	void updateAlignmentTurn();
	void updateReferenceShapeBaseFromBuffer();
	void resetBuffersAfterGrowth();
	void calculateZVecForTissueCoordAlignment(double* u);
	void calculateXVecForTissueCoordAlignment(double* u);
	void calculateReferenceVolume();

public:
	Tetrahedron(int* NodeIds,vector<Node*>& Nodes, int CurrId);
	~Tetrahedron();
	void  setElasticProperties(double EApical, double EBasal, double EMid,double v);
	//void  calculateBasalNormal(double * normal);
	//void  AlignReferenceBaseNormalToZ();
	void  calculateReferenceStiffnessMatrix();
	void  calculateForces(int RKid, double** SystemForces);
	void  checkHealth();
};

#endif /* TETRAHEDRON_H */
