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
private:
	void calculateReferenceStiffnessMatrix();
	void getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape);
	void setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer,double eta, double zeta, double nu);
	void setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack);
	void setCoeffMat();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, double eta, double zeta, double nu);


public:
	Prism(int* NodeIds,vector<double*>& Nodes, int CurrId);
	~Prism();
	void setElasticProperties(double E,double v);
	void calculateForces(int RKid, double** SystemForces);
};

#endif
