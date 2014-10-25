#ifndef PrismLateral_H
#define PrismLateral_H


#include "Prism.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class PrismLateral : public Prism{

private:
	void calculateNormals();
	void calculateNormalToBasalSide();
	void calculateReferenceNormalToBasalSide();
	void calculateNormalToApicalSide();
	void calculateReferenceNormalToApicalSide();
	void getCurrentAlignmentSides(double*, double*);
	void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	void calculateXVecForTissueCoordAlignment(double* u);
public:
	PrismLateral(int* NodeIds,vector<Node*>& Nodes, int CurrId);
	~PrismLateral();
	void  setElasticProperties(double E,double v);
	void  setViscosity(double ApicalVisc,double BasalVisc, vector <Node*>& Nodes);
};

#endif
