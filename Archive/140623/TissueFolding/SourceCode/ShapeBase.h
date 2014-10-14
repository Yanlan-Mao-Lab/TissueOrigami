#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "ReferenceShapeBase.h"

using namespace std;

class ShapeBase{
protected:
	int ShapeType;
	int Id;
	int nNodes;
	int* NodeIds;
	int nDim;
	ReferenceShapeBase* ReferenceShape;
	int* IdentifierColour;

	void setShapeType(string TypeName);
	void readNodeIds(int* tmpNodeIds);
	void setPositionMatrix(vector<double*>& Nodes);
	void setReferencePositionMatrix();
	void setIdentificationColour();

	bool InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse, double& det);
	int determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm);
	boost::numeric::ublas::matrix<double> D;
	boost::numeric::ublas::matrix<int> CoeffMat;
	boost::numeric::ublas::matrix<double> k;
	boost::numeric::ublas::matrix<double> B;
	boost::numeric::ublas::vector<double> Forces;
public:
	double** Positions;
	bool FixedLowerZ;
	int getId();
	string getName();
	int getShapeType();
	int getNodeNumber();
	int getDim();
	int* getIdentifierColour();
	double** getReferencePos();
	void displayName();
	void displayPositions();
	void displayIdentifierColour();
	void alignReference();
	void calculateForces(int RKid, double **SystemForces);
	void updatePositions(vector<double*>& Nodes);
	void growShape(float scale);
	void displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
};

#endif
