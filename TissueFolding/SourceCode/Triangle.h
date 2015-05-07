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
	int normalCrossOrder[2];	///< The ids of the nodes that will be used in calculation of normal pointing towards the apical surface.
	double slabHeight;			///< The slab height of the triangle, which is the amount of tissue in 3D, that the 2D element is be representing.
	void setTissueCoordsRotationsBuffers();  ///< The function to set the buffers for rotation matrices
	void getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape); ///< The function will return the positions of the reference shape of the element
	void setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer,double eta, double nu);	///< The function will set the shape function derivatives for given barycentric coordinates for stiffness matrix calculation
	void setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack);	///< The function will set the shape function derivatives stack from already calculated shape function derivatives. The matrix form is necessary in calculation of the stiffness matrix.
	void setCoeffMat();	///<The function will set the coefficient matrix of the elemetn necessary for stiffness matrix calculation
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double nu);	//This funciton will calculate the stiffness matrix for given barycentric coordinates, necessary for each step of the gaussian integration for calculation of the stiffness matrix
	void calculateReferenceVolume();	///< This function will calculate the the volume of the reference element, using reference shape positions and Triangle#slabHeight
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
	void AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, int RKId,  double ***SystemForces, double ***PackingForces, vector<Node*> &Nodes);
	bool IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode);
	void getApicalNodePos(double* posCorner);
	void getBasalNodePos(double* posCorner);
	void calculateNormalForPacking(int tissueplacement);
	bool IspointInsideTriangle(int tissueplacement, double x, double y,double z);
};

#endif /* TRIANGLE_H_ */
