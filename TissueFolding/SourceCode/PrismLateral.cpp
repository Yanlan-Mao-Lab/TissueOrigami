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
	//cout<<"constructed lateral prism"<<endl;
}

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
