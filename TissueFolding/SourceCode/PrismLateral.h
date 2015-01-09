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

public:
	PrismLateral(int* NodeIds,vector<Node*>& Nodes, int CurrId);
	~PrismLateral();
	void  setElasticProperties(double EApical, double EBasal, double EMid,double v);
	void  calculateBasalNormal(double * normal);
};

#endif
