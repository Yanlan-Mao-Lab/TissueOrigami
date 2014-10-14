#ifndef Simulation_H
#define Simulation_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include "ShapeBase.h"

using namespace std;

class Simulation {
private:
	int currElementId;
	void initiateSinglePrismNodes();
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void initiateElementsByRowAndColumn(int Row, int Column);
public:
	float time;
	vector <double*> Nodes;
	vector <ShapeBase*> Elements;
	double** SystemForces;
	double SystemCentre[3];
	Simulation();
	~Simulation();
	void initiateNodes();
	void initiateElements();
	void initiateSystemForces();
	void calculateSystemCentre();
	void RunOneStep();
	void growSystem(float scale);
};

#endif
