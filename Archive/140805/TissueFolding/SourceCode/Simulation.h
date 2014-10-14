#ifndef Simulation_H
#define Simulation_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include "ModelInputObject.h"
#include "ShapeBase.h"
#include "Node.h"

class ModelInputObject;
using namespace std;

class Simulation {
private:
	int currElementId;
	ModelInputObject* ModInp;

	void setDefaultParameters();
	void initiateSystemForces();
	void initiateMesh(int MeshType);
	void initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight);
	void initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight );
	void initiateMesh(int MeshType, string inputtype );
	void initiateSinglePrismNodes();
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void CalculateStiffnessMatrices();
	void fixAllD(int i);
	void fixZ(int i);
	void zeroForcesOnNode(int i);

public:
	float dt;
	float time;		//in seconds
	int SimLength;	//in time steps
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	double** SystemForces;
	double SystemCentre[3];
	Simulation();
	~Simulation();
	bool initiateSystem();

	void calculateSystemCentre();
	void RunOneStep();
	void growSystem(float scale);
};

#endif
