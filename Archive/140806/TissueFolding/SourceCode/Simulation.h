#ifndef Simulation_H
#define Simulation_H

#include <stdio.h>
#include <iostream>
#include <ostream>
#include <fstream>
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
	ofstream dataSaveFile;

	void setDefaultParameters();
	bool openFiles();
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
	void saveStep();
	void writeSaveFileStepHeader();
	void writeNodes();
	void writeElements();
	void writeSaveFileStepFooter();

public:
	float dt;
	int timestep;
	int SimLength;	//in time steps
	string saveDirectory;
	bool saveImages;
	bool saveData;
	string name_saveFile;
	int imageSaveInterval;
	int dataSaveInterval;
	double E, poisson;
	double ApicalVisc, BasalVisc;
	int noiseOnPysProp[4];

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
