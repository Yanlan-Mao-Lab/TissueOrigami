/*
 * ModelInputObject.h
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#ifndef MODELINPUTOBJECT_H_
#define MODELINPUTOBJECT_H_

#include <stdio.h>
#include <iostream>
#include "Simulation.h"

class Simulation;
using namespace std;

class ModelInputObject {
private:
	bool checkFileStatus(ifstream &file, string fileName);
	bool readPysicalProperties(ifstream &file);
	bool readSaveOptions(ifstream &file);
	bool readMeshParameters(ifstream& file);
	bool readPeripodiumParameters(ifstream& file);
	bool readTimeParameters(ifstream &file);
	bool readMeshType2(ifstream& file);
	bool readGrowthOptions(ifstream& file);
	bool readGrowthType1(ifstream& file);
	bool readGrowthType2(ifstream& file);
	bool readShapeChangeOptions(ifstream& file);
	bool readShapeChangeType1(ifstream& file);
public:
	Simulation* Sim;
	const char*  parameterFileName;
	string meshFileName;
	ModelInputObject();
	~ModelInputObject();

	bool readParameters();
};

#endif /* MODELINPUTOBJECT_H_ */
