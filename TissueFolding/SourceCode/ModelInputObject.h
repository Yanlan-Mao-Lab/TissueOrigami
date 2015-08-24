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

/*! ModelInputObject class */
class ModelInputObject {
private:
	bool checkFileStatus(ifstream &file, string fileName);
	bool readPysicalProperties(ifstream &file);
	bool readSaveOptions(ifstream &file);
	bool readMeshParameters(ifstream& file);
	bool readPeripodialMembraneParameters(ifstream& file);
	bool readNodeFixingParameters(ifstream& file);
	bool readTimeParameters(ifstream &file);
	bool readMeshType2(ifstream& file);
	bool readMeshType4(ifstream& file);
	bool readGrowthOptions(ifstream& file);
	bool readGrowthType1(ifstream& file);	//uniform growth
	bool readGrowthType2(ifstream& file);	//ring type growth
	bool readGrowthType3(ifstream& file);	//node based growth read from file;
	bool readShapeChangeOptions(ifstream& file);
	bool readShapeChangeType1(ifstream& file);
	bool readStretcherSetup(ifstream& file);
	bool readPipetteSetup(ifstream& file);
public:
    /** An enum type.
     *  The documentation block cannot be put after the enum!
     */
	Simulation* Sim;
	const char*  parameterFileName;
	string meshFileName;
	ModelInputObject();
	~ModelInputObject();

	bool readParameters();
};

#endif /* MODELINPUTOBJECT_H_ */
