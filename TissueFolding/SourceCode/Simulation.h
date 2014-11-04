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
	ofstream saveFileMesh;
	ofstream saveFileTensionCompression;
	ofstream saveFileVelocities;
	ofstream saveFileForces;
	ofstream saveFileSimulationSummary;
	ifstream saveFileToDisplayMesh;
	ifstream saveFileToDisplayTenComp;
	ifstream saveFileToDisplayForce;
	ifstream saveFileToDisplayVel;
	ifstream saveFileToDisplaySimSum;
	bool TensionCompressionSaved;
	bool ForcesSaved;
	bool VelocitiesSaved;
	vector <int> LateralNodeList;
	int	 nLateralNodes;
	int DVRight,DVLeft;

	bool readModeOfSim(int& i, int argc, char **argv);
	bool readParameters(int& i, int argc, char **argv);
	bool readOutputDirectory(int& i, int argc, char **argv);
	bool readSaveDirectoryToDisplay(int& i, int argc, char **argv);
	bool openFilesToDisplay();
	bool readSystemSummaryFromSave();
	void initiateNodesFromSave();
	void initiateElementsFromSave();
	void initiatePrismFromSave();
	void initiateLateralPrismFromSave();
	void initiatePrismFromSaveForUpdate(int k);
	void removeElementFromEndOfList();
	void updateNodeNumberFromSave();
	void updateNodePositionsFromSave();
	void updateElementStatesFromSave();
	void updateForcesFromSave();
	void updateTensionCompressionFromSave();
	void updateVelocitiesFromSave();
	void reInitiateSystemForces(int oldSize);
	bool checkInputConsistency();
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
	void fixApicalBasalNodes(vector<int> &NodesToFix);
	void fixLateralNodes();
	void GenerateLateralNodeList(vector<int> &NodesToFix, int nLastRow);
	void GenerateLateralNodes();
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void calculateStiffnessMatrices();
	void fixAllD(int i);
	void fixZ(int i);
	void zeroForcesOnNode(int RKId, int i);
	void updateDisplaySaveValuesFromRK();
	void saveStep();
	void writeSimulationSummary();
	void writeSaveFileStepHeader();
	void writeNodes();
	void writeElements();
	void writeSaveFileStepFooter();
	void writeTensionCompression();
	void writePeripodiumTensionCompression();
	void writeForces();
	void writeVelocities();
	void calculateGrowth();
	void cleanUpGrowthRates();
	void cleanUpShapeChangeRates();
	void calculateGrowthUniform(int currIndexForParameters);
	void calculateGrowthRing(int currIndexForParameters);
	void changeCellShapesInSystem();
	void changeCellShapeRing(int currIndexForParameters);
	double calculatePeripodiumArea(int RKId);
	double calculatePeripodiumResistance(int RKId);
	double calculatePeripodiumResistanceForce(int RKId);
	void addPeripodiumResistance(int RKId);
public:
	ofstream outputFile;
	bool displayIsOn;
	bool DisplaySave;
	bool reachedEndOfSaveFile;
	float dt;
	int timestep;
	int SimLength;	//in time steps
	string saveDirectory;
	string saveDirectoryToDisplayString;
	bool saveImages;
	bool saveData;
	string name_saveFile;
	int imageSaveInterval;
	int dataSaveInterval;
	double E, poisson;
	double ApicalVisc, BasalVisc;
	int noiseOnPysProp[4];
	int MeshType;
	int Row;
	int Column;
	float SideLength;
	float zHeight;
	bool ApicalNodeFix[2];
	bool BasalNodeFix[2];
	bool LateralNodeFix[2];
	vector <int> PeripodiumAnchorNodeList;
	double ReferencePeripodiumArea;
	double PeripodiumElasticity;
	int nGrowthFunctions;
	vector <int> GrowthFunctionTypes;
	vector <float> GrowthParameters;

	int nShapeChangeFunctions;
	vector <int> ShapeChangeFunctionTypes;
	vector <float> ShapeChangeParameters;
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	double*** SystemForces;
	double SystemCentre[3];
	double PeripodiumStrain;
	double RK1PeripodiumStrain;
	bool AddLateralNodes;
	bool AddPeripodialArea;
	Simulation();
	~Simulation();
	bool readExecutableInputs(int argc, char **argv);
	bool initiateSystem();
	void calculateSystemCentre();
	void cleanGrowthData();
	void cleanMatrixUpdateData();
	void resetForces();
	void runOneStep();
	void updateNodePositions(int RKId);
	void updateElementPositions(int RKId);
	bool initiateSavedSystem();
	void updateOneStepFromSave();
	void alignTissueDVToXPositive();
	void calculateDVDistance();


};

#endif
