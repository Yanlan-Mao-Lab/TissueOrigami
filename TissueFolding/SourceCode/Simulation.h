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
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"

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
	vector <int> ColumnarCircumferencialNodeList;
	vector <int> PeripodialMembraneCircumferencialNodeList;
	vector <int> ApicalColumnarCircumferencialNodeList;
	int	 nCircumferencialNodes;
	int DVRight,DVLeft;
	double StretchVelocity;
	double BoundingBoxSize[3];
	bool ContinueFromSave;


	bool readModeOfSim(int& i, int argc, char **argv);
	bool readParameters(int& i, int argc, char **argv);
	bool readOutputDirectory(int& i, int argc, char **argv);
	bool readSaveDirectoryToDisplay(int& i, int argc, char **argv);
	bool openFilesToDisplay();
	bool readSystemSummaryFromSave();
	void initiateNodesFromSave();
	bool readNodeDataToContinueFromSave();
	void initiateNodesFromMeshInput();
	void initiateElementsFromSave();
	bool readElementDataToContinueFromSave();
	void initiateElementsFromMeshInput();
	void initiatePrismFromSave();
	bool readShapeData(int i);
	void initiateTriangleFromSave(double height);
	void initiatePrismFromMeshInput();
	void initiateTetrahedraFromMeshInput();
	void initiateTriangleFromMeshInput();
	void initiatePrismFromSaveForUpdate(int k);
	void initiateTriangleFromSaveForUpdate(int k,double height);
	void removeElementFromEndOfList();
	void updateNodeNumberFromSave();
	void updateNodePositionsFromSave();
	void updateElementStatesFromSave();
	void updateForcesFromSave();
	void updateTensionCompressionFromSave();
	void readTensionCompressionToContinueFromSave();
	void updateVelocitiesFromSave();
	bool readFinalSimulationStep();
	void reInitiateSystemForces(int oldSize);
	bool checkInputConsistency();
	void setDefaultParameters();
	bool openFiles();
	bool reOpenOutputFile();
	void initiateSystemForces();
	bool initiateMesh(int MeshType,float zHeight);
	bool initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight);
	bool initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight );
	bool initiateMesh(int MeshType);
	bool addPeripodialMembraneToTissue();
	bool generateColumnarCircumferenceNodeList();
	void sortColumnarCircumferenceNodeList();
	void calculateCentreOfNodes(double* centre);
	void getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength);
	bool isColumnarLayer3D();
	bool CalculateTissueHeight();
	void addPeripodialMembraneNodes(vector <int*> &trianglecornerlist, double height, double d);
	void AddPeripodialMembraneCircumference(double height, int& index_begin, int &index_end);
	void AddHorizontalRowOfPeripodialMembraneNodes(vector <int*> &trianglecornerlist, double d, int &index_begin, int &index_end);
	void AddVerticalRowOfPeripodialMembraneNodes(int& layerCount, int nLayers, vector <int*> &trianglecornerlist, double height, int &index_begin, int &index_end);
	void AddPeripodialMembraneCapToMidAttached(int layerCount, vector <int*> &trianglecornerlist, double height, int index_begin, int index_end);
	void AddPeripodialMembraneCapToApicalAttached(int layerCount, vector <int*> &trianglecornerlist, double height, int index_begin, int index_end);
	void AddPeripodialMembraneCapToHoovering(int layerCount, vector <int*> &trianglecornerlist, double height, int index_begin, int index_end);
	void FillNodeAssociationDueToPeripodialMembrane();
	void assignMassWeightsDueToPeripodialMembrane();
	void addPeripodialMembraneElements(vector <int*> &trianglecornerlist, double height);
	void initiateSinglePrismNodes(float zHeight);
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void fixApicalBasalNodes(vector<int> &NodesToFix);
	//void GenerateColumnarCircumferencialNodeList(vector<int> &NodesToFix, int nLastRow);
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void calculateStiffnessMatrices();
	void assignNodeMasses();
	void assignConnectedElementsAndWeightsToNodes();
	void fixAllD(int i);
	void fixZ(int i);
	void zeroForcesOnNode(int RKId, int i);
	void updateDisplaySaveValuesFromRK();
	void saveStep();
	void writeSimulationSummary();
	void writeMeshFileSummary();
	void writeGrowthRatesSummary();
	void writeSaveFileStepHeader();
	void writeNodes();
	void writeElements();
	void writeSaveFileStepFooter();
	void writeTensionCompression();
	void writeForces();
	void writeVelocities();
	void calculateGrowth();
	void cleanUpGrowthRates();
	void cleanUpShapeChangeRates();
	void calculateGrowthUniform(GrowthFunctionBase* currGF);
	void calculateGrowthRing(GrowthFunctionBase* currGF);
	void calculateGrowthGridBased(GrowthFunctionBase* currGF);
	void calculatePeripodialGrowthGridBased(GrowthFunctionBase* currGF);
	void changeCellShapesInSystem();
	void changeCellShapeRing(int currIndexForParameters);
	void setStretch();
	void addStretchForces(int RKId);
	void setupPipetteExperiment();
	void addPipetteForces(int RKId);
	void LaserAblate(double OriginX, double OriginY, double Radius);
	void fillInNodeNeighbourhood();
	void updateElementVolumesAndTissuePlacements();
	void clearNodeMassLists();
	void clearLaserAblatedSites();
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
	string inputMeshFileName;
	string outputFileString;
	bool saveImages;
	bool saveData;
	string name_saveFile;
	int imageSaveInterval;
	int dataSaveInterval;
	double EApical,EBasal,EMid;
	double poisson;
	double ApicalVisc, BasalVisc;
	int noiseOnPysProp[4];
	int MeshType;
	int Row;
	int Column;
	float SideLength;
	float zHeight;
	bool ApicalNodeFix[2];
	bool BasalNodeFix[2];
	double PeripodialElasticity;
	double PeripodialViscosity;
	double PeripodialThicnessScale;
	double lumenHeight;
	double lumenHeightScale;
	int nGrowthFunctions;
	vector <int> GrowthFunctionTypes;
	vector <float> GrowthParameters;
	vector <double***> GrowthMatrices;

	int nShapeChangeFunctions;
	vector <int> ShapeChangeFunctionTypes;
	vector <float> ShapeChangeParameters;
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	double*** SystemForces;
	double*** PackingForces;
	double SystemCentre[3];
	bool AddPeripodialMembrane;
	int PeripodialMembraneType;
	bool stretcherAttached;
	bool PipetteSuction;
	bool ApicalSuction;
	int StretchInitialStep, StretchEndStep;
	int PipetteInitialStep, PipetteEndStep;
	double pipetteCentre[3];
	double pipetteDepth;
	double pipetteRadius;
	double pipetteRadiusSq;
	double effectLimitsInZ[2];

	double SuctionPressure[3];
	double StretchMin, StretchMax, StretchStrain;
	vector <int*> TrianglesToDraw;
	vector <double*> NodesToDraw;
	double TissueHeight;
	int TissueHeightDiscretisationLayers;
	double boundingBox[2][3];
	vector<GrowthFunctionBase*> GrowthFunctions;

	Simulation();
	~Simulation();
	bool readExecutableInputs(int argc, char **argv);
	bool initiateSystem();
	void calculateSystemCentre();
	void cleanGrowthData();
	void cleanMatrixUpdateData();
	void resetForces();
	void calculateColumnarLayerBoundingBox();
	void runOneStep();
	void calculatePacking(int RKId, double PeriThreshold, double ColThreshold);
	void checkPackingToPipette(bool& packsToPip, double* pos, double* pipF,double mass, int id);
	void getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner, bool& bothperipodial);
	void getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	inline void CapPackingForce(double& Fmag);
	void redistributePeripodialMembraneForces(int RKId);
	void updateNodePositions(int RKId);
	void updateNodePositionsForPeripodialMembraneCircumference(int RKId);
	void realignPositionsForMidAttachedPeripodialMembrane(int RKId);
	void realignPositionsForApicalAttachedPeripodialMembrane(int RKId);
	void realignPositionsForHooveringPeripodialMembrane(int RKId);
	void updateElementPositions(int RKId);
	void updateElementPositionsSingle(int RKId, int i );
	bool initiateSavedSystem();
	void updateOneStepFromSave();
	void alignTissueDVToXPositive();
	void calculateDVDistance();
	void TissueAxisPositionDisplay();
	void CoordinateDisplay();
	void calculatePersonalisedGrowthRates();
	void correctTiltedGrowthForGrowthFunction(GrowthFunctionBase* currGF);
	void calculatePersonalisedUniformOrRingGrowth(ShapeBase* currElement, GrowthFunctionBase* currGF, double* NewGrowth);
	void calculatePersonalisedGridBasedGrowth(ShapeBase* currElement, GrowthFunctionBase* currGF, double* NewGrowth);
	void calculatePersonalisedPeripodialGridBasedGrowth(ShapeBase* currElement, GrowthFunctionBase* currGF, double* NewGrowth);
	void calculateTiltedElementPositionOnBase(ShapeBase* currElement);
	void calculatePersonalisedSingleGrowthRate(ShapeBase* currElement, double DVGrowth, double APGrowth, double ABGrowth, double* NewGrowth);
	void calculateBaseElementsFinalPosition(int Id,double DVGrowth, double APGrowth, double ABGrowth);
	void calculateCurrentElementsFinalPosition(ShapeBase* currElement);

};

#endif
