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
    ofstream saveFileGrowth;
	ofstream saveFileVelocities;
	ofstream saveFileForces;
	ofstream saveFileSimulationSummary;
	ifstream saveFileToDisplayMesh;
	ifstream saveFileToDisplayTenComp;
    ifstream saveFileToDisplayGrowth;
	ifstream saveFileToDisplayForce;
	ifstream saveFileToDisplayVel;
	ifstream saveFileToDisplaySimSum;
	bool TensionCompressionSaved;
    bool GrowthSaved;
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
    int growthRotationUpdateFrequency;

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
    void updateGrowthFromSave();
	void readTensionCompressionToContinueFromSave();
    void readGrowthToContinueFromSave();
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
    void AddPeripodialMembraneNewScratch();
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
    void calculateShapeFunctionDerivatives();
    void assignNodeMasses();
    void updateNodeMasses();
    void updateElementToConnectedNodes(vector <Node*>& Nodes);
	void assignConnectedElementsAndWeightsToNodes();
	void fixAllD(int i);
	void fixZ(int i);
	void zeroForcesOnNode(int RKId, int i);
    void processDisplayDataAndSave();
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
    void writeGrowth();
	void writeForces();
	void writeVelocities();
	void calculateGrowth();
	void cleanUpGrowthRates();
    void updateGrowthRotationMatrices();
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
    void addPipetteForces(gsl_matrix *gExt);
	void LaserAblate(double OriginX, double OriginY, double Radius);
	void fillInNodeNeighbourhood();
	void updateElementVolumesAndTissuePlacements();
	void clearNodeMassLists();
	void clearLaserAblatedSites();
    void manualPerturbationToInitialSetup(bool deform, bool rotate);
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
    double pipetteThickness;
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
    void calculateZProjectedAreas();
    void correctzProjectedAreaForMidNodes();
    void clearProjectedAreas();
	void runOneStep();
    void updateStep4RK();
    void updateStepNR();
    void constructUnMatrix(gsl_matrix* un);
    void constructLumpedMassViscosityDtMatrix(gsl_matrix* mviscdt);
    void calculateElasticForcesForNR(gsl_matrix* ge);
    void calculateViscousForcesForNR(gsl_matrix* gv, gsl_matrix* mviscdt, gsl_matrix* uk, gsl_matrix* un);
    void calculateImplicitKElastic(gsl_matrix* K);
    void calculateImplucitKViscous(gsl_matrix* K, gsl_matrix*  mviscdt);
    void calculateImplucitKViscousNumerical(gsl_matrix*  mviscdt, gsl_matrix*  un, gsl_matrix* uk);
    void calculateImplucitKElasticNumerical(gsl_matrix* K,gsl_matrix* geNoPerturbation);
    void solveForDeltaU(gsl_matrix* K, gsl_vector* g, gsl_vector *deltaU);
    void constructiaForPardiso(gsl_matrix* K, int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec);
    void writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a);
    void writeginPardisoFormat(gsl_vector* g, double* b, const int n);
    int  solveWithPardiso(double* a, double*b, int* ia, int* ja, gsl_vector* x ,const int n_variables);
    void EliminateNode0Displacement(gsl_vector* deltaU);
    void updateUkInNR(gsl_matrix* uk, gsl_vector* deltaU);
    void updateElementPositionsinNR(gsl_matrix* uk);
    bool checkConvergenceViaDeltaU(gsl_vector* deltaU);
    bool checkConvergenceViaForce(gsl_vector* gSum);
    void updateNodePositionsNR(gsl_matrix* uk);
    void calcutateFixedK(gsl_matrix* K, gsl_vector* g);
    void smallStrainrunOneStep();
    void packToPipetteWall(gsl_matrix* gExt, gsl_vector* gSum);
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
	void correctTiltedGrowthForGrowthFunction(GrowthFunctionBase* currGF);
	void calculateTiltedElementPositionOnBase(ShapeBase* currElement);
	void calculateBaseElementsFinalPosition(int Id,double DVGrowth, double APGrowth, double ABGrowth);
	void calculateCurrentElementsFinalPosition(ShapeBase* currElement);

};

#endif
