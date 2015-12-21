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
#include "MyosinFunction.h"

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
	ofstream saveFileProteins;
	ofstream saveFileSimulationSummary;
	ifstream saveFileToDisplayMesh;
	ifstream saveFileToDisplayTenComp;
    ifstream saveFileToDisplayGrowth;
	ifstream saveFileToDisplayForce;
	ifstream saveFileToDisplayProteins;
	ifstream saveFileToDisplayVel;
	ifstream saveFileToDisplaySimSum;
	bool TensionCompressionSaved;
    bool GrowthSaved;
	bool ForcesSaved;
	bool VelocitiesSaved;
	bool ProteinsSaved;
	int	 nCircumferencialNodes;
	int dorsalTipIndex,ventralTipIndex,anteriorTipIndex,posteriorTipIndex;
	double StretchDistanceStep;
	bool recordForcesOnFixedNodes;
	double boundingBoxSize[3];
	//double columnarBoundingBoxSize[3];
	//double peripodialBoundingBoxSize[3];
	bool ContinueFromSave;
    int growthRotationUpdateFrequency;
    vector <Node*> symmetricYBoundaryNodes;

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
    void updateProteinsFromSave();
	void readTensionCompressionToContinueFromSave();
    void readGrowthToContinueFromSave();
    void readProteinsToContinueFromSave();
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
	bool generateColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void sortColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void clearCircumferenceDataFromSymmetricityLine();
	//void removeSymmetryBorderFromColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength);
	bool isColumnarLayer3D();
	bool calculateTissueHeight();
	bool addPeripodialMembraneToTissue();
	void calculateDiscretisationLayers(double &hColumnar, int& LumenHeightDiscretisationLayers, double &hLumen, double &peripodialHeight, int& peripodialHeightDiscretisationLayers, double& hPeripodial);
	void fillColumnarBasedNodeList(vector< vector<int> > &ColumnarBasedNodeArray, vector <int> &ColumnarCircumferencialNodeList);
	void calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double* pos, double sideThickness);
	void calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double* pos, double sideThickness);
	void addNodesForPeripodialOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray, double hColumnar, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void addNodesForPeripodialOnColumnarCircumference (vector< vector<int> > &ColumnarBasedNodeArray, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void addLateralPeripodialElements(int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray);
	void addNodesForPeripodialOnCap(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &PeripodialCapNodeArray, int TissueHeightDiscretisationLayers, int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void constructTriangleCornerListOfApicalSurface( vector< vector<int> > &TriangleList);
	void addCapPeripodialElements( vector< vector<int> > &TriangleList, vector< vector<int> > &PeripodialCapNodeArray, int peripodialHeightDiscretisationLayers);
	void correctCircumferentialNodeAssignment(vector< vector<int> > OuterNodeArray);
	void assignTips();

	void initiateSinglePrismNodes(float zHeight);
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void fixApicalBasalNodes(vector<int> &NodesToFix);
	void checkForNodeFixing();
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void checkForZeroViscosity();
	void calculateStiffnessMatrices();
    void calculateShapeFunctionDerivatives();
    void assignNodeMasses();
    void updateNodeMasses();
    void updateElementToConnectedNodes(vector <Node*>& Nodes);
	void assignConnectedElementsAndWeightsToNodes();
	void fixAllD(int i);
	void fixAllD(Node* currNode);
	void fixX(int i);
	void fixX(Node* currNode);
	void fixY(int i);
	void fixY(Node* currNode);
	void fixZ(int i);
	void fixZ(Node* currNode);
	void zeroForcesOnNode(int i);
    void processDisplayDataAndSave();
	//void updateDisplaySaveValuesFromRK();
	void saveStep();
	void writeSimulationSummary();
	void writeMeshFileSummary();
	void writeGrowthRatesSummary();
	void writeMyosinSummary();
	void writeSaveFileStepHeader();
	void writeNodes();
	void writeElements();
	void writeSaveFileStepFooter();
	void writeTensionCompression();
    void writeGrowth();
	void writeForces();
	void writeVelocities();
	void writeProteins();
	void calculateMyosinForces();
	void cleanUpMyosinForces();
	void checkForMyosinUpdates();
	void updateEquilibriumMyosinsFromInputSignal(MyosinFunction* currMF);
	void calculateGrowth();
	void calculateShapeChange();
	void cleanUpGrowthRates();
    void updateGrowthRotationMatrices();
	void cleanUpShapeChangeRates();
	void calculateGrowthUniform(GrowthFunctionBase* currGF);
	void calculateGrowthRing(GrowthFunctionBase* currGF);
	void calculateGrowthGridBased(GrowthFunctionBase* currGF);
	//void calculatePeripodialGrowthGridBased(GrowthFunctionBase* currGF);
	void calculateShapeChangeUniform (GrowthFunctionBase* currSCF);
	void changeCellShapesInSystem();
	void changeCellShapeRing(int currIndexForParameters);
	void setStretch();
	void setUpClampBorders(vector<int>& clampedNodeIds);
	void moveClampedNodesForStretcher();
	void recordForcesOnClampBorders();
	void setupPipetteExperiment();
    void addPipetteForces(gsl_matrix *gExt);
    void addMyosinForces(gsl_matrix* gExt);
	void laserAblate(double OriginX, double OriginY, double Radius);
	void fillInNodeNeighbourhood();
	void updateElementVolumesAndTissuePlacements();
	void clearNodeMassLists();
	void clearLaserAblatedSites();
    void manualPerturbationToInitialSetup(bool deform, bool rotate);
    void pokeElement(int elementId, double dx, double dy, double dz);
    void addCurvatureToColumnar(double h);
    void setupYsymmetricity();
    void ablateSpcific();
    //void setSymmetricNode(Node* currNode, double yLimPos);


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
	bool zeroViscosity;
	int noiseOnPysProp[4];
	int MeshType;
	int Row;
	int Column;
	float SideLength;
	float zHeight;
	bool ApicalNodeFix[3];
	bool BasalNodeFix[3];
	bool CircumferentialNodeFix[3][3]; //row 0: apical circumferece x,y,z ; row 1: basal circumference x,y,z; row 2: all circumference x,y,z
	double PeripodialElasticity;
	double PeripodialViscosity;
	double PeripodialThicnessScale;
	double PeripodialLateralThicnessScale; //the thickness of the side region linking two layers, as a fraction of tissue thickness
	double lumenHeight;
	double lumenHeightScale;
	int nGrowthFunctions;

	vector <double***> GrowthMatrices;
	vector<GrowthFunctionBase*> GrowthFunctions;

	int nShapeChangeFunctions;
	vector<GrowthFunctionBase*> ShapeChangeFunctions;

	int nMyosinFunctions;
	vector<MyosinFunction*> myosinFunctions;
	double kMyo;
	double forcePerMyoMolecule;

	bool addCurvatureToTissue;
	double tissueCurvatureDepth;
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	int nElements;
	int nNodes;
	double** SystemForces;
	double** PackingForces;
	double** PackingForcesPreviousStep;
	double** PackingForcesTwoStepsAgoStep;
	double** FixedNodeForces;
	double SystemCentre[3];
	bool AddPeripodialMembrane;
    bool symmetricY;
	bool stretcherAttached;
	vector <int> leftClampBorder;
	vector <int> rightClampBorder;
	double leftClampForces[3];
	double rightClampForces[3];
	bool DVClamp;
	int distanceIndex;	//the index of dimension for stretcher clamp position,
	bool PipetteSuction;
	bool ApicalSuction;
	vector <int> TransientZFixListForPipette;
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
	vector <int> pacingNodeCouples0;
	vector <int> pacingNodeCouples1;
	//double columnarBoundingBox[2][3];
	//double peripodialBoundingBox[2][3];

	Simulation();
	~Simulation();
	bool readExecutableInputs(int argc, char **argv);
	bool initiateSystem();
	void calculateSystemCentre();
	//void cleanGrowthData();
	void cleanMatrixUpdateData();
	void resetForces(bool resetPacking);
	void calculateApicalSize();
	void calculateBoundingBox();
	//void calculateColumnarLayerBoundingBox();
	//void calculatePeripodialBoundingBox();
    void calculateZProjectedAreas();
    void correctzProjectedAreaForMidNodes();
    void clearProjectedAreas();
    void checkForExperimentalSetupsBeforeIteration();
    void checkForExperimentalSetupsWithinIteration();
    void checkForExperimentalSetupsAfterIteration();
	//void bakeTissue();
    void runOneStep();
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
    void updateUkInNR(gsl_matrix* uk, gsl_vector* deltaU);
    void updateElementPositionsinNR(gsl_matrix* uk);
    bool checkConvergenceViaDeltaU(gsl_vector* deltaU);
    bool checkConvergenceViaForce(gsl_vector* gSum);
    void updateNodePositionsNR(gsl_matrix* uk);
    void calcutateFixedK(gsl_matrix* K, gsl_vector* g);
    void smallStrainrunOneStep();
    void packToPipetteWall();
    void calculatePacking2(double PeriThreshold, double ColThreshold);
    void addToEdgeList(Node* nodePointer, ShapeBase* elementPointer, vector <int> & edgeNodeData0, vector <int> & edgeNodeData1);
    bool checkIfPointIsOnEdge(int node0, int node1, double x, double y, double z, double& vx, double& vy, double& vz);
   	void detectPacingNodes();
    void calculatePacking();
    void calculatePackingImplicit();
    void calculatePackingK(gsl_matrix* K);
    void calculatePackingNumerical(gsl_matrix* K);
    void addPackingForces(gsl_matrix* gExt);
	void checkPackingToPipette(bool& packsToPip, double* pos, double* pipF,double mass, int id);
	void getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	inline void CapPackingForce(double& Fmag);
	void redistributePeripodialMembraneForces(int RKId);
	void updateElementPositions();
	void updateElementPositionsSingle(int i );
	bool initiateSavedSystem();
	void updateOneStepFromSave();
	void alignTissueDVToXPositive();
	void alignTissueAPToXYPlane();
	void wrapUpAtTheEndOfSimulation();
	void writeRelaxedMeshFromCurrentState();
	void calculateDVDistance();
	void TissueAxisPositionDisplay();
	void coordinateDisplay();
	void correctTiltedGrowthForGrowthFunction(GrowthFunctionBase* currGF);
	void calculateTiltedElementPositionOnBase(ShapeBase* currElement);
	void calculateBaseElementsFinalPosition(int Id,double DVGrowth, double APGrowth, double ABGrowth);
	void calculateCurrentElementsFinalPosition(ShapeBase* currElement);
	void fixNode0InPosition(double x, double y, double z);


};

#endif
