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
#include "NewtonRaphsonSolver.h"
#include "CellMigration.h"


#include <omp.h>
//test for rici pull

class ModelInputObject;
using namespace std;

class Simulation {
private:
	int currElementId;
	ModelInputObject* ModInp;
	ofstream saveFileMesh;
	ofstream saveFileTensionCompression;
    ofstream saveFileGrowth;
    ofstream saveFileGrowthRate;
	ofstream saveFileForces;
	ofstream saveFileProteins;
	ofstream saveFilePhysicalProp;
	ofstream saveFilePacking;
	ofstream saveFileSimulationSummary;
	ifstream saveFileToDisplayMesh;
	ifstream saveFileToDisplayTenComp;
    ifstream saveFileToDisplayGrowth;
    ifstream saveFileToDisplayGrowthRate;
	ifstream saveFileToDisplayForce;
	ifstream saveFileToDisplayProteins;
	ifstream saveFileToDisplayPhysicalProp;
	ifstream saveFileToDisplayPacking;
	ifstream saveFileToDisplayVel;
	ifstream saveFileToDisplaySimSum;

	bool TensionCompressionSaved;
    bool GrowthSaved;
    bool GrowthRateSaved;
	bool ForcesSaved;
	bool ProteinsSaved;
	bool physicalPropertiesSaved;
	bool PackingSaved;
	int	 nCircumferencialNodes;
	int dorsalTipIndex,ventralTipIndex,anteriorTipIndex,posteriorTipIndex;
	double StretchDistanceStep;
	bool recordForcesOnFixedNodes;
	double boundingBoxSize[3];
	//double columnarBoundingBoxSize[3];
	//double peripodialBoundingBoxSize[3];
	bool ContinueFromSave;
    int growthRotationUpdateFrequency;
    //vector <Node*> symmetricYBoundaryNodes;
    //vector <Node*> symmetricXBoundaryNodes;
    vector <int> AblatedNodes;

    NewtonRaphsonSolver *NRSolver;



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
    void updateGrowthRateFromSave();
    void updateProteinsFromSave();
    void updatePhysicalPropFromSave();
    void updatePackingFromSave();
	void readTensionCompressionToContinueFromSave();
    void readGrowthToContinueFromSave();
    void readGrowthRateToContinueFromSave();
    void readProteinsToContinueFromSave();
    void readPhysicalPropToContinueFromSave();
	bool readFinalSimulationStep();
	void reInitiateSystemForces(int oldSize);
	bool checkInputConsistency();
	void setDefaultParameters();
	bool openFiles();
	bool reOpenOutputFile();
	void initiateSystemForces();
	bool initiateMesh(int MeshType,float zHeight);
	bool initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight);
	//bool initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight );
	bool initiateMesh(int MeshType);
	void readInTissueWeights();
	bool checkIfTissueWeightsRecorded();
	bool generateColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void sortColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void clearCircumferenceDataFromSymmetricityLine();
	//void removeSymmetryBorderFromColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength);
	bool isColumnarLayer3D();
	bool checkIfThereIsPeripodialMembrane();
	void setLinkerCircumference();
	bool calculateTissueHeight();
	void assignInitialZPositions();
	bool addStraightPeripodialMembraneToTissue();
	bool addCurvedPeripodialMembraneToTissue();
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


	void initiateSinglePrismNodes(float zHeight);
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void fixApicalBasalNodes(vector<int> &NodesToFix);
	void checkForNodeFixing();
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void checkForZeroExternalViscosity();
	void calculateStiffnessMatrices();
    void calculateShapeFunctionDerivatives();
    void assignNodeMasses();
    void updateNodeMasses();
    void assignElementalSurfaceAreaIndices();
    void updateNodeViscositySurfaces();
    void updateElementToConnectedNodes(vector <Node*>& Nodes);
	void assignConnectedElementsAndWeightsToNodes();
	void fixAllD(int i, bool fixWithViscosity);
	void fixAllD(Node* currNode, bool fixWithViscosity);
	void fixX(int i, bool fixWithViscosity);
	void fixX(Node* currNode, bool fixWithViscosity);
	void fixY(int i, bool fixWithViscosity);
	void fixY(Node* currNode, bool fixWithViscosity);
	void fixZ(int i, bool fixWithViscosity);
	void fixZ(Node* currNode, bool fixWithViscosity);
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
	void writePacking();
	void writeProteins();
	void writePhysicalProp();
	void calculateMyosinForces();
	void cleanUpMyosinForces();
	void checkForMyosinUpdates();
	void updateEquilibriumMyosinsFromInputSignal(MyosinFunction* currMF);
	void calculateGrowth();
	void calculateShapeChange();
	void cleanUpGrowthRates();
	void assignIfElementsAreInsideEllipseBands();
	void checkForPinningPositionsUpdate();
	void updatePinningPositions();
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
	void laserAblateTissueType(int ablationType);
	void fillInNodeNeighbourhood();
	void fillInElementColumnLists();
	void updateElementVolumesAndTissuePlacements();
	void clearNodeMassLists();
	void clearLaserAblatedSites();
    void manualPerturbationToInitialSetup(bool deform, bool rotate);
    void pokeElement(int elementId, double dx, double dy, double dz);
    void addCurvatureToColumnar(double h);
    void addSoftPeriphery(double* fractions);
    void setupYsymmetricity();
    void setupXsymmetricity();
    void ablateSpcific();
    void setUpECMMimicingElements();
    void setUpActinMimicingElements();
    //void setSymmetricNode(Node* currNode, double yLimPos);


public:

	ofstream outputFile;
	bool displayIsOn;
	bool DisplaySave;
	bool reachedEndOfSaveFile;
	float dt;
	int timestep;
	double currSimTimeSec;
	double SimLength;	//in seconds
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
	double discProperApicalViscosity;
	double discProperBasalViscosity;
	double discProperMidlineViscosity;
	int noiseOnPysProp[4];
	bool zeroExternalViscosity[3]; //The boolean stating if there is zero external viscosity on any of the 3 dimensions
	bool extendExternalViscosityToInnerTissue;
	bool softenedECM;
	double externalViscosityDPApical;
	double externalViscosityDPBasal;
	double externalViscosityPMApical;
	double externalViscosityPMBasal;
	double externalViscosityLZApical;
	double externalViscosityLZBasal;
	int MeshType;
	int Row;
	int Column;
	float SideLength;
	float zHeight;
	bool fixWithExternalViscosity;
	double fixingExternalViscosity[3];
	bool ApicalNodeFix[3];
	bool BasalNodeFix[3];
	bool CircumferentialNodeFix[5][3]; //row 0: apical circumferece x,y,z ; row 1: basal circumference x,y,z; row 2: linker apical circumference x,y,z, row 3: linker basal circumference x,y,z, row 4: all circumference x,y,z
	double PeripodialElasticity;
	double peripodialApicalViscosity;
	double peripodialBasalViscosity;
	double peripodialMidlineViscosity;
	double currViscMidline;
	double PeripodialThicnessScale;
	double PeripodialLateralThicnessScale; //the thickness of the side region linking two layers, as a fraction of tissue thickness
	double lumenHeight;
	double lumenHeightScale;
	//linker zone parameters
	bool BaseLinkerZoneParametersOnPeripodialness;
	double LinkerZoneApicalElasticity;
	double LinkerZoneBasalYoungsModulus;
	double linkerZoneApicalViscosity;
	double linkerZoneBasalViscosity;
	double linkerZoneMidlineViscosity;

	int nGrowthFunctions;
	bool GridGrowthsPinnedOnInitialMesh;
	int nGrowthPinning; ///number of times growth pinning positions will be updated, the relative positions of elements in tissue.
	int* growthPinUpdateTime;
	bool* growthPinUpdateBools;
	int gridGrowthsInterpolationType; //0 = no interpolation, step function, 1 = linear interpolation (default = 1).
	vector <double***> GrowthMatrices;
	vector<GrowthFunctionBase*> GrowthFunctions;

	int nShapeChangeFunctions;
	vector<GrowthFunctionBase*> ShapeChangeFunctions;

	bool thereIsPlasticDeformation;
	bool plasticDeformationAppliedToPeripodial;
	bool plasticDeformationAppliedToColumnar;
	bool volumeConservedInPlasticDeformation;
	double plasticDeformationHalfLife;
	double zRemodellingLowerThreshold;
	double zRemodellingUpperThreshold;

	int nMyosinFunctions;
	vector<MyosinFunction*> myosinFunctions;
	double kMyo;
	double forcePerMyoMolecule;
	bool thereIsMyosinFeedback;
	double MyosinFeedbackCap;
	bool addCurvatureToTissue;
	double tissueCurvatureDepth;
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	int nElements;
	int nNodes;
	double** SystemForces;
	double** PackingForces;
	//double** PackingForcesPreviousStep;
	//double** PackingForcesTwoStepsAgoStep;
	double** FixedNodeForces;
	vector <double> randomForces;
	double randomForceMean;
	double randomForceVar;
	double SystemCentre[3];
	bool needPeripodialforInputConsistency;
	bool thereIsPeripodialMembrane;
	bool AddPeripodialMembrane;
    bool symmetricY;
    bool symmetricX;
    bool addingRandomForces;
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
	int StretchInitialTime, StretchEndTime;
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
	vector <int> pacingNodeSurfaceList0; // this is the id of the base node that is packing
	vector <int> pacingNodeSurfaceList1; // lists 1 to 3 are the edges of the triangle
	vector <int> pacingNodeSurfaceList2;
	vector <int> pacingNodeSurfaceList3;
	vector <int> pacingNodeEdgeList0; // this is the id of the base node that is packing
	vector <int> pacingNodeEdgeList1; // lists 1 to 2 are the nodes of the packing edge
	vector <int> pacingNodeEdgeList2;
	vector <int> pacingNodePointList0; // this is the id of the base node that is packing
	vector <int> pacingNodePointList1; // list 1 is the packing node

	vector <int> initialSignsSurfacex;
	vector <int> initialSignsSurfacey;
	vector <int> initialSignsSurfacez;
	vector <int> initialSignsEdgex;
	vector <int> initialSignsEdgey;
	vector <int> initialSignsEdgez;
	vector <int> initialSignsPointx;
	vector <int> initialSignsPointy;
	vector <int> initialSignsPointz;

	vector <double> initialWeightSurfacex;
	vector <double> initialWeightSurfacey;
	vector <double> initialWeightSurfacez;
	vector <double> initialWeightEdgex;
	vector <double> initialWeightEdgey;
	vector <double> initialWeightEdgez;
	vector <double> initialWeightPointx;
	vector <double> initialWeightPointy;
	vector <double> initialWeightPointz;

	int	nMarkerEllipseRanges;
	vector<double> markerEllipseBandXCentres;
	vector<double> markerEllipseBandR1Ranges;
	vector<double> markerEllipseBandR2Ranges;
    	bool 	thereIsECMSoftening;	
	int numberOfSoftenedEllipseBands;
	vector<int> ECMSofteningEllipseBandIds;
    //double 	ECMSofteningXRange[2];
	double 	softeningBeginTimeInSec;
	double  softeningEndTimeInSec;
	double 	ECMSofteningFraction;
	bool 	softenBasalECM;
	bool 	softenApicalECM;

	bool 	thereIsECMRemodellinbgWithDeforamtionRate;
	bool 	remodelBasalECM;
	bool 	remodelApicalECM;
	double	ECMRemodellingFraction;
	double 	remodellingThresholdFraction; //The displacement should be above the average with this fraction.

	double packingDetectionThreshold;
	double packingThreshold;
	//soft periphery parameters:
	bool 	softPeriphery;
	double 	softDepth;
	double 	softnessFraction;
	bool 	softPeripheryBooleans[4]; //  [applyToApical]  [applyToBasal]  [applyToColumnar]  [applyToPeripodial]

	bool implicitPacking;
	bool thereIsCellMigration;
	CellMigration* cellMigrationTool;
	vector <double> drawingPointsX, drawingPointsY, drawingPointsZ;

	bool thereIsExplicitECM;
	bool thereIsExplicitActin;
	double ECMRenawalHalfLife; //The half life for ECM renewal inside plastic deformation
	Simulation();
	~Simulation();
	void assignTips();
	bool readExecutableInputs(int argc, char **argv);
	bool initiateSystem();
	void calculateSystemCentre();
	void cleanMatrixUpdateData();
	void resetForces(bool resetPacking);
	void calculateApicalSize();
	void calculateBoundingBox();
    void calculateZProjectedAreas();
    void correctzProjectedAreaForMidNodes();
    void clearProjectedAreas();
    void checkForExperimentalSetupsBeforeIteration();
    void checkForExperimentalSetupsWithinIteration();
    void checkForExperimentalSetupsAfterIteration();
    void setLateralElementsRemodellingPlaneRotationMatrices();

    void checkECMSoftening();
    double calculateAverageDisplacement();
    void updateECMVisocityWithDeformationRate();
    bool runOneStep();
    void updatePlasticDeformation();
    void updateStepNR();
    void calculateNumericalJacobian(bool displayMatricesDuringNumericalCalculation);
    void updateElementPositionsinNR(gsl_matrix* uk);
    void updateNodePositionsNR(gsl_matrix* uk);
    void calculateRandomForces();
    void addRandomForces(gsl_matrix* gExt);
    void smallStrainrunOneStep();
    void packToPipetteWall();
    void calculatePacking2(double PeriThreshold, double ColThreshold);
    void addToEdgeList(Node* nodePointer, ShapeBase* elementPointer, vector <int> & edgeNodeData0, vector <int> & edgeNodeData1);
    bool checkIfPointIsOnEdge(int node0, int node1, double x, double y, double z, double& vx, double& vy, double& vz);
   	void detectPacingNodes();
   	void detectPacingCombinations();
   	void cleanUpPacingCombinations();
    void calculatePacking();
    void calculatePackingK(gsl_matrix* K);
    void calculatePackingNumerical(gsl_matrix* K);
    void calculatePackingForcesImplicit3D();
    void calculatePackingForcesExplicit3D();
    void calculatePackingJacobianNumerical3D(gsl_matrix* K);
    void calculatePackingImplicit3DnotWorking();
    void calculatePackingNumerical3DnotWorking(gsl_matrix* K);
    void addValueToMatrix(gsl_matrix* K, int i, int j, double value);
    void addPackingForces(gsl_matrix* gExt);
	void checkPackingToPipette(bool& packsToPip, double* pos, double* pipF,double mass);
	void getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	inline void CapPackingForce(double& Fmag);
	void bringMyosinStimuliUpToDate();
	void redistributePeripodialMembraneForces(int RKId);
	void updateElementPositions();
	void updateElementPositionsSingle(int i );
	bool initiateSavedSystem();
	void updateOneStepFromSave();
	void alignTissueDVToXPositive();
	void alignTissueAPToXYPlane();
	bool checkFlip();
	void flagElementsThatNeedRefinement();
	void refineElements();
	void addNodesForRefinement(ShapeBase* currElement, int* newNodeIdList);
	void addElementsForRefinement(int* elementsIdsOnColumn, int* newNodeIdList);
	void wrapUpAtTheEndOfSimulation();
	void writeRelaxedMeshFromCurrentState();
	void writeMeshRemovingAblatedRegions();
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
