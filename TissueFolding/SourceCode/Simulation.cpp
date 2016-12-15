
#include "Simulation.h"
#include "Prism.h"
#include "RandomGenerator.h"
#include <string.h>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


using namespace std;

Simulation::Simulation(){
	currElementId = 0;
	ModInp = new ModelInputObject();
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	TissueHeight = 0.0;
	TissueHeightDiscretisationLayers= 1;
	timestep = 0;
	currSimTimeSec = 0;
	reachedEndOfSaveFile = false;
	AddPeripodialMembrane = false;
	thereIsPeripodialMembrane = false;
	needPeripodialforInputConsistency = false;
	lumenHeight = -20;
	boundingBoxSize[0]=1000.0; boundingBoxSize[1]=1000.0; boundingBoxSize[2]=1000.0;
	//columnarBoundingBoxSize[0]=1000.0; columnarBoundingBoxSize[1]=1000.0; columnarBoundingBoxSize[2]=1000.0;
	//peripodialBoundingBoxSize[0]=1000.0; peripodialBoundingBoxSize[1]=1000.0; peripodialBoundingBoxSize[2]=1000.0;
	ContinueFromSave = false;
    growthRotationUpdateFrequency = 60.0/dt;
    nElements = 0;
    nNodes = 0;
    if (growthRotationUpdateFrequency<1) {growthRotationUpdateFrequency =1;}
	setDefaultParameters();
	implicitPacking = true;
}

Simulation::~Simulation(){
    cerr<<"destructor for simulation called"<<endl;
	delete ModInp;
	if (nGrowthPinning>0){
		delete [] growthPinUpdateTime;
		delete [] growthPinUpdateBools;
	}
	//n nodes
	for (int j=0;j<nNodes;++j){
		delete[] SystemForces[j];
		delete[] PackingForces[j];
		//delete[] PackingForcesPreviousStep[j];
		//delete[] PackingForcesTwoStepsAgoStep[j];
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
    delete[] PackingForces;
    //delete[] PackingForcesPreviousStep;
    //delete[] PackingForcesTwoStepsAgoStep;
    delete[] FixedNodeForces;
    cout<<"deleting elements"<<endl;
	while(!Elements.empty()){
		ShapeBase* tmp_pt;
		tmp_pt = Elements.back();
		Elements.pop_back();
		delete tmp_pt;
        //cerr<<"Element list size: "<<Elements.size()<<endl;
	}
    cout<<"deleting nodes"<<endl;
	while(!Nodes.empty()){
		Node* tmp_pt;
		tmp_pt = Nodes.back();
		Nodes.pop_back();
		delete tmp_pt;
	}
    cout<<"deleting GrowthFunctions"<<endl;
	while(!GrowthFunctions.empty()){
		GrowthFunctionBase* tmp_GF;
		tmp_GF = GrowthFunctions.back();
		GrowthFunctions.pop_back();
		delete tmp_GF;
        //cerr<<"GrowtFuncitons list size: "<<GrowthFunctions.size()<<endl;
	}
    cout<<"deleting NRSolver"<<endl;
    delete NRSolver;
    cout<<"deletion complete"<<endl;

}

void Simulation::setDefaultParameters(){
	dt = 0.01;				//sec
	SimLength = 10.0; 		//10 sec of simulation
	saveImages = false;		//do not save simulation images by default
	saveData = false;		//do not save simulation data by default
	imageSaveInterval = 60.0/dt;	//save images every minute
    dataSaveInterval  = 60.0/dt;	//save data every minute
	saveDirectory = "Not-Set";	//the directory to save the images and data points
	saveDirectoryToDisplayString  = "Not-Set"; //the file whcih will be read and displayed - no simulation
	EApical = 10.0;
	EBasal = 10.0;
	EMid = 10.0;
	PeripodialElasticity = 10.0;
	poisson = 0.3;
	discProperApicalViscosity = 0.0;
	discProperBasalViscosity = 0.0;
	for (int i=0; i<3; ++i){
		zeroExternalViscosity[i] = true;
	}
	memset(noiseOnPysProp,0,4*sizeof(int));
	// The default input is a calculated mesh of width 4 elements, each element being 2.0 unit high
	// and having 1.0 unit sides of the triangles.
	MeshType = 2;
	Row = 4;
	Column = Row-2;
	SideLength=1.0;
	zHeight = 2.0;
	fixWithExternalViscosity = false;
	for (int i=0; i<5; ++i){
		for (int j=0; j<3; j++){
			CircumferentialNodeFix[i][j] = false;
			ApicalNodeFix[j] = false;
			BasalNodeFix[j] = false;
			fixingExternalViscosity[j] = 0;
		}
	}
	nGrowthFunctions = 0;
	GridGrowthsPinnedOnInitialMesh = false;
	nGrowthPinning = 0;
	gridGrowthsInterpolationType = 1;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
    GrowthSaved = true;
    GrowthRateSaved = true;
	ForcesSaved = true;
	ProteinsSaved = true;
	PackingSaved = true;
	physicalPropertiesSaved = true;
	PeripodialElasticity = 0.0;
	peripodialApicalViscosity = discProperApicalViscosity;
	peripodialBasalViscosity  = discProperBasalViscosity;
	PeripodialThicnessScale = 1.0;
	PeripodialLateralThicnessScale = 0.3;
	lumenHeightScale = 0.3;
	dorsalTipIndex = 0;
	ventralTipIndex = 1;
	anteriorTipIndex = 0;
	posteriorTipIndex = 0;
	stretcherAttached = false;
	recordForcesOnFixedNodes = false;
	distanceIndex = false;
	StretchDistanceStep = 0.0;
	DVClamp = true;
	StretchInitialTime = -100;
	StretchEndTime = -100;
	PipetteSuction = false;
	PipetteInitialStep= -100;
	PipetteEndStep = 0.0;
	ApicalSuction = true;
	pipetteCentre[0] = 0.0;
	pipetteCentre[1] = 0.0;
	pipetteCentre[2] = 0.0;
	pipetteDepth = 0.0;
	pipetteRadius =0.0;
	SuctionPressure[0] = 0.0;
	SuctionPressure[1] = 0.0;
	SuctionPressure[2] = 0.0;
	pipetteRadiusSq = pipetteRadius*pipetteRadius;
    pipetteThickness = 11.7;
	effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
	effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;

	boundingBox[0][0] =  1000.0;	//left x
	boundingBox[0][1] =  1000.0;	//low y
	boundingBox[0][2] =  1000.0;	//bottom z
	boundingBox[1][0] = -1000.0;	//right x
	boundingBox[1][1] = -1000.0;	//high y
	boundingBox[1][2] = -1000.0;	//top z

	addCurvatureToTissue = false;
	tissueCurvatureDepth = 0.0;
	symmetricY = false;
	symmetricX = false;
	addingRandomForces = false;
	randomForceMean = 0.0;
	randomForceVar = 0.0;
	packingDetectionThreshold = 5.5;
	packingThreshold = 6.0;

	softPeriphery = false;
	softDepth = 0.0;
	softnessFraction = 0.0;
	softPeripheryBooleans[0] = false; //apical
	softPeripheryBooleans[1] = false; //basal
	softPeripheryBooleans[2] = false; //columnar
	softPeripheryBooleans[3] = false; //peripodial

	thereIsPlasticDeformation = false;
	plasticDeformationAppliedToPeripodial = false;
	plasticDeformationAppliedToColumnar = false;
	volumeConservedInPlasticDeformation = false;
	plasticDeformationHalfLife = 0.0;
	zRemodellingLowerThreshold = 0.5;
	zRemodellingUpperThreshold = 2.0;

	kMyo = 0.09873;
	forcePerMyoMolecule = 1.0;
	thereIsMyosinFeedback = false;
	MyosinFeedbackCap = 0.0;

	BaseLinkerZoneParametersOnPeripodialness = true;
	LinkerZoneApicalElasticity = 0.0;
	LinkerZoneBasalYoungsModulus = 0.0;
	linkerZoneApicalViscosity = discProperApicalViscosity;
	linkerZoneBasalViscosity = discProperBasalViscosity;

	extendExternalViscosityToInnerTissue = false;
	externalViscosityDPApical = 0.0;
	externalViscosityDPBasal  = 0.0;
	externalViscosityPMApical = 0.0;
	externalViscosityPMBasal  = 0.0;
	externalViscosityLZApical = 0.0;
	externalViscosityLZBasal  = 0.0;

	softenedECM = false;
    thereIsECMSoftening = false;
    numberOfSoftenedRanges = 0;
	//ECMSofteningXRange[0] = 0.0;
    //ECMSofteningXRange[1] = 1.0;
	softeningEndTimeInSec = -1;
	softeningEndTimeInSec = -1;
	ECMSofteningFraction = 0.0;
	softenBasalECM = false;
	softenApicalECM = false;
	thereIsECMRemodellinbgWithDeforamtionRate = false;
	remodellingThresholdFraction = 1.10;
	remodelBasalECM = false;
	remodelApicalECM = false;
	ECMRemodellingFraction = 0.0;

	thereIsCellMigration = false;
	thereIsExplicitECM = false;
	ECMRenawalHalfLife = 0.0;
	thereIsExplicitActin = false;
}

bool Simulation::readExecutableInputs(int argc, char **argv){
	int i = 1;
	bool Success = true;
	while(i<argc){
		const char *inptype = argv[i];
		if (string(inptype) == "-mode"){
			Success = readModeOfSim(i, argc, argv);
		}
		else if (string(inptype) == "-i"){
			Success = readParameters(i, argc, argv);
		}
		else if (string(inptype) == "-od"){
			Success = readOutputDirectory(i, argc, argv);
		}
		else if (string(inptype) == "-dInput"){
			Success = readSaveDirectoryToDisplay(i, argc, argv);
		}
		else {
			cerr<<"Please enter a valid option key: {-mode,-i, -od, -dInput}, current string: "<<inptype<<endl;
			return false;
		}
		i++;
		if (!Success){
			return Success;
		}
	}
	Success = checkInputConsistency();
	return Success;
}

bool Simulation::readModeOfSim(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the mode of simulation: {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	if (string(inpstring) == "DisplaySave"){
		DisplaySave = true;
		return true;
	}
	else if (string(inpstring) == "SimulationOnTheGo" || string(inpstring) == "Default"){
		DisplaySave = false;
		return true;
	}
	else if (string(inpstring) == "ContinueFromSave"){
		ContinueFromSave = true;
		DisplaySave = false;
		//DisplaySave = true;
		return true;
	}
	else{
		cerr<<"Please provide input mode: -mode {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}";
		return false;
	}
}

bool Simulation::readParameters(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the model input file"<<endl;
		return false;
	}
	//cerr<<"Reading parameter input file"<<endl;
	const char* inpstring = argv[i];
	ModInp->Sim=this;
	ModInp->parameterFileName =  inpstring;
	//cerr<<" Reading parameters from file: "<<ModInp->parameterFileName<<endl;
	bool Success = ModInp->readParameters();
	if (!Success){
		return Success;
	}
	return true;
}

bool Simulation::readSaveDirectoryToDisplay(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the save directory, contents of which will be displayed"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	saveDirectoryToDisplayString = string(inpstring);
	return true;
}

bool Simulation::readOutputDirectory(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the save directory"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	//This will set the save directory, but will not change the safe file boolean.
	//If your model input file states no saving, then the error and output files
	//will be directed into this directory, but the frame saving will not be toggled
	saveDirectory= string(inpstring);
	return true;

}

bool Simulation::readFinalSimulationStep(){
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
	//backing up data save interval and time step before reading the system summary:
	double timeStepCurrentSim = dt;
	int dataSaveIntervalCurrentSim = dataSaveInterval;
	//bool ZerothFrame = true;
	//reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	while(reachedEndOfSaveFile == false){
		//skipping the header:
		getline(saveFileToDisplayMesh,currline);
		if(saveFileToDisplayMesh.eof()){
			reachedEndOfSaveFile = true;
			break;
		}
		success = readNodeDataToContinueFromSave();
		cout<<"dt after  readNodeDataToContinueFromSave "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<" dataSaveIntervalCurrentSim: "<<dataSaveIntervalCurrentSim<<endl;
		if (!success){
			return false;
		}
		readElementDataToContinueFromSave();
		cout<<"dt after  readElementDataToContinueFromSave "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<" dataSaveIntervalCurrentSim: "<<dataSaveIntervalCurrentSim<<endl;
		//assignPhysicalParameters();

		if (TensionCompressionSaved){
			readTensionCompressionToContinueFromSave();
		}
        if (GrowthSaved){
            readGrowthToContinueFromSave();
        }
        if (GrowthRateSaved){
            readGrowthRateToContinueFromSave();
        }
        if (ProteinsSaved){
        	readProteinsToContinueFromSave();
		}
        if (physicalPropertiesSaved){
        	readPhysicalPropToContinueFromSave();
        }
		//if (ZerothFrame){
		//	ZerothFrame = false;
		//	cout<<"dt after  ZerothFrame if clause "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<endl;
		//}
		//else{
			timestep = timestep + dataSaveInterval;
			currSimTimeSec += dt*dataSaveInterval;
		//}
		cout<<"current time step "<<timestep<<" currSimTimeSec: "<<currSimTimeSec<<endl;

		//skipping the footer:
		getline(saveFileToDisplayMesh,currline);
		while (currline.empty() && !saveFileToDisplayMesh.eof()){
			//skipping empty line
			getline(saveFileToDisplayMesh,currline);
		}
	}
	updateElementVolumesAndTissuePlacements();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	calculateBoundingBox();
	bringMyosinStimuliUpToDate();

	//During a simulation, the data is saved after the step is run, and the time step is incremented after the save. Now I have read in the final step,
	//I need to increment my time step to continue the next time step from here.
	//currSimTimeSec += dt*dataSaveInterval;
	//timestep++;

	//bring the time step and data save time steps to the main modelinput:
	dataSaveInterval = dataSaveIntervalCurrentSim;
	dt = timeStepCurrentSim;

	return true;
}

bool Simulation::checkInputConsistency(){
	if (ContinueFromSave &&  saveDirectoryToDisplayString == "Not-Set"){
		cerr <<"The mode is set to continue from saved simulation, please provide an input directory containing the saved profile, using -dInput"<<endl;
		return false;
	}
	if (saveData || saveImages){
		if (saveDirectory == "Not-Set"){
			cerr <<"Modelinput file requires saving, please provide output directory, using -od tag"<<endl;
			return false;
		}
	}
	if (DisplaySave){
		if (saveDirectoryToDisplayString == "Not-Set"){
			cerr <<"The mode is set to display from save, please provide an input directory, using -dInput"<<endl;
			return false;
		}
	}
	if (AddPeripodialMembrane == false){
		for (int i=0; i<nGrowthFunctions; ++i){
			if(GrowthFunctions[i]->applyToPeripodialMembrane){
				cerr<<"There is no peripodial membrane, while growth function "<<i<<" is applicable to peropodial membrane, further checks needed"<<endl;
				needPeripodialforInputConsistency = true;
			}
		}
	}
	return true;
}

bool Simulation::initiateSystem(){
	bool Success = openFiles();
	if (!Success){
		return Success;
	}
	if (MeshType == 1){
		Success = initiateMesh(MeshType, zHeight); //zHeight
	}
	else if (MeshType == 2){
		Success = initiateMesh(MeshType, Row, Column,  SideLength,  zHeight);
	}
	else if(MeshType == 4){
		cout<<"mesh type : 4"<<endl;
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	if (symmetricY || symmetricX){
		 clearCircumferenceDataFromSymmetricityLine();
	}
	Success = checkIfThereIsPeripodialMembrane();
	Success = calculateTissueHeight(); //Calculating how many layers the columnar layer has, and what the actual height is.
	calculateBoundingBox();
	assignInitialZPositions();
	if (!Success){
		return Success;
	}
	if (AddPeripodialMembrane){
		if (thereIsPeripodialMembrane){
			Success = false;
			cerr<<"Error-there is already peripodial membrane added to the mesh, but modelinput file asks to add another"<<endl;
		}
		else{
			//Success = addCurvedPeripodialMembraneToTissue();
			Success = addStraightPeripodialMembraneToTissue();
			if (Success){
				thereIsPeripodialMembrane = true;
			}
		}
	}
	if (needPeripodialforInputConsistency){
		if (!thereIsPeripodialMembrane){
			cerr<<"There is no peripodial membrane but at least one growth function desires one"<<endl;
			Success = false;
		}
	}
	if (addCurvatureToTissue){
		addCurvatureToColumnar(tissueCurvatureDepth);
	}

	if (thereIsExplicitECM){
		setUpECMMimicingElements();
		//TO DO: NEED TO UPDATE THE PHYSICAL PROPERTIES AFTER THIS!
	}
	if (thereIsExplicitActin){
		setUpActinMimicingElements();
	}
	fillInNodeNeighbourhood();
	fillInElementColumnLists();
	checkForNodeFixing();
	assignTips();
	if (!Success){
		return Success;
	}
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
	checkForZeroExternalViscosity();
    //calculateStiffnessMatrices();
	calculateShapeFunctionDerivatives();
	assignNodeMasses();
	assignElementalSurfaceAreaIndices();
	assignConnectedElementsAndWeightsToNodes();
    alignTissueDVToXPositive();
    //alignTissueAPToXYPlane();
    calculateBoundingBox();
    calculateDVDistance();
	vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
    	(*itElement)->setInitialRelativePosInBoundingBox();
    }
	if (stretcherAttached){
		setStretch();
	}
    cout<<"setting the pipette"<<endl;
	if(PipetteSuction){
		setupPipetteExperiment();
	}

	if (symmetricY){
		setupYsymmetricity();
	}
	if (symmetricX){
		setupXsymmetricity();
	}
	if (ContinueFromSave){
		cout<<"Reading Final SimulationStep: "<<endl;
		Success = readFinalSimulationStep();
		if (!Success){
			return Success;
		}
	}
	if (saveData){
		cout<<"writing the summary current simulation parameters"<<endl;
		writeSimulationSummary();
	}
	nNodes = Nodes.size();
	nElements = Elements.size();
	//initiating the NR solver object, that generates the necessary matrices
	//for solving the tissue dynamics
    NRSolver = new NewtonRaphsonSolver(Nodes[0]->nDim,nNodes);

    if (thereIsCellMigration) {
    	cout<<"initiation of cell migration"<<endl;
    	double cellMigrationOriginAngle = M_PI/2.0;
    	cellMigrationTool = new CellMigration(nElements,0.5); //50% leaves the tissue region per hour
        cellMigrationTool->assignElementRadialVectors(Elements);
        cellMigrationTool->assignOriginOfMigration(Elements, cellMigrationOriginAngle);
        cellMigrationTool->assignElementConnectivity(Nodes,Elements);
        cellMigrationTool->generateListOfRateFractions(Elements);
    }
    if (thereIsPlasticDeformation){
    	setLateralElementsRemodellingPlaneRotationMatrices();
    }

	return Success;
}

void Simulation::setLateralElementsRemodellingPlaneRotationMatrices(){
	calculateSystemCentre();
	double cx = SystemCentre[0];
	double cy = SystemCentre[1];
	if (symmetricX){
		cx = 0;
	}
	if (symmetricY){
		cy = 0;
	}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->tissueType == 2){
			(*itElement)->setLateralElementsRemodellingPlaneRotationMatrix(cx,cy);
		}
	}
}

void Simulation::assignNodeMasses(){
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    (*itElement)->assignVolumesToNodes(Nodes);
		//(*itElement)->assignSurfaceAreaToNodes(Nodes);
	}
}

void Simulation::assignElementalSurfaceAreaIndices(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->assignExposedSurfaceAreaIndices(Nodes);
	}
}

void Simulation::assignConnectedElementsAndWeightsToNodes(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->assignElementToConnectedNodes(Nodes);
	}
}

bool Simulation::openFiles(){
	bool Success;
	if (saveData){
		//Mesh information at each step:
		string saveFileString = saveDirectory +"/Save_Frame";
		const char* name_saveFileMesh = saveFileString.c_str();
		saveFileMesh.open(name_saveFileMesh, ofstream::out);
		if (saveFileMesh.good() && saveFileMesh.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileMesh<<endl;
			Success = false;
		}

		saveFileString = saveDirectory +"Save_Summary";
		const char* name_saveFileSimulationSummary = saveFileString.c_str();
		saveFileSimulationSummary.open(name_saveFileSimulationSummary, ofstream::out);
		if (saveFileSimulationSummary.good() && saveFileSimulationSummary.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileSimulationSummary<<endl;
			Success = false;
		}

		//tension compression information at each step:
		saveFileString = saveDirectory +"/Save_TensionCompression";
		const char* name_saveFileTenComp = saveFileString.c_str();
		cout<<"opening the file" <<name_saveFileTenComp<<endl;
		saveFileTensionCompression.open(name_saveFileTenComp, ofstream::binary);
		if (saveFileTensionCompression.good() && saveFileTensionCompression.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileTenComp<<endl;
			Success = false;
		}

        //Growth information at each step (Fg):
        saveFileString = saveDirectory +"/Save_Growth";
        const char* name_saveFileGrowth = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowth<<endl;
        saveFileGrowth.open(name_saveFileGrowth, ofstream::binary);
        if (saveFileGrowth.good() && saveFileGrowth.is_open()){
            Success = true;
        }
        else{
            cerr<<"could not open file: "<<name_saveFileGrowth<<endl;
            Success = false;
        }

        //Growth rate information at each step (rx,ry,rz, (as in exp(rx*dt) )for display purposes only):
        saveFileString = saveDirectory +"/Save_GrowthRate";
        const char* name_saveFileGrowthRate = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowthRate<<endl;
        saveFileGrowthRate.open(name_saveFileGrowthRate, ofstream::binary);
        if (saveFileGrowthRate.good() && saveFileGrowthRate.is_open()){
            Success = true;
        }
        else{
            cerr<<"could not open file: "<<name_saveFileGrowthRate<<endl;
            Success = false;
        }

		//Force information at each step:
		saveFileString = saveDirectory +"/Save_Force";
		const char* name_saveFileForces = saveFileString.c_str();
		saveFileForces.open(name_saveFileForces, ofstream::binary);
		if (saveFileForces.good() && saveFileForces.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileForces<<endl;
			Success = false;
		}
		//opening the packing information file:
		saveFileString = saveDirectory +"/Save_Packing";
		const char* name_saveFilePacking = saveFileString.c_str();
		saveFilePacking.open(name_saveFilePacking, ofstream::binary);
		if (saveFilePacking.good() && saveFilePacking.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFilePacking<<endl;
			Success = false;
		}
		//opening the protein information file:
		saveFileString = saveDirectory +"/Save_Proteins";
		const char* name_saveFileProteins = saveFileString.c_str();
		saveFileProteins.open(name_saveFileProteins, ofstream::binary);
		if (saveFileProteins.good() && saveFileProteins.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileProteins<<endl;
			Success = false;
		}
		//opening the physical property information file:
		saveFileString = saveDirectory +"/Save_PhysicalProp";
		const char* name_saveFilePhysicalProp = saveFileString.c_str();
		saveFilePhysicalProp.open(name_saveFilePhysicalProp, ofstream::binary);
		if (saveFilePhysicalProp.good() && saveFilePhysicalProp.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFilePhysicalProp<<endl;
			Success = false;
		}
	}
	if (saveDirectory == "Not-Set"){
		cerr<<"Output directory is not set, outputting on Out file in current directory"<<endl;
		outputFileString = "./Out";
	}
	else {
		outputFileString = saveDirectory +"/Out";
	}
	const char* name_outputFile = outputFileString.c_str();
	outputFile.open(name_outputFile, ofstream::out);
	if (outputFile.good() && outputFile.is_open()){
		Success = true;
	}
	else{
		cerr<<"could not open file: "<<name_outputFile<<endl;
		Success = false;
	}
	return Success;
}

bool Simulation::reOpenOutputFile(){
	outputFile.close();
	const char* name_outputFile = outputFileString.c_str();
	outputFile.open(name_outputFile, ofstream::out);
	if (!(outputFile.good() && outputFile.is_open())){
		cerr<<"at step: "<<currSimTimeSec<<" could not open file: "<<name_outputFile<<endl;
		return false;
	}
	return true;
}

void Simulation::writeSimulationSummary(){
	saveFileSimulationSummary<<"TimeStep(sec):  ";
	saveFileSimulationSummary<<dt<<endl;
	saveFileSimulationSummary<<"DataSaveInterval(sec):  ";
	saveFileSimulationSummary<<dataSaveInterval*dt<<endl;
	saveFileSimulationSummary<<"ModelinputName:  ";
	saveFileSimulationSummary<<ModInp->parameterFileName<<endl;
	saveFileSimulationSummary<<"Mesh_Type:  ";
	saveFileSimulationSummary<<MeshType<<endl;
	saveFileSimulationSummary<<"Symmetricity-x: "<<symmetricX<<" Symmetricity-y: "<<symmetricY<<endl;
	writeMeshFileSummary();
	writeGrowthRatesSummary();
	writeMyosinSummary();
}

void Simulation::writeMeshFileSummary(){
	if ( MeshType == 2){
		saveFileSimulationSummary<<"	Row: ";
		saveFileSimulationSummary<<Row;
		saveFileSimulationSummary<<"	Column: ";
		saveFileSimulationSummary<<Column;
		saveFileSimulationSummary<<"	SideLength: ";
		saveFileSimulationSummary<<SideLength;
		saveFileSimulationSummary<<"	zHeight: ";
		saveFileSimulationSummary<<zHeight;
		saveFileSimulationSummary<<endl;
	}
	if ( MeshType == 4){
		saveFileSimulationSummary<<"	MeshFileName: ";
		saveFileSimulationSummary<<inputMeshFileName;
		saveFileSimulationSummary<<endl;
	}
}

void Simulation::writeGrowthRatesSummary(){
	for (int i=0; i<nGrowthFunctions; ++i){
		if(GrowthFunctions[i]->Type == 3){ //grid based function does not need dt in summary
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary);
		}
		else{
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary,dt);
		}
	}
}

void Simulation::writeMyosinSummary(){
	for (int i=0; i<nMyosinFunctions; ++i){
		myosinFunctions[i]->writeSummary(saveFileSimulationSummary);
	}
}

void Simulation::writeRelaxedMeshFromCurrentState(){
	string meshSaveString = saveDirectory +"/MeshFromEndPoint.mesh";
	const char* name_meshSaveString = meshSaveString.c_str();;
	ofstream file;
	file.open(name_meshSaveString, ofstream::out);
	file<<nNodes;
	file<<endl;
	for (int i=0; i<nNodes; ++i){
		file << Nodes[i]->Position[0];
		file<<" \t";
		file << Nodes[i]->Position[1];
		file<<" \t";
		file << Nodes[i]->Position[2];
		file<<" \t";
		file << Nodes[i]->tissuePlacement;
		file<<" \t";
		file << Nodes[i]->tissueType;
		file<<" \t";
		file << Nodes[i]->atCircumference;
		file << endl;
	}
	file<<nElements;
	file<<endl;
	for (int i=0; i<nElements; ++i){
		file<< Elements[i]->getShapeType();
		file<<" \t";
		const int n = Elements[i]->getNodeNumber();
		int* NodeIds;
		NodeIds = new int[n];
		NodeIds = Elements[i]->getNodeIds();
		for (int j=0; j<n; ++j){
			file<< NodeIds[j];
			file<<" \t";
		}
		for (int j=0; j<n; ++j){
			file << Nodes[NodeIds[j]]->Position[0];
			file<<" \t";
			file << Nodes[NodeIds[j]]->Position[1];
			file<<" \t";
			file << Nodes[NodeIds[j]]->Position[2];
			file<<" \t";
		}
		file << endl;
	}
	file.close();
}

bool Simulation::openFilesToDisplay(){
	string saveFileString = saveDirectoryToDisplayString +"/Save_Frame";
	const char* name_saveFileToDisplayMesh = saveFileString.c_str();;
	saveFileToDisplayMesh.open(name_saveFileToDisplayMesh, ifstream::in);
	if (!(saveFileToDisplayMesh.good() && saveFileToDisplayMesh.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayMesh<<endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Summary";
	const char* name_saveFileToDisplaySimSum = saveFileString.c_str();;
	saveFileToDisplaySimSum.open(name_saveFileToDisplaySimSum, ifstream::in);
	if (!(saveFileToDisplaySimSum.good() && saveFileToDisplaySimSum.is_open())){
		cerr<<"Cannot open the simulation summary: "<<name_saveFileToDisplaySimSum<<endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_TensionCompression";
	const char* name_saveFileToDisplayTenComp = saveFileString.c_str();
	saveFileToDisplayTenComp.open(name_saveFileToDisplayTenComp, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayTenComp<<endl;
		TensionCompressionSaved = false;
	}

    saveFileString = saveDirectoryToDisplayString +"/Save_Growth";
    const char* name_saveFileToDisplayGrowth = saveFileString.c_str();
    saveFileToDisplayGrowth.open(name_saveFileToDisplayGrowth, ifstream::in);
    if (!(saveFileToDisplayGrowth.good() && saveFileToDisplayGrowth.is_open())){
        cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowth<<endl;
        GrowthSaved = false;
    }

    saveFileString = saveDirectoryToDisplayString +"/Save_GrowthRate";
	const char* name_saveFileToDisplayGrowthRate = saveFileString.c_str();
	saveFileToDisplayGrowthRate.open(name_saveFileToDisplayGrowthRate, ifstream::in);
	if (!(saveFileToDisplayGrowthRate.good() && saveFileToDisplayGrowthRate.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowthRate<<endl;
		GrowthRateSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Force";
	const char* name_saveFileToDisplayForce = saveFileString.c_str();;
	saveFileToDisplayForce.open(name_saveFileToDisplayForce, ifstream::in);
	if (!(saveFileToDisplayForce.good() && saveFileToDisplayForce.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayForce<<endl;
		ForcesSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Proteins";
	const char* name_saveFileToDisplayProteins = saveFileString.c_str();;
	saveFileToDisplayProteins.open(name_saveFileToDisplayProteins, ifstream::in);
	if (!(saveFileToDisplayProteins.good() && saveFileToDisplayProteins.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayProteins<<endl;
		ProteinsSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Packing";
	const char* name_saveFileToDisplayPacking = saveFileString.c_str();;
	saveFileToDisplayPacking.open(name_saveFileToDisplayPacking, ifstream::in);
	if (!(saveFileToDisplayPacking.good() && saveFileToDisplayPacking.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayPacking<<endl;
		PackingSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_PhysicalProp";
	const char* name_saveFileToDisplayPhysicalProp = saveFileString.c_str();;
	saveFileToDisplayPhysicalProp.open(name_saveFileToDisplayPhysicalProp, ifstream::in);
	if (!(saveFileToDisplayPhysicalProp.good() && saveFileToDisplayPhysicalProp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayPhysicalProp<<endl;
		physicalPropertiesSaved = false;
	}

	return true;
}

bool Simulation::initiateSavedSystem(){
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
	//reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	initiateNodesFromSave();
	initiateElementsFromSave();
	assignPhysicalParameters();
	initiateSystemForces();
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}
    if (GrowthSaved){
        updateGrowthFromSave();
    }
    if (GrowthRateSaved){
        updateGrowthRateFromSave();
    }
	if (ForcesSaved){
		updateForcesFromSave();
	}
	updateElementVolumesAndTissuePlacements();
    //cleanMatrixUpdateData();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
    //calculateStiffnessMatrices();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	//skipping the footer:
	getline(saveFileToDisplayMesh,currline);
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		getline(saveFileToDisplayMesh,currline);
	}
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return true;
	}
	//cout<<"skipped footer: "<<currline<<endl;
	return true;
}

bool Simulation::readSystemSummaryFromSave(){
	string dummystring;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: TimeStep(sec):"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "TimeStep(sec):"
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: dt value"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dt;

	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: DataSaveInterval(sec):"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "DataSaveInterval(sec):"
	double dummydouble;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: save interval value"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummydouble;
	dataSaveInterval =  (int) ceil(dummydouble/dt);


	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: ModelinputName:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "ModelinputName: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting name of the model input file"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading the model input file
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Mesh_Type:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Mesh_Type: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: the mesh type(int)"<<endl;
		return false;
	}
	int dummyint;
	saveFileToDisplaySimSum >> dummyint; //reading the mesh type
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Symmetricity-x:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-x: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: symmetricitX boolean:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricX;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Symmetricity-y:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-y: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: symmetricitY boolean:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricY;
	cout<<"read the summary, symmetricity data: "<<symmetricX<<" "<<symmetricY<<endl;
	return true;
}

bool Simulation::readNodeDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	cout<<"number of nodes: "<<n<<endl;
	if(nNodes != n){
		cerr<<"The node number from save file("<<n<<") and model input("<<Nodes.size()<<") are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	for (int i=0; i<n; ++i){
		int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		if (Nodes[i]->tissuePlacement != tissuePlacement || Nodes[i]->atCircumference != atCircumference || Nodes[i]->tissueType != tissueType ){
			cerr<<"Node "<<i<<" properties are not consistent  - cannot continue simulation from save "<<endl;
			return false;
		}
	}
	return true;
}

void Simulation::initiateNodesFromSave(){
	saveFileToDisplayMesh >> nNodes;
	cout<<"number of nodes: "<<nNodes<<endl;
	Node* tmp_nd;
	for (int i=0; i<nNodes; ++i){
		double* pos = new double[3];
		int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> pos[0];
		saveFileToDisplayMesh >> pos[1];
		saveFileToDisplayMesh >> pos[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		tmp_nd = new Node(i, 3, pos,tissuePlacement, tissueType);
		tmp_nd-> atCircumference = atCircumference;
		Nodes.push_back(tmp_nd);
		delete[] pos;
	}
	cout<<"number of nodes: "<<nNodes<<endl;
}

void Simulation::initiateNodesFromMeshInput(){
	int n;
	saveFileToDisplayMesh >> n;
	Node* tmp_nd;
	for (int i=0; i<n; ++i){
		double* pos = new double[3];
		int tissuePos = -2;
		int tissueType = -2;
		int atCircumference;
		saveFileToDisplayMesh >> pos[0];
		saveFileToDisplayMesh >> pos[1];
		saveFileToDisplayMesh >> pos[2];
		saveFileToDisplayMesh >> tissuePos;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		tmp_nd = new Node(i, 3, pos,tissuePos, tissueType);
		tmp_nd->atCircumference = atCircumference;
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		delete[] pos;
	}
}

void Simulation::initiateElementsFromMeshInput(){
	int n;
	saveFileToDisplayMesh >> n;
	for (int i=0; i<n; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			initiatePrismFromMeshInput();
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}
	}
}


void Simulation::initiateElementsFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	cout<<"number of elements: "<<n<<endl;
	for (int i=0; i<n; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			initiatePrismFromSave();
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}
	}
	cout<<"number of elements: "<<nElements<<endl;
}

bool Simulation::readElementDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	if (nElements != n){
		cerr<<"The element number from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		cerr<<"n: "<<n<<" nElements: "<<endl;
		return false;
	}
	for (int i=0; i<nElements; ++i){
		//string line;
		//getline(saveFileToDisplayMesh, line);

		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (Elements[i]->getShapeType() != shapeType){
			cerr<<"The element type from save file and model input are not consistent - cannot continue simulation from save"<<endl;
			return false;
		}
		if (shapeType == 1){
			bool success = readShapeData(i);
			if (!success){
				cerr<<"Error reading shape data, element: "<<i<<endl;
				return false;
			}
		}
		else if (shapeType == 4){
			double height = 0.0;
			saveFileToDisplayMesh >> height;
			if (Elements[i]->ReferenceShape->height != height){
				cerr<<"The element height from save file and model input are not consistent - cannot continue simulation from save"<<endl;
				return false;
			}
			bool success = readShapeData(i);
			if (!success){
				cerr<<"Error reading shape data, element: "<<i<<endl;
				return false;
			}
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}

	}
	//string line;
	//getline(saveFileToDisplayMesh, line);

	return true;
}

void Simulation::initiatePrismFromSave(){
	//inserts a new prism at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		NodeIds[i] = 0;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
	//cout<<"Element: "<<PrismPnt01->Id<<endl;
	//Elements[Elements.size()-1]->displayPositions();
	//Elements[Elements.size()-1]->displayReferencePositions();
}

bool Simulation::readShapeData(int i){
	bool IsAblated;
	saveFileToDisplayMesh >> IsAblated;
	if ( Elements[i]->IsAblated != IsAblated){
		cerr<<"The element "<<i<<" ablation from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	bool success = Elements[i]->readNodeIdData(saveFileToDisplayMesh);
	if (!success){
		cerr<<"The element "<<i<<" node ids from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	success = Elements[i]->readReferencePositionData(saveFileToDisplayMesh);
	if (!success){
		cerr<<"The element "<<i<<" reference shape from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	return true;
}


void Simulation::initiatePrismFromMeshInput(){
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		int savedId;
		saveFileToDisplayMesh >> savedId;
		NodeIds[i] = savedId;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateReferencePositionMatrixFromMeshInput(saveFileToDisplayMesh);
	PrismPnt01->checkRotationConsistency3D();
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
	//Elements[Elements.size()-1]->displayReferencePositions();
}



void Simulation::reInitiateSystemForces(int oldSize){
	//deleting the old system forces:
	for (int j=0;j<oldSize;++j){
		delete[] SystemForces[j];
		delete[] PackingForces[j];
		//delete[] PackingForcesPreviousStep[j];
		//delete[] PackingForcesTwoStepsAgoStep[j];
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
	delete[] PackingForces;
	//delete[] PackingForcesPreviousStep;
	//delete[] PackingForcesTwoStepsAgoStep;
	delete[] FixedNodeForces;
	//reinitiating with the new size:
	const int n = nNodes;
	SystemForces = new double*[n];
	PackingForces = new double*[n];
	//PackingForcesPreviousStep = new double*[n];
	//PackingForcesTwoStepsAgoStep = new double*[n];
	FixedNodeForces = new double*[n];
	for (int j=0;j<n;++j){
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		//PackingForcesPreviousStep[j] = new double[3];
		//PackingForcesTwoStepsAgoStep[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
		//PackingForcesPreviousStep[j][0] = 0.0;
		//PackingForcesPreviousStep[j][1] = 0.0;
		//PackingForcesPreviousStep[j][2] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][0] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][1] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][2] = 0.0;
		FixedNodeForces[j][0] = 0.0;
		FixedNodeForces[j][1] = 0.0;
		FixedNodeForces[j][2] = 0.0;
	}
}

void Simulation::updateForcesFromSave(){
	for (int i=0;i<nNodes;++i){
		saveFileToDisplayForce.read((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][1], sizeof SystemForces[i][1]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][2], sizeof SystemForces[i][2]);
	}
	for (int i=0;i<nElements;++i){
		int n = Elements[i]->getNodeNumber();
		for (int j=0; j<n; ++j){
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][0], sizeof &Elements[i]->MyoForce[j][0]);
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][1], sizeof &Elements[i]->MyoForce[j][1]);
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][2], sizeof &Elements[i]->MyoForce[j][2]);
		}
	}
}

void Simulation::updateTensionCompressionFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		for (int j=0; j<6; ++j){
            double S = gsl_matrix_get((*itElement)->Strain,j,0);
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set((*itElement)->Strain,j,0,S);
        }
	}
}

void Simulation::updateGrowthFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		gsl_matrix* currFg = gsl_matrix_calloc(3,3);
        for (int j=0; j<3; ++j){
            for(int k=0; k<3; ++k){
                double Fgjk;
                saveFileToDisplayGrowth.read((char*) &Fgjk, sizeof Fgjk);
                gsl_matrix_set(currFg,j,k,Fgjk);
            }
        }
        (*itElement)->setFg(currFg);
        gsl_matrix_free(currFg);
    }
}

void Simulation::updateGrowthRateFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double rx=0.0,ry=0.0,rz=0.0;
		saveFileToDisplayGrowthRate.read((char*) &rx, sizeof rx);
		saveFileToDisplayGrowthRate.read((char*) &ry, sizeof ry);
		saveFileToDisplayGrowthRate.read((char*) &rz, sizeof rz);
		(*itElement)->setGrowthRateExpFromInput(rx,ry,rz);
	}
}

void Simulation::updateProteinsFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double c[4];
		double cEq[4];
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &c[j], sizeof c[j]);
		}
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &cEq[j], sizeof cEq[j]);
		}
		(*itElement)->setMyosinLevels(c[0], c[1], c[2], c[3]);
		(*itElement)->setEquilibriumMyosinLevels(cEq[0], cEq[1], cEq[2], cEq[3]);
	}
}

void Simulation::updatePhysicalPropFromSave(){
	readPhysicalPropToContinueFromSave();
}

void Simulation::updatePackingFromSave(){
	//reading the number of packing nodes:
	int n;
	saveFileToDisplayPacking.read((char*) &n, sizeof n);
	//emptying the packing node vectors:
	pacingNodeCouples0.empty();
	pacingNodeCouples1.empty();
	//filling in the packing node vectors for step
	for(int i=0; i<n; ++i){
		int elementId;
		saveFileToDisplayPacking.read((char*) &elementId, sizeof elementId);
		pacingNodeCouples0.push_back(elementId);
		saveFileToDisplayPacking.read((char*) &elementId, sizeof elementId);
		pacingNodeCouples1.push_back(elementId);
	}
	//reading the forces:
	for(int i=0; i<n; ++i){
		double Fx,Fy,Fz;
		saveFileToDisplayPacking.read((char*) &Fx, sizeof Fx);
		saveFileToDisplayPacking.read((char*) &Fy, sizeof Fy);
		saveFileToDisplayPacking.read((char*) &Fz, sizeof Fz);
		PackingForces[pacingNodeCouples0[i]][0] = Fx;
		PackingForces[pacingNodeCouples0[i]][1] = Fy;
		PackingForces[pacingNodeCouples0[i]][2] = Fz;
		saveFileToDisplayPacking.read((char*) &Fx, sizeof Fx);
		saveFileToDisplayPacking.read((char*) &Fy, sizeof Fy);
		saveFileToDisplayPacking.read((char*) &Fz, sizeof Fz);
		PackingForces[pacingNodeCouples1[i]][0] = Fx;
		PackingForces[pacingNodeCouples1[i]][1] = Fy;
		PackingForces[pacingNodeCouples1[i]][2] = Fz;
	}
}

void Simulation::readTensionCompressionToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		for (int j=0; j<6; ++j){
            double S;
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set((*itElement)->Strain,j,0,S);
		}
	}
}

void Simulation::readGrowthToContinueFromSave(){
	updateGrowthFromSave();
}

void Simulation::readGrowthRateToContinueFromSave(){
	updateGrowthRateFromSave();
}


void Simulation::readProteinsToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double c[4];
		double cEq[4];
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &c[j], sizeof c[j]);
		}
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &cEq[j], sizeof cEq[j]);
		}
		(*itElement)->setMyosinLevels(c[0], c[1], c[2], c[3]);
		(*itElement)->setEquilibriumMyosinLevels(cEq[0], cEq[1], cEq[2], cEq[3]);
	}
}

void Simulation::readPhysicalPropToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	double E;
	double internalViscposity;
	double externalViscosity[3];
	double zRemodelling;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		saveFileToDisplayPhysicalProp.read((char*) &E, sizeof E);
		saveFileToDisplayPhysicalProp.read((char*) &internalViscposity, sizeof internalViscposity);
		saveFileToDisplayPhysicalProp.read((char*) &zRemodelling, sizeof zRemodelling);
		(*itElement)->setYoungsModulus(E);
		(*itElement)->setViscosity(internalViscposity);
		(*itElement)->setZRemodellingSoFar(zRemodelling);
	}
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int i=0; i<3; ++i){
			saveFileToDisplayPhysicalProp.read((char*) &externalViscosity[i], sizeof externalViscosity[i]);
			(*itNode)->externalViscosity[i] = externalViscosity[i];
		}
	}
}


void Simulation::initiatePrismFromSaveForUpdate(int k){
	//inserts a new prism at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		//the node ids will be updated in function: updateShapeFromSave
		NodeIds[i] = 0;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	vector<ShapeBase*>::iterator it = Elements.begin();
	it += k;
	Elements.insert(it,PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
}


void Simulation::updateOneStepFromSave(){
	cout<<"updating step from save"<<endl;
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
	//cout<<"skipped header: "<<currline<<endl;
    if (implicitPacking){
        resetForces(true);	// reset packing forces
    }
    else{
        resetForces(false);	// do not reset packing forces
    }
    updateNodeNumberFromSave();
	updateNodePositionsFromSave();
	updateElementStatesFromSave();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		//This is updating positions from save.
		(*itElement)->updatePositions(Nodes);
	}
	if (TensionCompressionSaved){
        //cout<<"updating tension compression: "<<endl;
		updateTensionCompressionFromSave();
	}
    if (GrowthSaved){
        updateGrowthFromSave();
    }
    if (GrowthRateSaved){
    	updateGrowthRateFromSave();
    }
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if(ProteinsSaved){
		updateProteinsFromSave();
	}
	if(physicalPropertiesSaved){
		updatePhysicalPropFromSave();
	}
	//cout<<"trying to update packing from save, PackingSaved: "<<PackingSaved<<endl;
	if (PackingSaved){
		updatePackingFromSave();
		//cout<<"updated packing, the size of packing couples: "<<pacingNodeCouples0.size()<<endl;
	}
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
	calculateBoundingBox();
	//calculateColumnarLayerBoundingBox();
	//if (thereIsPeripodialMembrane){
	//	calculatePeripodialBoundingBox();
	//}
	//skipping the footer:
	string currline2;
	getline(saveFileToDisplayMesh,currline2);
	//cout<<"currline 1st reading: "<<currline2<<endl;
	getline(saveFileToDisplayMesh,currline2);
	//cout<<"currline 2nd reading: "<<currline2<<endl;
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		//cout<<"skipping empty line"<<endl;
		getline(saveFileToDisplayMesh,currline2);
	}
	//if(saveFileToDisplayMesh.eof()){
	//	reachedEndOfSaveFile = true;
	//	return;
	//}
	//cout<<"in step update, skipped footer: "<<currline2<<endl;
	timestep = timestep + dataSaveInterval;
	currSimTimeSec += dt*dataSaveInterval;
}

void  Simulation::updateNodeNumberFromSave(){
	//cout<<"Updating number of nodes from save"<<endl;
	//cout<<"Is save file open: "<<saveFileToDisplay.is_open()<<" is file good? "<<saveFileToDisplay.good()<<endl;
	int n;
	saveFileToDisplayMesh>> n;
	int currNodeNumber = Nodes.size();
	//cout<<"number of nodes from save file: "<<n <<" number of nodes on the vector: "<<currNodeNumber<<endl;
	if (n>currNodeNumber){
		Node* tmp_nd;
		for (int i = 0; i<(n-currNodeNumber); ++i){
			double* pos = new double[3];
			pos[0]=0.0;
			pos[1]=0.0;
			pos[2]=0.0;
			//the positions will be read and updated in function updateNodePositionsFromSave
			tmp_nd = new Node(i, 3, pos,-1, -1);
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			delete[] pos;
		}
	}
	else{
		for (int i = 0; i<(currNodeNumber-n); ++i){
			Node* tmp_nd;
			tmp_nd = Nodes.back();
			Nodes.pop_back();
			nNodes = Nodes.size();
			delete tmp_nd;
		}
	}
	n = Nodes.size();
	if ( n != nNodes){
		//the node number is change, I updated the node list, now I need to fix system forces:
		reInitiateSystemForces(nElements);
	}
	//cout<<"end of funciton, number of nodes from save file: "<<n <<" number of nodes on the vector: "<<Nodes.size()<<endl;
}

void  Simulation::updateElementStatesFromSave(){
	//cout<<"Updating element states from save"<<endl;
	int n;
	saveFileToDisplayMesh >> n;
	int currElementNumber = Elements.size();
	//The elements list is bigger than necessary, I am deleting elements form the end of the list
	//cout<<"number of elements from save file: "<<nElements <<" number of element on the Elements vector: "<<currElementNumber<<endl;
	while(currElementNumber>n){
		//cout<<"removing elements from list, current number vector length: "<<currElementNumber<<endl;
		removeElementFromEndOfList();
		currElementNumber = Elements.size();
	}
	for (int i=0; i<currElementNumber; ++i){
		//cout<<"reading element: "<<i<<endl;
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		int currShapeType = Elements[i]->getShapeType();
		if (shapeType == currShapeType || shapeType == 2 || shapeType ==4){
			//cout<<"shape type is correct, moving on to update"<<endl;
			//The current shape on the list, and the shape I am reading are of same type, I can read it:
			if (shapeType ==4){
				double height;
				saveFileToDisplayMesh >> height;
			}
			Elements[i]->updateShapeFromSave(saveFileToDisplayMesh);
		}
		else{
			//cout<<"shape type is wrong"<<endl;
			//the current element is a different type, I need to insert generate a new element, and insert it here:
			if (shapeType == 1){
				//the new shape is a prism, I will add it now;
				initiatePrismFromSaveForUpdate(i);
			}
			removeElementFromEndOfList();
			currElementNumber = Elements.size();
		}
	}
	while(n>currElementNumber){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			int i = Elements.size()-1;
			//this will initiate a prism at the current point in the Elements vector
			//then it will read the node ids and reference positions from the save file
			//it will update the node ids, and reference positions.
			//the normal positions will be updated using function updatePositions, called by updateOneStepFromSave
			initiatePrismFromSaveForUpdate(i);
			currElementNumber = Elements.size();
		}
	}
}

void Simulation::removeElementFromEndOfList(){
	ShapeBase* tmp_element;
	tmp_element = Elements.back();
	//int shapetype = tmp_element->getShapeType();
	Elements.pop_back();
	nElements = Elements.size();
	delete tmp_element;
}

void Simulation::updateNodePositionsFromSave(){
	//cout<<"Updating node positions from save"<<endl;
	//I have already read the current number of nodes in function "updateNodeNumberFromSave",
	//and I updated the number of nodes where necessary
	//the "cursor" in the file progressed, and is at the beginning of positions now
	for (int i=0; i<nNodes; ++i){
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> Nodes[i]->tissuePlacement;
		saveFileToDisplayMesh >> Nodes[i]->tissueType;
		saveFileToDisplayMesh >> Nodes[i]->atCircumference;
	}
}

bool Simulation::initiateMesh(int MeshType, float zHeight){
	if (MeshType == 1 ){
		initiateSinglePrismNodes(zHeight);
		initiateSinglePrismElement();
	}
	else if (MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 || MeshType ==4 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
		return false;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

bool Simulation::initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		//The necessary parameters:
		//InputMeshParameters:
		//  MeshInputMode(int-seeDocumentation): 2
		//  MeshRow(int): 3
		//  MeshColumn(int): 1
		//  SideLength: 1.0
		//  zHeight: 1.0
		//  ApicalCircumferenceFixZ: 0
		//  ApicalCircumferenceFixXY: 0
		//  BasalCircumferenceFixZ: 0
		//  BasalCircumferenceFixXY: 0
		initiateNodesByRowAndColumn(Row,Column,SideLength,zHeight);
		initiateElementsByRowAndColumn(Row,Column);
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments for mesh triangulation"<<endl;
		return false;
	}
	else if ( MeshType ==4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		return false;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

/*bool Simulation::initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight ){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 ){
		//this will be inputting circumference of the tissue, and the sidelength and z-height
		//generate mesh by triangulation
		return false;
	}
	else if ( MeshType == 4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		return false;
		//this will be reading full mesh data
		//read mesh data from file
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}*/

bool Simulation::initiateMesh(int MeshType){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
		return false;
	}
	else if ( MeshType ==4 ){
		//this will be reading full mesh data
		//read mesh data from file
		const char* name_inputMeshFile = inputMeshFileName.c_str();;
		saveFileToDisplayMesh.open(name_inputMeshFile, ifstream::in);
		if (!(saveFileToDisplayMesh.good() && saveFileToDisplayMesh.is_open())){
			cerr<<"Cannot open the save file to display: "<<name_inputMeshFile<<endl;
			return false;
		}
		cout<<"initiating nodes"<<endl;
		initiateNodesFromMeshInput();
		initiateElementsFromMeshInput();
		cout<<" nNodes: "<<nNodes<<" nEle: "<<nElements<<endl;
		bool areTissueWeightsRecorded = checkIfTissueWeightsRecorded();
		if (areTissueWeightsRecorded){
			readInTissueWeights();
			cout<<" read in tissue weights"<<endl;
		}

		//addCurvatureToColumnar(5.0);
		saveFileToDisplayMesh.close();
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

bool Simulation::checkIfTissueWeightsRecorded(){
	bool tissueWeightsRecorded;
	saveFileToDisplayMesh >> tissueWeightsRecorded;
	return tissueWeightsRecorded;
}

void Simulation::readInTissueWeights(){
	vector<ShapeBase*>::iterator itEle;
	double wPeri;
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		saveFileToDisplayMesh >> wPeri;
		(*itEle)->setGrowthWeightsViaTissuePlacement(wPeri); //peripodialness weight is recorded
	}
}

bool Simulation::generateColumnarCircumferenceNodeList(	vector <int> &ColumnarCircumferencialNodeList){
	//generating a list of nodes that are at the circumference and at the basal surface
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->atCircumference && Nodes[i]->tissuePlacement == 0){ // tissuePlacement = 0 -> basal node
			ColumnarCircumferencialNodeList.push_back(i);

		}
	}
	int n = ColumnarCircumferencialNodeList.size();
	if (n<=0){
		cerr<<"No circumferncial nodes indicated! Cannot generate PeripodialMembrane"<<endl;
		AddPeripodialMembrane = false;
		thereIsPeripodialMembrane = false;
		return false;
	}
	return true;
}

void Simulation::clearCircumferenceDataFromSymmetricityLine(){
	//I will take out anything that was at the border of symmetry, but I need the nodes at the tips. So find the x tips first:
	if (symmetricY){
		double xTipPos = -1000, xTipNeg = 1000;
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->Position[0] > xTipPos ){
				xTipPos = (*itNode)->Position[0];
			}
			if ((*itNode)->Position[0] < xTipNeg ){
				xTipNeg = (*itNode)->Position[0];
			}
		}
		double yLimPos = 0.1;
		double yLimNeg = (-1.0)*yLimPos;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			double x = (*itNode)->Position[0];
			double y = (*itNode)->Position[1];
			if ( y < yLimPos){
				if ( y  > yLimNeg){
					(*itNode)->atSymmetricityBorder = true;
					fixY((*itNode),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
					bool atTip = false;
					if ( x < xTipPos+yLimPos && x >xTipPos+yLimNeg){
						atTip = true;
					}
					if (!symmetricX){
						//This if clause is checking the x-tip at the negative end.
						//As described above, the node should still be a circumference node,
						//if it is at the line of symmetry, but at the tip as well.
						//BUT, if there is also xSymmetry, this "negative end tip" will
						//be the tip at the x symmetry line. (If I am modelling a quarter of a circle,
						//with x and y symmetry, this will point be the centre of the circle.)
						//Under these conditions, it is not actually at the tip. It should be
						//removed from circumference node list.
						//Therefore, I am carrying the atTip? check for the negative end only if
						//there is no x-symmetry.
						if ( x > xTipNeg+yLimNeg && x < xTipNeg+yLimPos){
							atTip = true;
						}
					}
					if (!atTip){
						//removing node from list:
						(*itNode)->atCircumference = false;
					}
				}
			}
		}
	}
	if (symmetricX){
		double yTipPos = -1000, yTipNeg = 1000;
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->Position[0] > yTipPos ){
				yTipPos = (*itNode)->Position[1];
			}
			if ((*itNode)->Position[0] < yTipNeg ){
				yTipNeg = (*itNode)->Position[1];
			}
		}
		double xLimPos = 0.1;
		double xLimNeg = (-1.0)*xLimPos;
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			double x = (*itNode)->Position[0];
			double y = (*itNode)->Position[1];
			if ( x < xLimPos){
				if ( x  > xLimNeg){
					(*itNode)->atSymmetricityBorder = true;
					fixX((*itNode),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
					bool atTip = false;
					if ( y < yTipPos+xLimPos && y >yTipPos+xLimNeg){
						atTip = true;
					}
					if ( y > yTipNeg+xLimNeg && y < yTipNeg+xLimPos){
						atTip = true;
					}
					if (!atTip){
						//removing node from list:
						(*itNode)->atCircumference = false;
					}
				}
			}
		}
	}
}
/*
void Simulation::removeSymmetryBorderFromColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList){
	//I will take out anything that was at the border of symmetry, but I need the nodes at the tips. So find the x tips first:
	double xTipPos = -1000, xTipNeg = 1000;
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->Position[0] > xTipPos ){
			xTipPos = (*itNode)->Position[0];
		}
		if ((*itNode)->Position[0] < xTipNeg ){
			xTipNeg = (*itNode)->Position[0];
		}
	}
	double yLimPos = 0.1;
	double yLimNeg = (-1.0)*yLimPos;
	int n = ColumnarCircumferencialNodeList.size();
	int i=0;
	while (i<n){
		double x = Nodes[ColumnarCircumferencialNodeList[i]]->Position[0];
		double y = Nodes[ColumnarCircumferencialNodeList[i]]->Position[1];
		if ( y < yLimPos){
			if ( y  > yLimNeg){
				Nodes[ColumnarCircumferencialNodeList[i]]->atSymmetricityBorder = true;
				fixY(Nodes[ColumnarCircumferencialNodeList[i]]);
				//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
				bool atTip = false;
				if ( x < xTipPos+yLimPos && x >xTipPos+yLimNeg){
					atTip = true;
				}
				if ( x > xTipNeg+yLimNeg && x < xTipNeg+yLimPos){
					atTip = true;
				}
				if (!atTip){
					//removing node from list:
					ColumnarCircumferencialNodeList.erase (ColumnarCircumferencialNodeList.begin()+i);
					Nodes[ColumnarCircumferencialNodeList[i]]->atCircumference = false;
					n--;
					i--;
				}
			}
		}
		i++;
	}
}
*/
void Simulation::sortColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList){
	//ordering the circumferencial nodes of the basal surface in clockwise rotation
	int n = ColumnarCircumferencialNodeList.size();
	vector <double> angles;
	for (int j =0 ; j<n; ++j){
		double x = Nodes[ColumnarCircumferencialNodeList[j]]->Position[0];
		double y = Nodes[ColumnarCircumferencialNodeList[j]]->Position[1];
		double tet = atan2(y,x);
		if (tet<0){tet += 2.0*3.14;}
		angles.push_back(tet);
	}

	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			if(angles[i]<angles[i-1]){
				int temp=ColumnarCircumferencialNodeList[i-1];
				ColumnarCircumferencialNodeList[i-1]=ColumnarCircumferencialNodeList[i];
				ColumnarCircumferencialNodeList[i]=temp;
				double td = angles[i-1];
				angles[i-1]=angles[i];
				angles[i]=td;
				swapped = true;
				//cout<<"swapped "<<i <<" with "<<i-1<<endl;
			}
		}
	}
}

void Simulation::getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength){
	double dsumPeri =0.0, dsumCol = 0.0;
	int colCounter =0, periCounter=0;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated){
			//do not count ablated elements
			if ((*itElement)->tissueType==0){ //element belongs to columnar layer
				dsumCol += (*itElement)->getApicalSideLengthAverage();
				colCounter++;
			}
			else{
				dsumPeri += (*itElement)->getApicalSideLengthAverage();
				periCounter++;
			}
		}
	}
	colAverageSideLength = dsumCol / (double)colCounter;
	if (periCounter>0){
		periAverageSideLength = dsumPeri / (double)periCounter;
	}
}

bool  Simulation::isColumnarLayer3D(){
	bool ColumnarLayer3D = false;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->ShapeDim == 3){
			ColumnarLayer3D  = true;
			break;
		}
	}
	return ColumnarLayer3D;
}

bool Simulation::calculateTissueHeight(){
	//check if the columnar layer is made of 3D elements:
	bool ColumnarLayer3D = isColumnarLayer3D();
	TissueHeight = 0;
	if (ColumnarLayer3D){
		//Find the first basal node on the Nodes List
		//Find the first element that has this node on the elements list
		//Move apically through the elements until you reach the apical surface - an apical node
		//Find a basal node:
		vector<Node*>::iterator itNode;
		bool foundNode = false;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if((*itNode)->tissueType == 0 && (*itNode)->tissuePlacement == 0){ //Columnar node is basal
				foundNode = true;
				break;
			}
		}
		if (!foundNode){
			return false;
		}
		//Find an element using the basal node, and move on the elements apically, until you reach the apical surface:
		int currNodeId = (*itNode)->Id;
		vector<ShapeBase*>::iterator itElement;
		bool foundElement = true;
		TissueHeightDiscretisationLayers = 0;
		while(Nodes[currNodeId]->tissuePlacement != 1 && foundElement){ //while the node is not apical, and I could find the next element
			foundElement = false;
			for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					foundElement = true;
					break;
				}
			}
			double currentH = (*itElement)->getElementHeight();
			TissueHeight += currentH;
			TissueHeightDiscretisationLayers++;
			currNodeId = (*itElement)->getCorrecpondingApical(currNodeId); //have the next node
		}
		if (!foundElement){
			return false;
		}
	}
	else{
		TissueHeight = Elements[0]->ReferenceShape->height;
		if (TissueHeight == 0){
			cout<<"The coulmanr layer is 2D, but the tissue height of the elements is not assigned properly, cannot obtain TissueHeight"<<endl;
			return false;
		}
	}
	if (thereIsPeripodialMembrane){
		double columnarTop = -1000.0, peripodialbottom = 1000.0;
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->tissueType == 0 ){
				//columnar node
				if ((*itNode)->tissuePlacement == 1){
					//apical node:
					if((*itNode)->Position[2]> columnarTop){
						columnarTop = (*itNode)->Position[2];
					}
				}
			}
			else if ((*itNode)->tissueType == 1 ){
				//peripodial node
				if ((*itNode)->tissuePlacement == 0){
					//basal node:
					if((*itNode)->Position[2]< peripodialbottom){
						peripodialbottom = (*itNode)->Position[2];
					}
				}
			}
		}
		lumenHeight = columnarTop - peripodialbottom;
	}
	return true;
}

void Simulation::assignInitialZPositions(){
	for(vector<ShapeBase*>::iterator itElement = Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setInitialZPosition(boundingBox[0][2], TissueHeight);  //minimum z of the tissue and the tissue height given as input
	}
}

void Simulation::calculateStiffnessMatrices(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<" setting up element :  "<<(*itElement)->Id<<" of "<<nElements<<endl;
		(*itElement)->calculateReferenceStiffnessMatrix();
	}
}

void Simulation::calculateShapeFunctionDerivatives(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<" setting up element :  "<<(*itElement)->Id<<" of "<<nElements<<endl;
		(*itElement)->calculateElementShapeFunctionDerivatives();
    }
}

void Simulation::fixAllD(Node* currNode, bool fixWithViscosity){
	for (int j =0 ; j<currNode->nDim; ++j){
		if(fixWithViscosity){
			currNode->externalViscosity[j] = fixingExternalViscosity[j];
			currNode->baseExternalViscosity[j] = currNode->externalViscosity[j];
			currNode->externalViscositySetInFixing[j] = true;
		}
		else{
			currNode->FixedPos[j]=true;
		}
	}
	//cout<<" fixAllD called on node: "<<currNode->Id<<endl;
}

void Simulation::fixAllD(int i, bool fixWithViscosity){
	for (int j =0 ; j<Nodes[i]->nDim; ++j){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[j] = fixingExternalViscosity[j];
			Nodes[i]->baseExternalViscosity[j] = Nodes[i]->externalViscosity[j];
			Nodes[i]->externalViscositySetInFixing[j] = true;
		}
		else{
			Nodes[i]->FixedPos[j]=true;
		}
	}
}

void Simulation::fixX(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>0){
		if(fixWithViscosity){
			currNode->externalViscosity[0] = fixingExternalViscosity[0];
			currNode->baseExternalViscosity[0] = currNode->externalViscosity[0];
			currNode->externalViscositySetInFixing[0] = true;
		}
		else{
			currNode->FixedPos[0]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have x-dimension"<<endl;
	}
}
void Simulation::fixX(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>0){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[0] = fixingExternalViscosity[0];
			Nodes[i]->baseExternalViscosity[0] = Nodes[i]->externalViscosity[0];
			Nodes[i]->externalViscositySetInFixing[0] = true;
		}
		else{
			Nodes[i]->FixedPos[0]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have x-dimension"<<endl;
	}

}

void Simulation::fixY(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>1){
		if(fixWithViscosity){
			currNode->externalViscosity[1] = fixingExternalViscosity[1];
			currNode->baseExternalViscosity[1] = currNode->externalViscosity[1];
			currNode->externalViscositySetInFixing[1] = true;
		}
		else{
			currNode->FixedPos[1]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have y-dimension"<<endl;
	}
}

void Simulation::fixY(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>1){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[1] = fixingExternalViscosity[1];
			Nodes[i]->baseExternalViscosity[1] = Nodes[i]->externalViscosity[1];
			Nodes[i]->externalViscositySetInFixing[1] = true;
		}
		else{
			Nodes[i]->FixedPos[1]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have y-dimension"<<endl;
	}
}

void Simulation::fixZ(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>2){
		if(fixWithViscosity){
			currNode->externalViscosity[2] = fixingExternalViscosity[2];
			currNode->baseExternalViscosity[2] = currNode->externalViscosity[2];
			currNode->externalViscositySetInFixing[2] = true;
		}
		else{
			currNode->FixedPos[2]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have z-dimension"<<endl;
	}
}

void Simulation::fixZ(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>2){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[2] = fixingExternalViscosity[2];
			Nodes[i]->baseExternalViscosity[2] = Nodes[i]->externalViscosity[2];
			Nodes[i]->externalViscositySetInFixing[2] = true;
		}
		else{
			Nodes[i]->FixedPos[2]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have z-dimension"<<endl;
	}
}

void Simulation::zeroForcesOnNode(int i){
	double ForceBalance[3];
	ForceBalance[0] = SystemForces[i][0];
	ForceBalance[1] = SystemForces[i][1];
	ForceBalance[2] = SystemForces[i][2];
	for (int i=0;i<nNodes;++i){
		SystemForces[i][0]-=ForceBalance[0];
		SystemForces[i][1]-=ForceBalance[1];
		SystemForces[i][2]-=ForceBalance[2];
	}
}

void Simulation::initiateSystemForces(){
	const int n = nNodes;
	//n nodes
	SystemForces = new double*[n];
	PackingForces = new double*[n];
	//PackingForcesPreviousStep = new double*[n];
	//PackingForcesTwoStepsAgoStep = new double*[n];
	FixedNodeForces = new double*[n];
	for (int j=0;j<n;++j){
		//3 dimensions
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		//PackingForcesPreviousStep[j] = new double[3];
		//PackingForcesTwoStepsAgoStep[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
		//PackingForcesPreviousStep[j][0]=0.0;
		//PackingForcesPreviousStep[j][1]=0.0;
		//PackingForcesPreviousStep[j][2]=0.0;
		//PackingForcesTwoStepsAgoStep[j][0]=0.0;
		//PackingForcesTwoStepsAgoStep[j][1]=0.0;
		//PackingForcesTwoStepsAgoStep[j][2]=0.0;
		FixedNodeForces[j][0] = 0.0;
		FixedNodeForces[j][1] = 0.0;
		FixedNodeForces[j][2] = 0.0;
		//cout<<"systemforces[i][j]: "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<endl;
	}
}

void Simulation::initiateSinglePrismNodes(float zHeight){
	double *pos = new double[3];
	Node* tmp_nd;
	pos[0]=0;pos[1]=1;pos[2]=0;
	tmp_nd = new Node(0, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=1;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(1, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(2, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=1;pos[2]=zHeight;
	tmp_nd = new Node(3, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=1;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(4, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(5, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	delete[] pos;
}

void Simulation::initiateSinglePrismElement(){
	int* NodeIds;
	NodeIds = new int[6];
	for (int i = 0; i < 6 ; i++){
		NodeIds[i]=i;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	fixZ(0, fixWithExternalViscosity);
	fixZ(1, fixWithExternalViscosity);
	fixZ(2, fixWithExternalViscosity);
}


void Simulation::initiateNodesByRowAndColumn(int Row, int Column, float SideLength, float zHeight){
	dorsalTipIndex = Row;
	ventralTipIndex = 0;
	//The height of the equilateral triangle with side length: SideLength
	double sqrt3 = 1.7321;
	float h = sqrt3/2*SideLength;
	vector <double> xPos, yPos;
	vector <int> NodesToFix;
	int toprowcounter = 0;	//number of nodes at the terminating end, the top row. I will need this to add the sideway prisms;
	for (int ColCount = 0; ColCount < Column+1; ++ColCount){
		double CurrY = ColCount*h;
		int CurrRowNum = Row + 1 - ColCount;
		double RowOffset = 0.5*SideLength*ColCount;
		for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
			double CurrX = RowOffset + RowCount * SideLength;
			xPos.push_back(CurrX);
			yPos.push_back(CurrY);
			if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
				NodesToFix.push_back(xPos.size()-1);
				if(ColCount == Column){
					toprowcounter++;
				}
			}
		}
		if (ColCount>0){
			CurrY = (-1.0)*CurrY;
			for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
				double CurrX = RowOffset + RowCount * SideLength;
				xPos.push_back(CurrX);
				yPos.push_back(CurrY);
				if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
					NodesToFix.push_back(xPos.size()-1);
				}
			}
		}
	}
	int n =  xPos.size();
	Node* tmp_nd;
	double* pos = new double[3];
	//Adding the basal level of nodes, all will form columnar elements:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = 0.0;
		tmp_nd = new Node(i, 3, pos,0,0);
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
	}
	//Adding the apical level, all will form columnar elements:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = zHeight;
		tmp_nd = new Node(n+i, 3, pos,1,0);
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		nNodes = Nodes.size();
	}
	delete[] pos;
}

void Simulation::setLinkerCircumference(){
	//First find the node that is further back on
	vector<Node*>::iterator itNode;
	double maxX = -100, ZofTip = -100;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 2){
			//The node is linker, is the x position higher than already recorded?
			if ((*itNode)->Position[0]> maxX){
				maxX  = (*itNode)->Position[0];
				ZofTip = (*itNode)->Position[2];
			}
		}
	}
	double thres = 0.2;
	//cout<<"Setting circumference, maxX: "<<maxX<<" Z of tip: "<<ZofTip<<endl;
	//Now I have the maxX, and the corresponding z height.
	//Declare all linkers at the basal/or apical surface are the circumferential nodes:
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 2){
			//The node is linker, if it is apical or basal, it can be circumferntial:
			if ( (*itNode)->tissuePlacement == 0 || (*itNode)->tissuePlacement == 1 ){
				//The node is linker, if it is in the range of the height of the tip then it is circumferential
				if ( (*itNode)->Position[2] < ZofTip+thres && (*itNode)->Position[2] > ZofTip-thres ){
					(*itNode)->atCircumference = true;
				}
			}

		}
	}
}

void Simulation::checkForNodeFixing(){
	//Are there any circumferential node fixing options enabled:
	bool thereIsCircumFix = false;
	if (!thereIsCircumFix){
		for (int i=0;i<5; ++i){
			for (int j=0;j<3; ++j){
				if (CircumferentialNodeFix[i][j] == true){
					thereIsCircumFix = true;
					break;
				}
			}
			if (thereIsCircumFix){
				break;
			}
		}
	}
	//If there is any circumferential Node fixing:
	if (thereIsCircumFix){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->atCircumference){
				//cout<<"Node "<<(*itNode)->Id<<" at circumference"<<endl;
				//The node is at circumference, I would like to fix only the columnar side,
				bool doNotFix = false;
				if ((*itNode)->tissueType == 1){ //the node is peripodial, no not fix this node
					doNotFix = true;
				}
				/*else if ((*itNode)->tissueType == 2){ //the node is linker (as will be the case for all circumference)
					//I will check if any of its neigs are peripodial:
					int nNeig = (*itNode)->immediateNeigs.size();
					for (int k = 0; k<nNeig; ++k){
						if (Nodes[(*itNode)->immediateNeigs[k]]->tissueType == 1){
							doNotFix = true;
							break;
						}
					}
				}*/
				if (!doNotFix){
					//cout<<"Node "<<(*itNode)->Id<<" in fix options"<<endl;
					//Node is at the circumference, now checking for all possibilities:
					// if i == 0 , I am checking for apical circumference
					// if i == 1 , I am checking for basal  circumference
					// if i == 2 , I am checking for the linker apical circumference
					// if i == 3 , I am checking for the linker basal circumference
					// if i == 4 , I am checking for all    circumference
					for (int i=0;i<5; ++i){
						if ( (i == 0 && (*itNode)->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
							 (i == 1 && (*itNode)->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
							 (i == 2 && (*itNode)->tissueType == 2 && (*itNode)->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
							 (i == 3 && (*itNode)->tissueType == 2 && (*itNode)->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
							 (i == 4)){										  //tissuePlacement is irrelevant, fixing all
							//The node is at circumference; if
							if (CircumferentialNodeFix[i][0]){
								fixX((*itNode),fixWithExternalViscosity);
							}
							if (CircumferentialNodeFix[i][1]){
								fixY((*itNode),fixWithExternalViscosity);
							}
							if (CircumferentialNodeFix[i][2]){
								fixZ((*itNode),fixWithExternalViscosity);
							}
						}
					}
				}
			}
		}
	}
	if (BasalNodeFix[0] || BasalNodeFix[1] || BasalNodeFix[2] ){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 0){
				if (BasalNodeFix[0]){
					fixX((*itNode), fixWithExternalViscosity);
				}
				if (BasalNodeFix[1]){
					fixY((*itNode), fixWithExternalViscosity);
				}
				if (BasalNodeFix[2]){
					fixZ((*itNode), fixWithExternalViscosity);
				}
			}
		}
	}
	if (ApicalNodeFix[0] || ApicalNodeFix[1]  || ApicalNodeFix[2]){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 1){
				if (ApicalNodeFix[0]){
					fixX((*itNode), fixWithExternalViscosity);
				}
				if (ApicalNodeFix[1]){
					fixY((*itNode), fixWithExternalViscosity);
				}
				if (ApicalNodeFix[2]){
					fixZ((*itNode), fixWithExternalViscosity);
				}
			}
		}
	}
}

void Simulation::initiateElementsByRowAndColumn(int Row, int Column){
	int xinit1 = 0;
	int xinit2 = xinit1+Row+1;
	int xinit3 = 0;
	int xinit4 = xinit1+2*(Row+1)-1;
	int n = Nodes.size() /2.0;
    //initialising the tissue elements:
	for (int ColCount = 0; ColCount < Column; ++ColCount){
		int CurrRowNum = Row + 1 - ColCount;
		for (int RowCount = 0; RowCount<CurrRowNum-1; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit1+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
		}
		for (int RowCount = 0; RowCount<CurrRowNum-2; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit2+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

		}
		xinit1 = xinit2;
		xinit2 = xinit4 + CurrRowNum-1;
		xinit3 = xinit4;
		xinit4 = xinit2 + CurrRowNum-2;
	}
	cout<<"finalised element initiation"<<endl;
}

void Simulation::calculateSystemCentre(){
	for (int i = 0; i< nNodes; ++i){
		for (int j =0; j<Nodes[i]->nDim; ++j){
			SystemCentre[j] += Nodes[i]->Position[j];
		}
	}
	SystemCentre[0]= SystemCentre[0]/nNodes;
	SystemCentre[1]= SystemCentre[1]/nNodes;
	SystemCentre[2]= SystemCentre[2]/nNodes;
}

void Simulation::addSoftPeriphery(double* fractions){
	double t2 = softDepth*softDepth;
	double midlineFraction = 0.0;
	double currSoftnessFraction = softnessFraction;
	double tissueTypeMultiplier = 1;
	if (softPeripheryBooleans[0] && softPeripheryBooleans[1]){
		midlineFraction = softnessFraction;
	}
	else if (softPeripheryBooleans[0] || softPeripheryBooleans[1]){
		midlineFraction = 1.0 - (1.0-softnessFraction)*0.5;
	}
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ( (*itNode)->atCircumference ){
			//this node is at circumference, I will calculate the distance of all elements to this node
			//if an element is close enough, update the fraction matrix set above:
			for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool applyToThisElement = true;
				if ( (*itElement)->tissuePlacement == 1 && !softPeripheryBooleans[0]){
					//element is apical, the softness is NOT applied to apical
					applyToThisElement = false;
				}
				if ((*itElement)->tissuePlacement == 0 && !softPeripheryBooleans[1]){
					//element is basal, the softness is NOT applied to apical
					applyToThisElement = false;
				}
				if ((*itElement)->tissueType == 0 && !softPeripheryBooleans[2]){
					//element is columnar, the softness is NOT applied to columnar
					applyToThisElement = false;
				}
				if ((*itElement)->tissueType == 1 && !softPeripheryBooleans[3]){
					//element is peripodial, the softness is NOT applied to peripodial
					applyToThisElement = false;
				}
				if ( applyToThisElement ){
					//scaling the softness fraction if necessary:
					currSoftnessFraction  = softnessFraction;
					if ((*itElement)->tissuePlacement == 2 ){ //element is on the midline
						currSoftnessFraction= midlineFraction;
					}
					//scaling the tissue type multiplier if necessary (important for linker elements:
					if ((*itElement)->tissueType == 2 ){ // element is on the linker zone
						if (softPeripheryBooleans[2] && softPeripheryBooleans[3]){
							//softness is applied to both peripodial and columnar zones, the linkers should not be scaling anything
							tissueTypeMultiplier = 1.0;
						}else if (softPeripheryBooleans[2] ){
							//softness only applied to columnar layer:
							tissueTypeMultiplier = (*itElement)->getColumnarness();
						}
						else if (softPeripheryBooleans[3] ){
							//softness only applied to peripodial layer:
							tissueTypeMultiplier = (*itElement)->getPeripodialness();
						}
					}
					else{ //element is not linker
						tissueTypeMultiplier = 1;
					}
					currSoftnessFraction  = 1.0 - (1.0-currSoftnessFraction)*tissueTypeMultiplier;
					double *c = (*itElement)->getCentre();
					double dx = c[0]-(*itNode)->Position[0];
					double dy = c[1]-(*itNode)->Position[1];
					//double dz = c[2]-(*itNode1)->Position[2];
					double d = dx*dx + dy*dy;
					if (d<t2){
						d = pow(d,0.5);
						double f = currSoftnessFraction + (1.0 - currSoftnessFraction)*d/softDepth;
						if ((softnessFraction< 1.0 && f < fractions[(*itElement)->Id]) ||(softnessFraction> 1.0 && f > fractions[(*itElement)->Id])){
							fractions[(*itElement)->Id] = f;
						}
					}
					delete[] c;
				}
			}
		}
	}
}

void Simulation::assignPhysicalParameters(){
	double* fractions;
	fractions = new double[(const int) nElements];
	for (int i=0; i<nElements; ++i){
		fractions[i] = 1.0;
	}
	if(softPeriphery){
		addSoftPeriphery(fractions);
	}
	//now I have softness fraction table, I can scale the element properties:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double r = (rand() % 200) / 100.0;	//random number between 0.00 and 2.00
		r = r - 1.0; 						//random number between -1.00 and 1.00
		float noise1 = r*noiseOnPysProp[0];	//percent noise on current element
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise2 = r*noiseOnPysProp[1];
		if ((*itElement)->tissueType == 0){ //Element is on the columnar layer
			double currEApical 	= fractions[(*itElement)->Id] * EApical*(1 + noise1/100.0);
			double currEBasal	= fractions[(*itElement)->Id] * EBasal*(1 + noise1/100.0);
			double currEMid		= fractions[(*itElement)->Id] * EMid*(1 + noise1/100.0);
			double currPoisson = poisson*(1 + noise2/100);
			if((*itElement)->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
			(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currPoisson);
			(*itElement)->setViscosity(discProperApicalViscosity,discProperBasalViscosity,discProperMidlineViscosity);
		}
		else if ((*itElement)->tissueType == 1){ //Element is on the peripodial membrane
			double currEApical = fractions[(*itElement)->Id] * PeripodialElasticity*(1 + noise1/100.0);
			double currEBasal = currEApical;
			double currEMid = currEApical;
			double currPoisson = poisson*(1 + noise2/100);
			if(thereIsExplicitECM){
				//The input file does not define apical and basal stiffness separately
				//for each element. If I have explicit ECM, I will change the basal ECM
				//stiffness such that it will be the same as basal of the columnar layer,
				//therefore the ECM.
				currEBasal = fractions[(*itElement)->Id] * EBasal*(1 + noise1/100.0);
				currEMid = fractions[(*itElement)->Id] * EMid*(1 + noise1/100.0);
			}
			if((*itElement)->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
			(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currPoisson);
			(*itElement)->setViscosity(peripodialApicalViscosity,peripodialBasalViscosity,peripodialMidlineViscosity);
		}
		else if ((*itElement)->tissueType == 2 ){ //Element is on the linker Zone,
			if (BaseLinkerZoneParametersOnPeripodialness){
				//I will weight the values:
				double currPeripodialE 	= fractions[(*itElement)->Id] * PeripodialElasticity * (1 + noise1/100.0);
				double currEApical 		= fractions[(*itElement)->Id] * EApical * (1 + noise1/100.0);
				double currEBasal 		= fractions[(*itElement)->Id] * EBasal * (1 + noise1/100.0);
				double currEMid 		= fractions[(*itElement)->Id] * EMid * (1 + noise1/100.0);
				double currPoisson      = poisson*(1 + noise2/100);
				double periWeight 		= (*itElement)->getPeripodialness();
				double colWeight = (*itElement)->getColumnarness();
				currEApical = colWeight * currEApical + periWeight * currPeripodialE;
				currEBasal  = colWeight * currEBasal  + periWeight * currPeripodialE;
				currEMid    = colWeight * currEMid    + periWeight * currPeripodialE;
				(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currPoisson);
				double currViscApical  = colWeight * discProperApicalViscosity  + periWeight * peripodialApicalViscosity;
				double currViscBasal   = colWeight * discProperBasalViscosity   + periWeight * peripodialBasalViscosity;
				(*itElement)->setViscosity(currViscApical,currViscBasal); //There is no midline in the linker zone.
			}
			else{
				//I have the inputs provided:
				double currEApical 	= fractions[(*itElement)->Id] * LinkerZoneApicalElasticity*(1 + noise1/100.0);
				double currEBasal	= fractions[(*itElement)->Id] * LinkerZoneBasalYoungsModulus*(1 + noise1/100.0);
				double currEMid		= fractions[(*itElement)->Id] * 0.5 * (LinkerZoneApicalElasticity+LinkerZoneBasalYoungsModulus)*(1 + noise1/100.0);
				double currPoisson  = poisson*(1 + noise2/100);
				if((*itElement)->isECMMimicing){
					//the Poisson ratio is zero so that the ECM layer will not thin!
					currPoisson = 0;
				}
				(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currPoisson);
				(*itElement)->setViscosity(linkerZoneApicalViscosity,linkerZoneBasalViscosity,linkerZoneMidlineViscosity);
			}
		}
	}
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		if ((*itNode)->tissueType == 2 ){ //linker, the viscosity can be averaged, or based on indivudual set of parameters:
			if (BaseLinkerZoneParametersOnPeripodialness){
				//take average of the two:
				double currExternalViscApical = 0.5*(externalViscosityDPApical + externalViscosityPMApical);
				double currExternalViscBasal  = 0.5*(externalViscosityDPBasal  + externalViscosityPMBasal);
				(*itNode)->setExternalViscosity(currExternalViscApical,currExternalViscBasal, extendExternalViscosityToInnerTissue);
				//(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3);
			}
			else{
				(*itNode)->setExternalViscosity(externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3, extendExternalViscosityToInnerTissue);
				//(*itNode)->setExternalViscosity(externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3, externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3);
			}
		}else if ((*itNode)->tissueType == 0){ //disc proper
			(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, extendExternalViscosityToInnerTissue);
			//(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3);
		}else if ((*itNode)->tissueType == 1){
			(*itNode)->setExternalViscosity(externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3, extendExternalViscosityToInnerTissue);
		}
	}
	delete[] fractions;
}

void Simulation::checkForZeroExternalViscosity(){
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int i=0; i<3; ++i){
			if ((*itNode)->externalViscosity[i] > 0){
				zeroExternalViscosity[i] = false;
			}
		}
	}
    if (zeroExternalViscosity[0] || zeroExternalViscosity[1] || zeroExternalViscosity[2]){
    	//At least one of the dimensions have zero viscosity
    	if(nNodes < 3) {
    		cerr<<"Cannot run zero viscosity simulation with less than 3 nodes, run simulation at your own risk!"<<endl;
    	}
        //fix x,y,z for 0
        /*Nodes[0]->FixedPos[0] = true;
        Nodes[0]->FixedPos[1] = true;
        Nodes[0]->FixedPos[2] = true;
        //fix x and z for 1
        Nodes[1]->FixedPos[0] = true;
        Nodes[1]->FixedPos[2] = true;
        //fix z for 2
        Nodes[2]->FixedPos[2] = true;

        //Nodes[4]->FixedPos[0] = true;
        //Nodes[4]->FixedPos[1] = true;
        //Nodes[4]->FixedPos[2] = true;

        */
        if (!symmetricY){
        	//TO DO!! In the node fixing options, there should be at least 2 different axes fixed.
        	//If there are any node fixes, on the surfaces in z, then I do not need to do fix z:
        	bool sufficientXFix = false;
        	bool sufficientYFix = false;
        	bool sufficientZFix = false;
        	for (int a=0; a<3; ++a){
        		//if zeroExternalViscosity[0] is false, it means there is already  sufficient x fix
        		//OR, if I spotted sufficientXFix in previous nodes, it persists,
        		//OR, there is circumferential fix on x dimension
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0];
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || CircumferentialNodeFix[a][1];
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
        	}
			//if there are no z fixe on the circumference, check the whole surfaces:
        	if (!sufficientXFix){
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
        	if (!sufficientYFix){
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || BasalNodeFix[1] ;
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || ApicalNodeFix[1] ;
        	}
        	if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}

        	if (!sufficientXFix){
        		Nodes[ventralTipIndex]->FixedPos[0] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[ventralTipIndex]->FixedPos[1] = true;
        	}
        	if (!sufficientZFix){
        		Nodes[ventralTipIndex]->FixedPos[2] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[dorsalTipIndex]->FixedPos[1] = true;
        	}
			if (!sufficientZFix){
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
			if (!sufficientZFix){
				//if there is symmetricity, then the mid-line nodes will be fixed in y, and I do not need to fix the third node.
				// in fact, fixing the position of the third node in z will cause problems.
				if (dorsalTipIndex != 1 && 	ventralTipIndex!= 1 ){
					Nodes[1]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 2 && ventralTipIndex!= 2 ){
					Nodes[2]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 3 && ventralTipIndex!= 3 ){
					Nodes[3]->FixedPos[2] = true;
				}
			}
        }
        else{
        	//There is symmetricity, then there is y fix, if there are any node fixes, on the surfaces in z, then I do not need to do fix z:
        	bool sufficientZFix = false;
        	for (int a=0; a<3; ++a){
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
			}
        	//if there are no z fixe on the circumferenece, check the whole surfaces:
			if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}
        	bool sufficientXFix = false;
			for (int a=0; a<3; ++a){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0] ;
			}
			//if there are no x fixes on the circumferenece, check the whole surfaces:
			if (!sufficientXFix){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
			//there is no circumference or surface fixing of suffieint nature, then I should fix the nodes:
			if (!sufficientXFix){
              	Nodes[ventralTipIndex]->FixedPos[0] = true;
			}
			Nodes[ventralTipIndex]->FixedPos[1] = true;
			if (!sufficientZFix){
				Nodes[ventralTipIndex]->FixedPos[2] = true;
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
			//cout<<"node fixing sufficiency, sufficientXFix: "<<sufficientXFix<<" sufficientZFix: "<<sufficientZFix<<endl;
/*
				vector<ShapeBase*>::iterator itElement;
				bool foundElement = false;
				int apicalId = 0;
				//find the apical node corresponding to the basal ventral tip:
				for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
					bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(ventralTipIndex);
					if (IsBasalOwner){
						foundElement = true;
						apicalId = (*itElement)->getCorrecpondingApical(ventralTipIndex); //have the next node
						break;
					}
				}
				if (foundElement){
					Nodes[apicalId]->FixedPos[2] = true;
				}
				else{
					cerr<<"Cannot run zero viscosity simulation, could not find the apical node corresponding to the basal ventral tip, run simulation at your own risk!"<<endl;
				}
*/
        }
    }
}

void Simulation::addCurvatureToColumnar(double h){
	//find the tips:
	double l1 = -1000, l2 = 1000, l3 = -1000;
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->Position[0] > l1 ){
			l1 = (*itNode)->Position[0]; //displacement from origin to positive x tip
		}
		if ((*itNode)->Position[0] < l2 ){
			l2 = (*itNode)->Position[0]; //displacement from origin to negative x tip
		}
		if ((*itNode)->Position[1] > l3 ){
			l3 = (*itNode)->Position[1];	//displacement from origin to positive y tip
		}
	}
	if (symmetricX){
		l2 = (-1.0)*l1;	//if there is symmetricity in x, then the negative tip will be at zero. This is not correct and the symmetricity should mean the negative displacement will be -l1
	}
	l2 *= -1.0;	//converting the displacement to distance, this is the only negative tip, the rest are already positive
	cout<<"l1: "<<l1 <<" l2: "<<l2<<" l3: "<<l3<<endl;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		double x = (*itNode)->Position[0];
		double y = (*itNode)->Position[1];
		double z = (*itNode)->Position[2];
		double a = l1;
		double c = l3;
		if ((*itNode)->Position[0] < 0){
			a = l2;
		}
		double value = (1 - x*x/a/a - y*y/c/c);
		if (value>0){
			double offset = pow((1 - x*x/a/a - y*y/c/c)*h*h,0.5);
			if (h<0){
				offset *= (-1.0);
			}
			if (thereIsPeripodialMembrane && (*itNode)->tissueType != 0){ //node is not columnar
				//there is peripodial membrane, any node above the mid-line of the lumen should be curved in the opposite direction of the columnar layer:
				//But I have a direction choice, if the columnar is curving down, the peripodial shoul curve up
				//If the columnar is curing up, the peripodial should still curve up!

				double heightThreshold = TissueHeight + lumenHeight/2.0;
				if ((*itNode)->Position[2] > heightThreshold) {
					if (offset > 0){
						offset *= (-1.0);
					}
				}
			}
			(*itNode)->Position[2] -= offset;
			if((*itNode)->Id ==8){
				cout<<" a: "<<a<<" c: "<<c<<" x "<<x<<" y "<< y<<" z "<<z<<" value: "<<value<<" offset: "<<offset<<endl;
				cout<<"Node[8] position: "<<(*itNode)->Position[0]<<" "<<(*itNode)->Position[1]<<" "<<(*itNode)->Position[2]<<endl;
			}
		}
	}

	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
		(*itElement)->updateReferencePositionsToCurentShape();
	}
}

void Simulation::fixNode0InPosition(double x, double y, double z){
	double dx = Nodes[0]->Position[0]-x;
	double dy = Nodes[0]->Position[1]-y;
	double dz = Nodes[0]->Position[2]-z;
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->Position[0] -=dx;
		(*itNode)->Position[1] -=dy;
		(*itNode)->Position[2] -=dz;
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
}

void Simulation::manualPerturbationToInitialSetup(bool deform, bool rotate){
	if(timestep==0){
		 for (vector<Node*>::iterator  itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			 if ((*itNode)->Id <3 ){
				//fixX((*itNode),0);
				//fixY((*itNode),0);
			 }
			/*if ((*itNode)->Id != 3 && (*itNode)->Id != 10){
				//fixAllD((*itNode),0);
			}
			else{
				//fixZ((*itNode),0);
				//fixY((*itNode),0);
				(*itNode)->Position[0] -= 5;
			}*/
			fixAllD((*itNode),0);
		}
		//laserAblateTissueType(1);
		//laserAblate(0.0, 0.0, 5.0);
		//deform = true;
        double scaleX = 2.0;
        double scaleY = 1.0;
        double scaleZ = 1.0;

        double PI = 3.14159265359;
        double tetX = 0 *PI/180.0;
        double tetY = 45 *PI/180.0;
        double tetZ = 0 *PI/180.0;
    	int axisToRotateOn = 1; //0: rotate around x axis, 1: around y-axis, and 2: around z-axis.
        if(deform){
        	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
        		(*itNode)->Position[0] *=scaleX;
        		(*itNode)->Position[1] *=scaleY;
                (*itNode)->Position[2] *=scaleZ;
            }
        }
        if (rotate){
        	cerr<<"Rotating system"<<endl;
            double R[3][3] =  {{1,0,0},{0,1,0},{0,0,1}};
            if (axisToRotateOn == 0){
            	double c = cos(tetX);
            	double s = sin(tetX);
            	//Rx = {{1,0,0},{0,c,-s},{0,s,c}};
            	R[1][1] = c;
				R[1][2] = -s;
				R[2][1] = s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 1){
            	double c = cos(tetY);
            	double s = sin(tetY);
            	//Ry = {{c,0,s},{0,1,0},{-s,0,c}};
            	R[0][0] = c;
				R[0][2] = s;
				R[2][0] = -s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 2){
            	double c = cos(tetZ);
            	double s = sin(tetZ);
            	//Rz = {{c,-s,0},{s,c,0},{0,0,1}};
            	R[0][0] = c;
				R[0][1] = -s;
				R[1][0] = s;
				R[1][1] = c;
            }
        	for (vector<Node*>::iterator  itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
                double x = (*itNode)->Position[0]*R[0][0] + (*itNode)->Position[1]*R[0][1] + (*itNode)->Position[2]*R[0][2];
                double y = (*itNode)->Position[0]*R[1][0] + (*itNode)->Position[1]*R[1][1] + (*itNode)->Position[2]*R[1][2];
                double z = (*itNode)->Position[0]*R[2][0] + (*itNode)->Position[1]*R[2][1] + (*itNode)->Position[2]*R[2][2];
                (*itNode)->Position[0]=x;
                (*itNode)->Position[1]=y;
                (*itNode)->Position[2]=z;
            }
        }
    	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    		(*itElement)->updatePositions(Nodes);
        }
    }
}

void Simulation::updateGrowthRotationMatrices(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
            (*itElement)->CalculateGrowthRotationByF();
        }
    }
}

void Simulation::calculateDiscretisationLayers(double &hColumnar, int& LumenHeightDiscretisationLayers, double &hLumen, double &peripodialHeight, int& peripodialHeightDiscretisationLayers, double& hPeripodial){
	//The average columnar layer element heigth would be TissueHeight divided by the number of element layers:
	hColumnar = TissueHeight/TissueHeightDiscretisationLayers;
	//Calculating the height of the lumen:
	lumenHeight = TissueHeight*lumenHeightScale;
	//calculating how many layer I need for this lumen:
	LumenHeightDiscretisationLayers = ceil(lumenHeight/hColumnar);
	//calculating the actual element height needed to produce the total height with the calculated layr number:
	hLumen = lumenHeight/LumenHeightDiscretisationLayers;
	//calculating the height of the peripodial membrane:
	peripodialHeight = PeripodialThicnessScale*TissueHeight;
	//calculating how many layers I need for this peripodial membrane:
	peripodialHeightDiscretisationLayers = ceil(peripodialHeight/hColumnar);
	//calculating the actual element height needed to produce the total height with the calculated layr number:
	hPeripodial = peripodialHeight/peripodialHeightDiscretisationLayers;
	//cout<<"Tissue height: "<<TissueHeight<<" columnar element height: "<<hColumnar<<" TissueHeightDiscretisationLayers: "<<TissueHeightDiscretisationLayers<<endl;
	//cout<<"Lumen height:  "<<lumenHeight<<" height of one lumen element: "<<hLumen<<" LumenHeightDiscretisationLayers: "<<LumenHeightDiscretisationLayers<<endl;
	//cout<<"Peripodial height:  "<<peripodialHeight<<" height of one peripodial element: "<<hPeripodial<<" LumenHeightDiscretisationLayers: "<<peripodialHeightDiscretisationLayers<<endl;
}

void Simulation::fillColumnarBasedNodeList(vector< vector<int> > &ColumnarBasedNodeArray, vector <int> &ColumnarCircumferencialNodeList){
	int nCircumference = ColumnarCircumferencialNodeList.size();
	for (int i =0; i<nCircumference; ++i){
    	int currNodeId = ColumnarBasedNodeArray[i][0];
		bool foundElement = true;
		bool finishedTissueThicness  = false;
		while(!finishedTissueThicness && foundElement){ //while the node is not apical, and I could find the next element
			foundElement = false;
			vector<ShapeBase*>::iterator itElement;
			for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					foundElement = true;
					break;
				}
			}
			//found the current node as basal node of an element
			//get the corresponding apical node, record in the list:
			currNodeId = (*itElement)->getCorrecpondingApical(currNodeId); //have the next node
			ColumnarBasedNodeArray[i].push_back(currNodeId);
			if (Nodes[currNodeId]->tissuePlacement == 1) {
				//Added the apical node now, I can stop:
				finishedTissueThicness = true;
			}
		}
    }
}

void Simulation::addNodesForPeripodialOnColumnarCircumference (vector< vector<int> > &ColumnarBasedNodeArray, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial ){
	//add nodes for the lumen and peripodial membrane on top of the columnar nodes:
	int nCircumference = ColumnarBasedNodeArray.size();
	for (int i=0; i<nCircumference; ++i){
		int nodeId0 = ColumnarBasedNodeArray[i][TissueHeightDiscretisationLayers];
		double* pos;
		pos = new double[3];
		for (int j=0; j<Nodes[nodeId0]->nDim; j++){
			pos[j] = Nodes[nodeId0]->Position[j];
		}
		//adding nodes for lumen:
		for (int j=1; j<LumenHeightDiscretisationLayers+1; ++j){
			pos[2] += hLumen;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			ColumnarBasedNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be apical, at it is the bottom of the peripodial membrane, looking into the lumen, change its placement:
		Nodes[Nodes.size()-1]->tissuePlacement = 1; //made tissueplacement apical
		//adding nodes for the peripodial membrane:
		for (int j=1; j<peripodialHeightDiscretisationLayers+1; ++j){
			pos[2] += hPeripodial;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			ColumnarBasedNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
		Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
		delete pos;
	}
}
void Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double* pos, double sideThickness){
	cout<<"nodeIDs: "<<nodeId0<<" "<<nodeId1<<" "<<nodeId2<<" sideThickness: "<<sideThickness<<endl;
	double* vec1;
	vec1 = new double[3];
	if (symmetricY && nodeId1 == -100){
		vec1[0] = -1.0;
		vec1[1] =  0.0;
		vec1[2] =  0.0;
	}
	else if (symmetricY && nodeId2 == -100){
		vec1[0] =  1.0;
		vec1[1] =  0.0;
		vec1[2] =  0.0;
	}
	else{
		double* vec2;
		vec2 = new double[3];
		for (int j=0; j<Nodes[nodeId0]->nDim; j++){
			vec1[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId1]->Position[j]);
			vec2[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId2]->Position[j]);
		}
		Elements[0]->normaliseVector3D(vec1);
		Elements[0]->normaliseVector3D(vec2);
		vec1[0] += vec2[0];
		vec1[1] += vec2[1];
		vec1[2] += vec2[2];
		double vec1Mag2 = vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2];
		if (vec1Mag2 < 1E-5){
			//the two nodes are linear, the resulting vector is of zero length.
			//I will rotate one of the vectors 90 degrees and it will be pointing in the correct orientation, the actual direction can be fixed in the below loop:
			//this vector is already normalised, I am skipping the normalisation step
			vec1[0] = -1.0*vec2[1];
			vec1[1] = vec2[0];
		}
		else{
			Elements[0]->normaliseVector3D(vec1);
		}
		//now I have the vector to point out from the base node 0. BUT, this will point outwards only if the tissue curvature is convex at all points
		//I need to check if it actually is pointing out, as the experimentally driven tissues can be concave at points.
		//the cross produc of vec[2] and the vector to the cell centre should have the opposite sign with the corss product of my orientation vector and vector 2.
		double* vecCentre;
		vecCentre = new double[3];
		vecCentre[0] =  SystemCentre[0] - Nodes[nodeId0]->Position[0];
		vecCentre[1] =  SystemCentre[1] - Nodes[nodeId0]->Position[1];
		vecCentre[2] =  SystemCentre[2] - Nodes[nodeId0]->Position[2];
		Elements[0]->normaliseVector3D(vecCentre);
		double* cross1;
		cross1 = new double[3];
		Elements[0]->crossProduct3D(vec2,vecCentre,cross1);
		Elements[0]->normaliseVector3D(cross1);
		double* cross2;
		cross2 = new double[3];
		Elements[0]->crossProduct3D(vec2,vec1,cross2);
		Elements[0]->normaliseVector3D(cross2);
		double dotp = Elements[0]->dotProduct3D(cross1,cross2);
		if (dotp >0 ){
			//the vectors are pointing to the same direction! Need to rotate vec1 180 degrees:
			vec1[0] *= -1.0;
			vec1[1] *= -1.0;
			vec1[2] *= -1.0;
		}
		delete[] vecCentre;
		delete[] cross1;
		delete[] cross2;
		delete[] vec2;
	}
	pos[0] = Nodes[nodeId0]->Position[0] + vec1[0]*sideThickness;
	pos[1] = Nodes[nodeId0]->Position[1] + vec1[1]*sideThickness;
	pos[2] = Nodes[nodeId0]->Position[2];
	delete[] vec1;

}


void Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double* pos, double sideThickness){
	double* vec0;
	vec0 = new double[3];
	double midpoint[3] = {0,0,0};
	for (int j=0; j<Nodes[nodeId0]->nDim; j++){
		vec0[j] = (Nodes[nodeId1]->Position[j] - Nodes[nodeId0]->Position[j])/2.0;
		midpoint[j] = Nodes[nodeId0]->Position[j] + vec0[j];
	}
	Elements[0]->normaliseVector3D(vec0);
	//The list is sorted counter-cock-wise, to point out, I will rotate normalised vector v0 -90 degrees on z axis:
	// (x,y,z) -> (y,-x,z);
	// then I will add this vector to the calculated mid point to gt the new node's position.
	pos[0] = midpoint[0] + vec0[1]*sideThickness;
	pos[1] = midpoint[1] - vec0[0]*sideThickness;
	pos[2] = midpoint[2];
	delete[] vec0;
}

void Simulation::addNodesForPeripodialOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray, double hColumnar, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial ){
	double peripodialSideConnectionThickness =PeripodialLateralThicnessScale*TissueHeight; //in microns
	double avrSide=0.0, dummy =0.0;
	getAverageSideLength(dummy,avrSide);	//first term will get you the average side length of the peripodial membrane elements, second is the columnar elements
	if (avrSide/peripodialSideConnectionThickness > 5 || avrSide/peripodialSideConnectionThickness< 0.2 ){
		cerr<<"WARNING, the lateral connection thickness between the peripodial membrane and the columnar layer is too different than average element size (more than 5 fold diference)"<<endl;
	}
	//Now I need the average side of an element, to add new nodes accordingly:
	int nCircumference = ColumnarBasedNodeArray.size();

	for (int i=0; i<nCircumference; ++i){
		//cout<<"at node in the list: "<<i<<endl;
		//adding 2 point based node:
/*		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		if( i == nCircumference - 1){
			nodeId1 = ColumnarBasedNodeArray[0][0];
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		double* pos = new double[3];
		calculateNewNodePosForPeripodialNodeAddition(nodeId0, nodeId1, pos, peripodialSideConnectionThickness);
		//cout<<" calculated pos : "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<endl;
		//Adding the array of new nodes:
*/
		//adding 3 point based node:
		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		int nodeId2;
		if( i == nCircumference - 1){
			if (symmetricY){
				//the node is the end tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to +x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId1 = -100;
			}
			else{
				nodeId1 = ColumnarBasedNodeArray[0][0];
			}
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		if ( i == 0 ){
			if (symmetricY){
				//the node is the beginning tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to -x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId2 = -100;
			}
			else{
				nodeId2 = ColumnarBasedNodeArray[nCircumference-1][0];
			}
		}
		else{
			nodeId2 = ColumnarBasedNodeArray[i-1][0];
		}
		double* pos = new double[3];
		calculateNewNodePosForPeripodialNodeAddition(nodeId0, nodeId1, nodeId2, pos, peripodialSideConnectionThickness);

		//cout<<" calculated pos : "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<endl;
		//Adding the array of new nodes:

		//adding the base:
		int newNodeId = Nodes.size();
		Node* tmp_nd = new Node(newNodeId, 3, pos, 0, 2); //Tissue placement basal (0), tissue type is linker zone (2)
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		OuterNodeArray[i].push_back(newNodeId);
		//adding the nodes for the columnar layer:
		for (int j=1; j<TissueHeightDiscretisationLayers+1; ++j){
			pos[2] += hColumnar;
			//cout<<" pos for columnar aligned new node: "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<" hColumnar: "<<hColumnar<<endl;

			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//adding nodes for lumen:
		for (int j=1; j<LumenHeightDiscretisationLayers+1; ++j){
			pos[2] += hLumen;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//adding nodes for the peripodial membrane:
		for (int j=1; j<peripodialHeightDiscretisationLayers+1; ++j){
			pos[2] += hPeripodial;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
		Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
		delete[] pos;
    }
}

void Simulation::addLateralPeripodialElements(int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray){
    // I need to add the elements:
    // Two elements are added for each element on the side:
    // First triangle base will be New node 0, oldNode1, oldnode0, top of the element will be read from the node id stacks
    // Second triangle base will be: New node 0, new node 1, old node 1
	int totalLayers = TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers+peripodialHeightDiscretisationLayers;
	double peripodialnessFractionStep = 1.0 / (double) (LumenHeightDiscretisationLayers+1.0); //this is the step increment in periopdiallness weight at each element of lumen side. If there is 1 layer, it should be 50% peripodial, 50% columnar etc, if there are two, it should be 33%, 2x33% etc.
	int nCircumference = ColumnarBasedNodeArray.size();
	int nLoop = nCircumference;
	if (symmetricY){
		//in symmetric setup, the last node is not connected to the first node, we simply ignore that step
		nLoop--;
	}
	for (int i=0;i<nLoop; ++i){
    	for (int j=0;j<totalLayers; ++j){
    		//calculating the fraction of peripodialness these two elements should have:
    		//counting which layer in the lumen they are:
    		double peripodialWeight = 1.0; //
    		if (j < TissueHeightDiscretisationLayers) {
    			//the currently added elements are still aligned with columnar, below the lumen.
    			peripodialWeight = 0.0;
    		}
    		else if (j<TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers){
    			//the currently added elements are at the lumen, below peripodial membrane:
    			peripodialWeight = (j-TissueHeightDiscretisationLayers + 1) * peripodialnessFractionStep;
    		}
    		// if (j>=TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers), then the added elements are above the lumen, should behave like peripodial, default was 1.0, no change.
    		//preparing the node lists:
    		int indiceTri0Corner0 = i;
			int indiceTri0Corner1 = i+1;
			int indiceTri0Corner2 = i;
			int indiceTri1Corner0 = indiceTri0Corner0;
			int indiceTri1Corner1 = i+1;
			int indiceTri1Corner2 = indiceTri0Corner1;
			if (indiceTri0Corner1 == nCircumference){
				indiceTri0Corner1 = 0;
				indiceTri1Corner2 = 0;
				indiceTri1Corner1 = 0;
			}
			int* NodeIds = new int[6];
			//adding the first element:
			NodeIds[0] = OuterNodeArray[indiceTri0Corner0][j];
			NodeIds[1] = ColumnarBasedNodeArray[indiceTri0Corner1][j];
			NodeIds[2] = ColumnarBasedNodeArray[indiceTri0Corner2][j];
			NodeIds[3] = OuterNodeArray[indiceTri0Corner0][j+1];
			NodeIds[4] = ColumnarBasedNodeArray[indiceTri0Corner1][j+1];
			NodeIds[5] = ColumnarBasedNodeArray[indiceTri0Corner2][j+1];
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
			//adding the second element:
			NodeIds[0] = OuterNodeArray[indiceTri1Corner0][j];
			NodeIds[1] = OuterNodeArray[indiceTri1Corner1][j];
			NodeIds[2] = ColumnarBasedNodeArray[indiceTri1Corner2][j];
			NodeIds[3] = OuterNodeArray[indiceTri1Corner0][j+1];
			NodeIds[4] = OuterNodeArray[indiceTri1Corner1][j+1];
			NodeIds[5] = ColumnarBasedNodeArray[indiceTri1Corner2][j+1];
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
    	}
    }
}

void Simulation::addNodesForPeripodialOnCap(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &PeripodialCapNodeArray, int TissueHeightDiscretisationLayers, int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, double hPeripodial){
    int ncurrNodes = nNodes;
    int nCircumference = ColumnarBasedNodeArray.size();
    for (int i = 0; i<ncurrNodes; ++i){
		if (Nodes[i]->tissuePlacement == 1 && Nodes[i]->tissueType == 0){ //Node is apical and on the columnar layer
			if (Nodes[i]->atCircumference){
				//copy form previous list,
				int lumenCapLayer =  TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers;
				for (int j=0; j<nCircumference; ++j){
					if(Nodes[i]->Id == ColumnarBasedNodeArray[j][TissueHeightDiscretisationLayers]){
						//found the current node on the existing list
						PeripodialCapNodeArray[j].push_back(ColumnarBasedNodeArray[j][TissueHeightDiscretisationLayers]);
						int n = ColumnarBasedNodeArray[j].size();
						for (int k=lumenCapLayer; k<n; ++k ){
							//copy all the previously added nodes to this list
							//the 0th layer will be the apical surface of the columnar layer, it will be used to read element structure
							//new elements will be added via the nodes on following layers
							PeripodialCapNodeArray[j].push_back(ColumnarBasedNodeArray[j][k]);
						}
						break;
					}
				}
			}
			else{
				//push back a new vector
				int n = PeripodialCapNodeArray.size();
				PeripodialCapNodeArray.push_back(vector<int>(0));
				PeripodialCapNodeArray[n].push_back(Nodes[i]->Id);
				//now I need to add the new nodes:
				//adding nodes for the peripodial membrane:
				double* pos;
				pos = new double[3];
				pos[0] = Nodes[i]->Position[0];
				pos[1] = Nodes[i]->Position[1];
				pos[2] = Nodes[i]->Position[2]+lumenHeight;
				for (int j=0; j<peripodialHeightDiscretisationLayers+1; ++j){
					int newNodeId = Nodes.size();
					int tissuePlacement = 2; //defalut is midlayer
					if (j==0){tissuePlacement = 1;} //The first node should be apical, as it is looking into the lumen, the basal node is corrected outside loop
					Node* tmp_nd = new Node(newNodeId, 3, pos, tissuePlacement, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
					Nodes.push_back(tmp_nd);
					nNodes = Nodes.size();
					PeripodialCapNodeArray[n].push_back(newNodeId);
					pos[2] += hPeripodial;
				}
				//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
				Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
				delete[] pos;
			}
		}
	}
 }

void Simulation::constructTriangleCornerListOfApicalSurface( vector< vector<int> > &TriangleList){
	int nTri = 0;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		vector <int> ApicalTriangles;
		(*itElement)->getApicalTriangles(ApicalTriangles);
		int nList = ApicalTriangles.size();
		for (int k=0; k<nList-2; k+=3){
			if (Nodes[ApicalTriangles[k]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+1]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+2]]->tissuePlacement == 1){
				TriangleList.push_back(vector<int>(3));
				TriangleList[nTri][0] = (Nodes[ApicalTriangles[k]]->Id);
				TriangleList[nTri][1] = (Nodes[ApicalTriangles[k+1]]->Id);
				TriangleList[nTri][2] = (Nodes[ApicalTriangles[k+2]]->Id);
				nTri++;
			}
		}
	}
}

void Simulation::addCapPeripodialElements( vector< vector<int> > &TriangleList, vector< vector<int> > &PeripodialCapNodeArray, int peripodialHeightDiscretisationLayers){
	int nTri = TriangleList.size();
	for (int i=0; i<nTri; ++i){
		int indiceTri0Corner0 =-1, indiceTri0Corner1=-1,indiceTri0Corner2=-1;
		int n = PeripodialCapNodeArray.size();
		//cout<<"size of the peripodial cap node array: "<<n<<endl;
		bool found0 = false, found1 = false, found2 = false;
		for (int j =0; j<n; ++j){
		   if (TriangleList[i][0] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner0 = j;
			   found0 = true;
		   }
		   else if (TriangleList[i][1] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner1 = j;
			   found1 = true;
		   }
		   else if (TriangleList[i][2] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner2 = j;
			   found2 = true;
		   }
		   if (found0 && found1 && found2 ){
			   break;
		   }
		}
		if(!(found0 && found1 && found2)){
			cerr<<"could not find all the corners in addCapPeripodialElements, will give error"<<endl;
			cerr.flush();
		}
		//now I have the indices of the corners, adding the elements, as many layers as necessary:
		int* NodeIds = new int[6];
		for (int j =1; j< peripodialHeightDiscretisationLayers+1; ++j){
			NodeIds[0] = PeripodialCapNodeArray[indiceTri0Corner0][j];
			NodeIds[1] = PeripodialCapNodeArray[indiceTri0Corner1][j];
			NodeIds[2] = PeripodialCapNodeArray[indiceTri0Corner2][j];
			NodeIds[3] = PeripodialCapNodeArray[indiceTri0Corner0][j+1];
			NodeIds[4] = PeripodialCapNodeArray[indiceTri0Corner1][j+1];
			NodeIds[5] = PeripodialCapNodeArray[indiceTri0Corner2][j+1];
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
		}
   }
}
void Simulation::correctCircumferentialNodeAssignment(vector< vector<int> > OuterNodeArray){
	//Now I added nodes and elements to the circumference of the tissue.
	//The Nodes that were at the circumference in the columnar layer are embedded in the tissue now,
	// while the new lateral nodes of the peripodial membrane are at the circumference.
	//I will correct this in node booleans:
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		//making all circumferential node flags of the columnar layer false
		if ( (*itNode)->tissueType == 0 && (*itNode)->atCircumference){
			(*itNode)->atCircumference = false;
		}
	}
	//I have a list of the Ids fot the OuterNodes I have added, changing their circumferential node flags to true:
	int n0 = OuterNodeArray.size();
	for (int i =0 ; i<n0; ++i){
		int n1 = OuterNodeArray[i].size();
		for (int j=0; j<n1; ++j){
			int currId = OuterNodeArray[i][j];
			//find the node with the curr id: (I know Ids and index order is the same now, but staying on the safe side for potential future changes
			for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
				//making all circumferential node flags of the columnar layer false
				if ( (*itNode)->Id == currId){
					(*itNode)->atCircumference = true;
					break;
				}
			}
		}
	}
}

bool Simulation::addCurvedPeripodialMembraneToTissue(){
	bool Success = true;
	//here I am calculating the height of the tissue and the discretisation layers used for columnar layer
	double hColumnar; //The average columnar layer element height
	int LumenHeightDiscretisationLayers;  //The number of elements that are used to discretised the lumen height
	double hLumen; //The lumen element height
	double peripodialHeight; 	//Height of the peropodial membrane
	int peripodialHeightDiscretisationLayers; //The number of elements that are used to discretised the peripodial membrane height
	double hPeripodial;  //The peripodial membrane element height
	calculateDiscretisationLayers(hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeight, peripodialHeightDiscretisationLayers, hPeripodial);
	//Now I have all the height information and the number of layers for each subsection of the tissue
	//I can start adding the nodes.
	//Initially, I am starting with the sides.
	//First I want the list of nodes at the basal circumference of the columnar layer:
	vector <int> ColumnarCircumferencialNodeList;
	Success = generateColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
	if (!Success){
		cerr<<"Error!! circumferential nodes not extracted properly"<<endl;
	}
	//Now I have the list, I need to sort in to rotate in one direction:
	calculateSystemCentre();
	//I will sort the list counter-clock-wise, while the (+)z is pointing at you
	sortColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
	//Now I will create the arrays of vectors that will keep my nodes organised:
	const int nCircumference = ColumnarCircumferencialNodeList.size();
	vector< vector<int> > ColumnarBasedNodeArray( nCircumference , vector<int>(0) );
	for (int i=0;i<nCircumference; ++i){
		ColumnarBasedNodeArray[i].push_back(ColumnarCircumferencialNodeList[i]);
	}
	//Fill the array of node IDs upwards, to cover the whole columnar layer:
	fillColumnarBasedNodeList(ColumnarBasedNodeArray, ColumnarCircumferencialNodeList);
	addNodesForPeripodialOnColumnarCircumference (ColumnarBasedNodeArray, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
	//now I will add the positions for outer surface:
	//starting from the base layer of the circumferencial node list, ending at the top of the array
	int divisionSteps = 6;
	vector< vector<int> > OuterNodeArray( nCircumference , vector<int>(0) );
	vector< vector<int> > InnerNodeArray( nCircumference , vector<int>(0) );
	//calculating angle from the centre of the circular edge to each new point:
	//I am thinking of a quarter circle, divided into divisionSteps/2 triangles. The central angle will be divided into segments of:
	double tet = (M_PI /2.0) / (divisionSteps / 2.0);
	double hMidPoint = TissueHeight + 0.5*hLumen;
	for (int i=0;i<nCircumference; ++i){
		int idBase = ColumnarBasedNodeArray[i][0];
		//int idTop = ColumnarBasedNodeArray[i][ColumnarBasedNodeArray[i].size()-1];
		//double posBase[3] = { Nodes[idBase]->Position[0],Nodes[idBase]->Position[1],Nodes[idBase]->Position[2]};
		//double posTop[3] = { Nodes[idTop]->Position[0],Nodes[idTop]->Position[1],Nodes[idTop]->Position[2]};

		//first find the vector pointing out from the base node:
		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		int nodeId2;
		if( i == nCircumference - 1){
			if (symmetricY){
				//the node is the end tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to +x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId1 = -100;
			}
			else{
				nodeId1 = ColumnarBasedNodeArray[0][0];
			}
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		if ( i == 0 ){
			if (symmetricY){
				//the node is the beginning tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to -x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId2 = -100;
			}
			else{
				nodeId2 = ColumnarBasedNodeArray[nCircumference-1][0];
			}
		}
		else{
			nodeId2 = ColumnarBasedNodeArray[i-1][0];
		}
		double* pos = new double[3];
		double* vec1;
		vec1 = new double[3];
		if (symmetricY && nodeId1 == -100){
			vec1[0] = -1.0;
			vec1[1] =  0.0;
			vec1[2] =  0.0;
		}
		else if (symmetricY && nodeId2 == -100){
			vec1[0] =  1.0;
			vec1[1] =  0.0;
			vec1[2] =  0.0;
		}
		else{
			double* vec2;
			vec2 = new double[3];
			for (int j=0; j<Nodes[nodeId0]->nDim; j++){
				vec1[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId1]->Position[j]);
				vec2[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId2]->Position[j]);
			}
			Elements[0]->normaliseVector3D(vec1);
			Elements[0]->normaliseVector3D(vec2);
			vec1[0] += vec2[0];
			vec1[1] += vec2[1];
			vec1[2] += vec2[2];
			Elements[0]->normaliseVector3D(vec1);
			//The list is sorted counter-cock-wise, to point out, I will rotate normalised vector v0 -90 degrees on z axis:
			// (x,y,z) -> (y,-x,z);
			// then I will add this vector to the calculated mid point to gt the new node's position.
			delete[] vec2;
		}
		delete[] pos;

		//first outer point will have two radia, its movement in z direction, and on the direction of the vector pointing out
		// from the base point.
		//similarly, first inner point will have two radia (r3, r4);
		for (int j=1; j<divisionSteps; ++j ){
			double c = cos(j*tet);
			double s = sin(j*tet);
			double r1 = hMidPoint *(1 -c);
			double r2 = hMidPoint * s;
			//outer node:
			drawingPointsX.push_back(Nodes[idBase]->Position[0] + r2 * vec1[0]);
			drawingPointsY.push_back(Nodes[idBase]->Position[1] + r2 * vec1[1]);
			drawingPointsZ.push_back(Nodes[idBase]->Position[2] + r1);
			cout<<"hMidPoint: "<<hMidPoint<<endl;
		}
		delete[] vec1;
	}



	for (int i=0;i<nCircumference; ++i){
		int nColumnarBasedNodeArray = ColumnarBasedNodeArray[i].size();
		for (int j =0; j< nColumnarBasedNodeArray; ++j){
			int id = ColumnarBasedNodeArray[i][j];
			drawingPointsX.push_back( Nodes[id]->Position[0]);
			drawingPointsY.push_back( Nodes[id]->Position[1]);
			drawingPointsZ.push_back( Nodes[id]->Position[2]);
		}
	}

	return Success;
}

bool Simulation::checkIfThereIsPeripodialMembrane(){
	bool Success = true;
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 1){
			thereIsPeripodialMembrane = true;
			//I have not added the peripodial membrane as yet,
			//if there is one, it came from the input mesh
			//Then I want to set up the circumferencial input properly.
			setLinkerCircumference();
			break;
		}
	}
	return Success;
}

bool Simulation::addStraightPeripodialMembraneToTissue(){
    //cout<<"adding peripodial membrane from scratch" <<endl;
	bool Success = true;
    //here I am calculating the height of the tissue and the discretisation layers used for columnar layer
    double hColumnar; //The average columnar layer element height
    int LumenHeightDiscretisationLayers;  //The number of elements that are used to discretised the lumen height
    double hLumen; //The lumen element height
    double peripodialHeight; 	//Height of the peropodial membrane
    int peripodialHeightDiscretisationLayers; //The number of elements that are used to discretised the peripodial membrane height
    double hPeripodial;  //The peripodial membrane element height
    calculateDiscretisationLayers(hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeight, peripodialHeightDiscretisationLayers, hPeripodial);
    //Now I have all the height information and the number of layers for each subsection of the tissue
    //I can start adding the nodes.
    //Initially, I am starting with the sides.
    //First I want the list of nodes at the basal circumference of the columnar layer:
    vector <int> ColumnarCircumferencialNodeList;
    Success = generateColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    if (!Success){
        cerr<<"Error!! circumferential nodes not extracted properly"<<endl;
    }
    //Now I have the list, I need to sort in to rotate in one direction:
    calculateSystemCentre();
    //I will sort the list counter-clock-wise, while the (+)z is pointing at you
    sortColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    //Now I will create the arrays of vectors that will keep my nodes organised:
    const int nCircumference = ColumnarCircumferencialNodeList.size();
    vector< vector<int> > ColumnarBasedNodeArray( nCircumference , vector<int>(0) );
    for (int i=0;i<nCircumference; ++i){
       	ColumnarBasedNodeArray[i].push_back(ColumnarCircumferencialNodeList[i]);
    }
    //Fill the array of node IDs upwards, to cover the whole columnar layer:
    fillColumnarBasedNodeList(ColumnarBasedNodeArray, ColumnarCircumferencialNodeList);
    addNodesForPeripodialOnColumnarCircumference (ColumnarBasedNodeArray, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
    //Now add nodes for the outer layer:
    vector< vector<int> > OuterNodeArray( nCircumference , vector<int>(0) );
    addNodesForPeripodialOnOuterCircumference (ColumnarBasedNodeArray, OuterNodeArray, hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
    //Now I need to add the elements:
    addLateralPeripodialElements(LumenHeightDiscretisationLayers,peripodialHeightDiscretisationLayers, ColumnarBasedNodeArray, OuterNodeArray);
    //now adding the nodes of the central region:
    // I am initialising the vector as the size of circumference node list, for any node at the circumference, I will copy from previuos lists,
    // for any node that is not at the circumference, I will add a new vector of integers to this list
    vector< vector<int> > PeripodialCapNodeArray( nCircumference , vector<int>(0) );
    addNodesForPeripodialOnCap(ColumnarBasedNodeArray, PeripodialCapNodeArray, TissueHeightDiscretisationLayers, LumenHeightDiscretisationLayers, peripodialHeightDiscretisationLayers, hPeripodial);
    //Now I will generate the triangle list from the apical surface of the columnar layer:
    vector< vector<int> > TriangleList( 0 , vector<int>(0) ); //The IDs of 3 nodes that are supposed to form a triangle
    constructTriangleCornerListOfApicalSurface(TriangleList);
    //Now I have the list of nodes forming triangles
    //For each triangle, find the map of indices in the 0th layer of PeripodialCapNodeArray, and add elements according to the node id stacks:
    addCapPeripodialElements( TriangleList, PeripodialCapNodeArray, peripodialHeightDiscretisationLayers);
    //Peripodial membrane is added, now correct position assignemnts at the circumference, and related node fixing:
    correctCircumferentialNodeAssignment(OuterNodeArray);
    return Success;
}

void Simulation::checkForExperimentalSetupsBeforeIteration(){
	if (stretcherAttached){
		moveClampedNodesForStretcher();
	}
}

void Simulation::checkForExperimentalSetupsWithinIteration(){
}

void Simulation::checkForExperimentalSetupsAfterIteration(){
	if (stretcherAttached){
		recordForcesOnClampBorders();
	}
}

void Simulation::checkECMSoftening(){
	//I have not carried out any softening as yet
	//I will if the time is after 32 hr
	if (currSimTimeSec >=softeningBeginTimeInSec && currSimTimeSec <softeningEndTimeInSec){
		if (softenedECM == false){
			softenedECM = true;
			//this is the first time step I am reducing the ECM viscosity.
			//I need to calculate rates first:
			for (int aa=0; aa<nNodes; ++aa){
				double timeDifferenceInHours = (softeningEndTimeInSec - softeningBeginTimeInSec)/3600;
				for (int i=0; i<3; ++i){
					Nodes[aa]->ECMViscosityReductionPerHour[i] = Nodes[aa]->externalViscosity[i]*(1-ECMSofteningFraction)/timeDifferenceInHours;
				}
			}
		}
		for (int aa=0; aa<nNodes; ++aa){
			if((Nodes[aa]->tissuePlacement == 0 && softenBasalECM ) || (Nodes[aa]->tissuePlacement == 1 && softenApicalECM )){
				double relativeX = (Nodes[aa]->Position[0] - boundingBox[0][0])/boundingBoxSize[0];
				for (int ECMReductionRangeCounter =0; ECMReductionRangeCounter<numberOfSoftenedRanges; ++ECMReductionRangeCounter){
					if (relativeX> ECMSofteningXRangeMins[ECMReductionRangeCounter] && relativeX<ECMSofteningXRangeMaxs[ECMReductionRangeCounter]){
						Nodes[aa]->externalViscosity[0] -= Nodes[aa]->ECMViscosityReductionPerHour[0]/3600*dt;
						Nodes[aa]->externalViscosity[1] -= Nodes[aa]->ECMViscosityReductionPerHour[1]/3600*dt;
						Nodes[aa]->externalViscosity[2] -= Nodes[aa]->ECMViscosityReductionPerHour[2]/3600*dt;
						Nodes[aa]->baseExternalViscosity[0] = Nodes[aa]->externalViscosity[0];
						Nodes[aa]->baseExternalViscosity[1] = Nodes[aa]->externalViscosity[1];
						Nodes[aa]->baseExternalViscosity[2] = Nodes[aa]->externalViscosity[2];
						//Nodes[aa]->externalViscosity[0] *= ECMSofteningFraction;
						//Nodes[aa]->externalViscosity[1] *= ECMSofteningFraction;
						//Nodes[aa]->externalViscosity[2] *= ECMSofteningFraction;
						//Nodes[aa]->baseExternalViscosity[0] = Nodes[aa]->externalViscosity[0];
						//Nodes[aa]->baseExternalViscosity[1] = Nodes[aa]->externalViscosity[1];
						//Nodes[aa]->baseExternalViscosity[2] = Nodes[aa]->externalViscosity[2];
					}
				}
			}
		}
	}
}

double Simulation::calculateAverageDisplacement(){
	double displacement = 0;
	int counter = 0;
	for (int aa=0; aa<nNodes; ++aa){
		if((Nodes[aa]->tissuePlacement == 0 && remodelBasalECM ) || (Nodes[aa]->tissuePlacement == 1 && remodelApicalECM )){
			double currDisplacement = Nodes[aa]->getDisplacement();
			displacement += currDisplacement;
			counter ++;
		}
	}
	if (counter >0) {
		displacement /= counter;
	}
	return displacement;
}

void Simulation::updateECMVisocityWithDeformationRate(){
	double averageDisplacement = calculateAverageDisplacement();
	double currECMRemodellingFraction = ECMRemodellingFraction /3600 * dt; //converting per hr tpo per time step
	for (int aa=0; aa<nNodes; ++aa){
		if((Nodes[aa]->tissuePlacement == 0 && remodelBasalECM ) || (Nodes[aa]->tissuePlacement == 1 && remodelApicalECM )){
			Nodes[aa]->updateECMVisocityWithDeformationRate(currECMRemodellingFraction, averageDisplacement*remodellingThresholdFraction);
		}
	}

}

bool Simulation::runOneStep(){
    bool Success = true;
	cout<<"entered run one step, time "<<currSimTimeSec<<endl;
	/*for (int i=0; i<nElements; ++i){
		if (i == 302 || i == 174 || i == 345 || i == 107){
			cout.precision(9);
			cout<<"Element: "<<Elements[i]->Id<<endl;
			Elements[i]->displayMatrix(Elements[i]->Fg,"FgAtTheBeginning");
		}
	}*/
	//ablateSpcific();
    if (currSimTimeSec == -16*3600) {
    	pokeElement(31,0,0,-0.1);pokeElement(34,0,0,-0.1);
    }
    manualPerturbationToInitialSetup(false,false); //bool deform, bool rotate
    resetForces(true); // reset the packing forces together with all the rest of the forces here
    int freq = 10.0/dt ;
    if (freq <1 ){freq =1;}
    if ((timestep - 1)% freq  == 0){
        cout<<"At time -- "<<currSimTimeSec<<" sec ("<<currSimTimeSec/3600<<" hours - "<<timestep<<" timesteps)"<<endl;
        alignTissueDVToXPositive();
        //alignTissueAPToXYPlane();
        calculateBoundingBox();
        calculateDVDistance();
        vector<ShapeBase*>::iterator itElement;
        for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        	(*itElement)->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
        }
	}
    //cout<<"after bounding box"<<endl;
    checkForExperimentalSetupsBeforeIteration();
    bool thereIsActinStrainFeedback = false;
    if (thereIsActinStrainFeedback){
    	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	        (*itElement)->calculateActinFeedback(dt);
    	        (*itElement)->updateElasticProperties();
    	}
    }
    if (thereIsECMSoftening) {
    	checkECMSoftening();
    }
    if (thereIsECMRemodellinbgWithDeforamtionRate){
    	updateECMVisocityWithDeformationRate();
    }
    if(nMyosinFunctions > 0){
    	checkForMyosinUpdates();
    }
    //cout<<"after checking for myoisn"<<endl;
    bool thereIsRefinement = false;
    if (thereIsRefinement){
    	int oldNodeSize = nNodes;
    	flagElementsThatNeedRefinement(); //this is done in parallel, the actual refinement cannot be parallel
    	refineElements();
    	cout<<"outside refine elements"<<endl;
    	if (nNodes != oldNodeSize){
    		cout<<"re-initiating system forces"<<endl;
    		reInitiateSystemForces(oldNodeSize);
    		NRSolver->reInitiateMatricesAfterRefinement(nNodes);
    	}
    }
    if(nGrowthFunctions>0 || nShapeChangeFunctions >0){
    	checkForPinningPositionsUpdate();
        //outputFile<<"calculating growth"<<endl;
		//if ((timestep - 1)% growthRotationUpdateFrequency  == 0){
			updateGrowthRotationMatrices();
		//}
        if(nGrowthFunctions>0){
        	calculateGrowth();
        }
        if(nShapeChangeFunctions>0){
        	calculateShapeChange();
        }
    }
    if(thereIsPlasticDeformation || thereIsExplicitECM){
    	updatePlasticDeformation();
    }
    if(thereIsCellMigration){
    	cellMigrationTool->updateMigratingElements(Elements);
    	cellMigrationTool->updateMigrationLists(Elements, dt);
    	cellMigrationTool->updateVolumesWithMigration(Elements);
    }
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	//(*itElement)->updateInternalViscosityTest();
        	(*itElement)->growShapeByFg();
        	//(*itElement)->defineFgByGrowthTemplate();
        	(*itElement)->changeShapeByFsc(dt);
        }
    }
    //cout<<"adding mass to nodes"<<endl;
    updateNodeMasses();
    updateNodeViscositySurfaces();
    //updateNodeSurfaces();
    //cout<<"updating connected node lists"<<endl;
    updateElementToConnectedNodes(Nodes);
    //cout<<"calculate Myosin Forces"<<endl;
    calculateMyosinForces();
    //cout<<"detect Pacing Nodes"<<endl;
    detectPacingNodes();
    if (implicitPacking == false){
    	calculatePackingForcesExplicit3D();
    }
    if (addingRandomForces){
    	calculateRandomForces();
    }
    //cout<<"starting NR"<<endl;
	/*for (int i =0 ; i<nElements; i++){
		if (i == 345 || i == 174 || i == 302 || i == 1041 || i == 107){
			cout.precision(6);
			cout<<"Before runOneStep: Time: "<<currSimTimeSec<<" Element : "<<Elements[i]->Id<<" grownVolume: "<<Elements[i]->GrownVolume<<" reference volume: "<<Elements[i]->ReferenceShape->Volume<<endl;
		}
	}*/
    updateStepNR();
	//Elements[2]->displayDebuggingMatrices();
	//Elements[115]->displayDebuggingMatrices();
    calculateBoundingBox();

    //calculateColumnarLayerBoundingBox();
	//if (thereIsPeripodialMembrane){
	//	calculatePeripodialBoundingBox();
	//}
    //cout<<" step: "<<timestep<<" Pressure: "<<SuctionPressure[2]<<" Pa, maximum z-height: "<<boundingBox[1][2]<<" L/a: "<<(boundingBox[1][2]-50)/(2.0*pipetteRadius)<<endl;
    Success = checkFlip();
    timestep++;
    currSimTimeSec += dt;
    if (Success){
    	processDisplayDataAndSave();
    }
    //if(nElements>201){
    //	Elements[99]->displayMatrix(Elements[99]->Strain,"Element99Strain");
    //	Elements[201]->displayMatrix(Elements[201]->Strain,"Element201Strain");
    //}
    return Success;
    //for (int i=0; i<nNodes; ++i){
    //	cout<<" Nodes["<<i<<"]->Position[0]="<<Nodes[i]->Position[0]<<"; Nodes["<<i<<"]->Position[1]="<<Nodes[i]->Position[1]<<";  Nodes["<<i<<"]->Position[2]="<<Nodes[i]->Position[2]<<"; "<<endl;
    //}
	//cout<<"Element: 38"<<endl;
	//Elements[38]->displayPositions();
	//cout<<"Element: 39"<<endl;
	//Elements[39]->displayPositions();
}

void Simulation::checkForPinningPositionsUpdate(){
	//cout<<"checking for pinning update, currTime: "<<currSimTimeSec<<" nGrowthPinning: "<<nGrowthPinning<<" GridGrowthsPinnedOnInitialMesh? "<<GridGrowthsPinnedOnInitialMesh<<endl;
	if (GridGrowthsPinnedOnInitialMesh){
		cout<<" grid is pinned "<<endl;
		for (int i=0; i<nGrowthPinning; ++i){
			//cout<<"growthPinUpdateTime["<<i<<"]: "<<growthPinUpdateTime[i]<<" growthPinUpdateBools["<<i<<"]: "<<growthPinUpdateBools[i]<<endl;
			if (currSimTimeSec >= growthPinUpdateTime[i] &&  !growthPinUpdateBools[i]){
				updatePinningPositions();
				growthPinUpdateBools[i] = true;
			}
		}
	}
}

void Simulation::updatePinningPositions(){
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setInitialRelativePosInBoundingBox();
	}
}

void Simulation::flagElementsThatNeedRefinement(){
	//cout<<"inside flag for refinement"<<endl;
	double thresholdSide = 5;
	double thresholdArea = thresholdSide/2.0 * thresholdSide/2.0 *pow(3,0.5);
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated && (*itElement)->tissueType ==0){//Refine only the columnar layer elements.
			if ((*itElement)->tissuePlacement == 1){ // element is apical, I should check with apical area
				(*itElement)->calculateApicalArea();
				(*itElement)->doesElementNeedRefinement(thresholdArea, 1); //checking refinement with apical surface;
			}
			else if ((*itElement)->tissuePlacement == 0){ // element is basal, I should check with basal area
				(*itElement)->calculateBasalArea();
				(*itElement)->doesElementNeedRefinement(thresholdArea, 0); //checking refinement with basal surface;
			}
		}
	}

}

void Simulation::refineElements(){
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		if (Elements[i]->willBeRefined){
			//I will refine all the elements in this column, they are
			//stored on elementsIdsOnSameColumn attribute of each element
			//This list starts from basal elements and moves to apical;
			//Clearing pu the flags that I may have on the other end of the tissue first:
			for (int j=0; j<TissueHeightDiscretisationLayers; ++j){
				Elements[Elements[i]->elementsIdsOnSameColumn[j]]->willBeRefined = false;
			}
			//I will add discretisation layers + 1 number of nodes, I want to keep the ids recorded for now:
			int* newNodeIdList = new int [TissueHeightDiscretisationLayers+1];
			addNodesForRefinement(Elements[i],newNodeIdList);
			addElementsForRefinement(Elements[i]->elementsIdsOnSameColumn,newNodeIdList);
			delete[] newNodeIdList;
		}
	}


}

void Simulation::addNodesForRefinement(ShapeBase* currElement, int* newNodeIdList){
	//cout<<" inside add nodes for refinement"<<endl;
	Node* tmp_nd;
	//now I will add nodes:
	//Adding the basal node first:
	int dividedElementId = currElement->elementsIdsOnSameColumn[0];
	double* pos = new double[3];
	Elements[dividedElementId]->getBasalCentre(pos);
	double* externalViscosity = Elements[dividedElementId]->getBasalMinViscosity(Nodes);
	//cout<<" basal centre: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
	int tissuePos = 0; //at basal position
	int tissueType = 0; //columnar layer tissue, I only refine columnar layer
	int atCircumference = 0; //newly added nodes will never be at circumferenece or at symmetry borders
	tmp_nd = new Node(nNodes, 3, pos,tissuePos, tissueType);
	tmp_nd->atCircumference = atCircumference;
	for (int j=0; j<3; ++j){
		tmp_nd->externalViscosity[j] = externalViscosity[j];
	}
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	newNodeIdList[0] = tmp_nd->Id;
	//Now adding the apical nodes:
	tissuePos = 2;//rest of the nodes are at mid-line, except for the last node to be added
	for (int i=0; i< TissueHeightDiscretisationLayers ;++i){
		if (i == TissueHeightDiscretisationLayers-1){
			tissuePos = 1; // the last node to be added is apical
		}
		dividedElementId = currElement->elementsIdsOnSameColumn[i];
		Elements[dividedElementId]->getApicalCentre(pos);
		//cout<<" apical centre: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
		double* externalViscosity = Elements[dividedElementId]->getApicalMinViscosity(Nodes);
		tmp_nd = new Node(nNodes, 3, pos,tissuePos, tissueType);
		tmp_nd->atCircumference = atCircumference;
		for (int j=0; j<3; ++j){
			tmp_nd->externalViscosity[j] = externalViscosity[j];
		}
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		newNodeIdList[i+1] = tmp_nd->Id;
		/*cout<<" Element: "<<(*itElement)->Id<<" column : ";
		for (int j=0; j<TissueHeightDiscretisationLayers; ++j){
			cout<<" "<<(*itElement)->elementsIdsOnSameColumn[j]<<" ";
		}
		cout<<endl;*/
	}
	delete[] pos;
}

void Simulation::addElementsForRefinement(int* elementsIdsOnColumn, int* newNodeIdList){
	//I will convert the existing element for each layer to the element using
	//corner 0-1-new basal node, 3 - 4 - new apical node. Then I will create two new
	//elements, that will use corners  1-2 and 2-0 (basally).

	for (int i=0; i< TissueHeightDiscretisationLayers ;++i){
		//cout<<"Started the layer: "<<i<<endl;
		ShapeBase* currElement = Elements[elementsIdsOnColumn[i]];
		const int dim = currElement->getDim();
		const int n = currElement->getNodeNumber();
		double** referenceOfShapeToBeRefined = currElement->getReferencePos();
		double * basalReferenceCentre = new double [dim];
		double * apicalReferenceCentre = new double [dim];
		currElement->getReferenceBasalCentre(basalReferenceCentre);
		currElement->getReferenceApicalCentre(apicalReferenceCentre);
		//cout<<" basal reference centre: "<<basalReferenceCentre[0]<<" "<<basalReferenceCentre[1]<<" "<<basalReferenceCentre[2]<<endl;
		//cout<<" apical reference centre: "<<apicalReferenceCentre[0]<<" "<<apicalReferenceCentre[1]<<" "<<apicalReferenceCentre[2]<<endl;

		//Removing the current element from the nodes it is attached to:
		for (int j=0; j<n; ++j){
			Nodes[currElement->NodeIds[j]]->removeFromConnectedElements(currElement->Id, currElement->VolumePerNode);
		}
		//adding the immediate neighbours before I change the node list of elemetn to be divided.
		//I am adding the new basal node to the lists of basal nodes of the element (0-2);
		Nodes[currElement->NodeIds[0]]->addToImmediateNeigs(newNodeIdList[i]);
		Nodes[currElement->NodeIds[1]]->addToImmediateNeigs(newNodeIdList[i]);
		Nodes[currElement->NodeIds[2]]->addToImmediateNeigs(newNodeIdList[i]);
		//I am adding the new apical node to the lists of apical nodes of the element (3-5);
		Nodes[currElement->NodeIds[3]]->addToImmediateNeigs(newNodeIdList[i+1]);
		Nodes[currElement->NodeIds[4]]->addToImmediateNeigs(newNodeIdList[i+1]);
		Nodes[currElement->NodeIds[5]]->addToImmediateNeigs(newNodeIdList[i+1]);
		//Now moving on to manipulating elements:
		int* NodeIds;
		NodeIds = new int[n];
		double** referencePos;
		referencePos = new double*[n];
		for (int j=0; j<n; ++j){
			referencePos[j] = new double[dim];
			for(int k=0; k<dim; ++k){
				referencePos[j][k] = 0.0;
			}
		}
		for (int k =0; k<3; ++k){
			int basalIdIndice0;
			int basalIdIndice1;
			int apicalIdIndice0;
			int apicalIdIndice1;
			if ( k == 0 ){
				basalIdIndice0 = 1;
				basalIdIndice1 = 2;
				apicalIdIndice0 = 4;
				apicalIdIndice1 = 5;
			}
			else if ( k == 1 ){
				basalIdIndice0 = 2;
				basalIdIndice1 = 0;
				apicalIdIndice0 = 5;
				apicalIdIndice1 = 3;
			}
			else if ( k == 2 ){
				//This is converting the main element into a smaller one!
				basalIdIndice0 = 0;
				basalIdIndice1 = 1;
				apicalIdIndice0 = 3;
				apicalIdIndice1 = 4;
			}
			NodeIds[0] = currElement->NodeIds[basalIdIndice0];
			NodeIds[1] = currElement->NodeIds[basalIdIndice1];
			NodeIds[2] = newNodeIdList[i];
			NodeIds[3] = currElement->NodeIds[apicalIdIndice0];
			NodeIds[4] = currElement->NodeIds[apicalIdIndice1];
			NodeIds[5] = newNodeIdList[i+1];
			for (int j=0; j<dim; ++j){
				referencePos[0][j] = referenceOfShapeToBeRefined[basalIdIndice0][j];
				referencePos[1][j] = referenceOfShapeToBeRefined[basalIdIndice1][j];
				referencePos[2][j] = basalReferenceCentre[j];
				referencePos[3][j] = referenceOfShapeToBeRefined[apicalIdIndice0][j];
				referencePos[4][j] = referenceOfShapeToBeRefined[apicalIdIndice1][j];
				referencePos[5][j] = apicalReferenceCentre[j];
			}
			if (k == 0 || k == 1){
				//cout<<" k is "<<k<<endl;
				Prism* PrismPnt01;
				PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
				PrismPnt01->updateReferencePositionMatrixFromInput(referencePos);
				PrismPnt01->checkRotationConsistency3D();
				PrismPnt01->copyElementInformationAfterRefinement(currElement,TissueHeightDiscretisationLayers,thereIsPlasticDeformation);
				PrismPnt01->calculateElementShapeFunctionDerivatives();
				//Add the weight of element to its nodes:
				for (int j = 0; j<n;++j){
					Nodes[NodeIds[j]]->addToConnectedElements(PrismPnt01->Id,PrismPnt01->VolumePerNode);
				}
				Elements.push_back(PrismPnt01);
				nElements = Elements.size();
				currElementId++;
			}
			else if (k == 2){
				//cout<<" k is two, manipulating old element"<<endl;
				//converting the current prism to a smaller version.
				currElement->updateNodeIdsForRefinement(NodeIds);
				currElement->updateReferencePositionMatrixFromInput(referencePos);
				currElement->checkRotationConsistency3D();
				currElement->updatePositions(Nodes);
				currElement->GrownVolume /=3.0;
				currElement->VolumePerNode /= 3.0;
				currElement->calculateElementShapeFunctionDerivatives();
				//Add the weight of element to its nodes, I have removed this element from its nodes above:
				for (int j = 0; j<n;++j){
					Nodes[currElement->NodeIds[j]]->addToConnectedElements(currElement->Id,currElement->VolumePerNode);
				}
				//cout<<"finalised manipulating old element: "<<currElement->Id<<endl;
			}
		}
		//clean up:
		delete[] NodeIds;
		for (int j=0; j<n; ++j){
			delete[] referencePos[j];
		}
		delete[] referencePos;
		delete[] apicalReferenceCentre;
		delete[] basalReferenceCentre;
	}
}

bool Simulation::checkFlip(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isFlipped){
			//there is a flipped element:
			outputFile<<"There is A flipped element: "<<(*itElement)->Id<<endl;
			return false;
		}
	}
	return true;
}
void Simulation::wrapUpAtTheEndOfSimulation(){
    alignTissueDVToXPositive();
    //alignTissueAPToXYPlane();
}


void Simulation::clearProjectedAreas(){
    for (int i=0;i<nNodes; ++i){
        Nodes[i]->zProjectedArea = 0.0;
    }
}

void Simulation::correctzProjectedAreaForMidNodes(){
    for (int i=0;i<nNodes; ++i){
        if (Nodes[i]->tissuePlacement == 2 ){ // the node is on midlayer
            //For the midline nodes, the area is added from apical and basal surfaces of elemetns on both sides.
            Nodes[i]->zProjectedArea /= 2.0;
        }
    }
}

void Simulation::calculateZProjectedAreas(){
    clearProjectedAreas();
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateZProjectedAreas();
        (*itElement)->assignZProjectedAreas(Nodes);
    }
    correctzProjectedAreaForMidNodes();
}

void Simulation::updatePlasticDeformation(){
	//double rate = plasticDeformationRate/3600.0*dt; //convert from per hour to per de(in sec)
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated ){
			if (thereIsExplicitECM &&  (*itElement)->isECMMimicing){
				//The simulation is defining an explicit ECM layer. All basal elements
				//should be subject to non-volum econserving plastic deformation.
				//The volume conservation should be false (first inout)
				//the half time used should be the ecm remodelling half time
				//there should not be thresholds for z remodelling, there is no
				//volume conservation, therefore these are not relevant anyway.
				(*itElement)->calculatePlasticDeformation3D(false,dt,ECMRenawalHalfLife, 0.1, 10.0);
			}
			else if (thereIsPlasticDeformation){
				//The ECM will always have plastic deformation, hence this function will be called.
				//Then for all other elements, I should first check if there is plastic deformation in the first place
				//Then I will check the tissue types
				//the lateral elements will have plastic deformation if either columnar or peripodial elemetns have the deformation.
				if(( ((*itElement)->tissueType == 0 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToColumnar) || ( ((*itElement)->tissueType == 1 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToPeripodial))
				{
					(*itElement)->calculatePlasticDeformation3D(volumeConservedInPlasticDeformation,dt,plasticDeformationHalfLife, zRemodellingLowerThreshold, zRemodellingUpperThreshold);
				}
			}
		}
		else{
			(*itElement)->setPlasticDeformationIncrement(1.0,1.0,1.0);
		}
	}
}

void Simulation::calculateNumericalJacobian(bool displayMatricesDuringNumericalCalculation){
	int dim = 3;
	gsl_matrix_set_zero(NRSolver->Knumerical);
	NRSolver->calculateDisplacementMatrix(dt);
	//PACKING SHOULD BE ADDED HERE If using this nuerical calculation!!!

	//Trying to see the manual values:
	resetForces(true); // reset the packing forces together with all the rest of the forces here
	//No perturbation:
	gsl_matrix* ge_noPerturb = gsl_matrix_calloc(dim*nNodes,1);
	gsl_matrix* gvInternal_noPerturb = gsl_matrix_calloc(dim*nNodes,1);
	NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
	NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
	gsl_matrix_memcpy(ge_noPerturb, NRSolver->ge);
	gsl_matrix_memcpy(gvInternal_noPerturb, NRSolver->gvInternal);
	//perturbation loop:
	gsl_matrix* uk_original = gsl_matrix_calloc(dim*nNodes,1);
	gsl_matrix_memcpy(uk_original,NRSolver->uk);
	for (int i=0; i<nNodes; ++i){
		for (int j=0; j<3; ++j){
	        if (implicitPacking){
	            resetForces(true);	// reset packing forces
	        }
	        else{
	            resetForces(false);	// do not reset packing forces
	        }
	        gsl_matrix_set(NRSolver->uk,i*3+j,0,gsl_matrix_get(NRSolver->uk,i*3+j,0)+1E-6);
			NRSolver->calculateDisplacementMatrix(dt);
			Nodes[i]->Position[j] += 1E-6;
			updateElementPositions();
			NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
			gsl_matrix* ge_withPerturb = gsl_matrix_calloc(dim*nNodes,1);
			gsl_matrix* gvInternal_withPerturb = gsl_matrix_calloc(dim*nNodes,1);
			NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
			gsl_matrix_memcpy(ge_withPerturb, NRSolver->ge);
			gsl_matrix_memcpy(gvInternal_withPerturb, NRSolver->gvInternal);
			//Calculate dg/dx:
			gsl_matrix_sub(ge_withPerturb,ge_noPerturb);
			gsl_matrix_sub(gvInternal_withPerturb,gvInternal_noPerturb);
			gsl_matrix_scale(ge_withPerturb,1.0/1E-6);
			gsl_matrix_scale(gvInternal_withPerturb,1.0/1E-6);
			for (int k=0; k<nNodes*3; ++k){
				double valueElastic =   0;//gsl_matrix_get(ge_withPerturb,k,0);
				double valueViscous = 	gsl_matrix_get(gvInternal_withPerturb,k,0);
				double value = valueElastic + valueViscous;
				value *= -1.0;
				gsl_matrix_set(NRSolver->K,i*3+j,k,value);
			}
			gsl_matrix_memcpy(NRSolver->uk,uk_original);
			//gsl_matrix_set(uk,i*3+j,0,gsl_matrix_get(uk,i*3+j,0)-1E-6);
			Nodes[i]->Position[j] -= 1E-6;
			updateElementPositions();
			gsl_matrix_free(ge_withPerturb);
			gsl_matrix_free(gvInternal_withPerturb);
		}
	}
	NRSolver->calcutateFixedK(Nodes);
	if (displayMatricesDuringNumericalCalculation){
		Elements[0]->displayMatrix(NRSolver->K,"numericalK");
	}
	gsl_matrix_memcpy(NRSolver->Knumerical,NRSolver->K);
	NRSolver->setMatricesToZeroInsideIteration();
	gsl_matrix_free(ge_noPerturb);
	gsl_matrix_free(gvInternal_noPerturb);
	gsl_matrix_free(uk_original);
}

void Simulation::updateStepNR(){
    int iteratorK = 0;
    bool converged = false;

    bool numericalCalculation = false;
    bool displayMatricesDuringNumericalCalculation = false;
    bool useNumericalKIncalculation = false;

    //double clock0 = ( std::clock() - simulationStartClock ) / (double) CLOCKS_PER_SEC;
    NRSolver->setMatricesToZeroAtTheBeginningOfIteration(numericalCalculation);
    NRSolver->constructUnMatrix(Nodes);
    //NRSolver->constructLumpedMassExternalViscosityMatrix(Nodes);
    NRSolver->initialteUkMatrix();
    while (!converged){
        cout<<"iteration: "<<iteratorK<<endl;
        if (implicitPacking){
            resetForces(true);	// reset packing forces
        }
        else{
            resetForces(false);	// do not reset packing forces
        }
        NRSolver->setMatricesToZeroInsideIteration();
        if (numericalCalculation){
        	calculateNumericalJacobian(displayMatricesDuringNumericalCalculation);
        }
        NRSolver->calculateDisplacementMatrix(dt);
        NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
	    //Writing elastic Forces and elastic Ke:
		NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
	    NRSolver->writeImplicitElementalKToJacobian(Elements);
	    if (numericalCalculation){
			NRSolver->calculateDifferenceBetweenNumericalAndAnalyticalJacobian(Nodes, displayMatricesDuringNumericalCalculation);
			if(useNumericalKIncalculation){
				NRSolver->useNumericalJacobianInIteration();
			}
		}
		NRSolver->calculateExternalViscousForcesForNR(Nodes);
	    NRSolver->addImplicitKViscousExternalToJacobian(Nodes,dt);
	    //These are the calculation of packing forces that would work if I wanted implicit packing
	    if (implicitPacking){
	    	calculatePackingForcesImplicit3D();
	    	calculatePackingJacobianNumerical3D(NRSolver->K);
	    }
	    //End of packing forces.
        //Now I will check if there are any nodes with zero mass, then I will be able to fill in the zero K matrix with identity if necessary.
        NRSolver->checkJacobianForAblatedNodes(AblatedNodes);
        NRSolver->calculateSumOfInternalForces();
        if (PipetteSuction && timestep >= PipetteInitialStep && timestep<PipetteEndStep){
			packToPipetteWall();
			calculateZProjectedAreas();
			addPipetteForces(NRSolver->gExt);
		}
		addMyosinForces(NRSolver->gExt);
		addPackingForces(NRSolver->gExt);
		if (addingRandomForces){
			addRandomForces(NRSolver->gExt);
		}
        NRSolver->addExernalForces();
        checkForExperimentalSetupsWithinIteration();
        NRSolver->calcutateFixedK(Nodes);
        //cout<<"displaying the jacobian after all additions"<<endl;
        //Elements[0]->displayMatrix(NRSolver->K,"theJacobian");
        //cout<<"checking convergence with forces"<<endl;
        //converged = NRSolver->checkConvergenceViaForce();
        if (converged){
            break;
        }
        //cout<<"solving for deltaU"<<endl;
        //Elements[0]->displayMatrix(NRSolver->K,"K");
        NRSolver->solveForDeltaU();
        //cout<<"checking convergence"<<endl;
        converged = NRSolver->checkConvergenceViaDeltaU();
        NRSolver->updateUkInIteration();
        updateElementPositionsinNR(NRSolver->uk);
        updateNodePositionsNR(NRSolver->uk);

        //int nNodeList = 6;
        //int nodeList[6] = {0,1,2,3,4,5};
        //for (int nodeListIterator=0; nodeListIterator<nNodeList;nodeListIterator++){
		//	int currNodeOfInterest = nodeList[nodeListIterator];
		//	if (nNodes>currNodeOfInterest){
		//		double currNodeOfInterestX = gsl_vector_get(NRSolver->deltaU,currNodeOfInterest*3);
		//		double currNodeOfInterestY = gsl_vector_get(NRSolver->deltaU,currNodeOfInterest*3+1);
		//		double currNodeOfInterestZ = gsl_vector_get(NRSolver->deltaU,currNodeOfInterest*3+2);
		//		cout<<" displacement for Node "<<currNodeOfInterest<<":  "<<currNodeOfInterestX<<" "<<currNodeOfInterestY<<" "<<currNodeOfInterestZ<<endl;
		//	}
        //}
        iteratorK ++;
        if (!converged && iteratorK > 20){
            cerr<<"Error: did not converge!!!"<<endl;
            converged = true;
        }
    }
    checkForExperimentalSetupsAfterIteration();
    //Now the calculation is converged, I update the node positions with the latest positions uk:
    updateNodePositionsNR(NRSolver->uk);
    //Element positions are already up to date.
    cout<<"finished run one step"<<endl;
}

void Simulation::calculateRandomForces(){
	randomForces.empty();
	randomForces=RandomGenerator::Obj().getNormRV( randomForceMean,randomForceVar, 3*nNodes );
	//making the sum of forces zero:
	double sumRandomForceX = 0.0;
	double sumRandomForceY = 0.0;
	double sumRandomForceZ = 0.0;

	vector<double>::iterator itDouble;
	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		sumRandomForceX+= (*itDouble);
		sumRandomForceY+= (*(itDouble+1));
		sumRandomForceZ+= (*(itDouble+2));
	}
	sumRandomForceX /= nNodes;
	sumRandomForceY /= nNodes;
	sumRandomForceZ /= nNodes;

	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		(*itDouble) -= sumRandomForceX;
		(*(itDouble+1)) -= sumRandomForceY;
		(*(itDouble+2)) -= sumRandomForceZ;
	}
}

void Simulation::addRandomForces(gsl_matrix* gExt){
		for (int j=0; j<3*nNodes; ++j){
			double F = randomForces[j];
			F += gsl_matrix_get(gExt,j,0);
			gsl_matrix_set(gExt,j,0,F);
		}
}

void Simulation::updateNodePositionsNR(gsl_matrix* uk){
    int dim = 3;
    for (int i = 0; i<nNodes; ++i){
    	Nodes[i]->updatePreviousPosition();
    	for (int j=0; j<dim; ++j){
            Nodes[i]->Position[j]=gsl_matrix_get(uk,dim*i+j,0);
        }
    }
    //cout<<"finised node pos update"<<endl;
}

void Simulation::updateElementPositionsinNR(gsl_matrix* uk){
    int dim = 3;
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        int* nodeIds = (*itElement)->getNodeIds();
        int nNodes= (*itElement)->getNodeNumber();
        for (int j=0; j<nNodes; ++j){
            double x = gsl_matrix_get(uk,dim*nodeIds[j],0);
            double y = gsl_matrix_get(uk,dim*nodeIds[j]+1,0);
            double z = gsl_matrix_get(uk,dim*nodeIds[j]+2,0);
            (*itElement)->Positions[j][0] = x;
            (*itElement)->Positions[j][1] = y;
            (*itElement)->Positions[j][2] = z;
        }
    }
}

void Simulation::calculatePackingForcesExplicit3D(){
	//cout<<"inside calculatePackingForcesExplicit3D, size of packing node couples: "<<pacingNodeCouples0.size()<<endl;
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		//if (id0 == 3445 || id1 == 3445){
		//	cout<<"calculating packing for nodes : "<<id0<<" "<<id1<<endl;
		//}

		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double d = pow((dx*dx + dy*dy + dz*dz),0.5);

		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = 5;

		double F = 5000.0 * averageMass / (1 + exp(sigmoidSaturation* d / packingThreshold));

		double Fx = F * initialWeightPointx[i];
		double Fy = F * initialWeightPointy[i];
		double Fz = F * initialWeightPointz[i];
		PackingForces[id0][0] += Fx;
		PackingForces[id0][1] += Fy;
		PackingForces[id0][2] += Fz;
		PackingForces[id1][0] -= Fx;
		PackingForces[id1][1] -= Fy;
		PackingForces[id1][2] -= Fz;
		if (id0 == 3444 || id0 == 3444 || id1 == 3444 || id1 == 3444){
			cout<<"id0: "<<id0<<" id1: "<<id1<<endl;
			cout<<" dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" d: "<<d<<endl;
			cout<<" Fz: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<" F: "<<F<<endl;
			cout<<" initialWeightPoins: "<<initialWeightPointx[i]<<" "<<initialWeightPointy[i]<<" "<<initialWeightPointz[i]<<endl;
			cout<<" PackingForces["<<id0<<"]: "<<PackingForces[id0][0]<<" "<<PackingForces[id0][1]<<" "<<PackingForces[id0][2]<<endl;
			cout<<" PackingForces["<<id1<<"]: "<<PackingForces[id1][0]<<" "<<PackingForces[id1][1]<<" "<<PackingForces[id1][2]<<endl;
		}
	}
}

void Simulation::calculatePackingForcesImplicit3D(){
	//cout<<"inside calculatePackingForcesImplicit3D, size of packing node couples: "<<pacingNodeCouples0.size()<<endl;
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double multiplier = 1000;
		//if (id0 == 3445 || id1 == 3445){
		//	cout<<"calculating packing for nodes : "<<id0<<" "<<id1<<endl;
		//}

/*		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		if (initialWeightPointx[i] < 0 ){
			dx *= -1.0;
		}
		double Fx = 100.0*dx;
		Fx *=  initialWeightPointx[i];
		PackingForces[id0][0] += Fx;
		PackingForces[id1][0] -= Fx;
*/


		//sigmoid test:
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = 5;

		if (initialWeightPointx[i]>0){
			dx *= -1.0;
		}
		if (initialWeightPointy[i]>0){
			dy *= -1.0;
		}
		if (initialWeightPointz[i]>0){
			dz *= -1.0;
		}
		double Fx = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dx)));
		Fx *= initialWeightPointx[i];

		double Fy = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dy)));
		Fy *= initialWeightPointy[i];

		double Fz = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dz)));
		Fz *= initialWeightPointz[i];

		PackingForces[id0][0] += Fx;
		PackingForces[id0][1] += Fy;
		PackingForces[id0][2] += Fz;
		PackingForces[id1][0] -= Fx;
		PackingForces[id1][1] -= Fy;
		PackingForces[id1][2] -= Fz;

		//if ((id0 == 3444 && id1 == 3444 )|| (id0 == 3190 || id1 == 3444)){
		//	cout<<"id0: "<<id0<<" id1: "<<id1<<endl;
		//	cout<<" dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" d: "<<d<<endl;
		//	cout<<" Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<" F: "<<F<<endl;
		//	cout<<" position ["<<id0<<"][0]: "<<Nodes[id0]->Position[0]<<"  position ["<<id1<<"][0]: "<<Nodes[id1]->Position[0]<<endl;
		//	cout<<" initialWeightPoins: "<<initialWeightPointx[i]<<" "<<initialWeightPointy[i]<<" "<<initialWeightPointz[i]<<endl;
		//	cout<<" PackingForces["<<id0<<"]: "<<PackingForces[id0][0]<<" "<<PackingForces[id0][1]<<" "<<PackingForces[id0][2]<<endl;
		//	cout<<" PackingForces["<<id1<<"]: "<<PackingForces[id1][0]<<" "<<PackingForces[id1][1]<<" "<<PackingForces[id1][2]<<endl;
		//}

		/*double multiplier = 1;
		double p = 2; //the power of the division (d/t);

		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double c = multiplier * averageMass;
		double Fx = c * (1 - pow(dx/packingThreshold,p));
		double Fy = c * (1 - pow(dy/packingThreshold,p));
		double Fz = c * (1 - pow(dz/packingThreshold,p));

		Fx *=  initialWeightPointx[i];
		Fy *=  initialWeightPointy[i];
		Fz *=  initialWeightPointz[i];
		PackingForces[id0][0] += Fx;
		PackingForces[id0][1] += Fy;
		PackingForces[id0][2] += Fz;
		PackingForces[id1][0] -= Fx;
		PackingForces[id1][1] -= Fy;
		PackingForces[id1][2] -= Fz;
		*/

		//if (id0 == 3445 || id1 == 3445){
		//	cout<<" packing force : "<<Fx <<" "<<Fy<<" "<<Fz<<endl;
		//	cout<<" node0("<<id0<<") pos: "<<Nodes[id0]->Position[0]<<" "<<Nodes[id0]->Position[1]<<" "<<Nodes[id0]->Position[2]<<endl;
		//	cout<<" node1("<<id1<<") pos: "<<Nodes[id1]->Position[0]<<" "<<Nodes[id1]->Position[1]<<" "<<Nodes[id1]->Position[2]<<endl;
		//	cout<<" dx, dy, dz: "<<dx<<" "<<dy<<" "<<dz<<endl;
		//	cout<<" average mass: "<<averageMass<<" dz: "<<dz<<endl;
		//}
	}
}
void Simulation::calculatePackingImplicit3DnotWorking(){
	double multiplier = 1;
	double p = 1; //the power of the division (d/t);
	double threshold = 7; //packing threshold;
	//surfaces:
	int n = pacingNodeSurfaceList0.size();
	for(int i = 0 ; i<n; ++i){
		//cout<<"checking surfaces, "<<i <<" of "<<n<<endl;
		int id0 = pacingNodeSurfaceList0[i];
		int id1 = pacingNodeSurfaceList1[i];
		int id2 = pacingNodeSurfaceList2[i];
		int id3 = pacingNodeSurfaceList3[i];
		double dx = Nodes[id0]->Position[0] - (Nodes[id1]->Position[0]+Nodes[id2]->Position[0]+Nodes[id3]->Position[0])/3.0;
		double dy = Nodes[id0]->Position[1] - (Nodes[id1]->Position[1]+Nodes[id2]->Position[1]+Nodes[id3]->Position[1])/3.0;
		double dz = Nodes[id0]->Position[2] - (Nodes[id1]->Position[2]+Nodes[id2]->Position[2]+Nodes[id3]->Position[2])/3.0;
		//if dz was initially negative, I need to convert it to positive, I do not care if it became positive now. If dx was initially positive, multiplier is 1 anyway
		//dx *= initialSignsSurfacex[i];
		//dy *= initialSignsSurfacey[i];
		//dz *= initialSignsSurfacez[i];
		double averageMass = 0.25 *( Nodes[id0]->mass + Nodes[id1]->mass +  Nodes[id2]->mass  + Nodes[id3]->mass);
		double c = multiplier * averageMass;
		//double Fmagx = initialSignsSurfacex[i] * c * (1 - pow(dx/threshold,p));
		//double Fmagy = initialSignsSurfacey[i] * c * (1 - pow(dy/threshold,p));
		//double Fmagz = initialSignsSurfacez[i] * c * (1 - pow(dz/threshold,p));


		double Fmag = c * (1 - pow((dx*dx + dy*dy + dz*dz)/threshold/threshold,p));
		double Fmagx = initialWeightSurfacex[i]*Fmag;
		double Fmagy = initialWeightSurfacey[i]*Fmag;
		double Fmagz = initialWeightSurfacez[i]*Fmag;

		PackingForces[id0][0] += Fmagx ;
		PackingForces[id0][1] += Fmagy ;
		PackingForces[id0][2] += Fmagz ;
		Fmagx /= 3.0;
		Fmagy /= 3.0;
		Fmagz /= 3.0;

		PackingForces[id1][0] -= Fmagx ;
		PackingForces[id1][1] -= Fmagy ;
		PackingForces[id1][2] -= Fmagz ;
		PackingForces[id2][0] -= Fmagx ;
		PackingForces[id2][1] -= Fmagy ;
		PackingForces[id2][2] -= Fmagz ;
		PackingForces[id3][0] -= Fmagx ;
		PackingForces[id3][1] -= Fmagy ;
		PackingForces[id3][2] -= Fmagz ;
		cout<<"ids: "<<id0<<" "<<id1<<" "<<id2<<" "<<id3<<" d: "<<dx<<" "<<dy<<" "<<dz<<" F: "<<Fmagx<<" "<<Fmagy<<" "<<Fmagz<<endl;
	}
	//edges:
	n = pacingNodeEdgeList0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeEdgeList0[i];
		int id1 = pacingNodeEdgeList1[i];
		int id2 = pacingNodeEdgeList2[i];
		//cout<<"checking edges, "<<i <<" of "<<n<<" node ids: "<<id0<<" "<<id1<<" "<<id2<<endl;
		double dx = Nodes[id0]->Position[0] - (Nodes[id1]->Position[0]+Nodes[id2]->Position[0])/2.0;
		double dy = Nodes[id0]->Position[1] - (Nodes[id1]->Position[1]+Nodes[id2]->Position[1])/2.0;
		double dz = Nodes[id0]->Position[2] - (Nodes[id1]->Position[2]+Nodes[id2]->Position[2])/2.0;
		//if dz was initially negative, I need to convert it to positive, I do not care if it became positive now. If dx was initially positive, multiplier is 1 anyway
		//dx *= initialSignsEdgex[i];
		//dy *= initialSignsEdgey[i];
		//dz *= initialSignsEdgez[i];
		double averageMass = 1.0/3.0 *( Nodes[id0]->mass + Nodes[id1]->mass +  Nodes[id2]->mass );
		double c = multiplier * averageMass;
		//double Fmagx = initialSignsEdgex[i] * c * (1 - pow(dx/threshold,p));
		//double Fmagy = initialSignsEdgey[i] * c * (1 - pow(dy/threshold,p));
		//double Fmagz = initialSignsEdgez[i] * c * (1 - pow(dz/threshold,p));

		double Fmag = c * (1 - pow((dx*dx + dy*dy + dz*dz)/threshold/threshold,p));
		double Fmagx = initialWeightEdgex[i]*Fmag;
		double Fmagy = initialWeightEdgey[i]*Fmag;
		double Fmagz = initialWeightEdgez[i]*Fmag;

		PackingForces[id0][0] += Fmagx ;
		PackingForces[id0][1] += Fmagy ;
		PackingForces[id0][2] += Fmagz ;
		Fmagx /= 2.0;
		Fmagy /= 2.0;
		Fmagz /= 2.0;

		PackingForces[id1][0] -= Fmagx ;
		PackingForces[id1][1] -= Fmagy ;
		PackingForces[id1][2] -= Fmagz ;
		PackingForces[id2][0] -= Fmagx ;
		PackingForces[id2][1] -= Fmagy ;
		PackingForces[id2][2] -= Fmagz ;
		cout<<"ids: "<<id0<<" "<<id1<<" "<<id2<<" d: "<<dx<<" "<<dy<<" "<<dz<<" F: "<<Fmagx<<" "<<Fmagy<<" "<<Fmagz<<endl;
	}
	//nodes:
	n = pacingNodePointList0.size();
	for(int i = 0 ; i<n; ++i){
		//cout<<"checking nodes, "<<i <<" of "<<n<<endl;
		int id0 = pacingNodePointList0[i];
		int id1 = pacingNodePointList1[i];
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		//if dz was initially negative, I need to convert it to positive, I do not care if it became positive now. If dx was initially positive, multiplier is 1 anyway
		//dx *= initialSignsPointx[i];
		//dy *= initialSignsPointy[i];
		//dz *= initialSignsPointz[i];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double c = multiplier * averageMass;
		//double Fmagx = initialSignsPointx[i] * c * (1 - pow(dx/threshold,p));
		//double Fmagy = initialSignsPointy[i] * c * (1 - pow(dy/threshold,p));
		//double Fmagz = initialSignsPointz[i] * c * (1 - pow(dz/threshold,p));

		double Fmag = c * (1 - pow((dx*dx + dy*dy + dz*dz)/threshold/threshold,p));
		double Fmagx = initialWeightPointx[i]*Fmag;
		double Fmagy = initialWeightPointy[i]*Fmag;
		double Fmagz = initialWeightPointz[i]*Fmag;


		//PackingForces[id0][0] += Fmagx ;
		//PackingForces[id0][1] += Fmagy ;
		PackingForces[id0][2] += Fmagz ;

		//PackingForces[id1][0] -= Fmagx ;
		//PackingForces[id1][1] -= Fmagy ;
		PackingForces[id1][2] -= Fmagz ;
		cout<<"ids: "<<id0<<" "<<id1<<" d: "<<dx<<" "<<dy<<" "<<dz<<" F: "<<Fmagx<<" "<<Fmagy<<" "<<Fmagz<<endl;
	}
}

void Simulation::calculatePackingJacobianNumerical3D(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double multiplier = 1000;
		/*double multiplier = 1;
		double p = 2; //the power of the division (d/t);
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );

		double c = multiplier * averageMass;
		double perturbation = 0.0000001;

		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double Fx = c * (1 - pow(dx/packingThreshold,p));
		double Fy = c * (1 - pow(dy/packingThreshold,p));
		double Fz = c * (1 - pow(dz/packingThreshold,p));

		Fx *=  initialWeightPointx[i];
		Fy *=  initialWeightPointy[i];
		Fz *=  initialWeightPointz[i];
		double F0[3] = {Fx, Fy,  Fz};

		dx += perturbation;
		double Fxp =  initialWeightPointx[i] * c * (1 - pow(dx/packingThreshold,p));
		double F1[3] = {Fxp, Fy,  Fz};

		dy += perturbation;
		double Fyp =  initialWeightPointy[i] * c * (1 - pow(dy/packingThreshold,p));;
		double F2[3] = {Fx, Fyp,  Fz};

		dz += perturbation;
		double Fzp =  initialWeightPointz[i] * c * (1 - pow(dz/packingThreshold,p));;
		double F3[3] = {Fx, Fy,  Fzp};

		double dgcxx = (F1[0]-F0[0])/ perturbation;
		double dgcyy = (F2[1]-F0[1])/ perturbation;
		double dgczz = (F3[2]-F0[2])/ perturbation;


		double value = gsl_matrix_get(K,3*id0,3*id0);
		value += dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		value = gsl_matrix_get(K,3*id0,3*id1);
		value -= dgcxx;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id1);
		value += dgcxx;
		gsl_matrix_set(K,3*id1,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id1,3*id0,value);

		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value += dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value += dgcyy;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);

		//cout<<"dgczz: "<<dgczz<<" dgcxy: "<<dgcxy<<" dgcyx: "<<dgcyx<<" c: "<<c<<endl;
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);

		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value += dgczz;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);
		*/


/*
		double dFxdx0 = 100.0*initialWeightPointx[i];
		double dFxdx1 = -100.0*initialWeightPointx[i];

		if (initialWeightPointx[i] < 0 ){
			dFxdx0 *= -1.0;
			dFxdx1 *= -1.0;
		}
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id0,3*id0,value);

		value = gsl_matrix_get(K,3*id0,3*id1);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id0,3*id1,value);

		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id1,3*id1,value);

		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id1,3*id0,value);
*/


		//sigmoid test:
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = 5;
		if (initialWeightPointx[i]>0){
			dx *= -1.0;
		}
		if (initialWeightPointy[i]>0){
			dy *= -1.0;
		}
		if (initialWeightPointz[i]>0){
			dz *= -1.0;
		}
		double sigmoidx =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dx) ));
		double dFxdx0 = sigmoidx * (1 - sigmoidx) * multiplier * averageMass * initialWeightPointx[i] * (sigmoidSaturation/packingThreshold);
		double dFxdx1 = -1.0*dFxdx0;
		if (initialWeightPointx[i]>0){
			dFxdx0 *= -1.0;
			dFxdx1 *= -1.0;
		}

		double sigmoidy =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dy) ));
		double dFydy0 = sigmoidy * (1 - sigmoidy) * multiplier * averageMass * initialWeightPointy[i] * (sigmoidSaturation/packingThreshold);
		double dFydy1 = -1.0*dFydy0;
		if (initialWeightPointy[i]>0){
			dFydy0 *= -1.0;
			dFydy1 *= -1.0;
		}

		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * averageMass * initialWeightPointz[i] * (sigmoidSaturation/packingThreshold);
		double dFzdz1 = -1.0*dFzdz0;
		if (initialWeightPointz[i]>0){
			dFzdz0 *= -1.0;
			dFzdz1 *= -1.0;
		}
		//x values:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id0,3*id0,value);
		value = gsl_matrix_get(K,3*id0,3*id1);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id1,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id1,3*id0,value);

		//y values:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value -= dFydy0;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value -= dFydy1;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value -= dFydy0;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dFydy1;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);

		//z values:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -= dFzdz1;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dFzdz1;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);

	}
}

void Simulation::calculatePackingNumerical3DnotWorking(gsl_matrix* K){
	double multiplier = 1;
	double p = 1; //the power of the division (d/t);
	double threshold = 6.0; //packing threshold
/*
	//surfaces:
	int n = pacingNodeSurfaceList0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeSurfaceList0[i];
		int id1 = pacingNodeSurfaceList1[i];
		int id2 = pacingNodeSurfaceList2[i];
		int id3 = pacingNodeSurfaceList3[i];
		double averageMass = 0.25 *( Nodes[id0]->mass + Nodes[id1]->mass +  Nodes[id2]->mass  + Nodes[id3]->mass);
		//calculating the base force (F0):
		double dx = Nodes[id0]->Position[0] - (Nodes[id1]->Position[0]+Nodes[id2]->Position[0]+Nodes[id3]->Position[0])/3.0;
		double dy = Nodes[id0]->Position[1] - (Nodes[id1]->Position[1]+Nodes[id2]->Position[1]+Nodes[id3]->Position[1])/3.0;
		double dz = Nodes[id0]->Position[2] - (Nodes[id1]->Position[2]+Nodes[id2]->Position[2]+Nodes[id3]->Position[2])/3.0;
		double c = multiplier * averageMass;
		double F0[3] = {0.0,0.0,0.0};
		F0[0] = c * (1 - pow(dx/threshold,p));
		F0[1] = c * (1 - pow(dy/threshold,p));
		F0[2] = c * (1 - pow(dz/threshold,p));
		//a perturbation to x position of base node (id0)    :  dx -> dx + (perturbation)
		//any perturbation to any of the slave nodes (id1-3) :  dx -> dx - (perturbation/3.0)
		double perturbation = 0.0000001;
		double F1[3] = {0.0,0.0,0.0};
		F1[0] = c * (1 - pow((dx+perturbation)/threshold,p));
		F1[1] = c * (1 - pow((dy+perturbation)/threshold,p));
		F1[2] = c * (1 - pow((dz+perturbation)/threshold,p));
		double dgcxx = (F1[0]-F0[0])/ perturbation;
		double dgcyy = (F1[1]-F0[1])/ perturbation;
		double dgczz = (F1[2]-F0[2])/ perturbation;
		double dgcxxSlave = -dgcxx/3.0;
		double dgcyySlave = -dgcyy/3.0;
		double dgczzSlave = -dgczz/3.0;
		//add to master_x - master_x:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value += dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		//add to master_y - master_y:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value += dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		//add to master_z - master_z:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);

		//add to master - slaves: ( how much master force changes upon slave movement)
		value = gsl_matrix_get(K,3*id0,3*id1);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id0,3*id2);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id2,value);
		value = gsl_matrix_get(K,3*id0,3*id3);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id3,value);

		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id2+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id2+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id3+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id3+1,value);

		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id2+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id2+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id3+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id3+2,value);
		//add to slave - slave: ( how much slave force changes upon slave movement)
		//slave to slave xx:
		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id2);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id2,value);
		value = gsl_matrix_get(K,3*id1,3*id3);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id3,value);

		value = gsl_matrix_get(K,3*id2,3*id2);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id2,3*id2,value);
		value = gsl_matrix_get(K,3*id2,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id2,3*id1,value);
		value = gsl_matrix_get(K,3*id2,3*id3);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id2,3*id3,value);

		value = gsl_matrix_get(K,3*id3,3*id3);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id3,3*id3,value);
		value = gsl_matrix_get(K,3*id3,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id3,3*id1,value);
		value = gsl_matrix_get(K,3*id3,3*id2);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id3,3*id2,value);

		//slave to slave yy:
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id2+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id2+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id3+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id3+1,value);

		value = gsl_matrix_get(K,3*id2+1,3*id2+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id2+1,3*id2+1,value);
		value = gsl_matrix_get(K,3*id2+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id2+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id2+1,3*id3+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id2+1,3*id3+1,value);

		value = gsl_matrix_get(K,3*id3+1,3*id3+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id3+1,3*id3+1,value);
		value = gsl_matrix_get(K,3*id3+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id3+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id3+1,3*id2+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id3+1,3*id2+1,value);

		//slave to slave zz:
		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id2+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id2+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id3+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id3+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id2+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id2+2,3*id2+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id2+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id3+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id2+2,3*id3+2,value);
		value = gsl_matrix_get(K,3*id3+2,3*id3+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id3+2,3*id3+2,value);
		value = gsl_matrix_get(K,3*id3+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id3+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id3+2,3*id2+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id3+2,3*id2+2,value);

		//add to slave - master: ( how much slave force changes upon master movement)
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id1,3*id0,value);
		value = gsl_matrix_get(K,3*id2,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id2,3*id0,value);
		value = gsl_matrix_get(K,3*id3,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id3,3*id0,value);

		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id2+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id2+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id3+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id3+1,3*id0+1,value);

		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id2+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id3+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id3+2,3*id0+2,value);
	}
	n = pacingNodeEdgeList0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeEdgeList0[i];
		int id1 = pacingNodeEdgeList1[i];
		int id2 = pacingNodeEdgeList2[i];
		double averageMass = 1.0/3.0 *( Nodes[id0]->mass + Nodes[id1]->mass +  Nodes[id2]->mass );
		//calculating the base force (F0):
		double dx = Nodes[id0]->Position[0] - (Nodes[id1]->Position[0]+Nodes[id2]->Position[0])/2.0;
		double dy = Nodes[id0]->Position[1] - (Nodes[id1]->Position[1]+Nodes[id2]->Position[1])/2.0;
		double dz = Nodes[id0]->Position[2] - (Nodes[id1]->Position[2]+Nodes[id2]->Position[2])/2.0;
		double c = multiplier * averageMass;
		double F0[3] = {0.0,0.0,0.0};
		F0[0] = c * (1 - pow(dx/threshold,p));
		F0[1] = c * (1 - pow(dy/threshold,p));
		F0[2] = c * (1 - pow(dz/threshold,p));
		//a perturbation to x position of base node (id0)    :  dx -> dx + (perturbation)
		//any perturbation to any of the slave nodes (id1-3) :  dx -> dx - (perturbation/2.0)
		double perturbation = 0.0000001;
		double F1[3] = {0.0,0.0,0.0};
		F1[0] = c * (1 - pow((dx+perturbation)/threshold,p));
		F1[1] = c * (1 - pow((dy+perturbation)/threshold,p));
		F1[2] = c * (1 - pow((dz+perturbation)/threshold,p));
		double dgcxx = (F1[0]-F0[0])/ perturbation;
		double dgcyy = (F1[1]-F0[1])/ perturbation;
		double dgczz = (F1[2]-F0[2])/ perturbation;
		double dgcxxSlave = -dgcxx/2.0;
		double dgcyySlave = -dgcyy/2.0;
		double dgczzSlave = -dgczz/2.0;
		//add to master_x - master_x:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value += dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		//add to master_y - master_y:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value += dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		//add to master_z - master_z:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);

		//add to master - slaves: ( how much master force changes upon slave movement)
		value = gsl_matrix_get(K,3*id0,3*id1);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id0,3*id2);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id2,value);


		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id2+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id2+1,value);

		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id2+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id2+2,value);
		//add to slave - slave: ( how much slave force changes upon slave movement)
		//slave to slave xx:
		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id2);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id2,value);

		value = gsl_matrix_get(K,3*id2,3*id2);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id2,3*id2,value);
		value = gsl_matrix_get(K,3*id2,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id2,3*id1,value);

		//slave to slave yy:
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id2+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id2+1,value);

		value = gsl_matrix_get(K,3*id2+1,3*id2+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id2+1,3*id2+1,value);
		value = gsl_matrix_get(K,3*id2+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id2+1,3*id1+1,value);

		//slave to slave zz:
		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id2+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id2+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id2+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id2+2,3*id2+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id2+2,3*id1+2,value);

		//add to slave - master: ( how much slave force changes upon master movement)
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id1,3*id0,value);
		value = gsl_matrix_get(K,3*id2,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id2,3*id0,value);

		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id2+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id2+1,3*id0+1,value);

		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id2+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id2+2,3*id0+2,value);
	}*/
	//KcalcForPoint
	int n = pacingNodePointList0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodePointList0[i];
		int id1 = pacingNodePointList1[i];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass);
		//calculating the base force (F0):
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		double c = multiplier * averageMass;

		double Fmag = c * (1 - pow((dx*dx + dy*dy + dz*dz)/threshold/threshold,p));
		double F0[3] = {0.0,0.0,0.0};
		F0[0] = initialWeightPointx[i]*Fmag;
		F0[1] = initialWeightPointy[i]*Fmag;
		F0[2] = initialWeightPointz[i]*Fmag;

		//perturbation to dx, as in moveing master in +perturbation dir:
		double perturbation = 0.0000001;
		double dxp = dx + perturbation;
		Fmag = c * (1 - pow((dxp*dxp + dy*dy + dz*dz)/threshold/threshold,p));
		double F1[3] = {0.0,0.0,0.0};
		F1[0] = initialWeightPointx[i]*Fmag;
		F1[1] = initialWeightPointy[i]*Fmag;
		F1[2] = initialWeightPointz[i]*Fmag;

		//perturbation to dy, as in moving master in +perturbation dir:
		double dyp = dy + perturbation;
		Fmag = c * (1 - pow((dx*dx + dyp*dyp + dz*dz)/threshold/threshold,p));
		double F2[3] = {0.0,0.0,0.0};
		F2[0] = initialWeightPointx[i]*Fmag;
		F2[1] = initialWeightPointy[i]*Fmag;
		F2[2] = initialWeightPointz[i]*Fmag;

		//perturbation to dz, as in moving master in +perturbation dir:
		double dzp = dz + perturbation;
		Fmag = c * (1 - pow((dx*dx + dy*dy + dzp*dzp)/threshold/threshold,p));
		double F3[3] = {0.0,0.0,0.0};
		F3[0] = initialWeightPointx[i]*Fmag;
		F3[1] = initialWeightPointy[i]*Fmag;
		F3[2] = initialWeightPointz[i]*Fmag;

		double dgcxx = (F1[0]-F0[0])/ perturbation;
		double dgcxy = (F1[1]-F0[1])/ perturbation;
		double dgcxz = (F1[2]-F0[2])/ perturbation;
		double dgcyy = (F2[1]-F0[1])/ perturbation;
		double dgcyx = (F2[0]-F0[0])/ perturbation;
		double dgcyz = (F2[2]-F0[2])/ perturbation;
		double dgczz = (F3[2]-F0[2])/ perturbation;
		double dgczx = (F3[0]-F0[0])/ perturbation;
		double dgczy = (F3[1]-F0[1])/ perturbation;
		dgcxx *= -1.0;
		dgcyy *= -1.0;
		dgczz *= -1.0;
		dgcxy *= -1.0;
		dgcxz *= -1.0;
		dgcyx *= -1.0;
		dgcyz *= -1.0;
		dgczx *= -1.0;
		dgczy *= -1.0;
		//dgcxx = 0.0;
		//dgcyy = 0.0;
		//dgczz = 0.0;
		//dgcxy = 0.0;
		//dgcxz = 0.0;
		//dgcyx = 0.0;
		//dgcyz = 0.0;
		//dgczx = 0.0;
		//dgczy = 0.0;

		double dgcxxSlave = -dgcxx;
		double dgcxySlave = -dgcxy;
		double dgcxzSlave = -dgcxz;
		double dgcyySlave = -dgcyy;
		double dgcyxSlave = -dgcyx;
		double dgcyzSlave = -dgcyz;
		double dgczzSlave = -dgczz;
		double dgczxSlave = -dgczx;
		double dgczySlave = -dgczy;


		//add to master_x - master_x:
		addValueToMatrix(K,3*id0,3*id0,dgcxx);
		//add to master_x - master_y:
		addValueToMatrix(K,3*id0,3*id0+1,dgcxy);
		//add to master_x - master_z:
		addValueToMatrix(K,3*id0,3*id0+2,dgcxz);
		//add to master_y - master_y:
		addValueToMatrix(K,3*id0+1,3*id0+1,dgcyy);
		//add to master_y - master_x:
		addValueToMatrix(K,3*id0+1,3*id0,dgcyx);
		//add to master_y - master_z:
		addValueToMatrix(K,3*id0+1,3*id0+2,dgcyz);
		//add to master_z - master_z:
		addValueToMatrix(K,3*id0+2,3*id0+2,dgczz);
		//add to master_z - master_x:
		addValueToMatrix(K,3*id0+2,3*id0,dgczx);
		//add to master_z - master_y:
		addValueToMatrix(K,3*id0+2,3*id0+1,dgczy);

		//add values for master To slave:
		addValueToMatrix(K,3*id0,3*id1,dgcxxSlave);
		addValueToMatrix(K,3*id0,3*id1+1,dgcxySlave);
		addValueToMatrix(K,3*id0,3*id1+2,dgcxzSlave);
		addValueToMatrix(K,3*id0+1,3*id1+1,dgcyySlave);
		addValueToMatrix(K,3*id0+1,3*id1,dgcyxSlave);
		addValueToMatrix(K,3*id0+1,3*id1+2,dgcyzSlave);
		addValueToMatrix(K,3*id0+2,3*id1+2,dgczzSlave);
		addValueToMatrix(K,3*id0+2,3*id1,dgczxSlave);
		addValueToMatrix(K,3*id0+2,3*id1+1,dgczySlave);

		//add values for slave To slave:
		addValueToMatrix(K,3*id1,3*id1,-dgcxxSlave);
		addValueToMatrix(K,3*id1,3*id1+1,-dgcxySlave);
		addValueToMatrix(K,3*id1,3*id1+2,-dgcxzSlave);
		addValueToMatrix(K,3*id1+1,3*id1+1,-dgcyySlave);
		addValueToMatrix(K,3*id1+1,3*id1,-dgcyxSlave);
		addValueToMatrix(K,3*id1+1,3*id1+2,-dgcyzSlave);
		addValueToMatrix(K,3*id1+2,3*id1+2,-dgczzSlave);
		addValueToMatrix(K,3*id1+2,3*id1,-dgczxSlave);
		addValueToMatrix(K,3*id1+2,3*id1+1,-dgczySlave);

		//add values for slave To master:
		addValueToMatrix(K,3*id1,3*id0,-dgcxx);
		addValueToMatrix(K,3*id1,3*id0+1,-dgcxy);
		addValueToMatrix(K,3*id1,3*id0+2,-dgcxz);
		addValueToMatrix(K,3*id1+1,3*id0+1,-dgcyy);
		addValueToMatrix(K,3*id1+1,3*id0,-dgcyx);
		addValueToMatrix(K,3*id1+1,3*id0+2,-dgcyz);
		addValueToMatrix(K,3*id1+2,3*id0+2,-dgczz);
		addValueToMatrix(K,3*id1+2,3*id0,-dgczx);
		addValueToMatrix(K,3*id1+2,3*id0+1,-dgczy);




	/*	//a perturbation to x position of base node (id0)    :  dx -> dx + (perturbation)
		//any perturbation to any of the slave nodes (id1-3) :  dx -> dx - (perturbation/2.0)
		double perturbation = 0.0000001;
		double F1[3] = {0.0,0.0,0.0};
		F1[0] = c * (1 - pow((dx+perturbation)/threshold,p));
		F1[1] = c * (1 - pow((dy+perturbation)/threshold,p));
		F1[2] = c * (1 - pow((dz+perturbation)/threshold,p));
		double dgcxx = (F1[0]-F0[0])/ perturbation;
		double dgcyy = (F1[1]-F0[1])/ perturbation;
		double dgczz = (F1[2]-F0[2])/ perturbation;
		double dgcxxSlave = -dgcxx;
		double dgcyySlave = -dgcyy;
		double dgczzSlave = -dgczz;

		//add to master_x - master_x:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value += dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		//add to master_y - master_y:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value += dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		//add to master_z - master_z:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);

		//add to master - slaves: ( how much master force changes upon slave movement)
		value = gsl_matrix_get(K,3*id0,3*id1);
		value += dgcxxSlave;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value += dgcyySlave;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value += dgczzSlave;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);

		//add to slave - slave: ( how much slave force changes upon slave movement)
		//slave to slave xx:
		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dgcxxSlave;
		gsl_matrix_set(K,3*id1,3*id1,value);

		//slave to slave yy:
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value -= dgcyySlave;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);

		//slave to slave zz:
		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value -= dgczzSlave;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);

		//add to slave - master: ( how much slave force changes upon master movement)
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dgcxx;
		gsl_matrix_set(K,3*id1,3*id0,value);

		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dgcyy;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);

		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);*/

	}
}

void Simulation::addValueToMatrix(gsl_matrix* K, int i, int j , double value){
	value += gsl_matrix_get(K,i,j);
	gsl_matrix_set(K,i,j,value);
}

void Simulation::calculatePackingNumerical(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double multiplier = 1;
		double c = multiplier * averageMass;
		double p = 2; //the power of the division (d/t);
		double threshold = 6.0; //packing threshold;

		double perturbation = 0.0000001;
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		double Fmag = c * (1 - pow(dz/threshold,p));
		double F0[3] = {0, 0, Fmag};
		
		dz += perturbation;
		Fmag = c * (1 - pow(dz/threshold,p));
		double F3[3] = {0 , 0 , Fmag};

		double dgczz = (F3[2]-F0[2])/ perturbation;

		//cout<<"dgczz: "<<dgczz<<endl;//" dgcxy: "<<dgcxy<<" dgcyx: "<<dgcyx<<" c: "<<c<<endl;
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id0+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id0+2)<<endl;
		double value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id0+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id0+2)<<endl;
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id1+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id1+2)<<endl;
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id1+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id1+2)<<endl;

		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value += dgczz;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);

		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);

		/*
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value +=dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		//y_i y_i
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value +=dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		//z_i z_i
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value +=dgczz;
		cout<<"K matrix 92,92: "<<gsl_matrix_get(K,92,92)<<endl;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		cout<<"K matrix 92,92: "<<gsl_matrix_get(K,92,92)<<endl;
		//x_i y_i  
		value = gsl_matrix_get(K,3*id0,3*id0+1);
		value +=dgcxy;
		gsl_matrix_set(K,3*id0,3*id0+1,value);
		//y_i x_i
		value = gsl_matrix_get(K,3*id0+1,3*id0);
		value +=dgcyx;
		gsl_matrix_set(K,3*id0+1,3*id0,value);
		//x_i z_i  
		value = gsl_matrix_get(K,3*id0,3*id0+2);
		value +=dgcxz;
		gsl_matrix_set(K,3*id0,3*id0+2,value);
		//z_i x_i
		value = gsl_matrix_get(K,3*id0+2,3*id0);
		value +=dgczx;
		gsl_matrix_set(K,3*id0+2,3*id0,value);
		//y_i z_i  
		value = gsl_matrix_get(K,3*id0+1,3*id0+2);
		value +=dgcyz;
		gsl_matrix_set(K,3*id0+1,3*id0+2,value);
		//z_i y_i
		value = gsl_matrix_get(K,3*id0+2,3*id0+1);
		value +=dgczy;
		gsl_matrix_set(K,3*id0+2,3*id0+1,value);

		cout<<"K matrix 92,92: before adding to slave "<<gsl_matrix_get(K,92,92)<<endl;
		//add the negative values to the slave:
		//x_i x_slave
		value = gsl_matrix_get(K,3*id0,3*id1);
		value -=dgcxx;
		gsl_matrix_set(K,3*id0,3*id1,value);
		//y_i y_slave
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value -=dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		//z_i z_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -=dgczz;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		//x_i y_slave  
		value = gsl_matrix_get(K,3*id0,3*id1+1);
		value -=dgcxy;
		gsl_matrix_set(K,3*id0,3*id1+1,value);
		//y_i x_slave
		value = gsl_matrix_get(K,3*id0+1,3*id1);
		value -=dgcyx;
		gsl_matrix_set(K,3*id0+1,3*id1,value);
		//x_i z_slave  
		value = gsl_matrix_get(K,3*id0,3*id1+2);
		value -=dgcxz;
		gsl_matrix_set(K,3*id0,3*id1+2,value);
		//z_i x_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1);
		value -=dgczx;
		gsl_matrix_set(K,3*id0+2,3*id1,value);
		//y_i z_slave  
		value = gsl_matrix_get(K,3*id0+1,3*id1+2);
		value -=dgcyz;
		gsl_matrix_set(K,3*id0+1,3*id1+2,value);
		//z_i y_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1+1);
		value -=dgczy;
		gsl_matrix_set(K,3*id0+2,3*id1+1,value);
*/
	}
}

void Simulation::calculatePackingK(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		 int id0 = pacingNodeCouples0[i];
		 int id1 = pacingNodeCouples1[i];
		 double multiplier = 1.0;
		 double p = 0.3; //the power of the division (d/t);
		 double threshold = 10.0; //packing threshold;
		 double tp = pow(threshold,p);

		 double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		 double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		 double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		 double d = dx*dx + dy*dy + dz*dz;  //distance between the two nodes at current itertion
		 d = pow(d,0.5);
		 double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		 double c = multiplier * averageMass;
		 //cout<<"C: "<<c<<endl;
		 double dp = pow(d,p);

		 double ddx = dx / d;
		 double ddy = dy / d;
		 double ddz = dz / d;

		 //double dgcxx = ( -1.0 * c /tp *p *dp /d * ddx) * dx/d + (c * (1 - dp/tp)) * ( 1.0 / d - dx /d/d * ddx );
		 //double dgcxy = ( -1.0 * c /tp *p *dp /d * ddy) * dx/d - (c * (1 - dp/tp)) * dx /d/d * ddy;
		 //double dgcxz = ( -1.0 * c /tp *p *dp /d * ddz) * dx/d - (c * (1 - dp/tp)) * dx /d/d * ddz;
		 //double dgcyy = ( -1.0 * c /tp *p *dp /d * ddy) * dy/d + (c * (1 - dp/tp)) * (1.0 / d - dy /d/d * ddy);
		 //double dgcyz = ( -1.0 * c /tp *p *dp /d * ddz) * dy/d - (c * (1 - dp/tp)) * dy /d/d * ddz;
		 //double dgczz = ( -1.0 * c /tp *p *dp /d * ddz) * dz/d + (c * (1 - dp/tp)) * (1.0 / d - dz /d/d * ddz);

		 double dgcxx = (-1.0 * c/d/d * ddx - (p-1)/tp*dp/d/d*ddx) * dx + c/d*(1.0 + dp/tp);
		 double dgcxy = (-1.0 * c/d/d * ddy - (p-1)/tp*dp/d/d*ddy) * dx;
		 double dgcxz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dx;
		 double dgcyy = (-1.0 * c/d/d * ddy - (p-1)/tp*dp/d/d*ddy) * dy + c/d*(1.0 + dp/tp);
		 double dgcyz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dy;
		 double dgczz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dz + c/d*(1.0 + dp/tp);
		 // cout<<"dgcxx: "<<dgcxx<<" dgcyy: "<<dgcyy<<" dgcxy: "<<dgcxy<<endl;
		 //x_i x_i
		 double value = gsl_matrix_get(K,3*id0,3*id0);
		 value +=dgcxx;
		 gsl_matrix_set(K,3*id0,3*id0,value);
		 //y_i y_i
		 value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		 value +=dgcyy;
		 gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		 //z_i z_i
		 value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		 value +=dgczz;
		 gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		 //x_i y_i  && y_i x_i
		 value = gsl_matrix_get(K,3*id0,3*id0+1);
		 value +=dgcxy;
		 gsl_matrix_set(K,3*id0,3*id0+1,value);
		 value = gsl_matrix_get(K,3*id0+1,3*id0);
		 value +=dgcxy;
		 gsl_matrix_set(K,3*id0+1,3*id0,value);
		 //x_i z_i  && z_i x_i
		 value = gsl_matrix_get(K,3*id0,3*id0+2);
		 value +=dgcxz;
		 gsl_matrix_set(K,3*id0,3*id0+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id0);
		 value +=dgcxz;
		 gsl_matrix_set(K,3*id0+2,3*id0,value);
		 //y_i z_i  && z_i y_i
		 value = gsl_matrix_get(K,3*id0+1,3*id0+2);
		 value +=dgcyz;
		 gsl_matrix_set(K,3*id0+1,3*id0+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id0+1);
		 value +=dgcyz;
		 gsl_matrix_set(K,3*id0+2,3*id0+1,value);

		 //add the negative values to the slave:
		 //x_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1);
		 value -=dgcxx;
		 gsl_matrix_set(K,3*id0,3*id1,value);
		 //y_i y_slave
		 value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		 value -=dgcyy;
		 gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		 //z_i z_slave
		 value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		 value -=dgczz;
		 gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		 //x_i y_slave  && y_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1+1);
		 value -=dgcxy;
		 gsl_matrix_set(K,3*id0,3*id1+1,value);
		 value = gsl_matrix_get(K,3*id0+1,3*id1);
		 value -=dgcxy;
		 gsl_matrix_set(K,3*id0+1,3*id1,value);
		 //x_i z_slave  && z_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1+2);
		 value -=dgcxz;
		 gsl_matrix_set(K,3*id0,3*id1+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id1);
		 value -=dgcxz;
		 gsl_matrix_set(K,3*id0+2,3*id1,value);
		 //y_i z_slave  && z_i y_slave
		 value = gsl_matrix_get(K,3*id0+1,3*id1+2);
		 value -=dgcyz;
		 gsl_matrix_set(K,3*id0+1,3*id1+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id1+1);
		 value -=dgcyz;
		 gsl_matrix_set(K,3*id0+2,3*id1+1,value);
	}
}

void Simulation::processDisplayDataAndSave(){
    //if (displayIsOn && !DisplaySave){
        //The simulation is not displaying a saved setup, it is running and displaying
        //I need to correct the values to be displayed, and store averages, otherwise
        //the displayed values will be from artificial setups of different RK steps. (RK1 or RK4 depending on parameter)
        //updateDisplaySaveValuesFromRK();
    //}
	cout<<"timestep: "<<timestep<<" dataSaveInterval "<<dataSaveInterval<<" currSimTimeSec: "<<currSimTimeSec<<endl;

	if (saveData && ((int)currSimTimeSec/ (int) dt) % dataSaveInterval == 0){ //timestep % dataSaveInterval == 0){
        //updateDisplaySaveValuesFromRK();
        saveStep();
    }
	cout<<"finished processDisplayDataAndSave"<<endl;
}

void Simulation::updateNodeMasses(){
    vector<Node*>::iterator itNode;
    for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
    	(*itNode)->mass = 0;
    }
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->assignVolumesToNodes(Nodes);
        }
    }
}

void Simulation::updateNodeViscositySurfaces(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated){
			(*itElement)->calculateViscositySurfaces();
		}
	}
    //calculated the areas, now assigning them
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->viscositySurface = 0;
	}
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	if (!(*itElement)->IsAblated){
    		(*itElement)->assignViscositySurfaceAreaToNodes(Nodes);
    	}
    }
    //for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
    //	cout<<" Viscosity surface of node: "<<(*itNode)->Id<<" is "<<(*itNode)->viscositySurface<<endl;
    //}
}

void 	Simulation:: updateElementToConnectedNodes(vector <Node*>& Nodes){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
	    if ((*itNode)->mass > 0){//an ablated node will have this as zero
			int n = (*itNode)->connectedElementIds.size();
			for (int i=0; i<n; ++i){
				//if(!Elements[Nodes[j]->connectedElementIds[i]] -> IsAblated){
				(*itNode)->connectedElementWeights[i] = Elements[(*itNode)->connectedElementIds[i]]->VolumePerNode/(*itNode)->mass;
				//}
			}
	    }
    }
}


void Simulation::fillInNodeNeighbourhood(){
	vector<ShapeBase*>::iterator itEle;
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->fillNodeNeighbourhood(Nodes);
	}
}
void Simulation::fillInElementColumnLists(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated && (*itElement)->tissueType ==0 ){//Check only the columnar layer elements.
			if ((*itElement)->tissuePlacement == 0){
				//start from the basal element and move up
				//then you can copy the element numbers to the other elements on the list
				(*itElement)->constructElementStackList(TissueHeightDiscretisationLayers, Elements);
			}
		}
	}
}

void Simulation::getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (!ElementPointer->ApicalNormalForPackingUpToDate){
		ElementPointer->calculateNormalForPacking(1);
	}
	normalForPacking[0] = ElementPointer->ApicalNormalForPacking[0];
	normalForPacking[1] = ElementPointer->ApicalNormalForPacking[1];
	normalForPacking[2] = ElementPointer->ApicalNormalForPacking[2];
	ElementPointer->getApicalNodePos(posCorner);
}

void Simulation::getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (!ElementPointer->BasalNormalForPackingUpToDate){
		ElementPointer->calculateNormalForPacking(0);
	}
	normalForPacking[0] = ElementPointer->BasalNormalForPacking[0];
	normalForPacking[1] = ElementPointer->BasalNormalForPacking[1];
	normalForPacking[2] = ElementPointer->BasalNormalForPacking[2];
	ElementPointer->getBasalNodePos(posCorner);
}

void Simulation::getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (NodePointer->tissuePlacement == 1){
		//node is apical, should pack to apical nodes only
		getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
	}
	else if (NodePointer->tissuePlacement == 0){
		//node is basal, should pack to basal nodes only
		getBasalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
	}
}

inline void Simulation::CapPackingForce(double& Fmag){
	double Fcap = 1e2;
	if (Fmag > Fcap){
		Fmag = Fcap;
	}
}

void Simulation::checkPackingToPipette(bool& packsToPip, double* pos, double* pipF, double mass){
	double pipThickness = 2; //microns
	double multiplier = 0.05; //just a term to scale the forces down
	double threshold = 1.0;	 //pipette packing forces start at 2 microns
	double t2 = threshold*threshold;
	double pipRange2[2] = {pipetteRadius-threshold, pipetteRadius+pipThickness+threshold};
	pipRange2[0] *= pipRange2[0];
	pipRange2[1] *= pipRange2[1];
	double dist[3] = {pos[0]-pipetteCentre[0], pos[1]-pipetteCentre[1], pos[2]-pipetteCentre[2]};
	double dist2[3] = {dist[0]*dist[0], dist[1]*dist[1],dist[2]*dist[2]};
	double xydist2 = dist2[0]+dist2[1];
	if (xydist2 > pipRange2[0] && xydist2<pipRange2[1]){
		//cout<<"Node "<<id<<" is within range - pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
		//the node is within the packing ring range
		//is it close enough to bottom?
		if (dist2[2]<t2){
			//yes the node should be pushe (-ve)z for apical suciton and (+)ve z for basal suction:
			if (dist2[2]<1e-2){dist2[2] = 1e-2;}
			double Fmag = mass * multiplier * (1.0/dist2[2] - 1.0/t2);
			pipF[0] = 0.0;
			pipF[1] = 0.0;
			pipF[2] = Fmag;
			if (ApicalSuction){
				pipF[2] *= -1.0;
			}
			packsToPip = true;
			//cout<<"Node is pushed from bottom -  packing force: "<<pipF[0]<<" "<<pipF[1]<<" "<<pipF[2]<<endl;
			return;
		}
		//the node is not close to the pipette tip in x&y but not in range to be pushed by tip bottom, it may be too far away
		//now I need to check the distance in z and then the walls:
		if( (ApicalSuction && dist[2]>0) || (!ApicalSuction && dist[2]<0) ) {
			double midWallDist2 =  pipetteRadius + 0.5*pipThickness;
			midWallDist2 *= midWallDist2;
			double d2 = 0.0;
			if (midWallDist2>xydist2){
				//the node should be pushed inside:
				multiplier *= -1.0;
				if(xydist2>pipetteRadius*pipetteRadius){
					//the node is inside the pipette: maximum force, d2 is zero:
					d2 =1e-2;
				}
				else{
					double dx = pos[0]-pipetteRadius;
					double dy = pos[1]-pipetteRadius;
					d2 = dx*dx+dy*dy;
				}
			}
			else{
				if(xydist2<(pipetteRadius+pipThickness)*(pipetteRadius+pipThickness)){
					//the node is inside the pipette: maximum force, d2 is zero:
					d2 =1e-2;
				}
				else{
					double dx = pos[0]-pipetteRadius-pipThickness;
					double dy = pos[1]-pipetteRadius-pipThickness;
					d2 = dx*dx+dy*dy;
				}
			}
			if (d2<1e-2){d2 = 1e-2;}
			double Fmag = mass * multiplier * (1.0/d2 - 1.0/t2);
			double xydist = pow(xydist2,0.5);
			Fmag /= xydist;
			pipF[0] = Fmag*dist[0];
			pipF[1] = Fmag*dist[1];
			pipF[2] = 0.0;
			//cout<<"Node "<<id<<" is pushed from side wall - packing force: "<<pipF[0]<<" "<<pipF[1]<<" "<<pipF[2]<<endl;
			packsToPip = true;
			return;
		}
	}
}
void Simulation::cleanUpPacingCombinations(){
	//each edge should be recorded once for a base node:
	//int m = pacingNodeEdgeList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"edgelist: "<<pacingNodeEdgeList0[a]<<" "<<pacingNodeEdgeList1[a]<<" "<<pacingNodeEdgeList2[a]<<endl;
	//}
	//m = pacingNodePointList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"nodelist: "<<pacingNodePointList0[a]<<" "<<pacingNodePointList1[a]<<endl;
	//}
	int n = pacingNodeEdgeList0.size();
	int i =0;
	while ( i<n){
		int currBase = pacingNodeEdgeList0[i];
		int edge1 = pacingNodeEdgeList1[i];
		int edge2 = pacingNodeEdgeList2[i];
		int j = i+1;
		while(j<n){
			if (pacingNodeEdgeList0[j] == currBase){
				if ( (edge1 == pacingNodeEdgeList1[j] && edge2 == pacingNodeEdgeList2[j]) || (edge1 == pacingNodeEdgeList2[j] && edge2 == pacingNodeEdgeList1[j])){
					//the edge is recorded again, remove these points:
					pacingNodeEdgeList0.erase(pacingNodeEdgeList0.begin()+j);
					pacingNodeEdgeList1.erase(pacingNodeEdgeList1.begin()+j);
					pacingNodeEdgeList2.erase(pacingNodeEdgeList2.begin()+j);
					initialSignsEdgex.erase(initialSignsEdgex.begin()+j);
					initialSignsEdgey.erase(initialSignsEdgey.begin()+j);
					initialSignsEdgez.erase(initialSignsEdgez.begin()+j);
					initialWeightEdgex.erase(initialWeightEdgex.begin()+j);
					initialWeightEdgey.erase(initialWeightEdgey.begin()+j);
					initialWeightEdgez.erase(initialWeightEdgez.begin()+j);
					j--;
					n--;
				}
			}
			j++;
		}
		i++;
	}
	//cleaning up nodes, each node should back to another node once:
	n = pacingNodePointList0.size();
	i =0;
	while (i<n){
		int node0 = pacingNodePointList0[i];
		int node1 = pacingNodePointList1[i];
		int j = i+1;
		while(j<n){
			if ((pacingNodePointList0[j] == node0 && pacingNodePointList1[j] == node1) || (pacingNodePointList1[j] == node0 && pacingNodePointList0[j] == node1)) {
				pacingNodePointList0.erase(pacingNodePointList0.begin()+j);
				pacingNodePointList1.erase(pacingNodePointList1.begin()+j);
				initialSignsPointx.erase(initialSignsPointx.begin()+j);
				initialSignsPointy.erase(initialSignsPointy.begin()+j);
				initialSignsPointz.erase(initialSignsPointz.begin()+j);
				initialWeightPointx.erase(initialWeightPointx.begin()+j);
				initialWeightPointy.erase(initialWeightPointy.begin()+j);
				initialWeightPointz.erase(initialWeightPointz.begin()+j);
				j--;
				n--;
			}
			j++;
		}
		i++;
	}
	//m = pacingNodeEdgeList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"edgelist: "<<pacingNodeEdgeList0[a]<<" "<<pacingNodeEdgeList1[a]<<" "<<pacingNodeEdgeList2[a]<<endl;
	//}
	//m = pacingNodePointList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"nodelist: "<<pacingNodePointList0[a]<<" "<<pacingNodePointList1[a]<<endl;
	//}
}

void Simulation::detectPacingCombinations(){
	double threshold = 7;	 //packing forces start at 4 microns - keep this lower than the force threshold!!
	double t2 = threshold*threshold;	//threshold square for rapid calculation
	//TO DO: make this the function  emptyPackingVectors();
	pacingNodeSurfaceList0.empty();
	pacingNodeSurfaceList1.empty();
	pacingNodeSurfaceList2.empty();
	pacingNodeSurfaceList3.empty();
	initialSignsSurfacex.empty();
	initialSignsSurfacey.empty();
	initialSignsSurfacez.empty();
	pacingNodeEdgeList0.empty();
	pacingNodeEdgeList1.empty();
	pacingNodeEdgeList2.empty();
	initialSignsEdgex.empty();
	initialSignsEdgey.empty();
	initialSignsEdgez.empty();
	pacingNodePointList0.empty();
	pacingNodePointList1.empty();
	initialSignsPointx.empty();
	initialSignsPointy.empty();
	initialSignsPointz.empty();
	initialWeightSurfacex.empty();
	initialWeightSurfacey.empty();
	initialWeightSurfacez.empty();
	initialWeightEdgex.empty();
	initialWeightEdgey.empty();
	initialWeightEdgez.empty();
	initialWeightPointx.empty();
	initialWeightPointy.empty();
	initialWeightPointz.empty();

	//
	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	// end of function
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
		if (NodeHasPacking){
			double* pos;
			pos = new double[3];
			(*itNode)->getCurrentPosition(pos);
			bool packedToSurface = false;
			bool packedToEdge = false;
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this node
				bool PackingToThisElement =  (*itEle)->checkPackingToThisNodeViaState(TissueHeightDiscretisationLayers, (*itNode));
				if (PackingToThisElement){
					int id1 = -1, id2 = -1, id3 = -1;
					(*itEle)->getRelevantNodesForPacking((*itNode)->tissuePlacement, id1, id2, id3);
					//make a separate function to detect packing to surface, filling the vecotrs if necessary and returning packedToSurface bool;
					//get mid point:
					/*
					//eliminating surf packing
					double midPos[3] = {
						(Nodes[id1]->Position[0] + Nodes[id2]->Position[0] + Nodes[id3]->Position[0])/3.0,
						(Nodes[id1]->Position[1] + Nodes[id2]->Position[1] + Nodes[id3]->Position[1])/3.0,
						(Nodes[id1]->Position[2] + Nodes[id2]->Position[2] + Nodes[id3]->Position[2])/3.0
						};
					double dx = pos[0] - midPos[0];
					double dy = pos[1] - midPos[1];
					double dz = pos[2] - midPos[2];
					double d2 = dx*dx + dy*dy + dz*dz;
					if (d2<t2){
						//close enough for packing , add to list:
						pacingNodeSurfaceList0.push_back((*itNode)->Id);
						pacingNodeSurfaceList1.push_back(id1);
						pacingNodeSurfaceList2.push_back(id2);
						pacingNodeSurfaceList3.push_back(id3);
						packedToSurface = true;
						//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id1<<" "<<id2<<" "<<id3<<" from element: "<<(*itEle)->Id<<endl;
						if (dx >0) {initialSignsSurfacex.push_back(1);}
						else {initialSignsSurfacex.push_back(-1);}
						if (dy >0) {initialSignsSurfacey.push_back(1);}
						else {initialSignsSurfacey.push_back(-1);}
						if (dz >0) {initialSignsSurfacez.push_back(1);}
						else {initialSignsSurfacez.push_back(-1);}
						double d = pow(d2,0.5);
						initialWeightSurfacex.push_back(dx/d);
						initialWeightSurfacey.push_back(dy/d);
						initialWeightSurfacez.push_back(dz/d);
					}
					//end of surface detection function
					 */
					/*
					//eliminatign edge packing
					if (!packedToSurface){
						//node did not pack to surface, check edges:
						//get first mid point:
						double midPos[3];
						midPos[0] =	(Nodes[id1]->Position[0] + Nodes[id2]->Position[0])/2.0;
						midPos[1] =	(Nodes[id1]->Position[1] + Nodes[id2]->Position[1])/2.0;
						midPos[2] =	(Nodes[id1]->Position[2] + Nodes[id2]->Position[2])/2.0;
						double dx = pos[0] - midPos[0];
						double dy = pos[1] - midPos[1];
						double dz = pos[2] - midPos[2];
						double d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodeEdgeList0.push_back((*itNode)->Id);
							pacingNodeEdgeList1.push_back(id1);
							pacingNodeEdgeList2.push_back(id2);
							packedToEdge = true;
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id1<<" "<<id2<<" from element: "<<(*itEle)->Id<<" midpos: "<<midPos[0]<<" "<<midPos[1]<<" "<<midPos[2]<<" d2: "<<d2<<endl;
							//cout<<"Node: "<<id1<<" pos: "<<Nodes[id1]->Position[0]<<" "<<Nodes[id1]->Position[1]<<" "<<Nodes[id1]->Position[2]<<endl;
							//cout<<"Node: "<<id2<<" pos: "<<Nodes[id2]->Position[0]<<" "<<Nodes[id2]->Position[1]<<" "<<Nodes[id2]->Position[2]<<endl;
							if (dx >0) {initialSignsEdgex.push_back(1);}
							else {initialSignsEdgex.push_back(-1);}
							if (dy >0) {initialSignsEdgey.push_back(1);}
							else {initialSignsEdgey.push_back(-1);}
							if (dz >0) {initialSignsEdgez.push_back(1);}
							else {initialSignsEdgez.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightEdgex.push_back(dx/d);
							initialWeightEdgey.push_back(dy/d);
							initialWeightEdgez.push_back(dz/d);
						}
						//get second mid point:
						midPos[0] =	(Nodes[id1]->Position[0] + Nodes[id3]->Position[0])/2.0;
						midPos[1] =	(Nodes[id1]->Position[1] + Nodes[id3]->Position[1])/2.0;
						midPos[2] =	(Nodes[id1]->Position[2] + Nodes[id3]->Position[2])/2.0;
						dx = pos[0] - midPos[0];
						dy = pos[1] - midPos[1];
						dz = pos[2] - midPos[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodeEdgeList0.push_back((*itNode)->Id);
							pacingNodeEdgeList1.push_back(id1);
							pacingNodeEdgeList2.push_back(id3);
							packedToEdge = true;
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id1<<" "<<id3<<" from element: "<<(*itEle)->Id<<" midpos: "<<midPos[0]<<" "<<midPos[1]<<" "<<midPos[2]<<" d2: "<<d2<<endl;
							//cout<<"Node: "<<id1<<" pos: "<<Nodes[id1]->Position[0]<<" "<<Nodes[id1]->Position[1]<<" "<<Nodes[id1]->Position[2]<<endl;
							//cout<<"Node: "<<id3<<" pos: "<<Nodes[id3]->Position[0]<<" "<<Nodes[id3]->Position[1]<<" "<<Nodes[id3]->Position[2]<<endl;
							if (dx >0) {initialSignsEdgex.push_back(1);}
							else {initialSignsEdgex.push_back(-1);}
							if (dy >0) {initialSignsEdgey.push_back(1);}
							else {initialSignsEdgey.push_back(-1);}
							if (dz >0) {initialSignsEdgez.push_back(1);}
							else {initialSignsEdgez.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightEdgex.push_back(dx/d);
							initialWeightEdgey.push_back(dy/d);
							initialWeightEdgez.push_back(dz/d);
						}
						//get third mid point:
						midPos[0] =	(Nodes[id2]->Position[0] + Nodes[id3]->Position[0])/2.0;
						midPos[1] =	(Nodes[id2]->Position[1] + Nodes[id3]->Position[1])/2.0;
						midPos[2] =	(Nodes[id2]->Position[2] + Nodes[id3]->Position[2])/2.0;
						dx = pos[0] - midPos[0];
						dy = pos[1] - midPos[1];
						dz = pos[2] - midPos[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodeEdgeList0.push_back((*itNode)->Id);
							pacingNodeEdgeList1.push_back(id2);
							pacingNodeEdgeList2.push_back(id3);
							packedToEdge = true;
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id2<<" "<<id3<<" from element: "<<(*itEle)->Id<<" midpos: "<<midPos[0]<<" "<<midPos[1]<<" "<<midPos[2]<<" d2: "<<d2<<endl;
							//cout<<"Node: "<<id2<<" pos: "<<Nodes[id2]->Position[0]<<" "<<Nodes[id2]->Position[1]<<" "<<Nodes[id2]->Position[2]<<endl;
							//cout<<"Node: "<<id3<<" pos: "<<Nodes[id3]->Position[0]<<" "<<Nodes[id3]->Position[1]<<" "<<Nodes[id3]->Position[2]<<endl;
							if (dx >0) {initialSignsEdgex.push_back(1);}
							else {initialSignsEdgex.push_back(-1);}
							if (dy >0) {initialSignsEdgey.push_back(1);}
							else {initialSignsEdgey.push_back(-1);}
							if (dz >0) {initialSignsEdgez.push_back(1);}
							else {initialSignsEdgez.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightEdgex.push_back(dx/d);
							initialWeightEdgey.push_back(dy/d);
							initialWeightEdgez.push_back(dz/d);
						}
					}
					//end of edge detection
					 */
					if (!packedToSurface && !packedToEdge){
						//did not pack to surface or edge, packing to nodes:
						double dx = pos[0] - Nodes[id1]->Position[0];
						double dy = pos[1] - Nodes[id1]->Position[1];
						double dz = pos[2] - Nodes[id1]->Position[2];
						double d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id1);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id1<<" from element: "<<(*itEle)->Id<<endl;
						}
						dx = pos[0] - Nodes[id2]->Position[0];
						dy = pos[1] - Nodes[id2]->Position[1];
						dz = pos[2] - Nodes[id2]->Position[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id2);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id2<<" from element: "<<(*itEle)->Id<<endl;
						}
						dx = pos[0] - Nodes[id3]->Position[0];
						dy = pos[1] - Nodes[id3]->Position[1];
						dz = pos[2] - Nodes[id3]->Position[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id3);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id3<<" from element: "<<(*itEle)->Id<<endl;
						}
					}

				}
			}
			delete[] pos;
		}
	}
	cleanUpPacingCombinations();
}

void Simulation::detectPacingNodes(){
	double periAverageSideLength = 0,colAverageSideLength = 0;
	getAverageSideLength(periAverageSideLength, colAverageSideLength);

	if (thereIsPeripodialMembrane){
		colAverageSideLength = (periAverageSideLength+colAverageSideLength)/2.0;
	}
	packingThreshold = 0.4*colAverageSideLength;  //0.6 was old value
	packingDetectionThreshold = 1.2 * packingThreshold; //0.9 * packingThreshold;
	double t2 = packingDetectionThreshold*packingDetectionThreshold;	//threshold square for rapid calculation
	pacingNodeCouples0.empty();
	pacingNodeCouples1.empty();
	initialWeightPointx.empty();
	initialWeightPointy.empty();
	initialWeightPointz.empty();
	//Added for parallelisation:
	//node based only:
	const int nArray = 16; //the size of array that I will divide the elements into.
	int parallellisationSegmentSize = ceil(nNodes/nArray);
	if (parallellisationSegmentSize<1){
		parallellisationSegmentSize =1;
	}
	const int nSegments = ceil(nNodes/parallellisationSegmentSize); //this does not need to be the same as desired size
	//if there are 5 nodes and I want 16 segments, the actual nSegments will be 5.
	const int segmentBoundariesArraySize = nSegments+1;
	int parallellisationSegmentBoundaries[segmentBoundariesArraySize];
	parallellisationSegmentBoundaries[0] = 0;
	parallellisationSegmentBoundaries[segmentBoundariesArraySize-1] = nNodes;
	for (int i=1; i<segmentBoundariesArraySize-1; ++i){
		parallellisationSegmentBoundaries[i] = i*parallellisationSegmentSize;
	}
	vector <vector<int> > arrayForParallelisationPacingNodeCouples0(nSegments, vector <int>(0));
	vector <vector<int> > arrayForParallelisationPacingNodeCouples1(nSegments, vector <int>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointx(nSegments, vector <double>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointy(nSegments, vector <double>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointz(nSegments, vector <double>(0));
	//parallelise this loop:
	//cout<<" at parallelisation loop for packing, within detectPacingNodes"<<endl;
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	/*cout<<"parallellisationSegmentBoundaries: "<<endl;
	for (int i = 0; i<segmentBoundariesArraySize; ++i){
		cout<<" "<<parallellisationSegmentBoundaries[i]<<" ";
	}
	cout<<endl;*/
	//#pragma omp parallel for //private(arrayForParallelisationPacingNodeCouples0, arrayForParallelisationPacingNodeCouples1, arrayForParallelisationInitialWeightPointx, arrayForParallelisationInitialWeightPointy, arrayForParallelisationInitialWeightPointz)
	for (int a =0; a<nSegments; ++a){
		int initialpoint = parallellisationSegmentBoundaries[a];
		int breakpoint =  parallellisationSegmentBoundaries[a+1];
		//cout<<" a = "<<a<<" initialpoint: "<<initialpoint<<" breakpoint "<<breakpoint<<" nNodes: "<<nNodes<<endl;
		for (vector<Node*>::iterator itNode=Nodes.begin()+ initialpoint; itNode<Nodes.begin()+breakpoint; ++itNode){
			//cout<<" inside loop, for node : "<<(*itNode)->Id<<endl;
			bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
			if (NodeHasPacking){
				double* pos;
				pos = new double[3];
				(*itNode)->getCurrentPosition(pos);
				for (vector<Node*>::iterator itNodeSlave=itNode+1; itNodeSlave<Nodes.end(); ++itNodeSlave){
					bool SlaveNodeHasPacking = (*itNodeSlave)->checkIfNodeHasPacking();
					//bool nodesOnSeperateSurfaces = (*itNodeSlave)->tissueType != (*itNode)->tissueType;
					//if (SlaveNodeHasPacking && nodesOnSeperateSurfaces){
					if (SlaveNodeHasPacking){
						if ((*itNode)->tissuePlacement == (*itNodeSlave)->tissuePlacement){
							//nodes can pack, are they connected?
							bool neigbours = false;
							int n = (*itNode)->immediateNeigs.size();
							for (int i= 0; i<n; ++i){
								if ((*itNode)->immediateNeigs[i] ==(*itNodeSlave)->Id ){
									neigbours = true;
									break;
								}
							}
							if (!neigbours){
								//the nodes can potentially pack, are they close enough?
								double* posSlave;
								posSlave = new double[3];
								(*itNodeSlave)->getCurrentPosition(posSlave);
								double dx = pos[0] - posSlave[0];
								double dy = pos[1] - posSlave[1];
								double dz = pos[2] - posSlave[2];
								double d2 = dx*dx + dy*dy + dz*dz;
								if (d2<t2){
									//close enough for packing , add to list:
									arrayForParallelisationPacingNodeCouples0[a].push_back((*itNode)->Id);
									arrayForParallelisationPacingNodeCouples1[a].push_back((*itNodeSlave)->Id);
									double d = pow (d2,0.5);
									arrayForParallelisationInitialWeightPointx[a].push_back(dx/d);
									arrayForParallelisationInitialWeightPointy[a].push_back(dy/d);
									arrayForParallelisationInitialWeightPointz[a].push_back(dz/d);
								}
								delete[] posSlave;
							}
						}
					}
				}
				delete[] pos;
			}
		}
		//calculate packing here
	}
	//outside parallellisation, need to combine the nodes:
	for (int a =0; a<nSegments-1; ++a){
		int n = arrayForParallelisationPacingNodeCouples0[a].size();
		for (int i=0; i<n; ++i){
			pacingNodeCouples0.push_back(arrayForParallelisationPacingNodeCouples0[a][i]);
			pacingNodeCouples1.push_back(arrayForParallelisationPacingNodeCouples1[a][i]);
			initialWeightPointx.push_back(arrayForParallelisationInitialWeightPointx[a][i]);
			initialWeightPointy.push_back(arrayForParallelisationInitialWeightPointy[a][i]);
			initialWeightPointz.push_back(arrayForParallelisationInitialWeightPointz[a][i]);
		}
	}
}

void Simulation::calculatePacking(){
	//cout<<"in calculate packing: "<<endl;
	double multiplier = 0.3; //just a term to scale the forces down
	double thresholdForNodeDetection = 20; //check every node within 5 microns of distance;
	double threshold = 10;	 //packing forces start at 2 microns
	double thresholdNeg = (-1.0)*threshold;
	double t2 = threshold*threshold;	//threshold square for rapid calculation
	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	//resetting all the normal update flags
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->ApicalNormalForPackingUpToDate = false;
		(*itEle)->BasalNormalForPackingUpToDate = false;
	}
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		/*double sumPack[3] = {0,0,0};
		for (int pp=0; pp<nNodes; ++pp){
			sumPack[0] += PackingForces[pp][0];
			sumPack[1] += PackingForces[pp][1];
			sumPack[2] += PackingForces[pp][2];
		}
		cout<<"checking node: "<<(*itNode)->Id<<" sum of Packing forces: "<<sumPack[0]<<" "<<sumPack[1]<<" "<<sumPack[2]<<endl;
		*/
		//if ((*itNode)->Id == 60) {cout<<" checking node : "<<(*itNode)->Id<<endl;}
		bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
		if (NodeHasPacking){
			//cout<<"	node has packing"<<endl;
			double* pos;
			pos = new double[3];
			(*itNode)->getCurrentPosition(pos);
			vector <int> edgeNodeData0, edgeNodeData1;
			vector <double> edgeNormalDataX, edgeNormalDataY, edgeNormalDataZ;
			bool pushedBySurface = false;
			bool pushedByEdge = false;
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this element
				//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking 1st to element  : "<<(*itEle)->Id<<endl;}
				bool PackingToThisElement =  (*itEle)->checkPackingToThisNodeViaState(TissueHeightDiscretisationLayers, (*itNode));
				if (PackingToThisElement){
					//cout<<"	node may pack to element: "<<(*itEle)->Id<<endl;
					//the node does not belong to the element, and placed correctly. Lets have a preliminary distance check:
					//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking 2nd to element  : "<<(*itEle)->Id<<endl;}
					PackingToThisElement = (*itEle)->IsPointCloseEnoughForPacking((*itNode)->Position, thresholdForNodeDetection, (*itNode)->tissuePlacement);
				}
				if (PackingToThisElement){
					//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" continuing packing to element  : "<<(*itEle)->Id<<endl;}
					//cout<<"	node can pack to element: "<<(*itEle)->Id<<endl;
					//all positions and stated are correct, this node can pack to this element
					double* normalForPacking;
					normalForPacking = new double[3];
					normalForPacking[0]=0.0;normalForPacking[1]=0.0;normalForPacking[2]=0.0;
					double* posCorner;
					posCorner = new double[3];
					//getting the normal to the element surface in the direction of the pushing force, and one corner of the surface
					getNormalAndCornerPosForPacking((*itNode),(*itEle),normalForPacking,posCorner);
					double dInNormalDir = 0.0;
					//if ((*itNode)->Id == 78 ) {cout<<" checking distance to surface for  : "<<(*itEle)->Id<<endl;}
					for (int i=0; i<3; ++i){
						dInNormalDir += (pos[i] - posCorner[i])*normalForPacking[i];
					}
					//cout<<"	distance to surface in direction of normal: "<<dInNormalDir<<endl;
					if ((dInNormalDir > 0 && dInNormalDir < threshold) || (dInNormalDir < 0 && dInNormalDir > thresholdNeg)){
						//if the edge is not composed of any of my neigs, then add to list:
						addToEdgeList((*itNode), (*itEle), edgeNodeData0, edgeNodeData1);
						edgeNormalDataX.push_back(normalForPacking[0]);
						edgeNormalDataY.push_back(normalForPacking[1]);
						edgeNormalDataZ.push_back(normalForPacking[2]);
						//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking projection to surface for  : "<<(*itEle)->Id<<endl;}
						//the node is close enough in the direction of the normal, now I should check if the projection falls within hte triangle:
						double projectedPoint[3];
						for (int i=0; i<3; ++i){
							projectedPoint[i] = pos[i] - dInNormalDir*normalForPacking[i];
						}
						//cout<<"	pos: " <<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
						//cout<<"	normalForPacking: "<<normalForPacking[0]<<" "<<normalForPacking[1]<<" "<<normalForPacking[2]<<endl;
						//cout<<"	projected point: "<<projectedPoint[0]<<" "<<projectedPoint[1]<<" "<<projectedPoint[2]<<endl;
						bool pointInsideTriangle = (*itEle)->IspointInsideTriangle((*itNode)->tissuePlacement,projectedPoint[0],projectedPoint[1],projectedPoint[2]);
						//cout<<"	pointInsideTriangle: "<<pointInsideTriangle<<endl;
						if (pointInsideTriangle){
							//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" calculating force for  : "<<(*itEle)->Id<<endl;}
							//the point is within packing distance, calculate packing force now:
							//if dInNormalDir is negative, then the node penetrated the surface
							// if this is the case, I do not calculate force via distance, I set distance to zero:
							float d2;
							if (dInNormalDir<0){
								d2=0;
							}
							else{
								d2 = dInNormalDir*dInNormalDir;
							}
							double averageMass = 0.5*((*itEle)->VolumePerNode + (*itNode)->mass);
							//double Fmag = averageMass* multiplier * (1.0/d2 - 1.0/t2);
							//double Fmag = averageMass* multiplier * (1.0 - (d2/t2));
							double Fmag = averageMass* multiplier * (1.0 - pow((d2/t2),0.5));
							//cout<<"	Fmag before capping: "<<Fmag<<endl;
							CapPackingForce(Fmag);
							//cout<<"	Fmag after capping: "<<Fmag<<endl;
							double F[3] = {Fmag*normalForPacking[0],  Fmag*normalForPacking[1], Fmag*normalForPacking[2]};
							//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
							//if (PackingForces[126][2]> 0 ){
							//	cout<<" force on 126 became non-zero via node, node Id: "<<(*itNode)->Id<<" prisim: "<<(*itEle)->Id<<endl;
							//}
							for(int i=0; i<3; ++i){
								if ((*itNode)->FixedPos[i]){
									F[i] = 0.0;
									//if the node is fixed on  a certain axis, then do not reflect the force on the surface
								}
							}

							bool allCornersFixed[3] = {false,false,false};
							(*itEle)->AddPackingToSurface((*itNode)->tissuePlacement, F[0],F[1],F[2], PackingForces, Nodes, allCornersFixed[0], allCornersFixed[1], allCornersFixed[2]);
							for(int i=0; i<3; ++i){
								if (!allCornersFixed[i] && !(*itNode)->FixedPos[i]){
									//only adding the force if I could reflect it on a surface and the node is not fixed (already made it zero above, keeping here as a reminder
									PackingForces[(*itNode)->Id][i] += F[i];
								}
							}

							//if (PackingForces[126][2]> 0 ){
							///	cout<<" force on 126 is non-zero after surface addition, node Id: "<<(*itNode)->Id<<" node tissue placement: "<<(*itNode)->tissuePlacement<<" has lateral owner? "<<(*itNode)->hasLateralElementOwner<<" prisim: "<<(*itEle)->Id<<" element tissue type: "<<(*itEle)->tissueType<<endl;
							//}
							pushedBySurface = true;
						}
					}
					delete[] normalForPacking;
					delete[] posCorner;
				}
			}
			int closestNode  = -1000;
			double currMinDist2 = 1000;
			double closestNodeVec[3] = {0.0,0.0,0.0};
			if (!pushedBySurface){
				//find closest node, if any edge is closer than the node, it will pack:
				int n = edgeNodeData0.size();
				bool isNodeChecked[(const int) n];
				bool isEdgeChecked[(const int) n];
				for( int i = 0; i<n; ++i){
					isNodeChecked[i] = false;
					isEdgeChecked[i] = false;
				}
				for( int i = 0; i<n; ++i){
					bool isNeig = (*itNode)->checkIfNeighbour(edgeNodeData0[i]);
					if (isNeig){
						isNodeChecked[i] = true;
						isEdgeChecked[i] = true;
					}
					else{
						//check the other end of the edge:
						isNeig = (*itNode)->checkIfNeighbour(edgeNodeData1[i]);
						if (isNeig){
								isEdgeChecked[i] = true;
						}
					}
				}
				//clear all the nodes and corresponding edges that are my immeditate neigs:

				for (int i=0; i<n; ++i){
					//looping over the nodes
					if (!isNodeChecked[i]){ //do not double check the node if it has been checked previously
						int node0 = edgeNodeData0[i];
						//go over all nodes in the list (after this one) and mark all same id as checked:
						for (int j=i; j<n; ++j){
							if (edgeNodeData0[i] == edgeNodeData0[j]){
								isNodeChecked[j] = true;
							}
						}
						//calculate distace:
						double dx = (*itNode)->Position[0] - Nodes[node0]->Position[0];
						double dy = (*itNode)->Position[1] - Nodes[node0]->Position[1];
						double dz = (*itNode)->Position[2] - Nodes[node0]->Position[2];
						double d2 = dx*dx +  dy*dy  + dz*dz ;
						if (d2 < currMinDist2){
							currMinDist2 = d2;
							closestNode = node0;
							closestNodeVec[0] = dx;
							closestNodeVec[1] = dy;
							closestNodeVec[2] = dz;
						}
					}
				}
				//found closest node
				//if ((*itNode)->Id == 78 ) {cout<<" checking edge packing for node : "<<(*itNode)->Id<<endl;}
				//the node have not been pushed by nay of the surfaces, I need to check if it is being pushed by edges:
				n = edgeNormalDataX.size();
				for (int i=0; i<n; ++i){
					//looping over the normals vectors
					for (int j =0; j<3; ++j){
						//three node couples recorded for each node data point
						if (!isEdgeChecked[3*i+j]){ //do not double check the edge if it has been checked previously
							double normal1[3] = {0.0,0.0,1.0};
							int node0 = edgeNodeData0[3*i+j];
							int node1 = edgeNodeData1[3*i+j];
							double normal0[3] = {edgeNormalDataX[i],edgeNormalDataY[i],edgeNormalDataZ[i]};
							bool secondNormalFound = false;
							isEdgeChecked[3*i+j] = true;
							//I need to start looing at the nodes that are after this triplet, to find the second occurace of this edge:
							for (int k=i+1; k<n; ++k){
								if (!isEdgeChecked[3*k] && (edgeNodeData0[3*k] == node0 || edgeNodeData0[3*k] == node1)){
									if (edgeNodeData1[3*k] == node0 || edgeNodeData1[3*k] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k] = true;
										secondNormalFound = true;
										break;
									}
								}
								if (!isEdgeChecked[3*k+1] && (edgeNodeData0[3*k+1] == node0 || edgeNodeData0[3*k+1] == node1)){
									if (edgeNodeData1[3*k+1] == node0 || edgeNodeData1[3*k+1] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k+1] = true;
										secondNormalFound = true;
										break;
									}
								}
								if (!isEdgeChecked[3*k+2] && (edgeNodeData0[3*k+2] == node0 || edgeNodeData0[3*k+2] == node1)){
									if (edgeNodeData1[3*k+2] == node0 || edgeNodeData1[3*k+2] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k+2] = true;
										secondNormalFound = true;
										break;
									}
								}
							}
							if (!secondNormalFound){
								//the edge is recorded only once, the second time is not recorded.
								//the edge is either on the edge of the tissue, or it is at the symmetricity border
								normal1[0] = normal0[0];
								normal1[2] = normal0[2];
								if (Nodes[node0]->atSymmetricityBorder && Nodes[node1]->atSymmetricityBorder){
									normal1[1] = (-1.0)*normal0[1];
								}
								else{
									normal1[1] = normal0[1];
								}
							}
							//if ((*itNode)->Id == 78 ) {
							//	cout<<" checking edge couple : "<<node0<<" "<<node1<<endl;
							//	cout<<"normals: "<<normal0[0]<<" "<<normal0[1]<<" "<<normal0[2]<<" / "<<normal1[0]<<" "<<normal1[1]<<" "<<normal1[2]<<endl;
							//}
							//now I have normals from both sides for this edge
							double v[3]; //vector from the edge to the point;
							bool isPointOnEdge = checkIfPointIsOnEdge(node0, node1, pos[0], pos[1], pos[2], v[0], v[1], v[2]);
							//if ((*itNode)->Id == 78 ) {
							//	cout<<" is point on the edge? : "<<isPointOnEdge<<endl;
							//}
							if (isPointOnEdge){
								double d2 = v[0]*v[0] +  v[1]*v[1] +  v[2]*v[2];
								//if ((*itNode)->Id == 78 ) {
								//	cout<<" point on the edge, distance squared is : "<<d2<<endl;
								//}
								if (d2 < t2 && d2 <currMinDist2){ // the edge is close enough to pack, and is also closer than the closest node:
									//if ((*itNode)->Id == 78 ) {
									//	cout<<"distance was close enough, nodes were: "<<node0<<" "<<node1<<endl;
									//}
									double d = pow(d2,0.5);
									v[0] /= d;
									v[1] /= d;
									v[2] /= d;
									//check if the direction is correct:
									double dot = (normal0[0] + normal1[0])*v[0] + (normal0[1] + normal1[1])*v[1] + (normal0[2] + normal1[2])*v[2];
									if (dot < 0 ){
										//point penetrated into the elemetns, the pushing force should be inverted:
										v[0] *= (-1.0);
										v[1] *= (-1.0);
										v[2] *= (-1.0);
									}
									//point close enough for packing:
									double averageMass = 1.0/3.0 *( (*itNode)->mass + Nodes[node0]->mass + Nodes[node1]->mass );
									//double Fmag = averageMass* multiplier * (1.0 - (d2/t2));
									double Fmag = averageMass* multiplier * (1.0 - pow((d2/t2),0.5));
									CapPackingForce(Fmag);
									double F[3] = {Fmag*v[0],  Fmag*v[1], Fmag*v[2]};
									//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
									for(int i=0; i<3; ++i){
										/*if (!(*itNode)->FixedPos[i]){
											PackingForces[(*itNode)->Id][i] += F[i];
										}
										//add to the edges
										if (!Nodes[node0]->FixedPos[i]){
											PackingForces[node0][i] -= 0.5*F[i];
										}
										if (!Nodes[node1]->FixedPos[i]){
											PackingForces[node1][i] -= 0.5*F[i];
										}*/
										if (!(*itNode)->FixedPos[i]){
											//the node being checked for packing is not fixed on this axis, may add forces:
											bool edgesFixed[2] = {Nodes[node0]->FixedPos[i], Nodes[node1]->FixedPos[i]};
											if (!edgesFixed[0] && !edgesFixed[1]){
												//both edges  are not fixed, add forces normally:
												PackingForces[(*itNode)->Id][i] += F[i];
												PackingForces[node0][i] -= 0.5*F[i];
												PackingForces[node1][i] -= 0.5*F[i];
											}
											else if(!edgesFixed[0] || !edgesFixed[1]){
												//at least one of the edges are not fixed, find which one, and add all the force on it:
												PackingForces[(*itNode)->Id][i] += F[i];
												if (!edgesFixed[0]){
													PackingForces[node0][i] -= F[i];
												}else {
													PackingForces[node1][i] -= F[i];
												}
											}
											//the other condition is that both edges are fixed, and the node should not be applying force in that direction
										}
									}
									pushedByEdge = true;
								}
							}
						}
					}
				}
			}
			//if ((*itNode)->Id == 78 ) {
			//	cout<<" pushedByEdge? : "<<pushedByEdge<<endl;
			//}
			if(!pushedBySurface && !pushedByEdge){
				//if ((*itNode)->Id == 78 ) {cout<<" checking node based packing for node : "<<(*itNode)->Id<<endl;}
				//checking point to point:
				if (currMinDist2 < t2 ){
					//if ((*itNode)->Id == 211 ) {
					//	cout<<" node 211 is packing to node : "<<closestNode<<endl;
					//}
					double d = pow(currMinDist2,0.5);
					closestNodeVec[0] /= d;
					closestNodeVec[1] /= d;
					closestNodeVec[2] /= d;
					//point close enough for packing:
					double averageMass = 1.0/2.0 *( (*itNode)->mass + Nodes[closestNode]->mass );
					//double Fmag = averageMass* multiplier * (1.0 - (currMinDist2/t2));
					double Fmag = averageMass* multiplier * (1.0 - pow((currMinDist2/t2),0.5));
					CapPackingForce(Fmag);
					double F[3] = {Fmag*closestNodeVec[0],  Fmag*closestNodeVec[1], Fmag*closestNodeVec[2]};
							//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
					for(int i=0; i<3; ++i){
						//if (!(*itNode)->FixedPos[i]){
						//	PackingForces[(*itNode)->Id][i] += F[i];
						//}
						//add to the edges
						//if (!Nodes[closestNode]->FixedPos[i]){
						//	PackingForces[closestNode][i] -= F[i];
						//}
						//if both the packing node and the corresponding node are not fixed:
						if (!(*itNode)->FixedPos[i] && !Nodes[closestNode]->FixedPos[i]){
							PackingForces[closestNode][i] -= F[i];
						}
					}
				}
			}
			//if ((*itNode)->Id == 211 ) {cout<<"node 211: packing to surf: "<<pushedBySurface<<" packing to edge: "<<pushedByEdge<<endl;}
			delete[] pos;
		}
	}
}

bool Simulation::checkIfPointIsOnEdge(int node0, int node1, double x, double y, double z, double& vx, double& vy, double& vz){
	double vec0[3] = {Nodes[node1]->Position[0] - Nodes[node0]->Position[0], Nodes[node1]->Position[1] - Nodes[node0]->Position[1], Nodes[node1]->Position[2]- Nodes[node0]->Position[2]};
	double vec1[3] = {x - Nodes[node0]->Position[0], y - Nodes[node0]->Position[1], z- Nodes[node0]->Position[2]};
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"vec of edge: "<<vec0[0]<<" "<<vec0[1]<<" "<<vec0[2]<<" vec from "<<node0<<" to curr node: " <<vec1[0]<<" "<<vec1[1]<<" "<<vec1[2]<<endl;
	//}
	//get length of edge and normalise edge vector
	double v0Mag = pow(vec0[0]*vec0[0] + vec0[1]*vec0[1]+vec0[2]*vec0[2],0.5);
	vec0[0] /= v0Mag;
	vec0[1] /= v0Mag;
	vec0[2] /= v0Mag;
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"vec of edge ("<<node0<<" - "<<node1<<") after normalisaiton: "<<vec0[0]<<" "<<vec0[1]<<" "<<vec0[1]<<" mag was:  "<<v0Mag<<endl;
	//}
	//dot product of normalised vec0 and vec1 should be positive, and should be smaller than the length of vec0:
	double dot = vec0[0]*vec1[0] + vec0[1]*vec1[1] + vec0[2]*vec1[2];
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"dotp: "<<dot<<endl;
	//}
	if (dot > 0 && dot < v0Mag){
		//the point falls within the edge,
		//now record the normal vector from edge to point:
		vx = vec1[0] - vec0[0]*dot;
		vy = vec1[1] - vec0[1]*dot;
		vz = vec1[2] - vec0[2]*dot;
		return true;
	}
	vx = 0.0;
	vy = 0.0;
	vz = 0.0;
	return false;
}

void Simulation::addToEdgeList(Node* nodePointer, ShapeBase* elementPointer, vector <int> & edgeNodeData0, vector <int> & edgeNodeData1 ){
	int  E0Index = -1, E1Index = -1, E2Index = -1;
	if (nodePointer->tissuePlacement == 0 ){ //checking basal surface,
		if (elementPointer->tissueType == 0){ //element is columnar
			E0Index = 0;
			E1Index = 1;
			E2Index = 2;
		}
		else{//element is periodial:
			E0Index = 3;
			E1Index = 4;
			E2Index = 5;
		}
	}
	else if (nodePointer->tissuePlacement  == 1 ){ //checking apical surface
		if (elementPointer->tissueType == 0){ //element is columnar
			E0Index = 3;
			E1Index = 4;
			E2Index = 5;
		}
		else{//element is periodial:
			E0Index = 0;
			E1Index = 1;
			E2Index = 2;
		}
	}

	edgeNodeData0.push_back(elementPointer->getNodeId(E0Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E1Index));
	edgeNodeData0.push_back(elementPointer->getNodeId(E0Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E2Index));
	edgeNodeData0.push_back(elementPointer->getNodeId(E1Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E2Index));
}

void Simulation::addPackingForces(gsl_matrix* gExt){
	double sumPack[3] = {0.0,0.0,0.0};
	//double sumPackPre[3] = {0.0,0.0,0.0};
	//double sumAv[3] = {0.0,0.0,0.0};
	//double sumgExt[3] = {0.0,0.0,0.0};
	//cout<<"in add packing forces"<<endl;
	for (int j=0; j<nNodes; ++j){
		//cout<<"calculating node: "<<j<<" of "<<nNodes<<endl;
		//double Fx = 0.333* ( PackingForces[j][0] + PackingForcesPreviousStep[j][0] + PackingForcesTwoStepsAgoStep[j][0] );
		//double Fy = 0.333* ( PackingForces[j][1] + PackingForcesPreviousStep[j][1] + PackingForcesTwoStepsAgoStep[j][1] );
		//double Fz = 0.333* ( PackingForces[j][2] + PackingForcesPreviousStep[j][2] + PackingForcesTwoStepsAgoStep[j][2] );
		double Fx = PackingForces[j][0];
		double Fy = PackingForces[j][1];
		double Fz = PackingForces[j][2];
		sumPack[0] += PackingForces[j][0];
		//sumPackPre[0] += PackingForcesPreviousStep[j][0];
		//sumAv[0] += Fx;
		sumPack[1] += PackingForces[j][1];
		//sumPackPre[1] += PackingForcesPreviousStep[j][1];
		//sumAv[1] += Fy;
		sumPack[2] += PackingForces[j][2];
		//sumPackPre[2] += PackingForcesPreviousStep[j][2];
		//sumAv[2] += Fz;
		//cout<<"node: "<<j<<" after summation"<<endl;
		int indice = j*3;
		Fx += gsl_matrix_get(gExt,indice,0);
		gsl_matrix_set(gExt,indice,0,Fx);
		Fy += gsl_matrix_get(gExt,indice+1,0);
		gsl_matrix_set(gExt,indice+1,0,Fy);
		Fz += gsl_matrix_get(gExt,indice+2,0);
		gsl_matrix_set(gExt,indice+2,0,Fz);
		//sumgExt[0] += gsl_matrix_get(gExt,indice,0);
		//sumgExt[1] += gsl_matrix_get(gExt,indice+1,0);
		//sumgExt[2] += gsl_matrix_get(gExt,indice+2,0);
	}
	//cout<<" sum pack: "<<sumPack[0]<<" "<<sumPack[1]<<" "<<sumPack[2]<<endl;
	//cout<<" sum gExt: "<<sumgExt[0]<<" "<<sumgExt[1]<<" "<<sumgExt[2]<<endl;
}

void Simulation::updateElementPositions(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
}

void Simulation::updateElementPositionsSingle(int i ){
	Elements[i]->updatePositions(Nodes);
}

void Simulation::assignTips(){
	double xTips[2] ={ 10000, -10000};
	double yTips[2] ={ 10000, -10000};
	for (int i =0; i<nNodes; ++i){
		if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
			if (Nodes[i]->Position[0] < xTips[0] ){
				dorsalTipIndex = i;
				xTips[0] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] < yTips[0] ){
				anteriorTipIndex = i;
				yTips[0] = Nodes[i]->Position[1];
			}
			if (Nodes[i]->Position[0] > xTips[1] ){
				ventralTipIndex = i;
				xTips[1] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] > yTips[1] ){
				posteriorTipIndex = i;
				yTips[1] = Nodes[i]->Position[1];
			}
		}
		//cout<<"DV Tip node indexes: "<<dorsalTipIndex<<" "<<ventralTipIndex<<endl;
		//cout<<"AP Tip node indexes: "<<anteriorTipIndex<<" "<<posteriorTipIndex<<endl;
	}
	//I have identified the tips assuming the tissue is represented in full.
	//If there is symmetry in the input, then I will need to correct this.
	//In the case of x symmetry, the minimum x will be correct in DV, but not the maximum x.
	//In the case of y symmetry, the maximum y will be correct in AP. but not the minimum.
	if (symmetricX){
		double delta = 0.02;
		//the ventralTipIndex is correct, it is at the minimum of x axis. The dorsalTipIndex is at x=0. but there are a series of nodes at x=0,
		//and the selected one should be at y=0 too.
		for (int i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
							dorsalTipIndex = i;
							break;
						}
					}
				}
		}
	}
	if (symmetricY){
		double delta = 0.02;
		//the posteriorTipIndex is correct, it is at the maximum of y axis. The anteriorTipIndex is at y=0. but there are a series of nodes at y=0,
		//and the selected one should be at x=0 too.
		for (int i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
							anteriorTipIndex = i;
							break;
						}
					}
				}
		}
	}
	cout<<"Tip node indexes - dorsalTipIndex: "<<dorsalTipIndex<<" ventralTipIndex: "<<ventralTipIndex<<endl;
	cout<<"Tip node indexes - anteriorTipIndex: "<<anteriorTipIndex<<" posteriorTipIndex: "<<posteriorTipIndex<<endl;
}

void Simulation::alignTissueDVToXPositive(){
	if (!symmetricX){
		double* u = new double[3];
		double* v = new double[3];
		//For simulations with no external viscosity, the position of Dorsal tip is fixed, Ventral tip is fixed in y and z, another node is fixed in z
		for (int i=0;i<3;++i){
			u[i] = Nodes[ventralTipIndex]->Position[i] - Nodes[dorsalTipIndex]->Position[i];
		}
		Elements[0]->normaliseVector3D(u);
		v[0]=1;v[1]=0;v[2]=0;
		double c, s;
		Elements[0]->calculateRotationAngleSinCos(u,v,c,s);
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		Elements[0]->calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
		Elements[0]->constructRotationMatrix(c,s,rotAx,rotMat);
		for(int i=1;i<nNodes;++i){
			for(int j = 0; j< Nodes[i]->nDim; ++j){
				u[j] = Nodes[i]->Position[j] - Nodes[0]->Position[j];
			}
			Elements[0]->rotateVectorByRotationMatrix(u,rotMat);

			for(int j = 0; j< Nodes[i]->nDim; ++j){
				Nodes[i]->Position[j] = u[j] + Nodes[0]->Position[j];
			}
		}
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			(*itElement)->updatePositions(Nodes);
		}
		delete[] rotAx;
		delete[] rotMat;
		delete[] u;
		delete[] v;
	}
}



void Simulation::alignTissueAPToXYPlane(){
	double* u = new double[3];
	double* v = new double[3];
	for (int i=0;i<3;++i){
		u[i] = Nodes[anteriorTipIndex]->Position[i] - Nodes[posteriorTipIndex]->Position[i];
	}
	Elements[0]->normaliseVector3D(u);
	v[0]=u[0];v[1]=u[1];v[2]=0;
	Elements[0]->normaliseVector3D(v);
	double c, s;
	Elements[0]->calculateRotationAngleSinCos(u,v,c,s);
	double *rotAx;
	rotAx = new double[3];
	double *rotMat;
	rotMat = new double[9]; //matrix is written in one row
	Elements[0]->calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
	Elements[0]->constructRotationMatrix(c,s,rotAx,rotMat);
	for(int i=0;i<nNodes;++i){
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			u[j] = Nodes[i]->Position[j];
		}
		Elements[0]->rotateVectorByRotationMatrix(u,rotMat);
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			Nodes[i]->Position[j] = u[j];
		}
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
	delete[] rotAx;
	delete[] rotMat;
	delete[] u;
	delete[] v;
}

void Simulation::calculateDVDistance(){
	double d[3];
	if (symmetricX){
		for (int i=0;i<3;++i){
			d[i] = Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[dorsalTipIndex]->Position[i] - Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	double dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	outputFile<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
	cout<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
	if (symmetricY){
		for (int i=0;i<3;++i){
			d[i] = Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[anteriorTipIndex]->Position[i] - Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	outputFile<<" AP distance is: "<<dmag<<endl;
	cout<<" AP distance is: "<<dmag<<endl;
	//cout<<"time: "<<currSimTimeSec<<" Node 0 position: "<<Nodes[0]->Position[0]<<" "<<Nodes[0]->Position[1]<<" "<<Nodes[0]->Position[2]<<endl;
	//cout<<"TIPS: "<<dorsalTipIndex<<" "<<ventralTipIndex<<" "<<anteriorTipIndex<<" "<<posteriorTipIndex<<endl;
}

void Simulation::resetForces(bool resetPacking){
	int dim = 3;
	//n nodes
	for (int j=0;j<nNodes;++j){
		//3 dimensions
		for (int k=0;k<dim;++k){
			SystemForces[j][k]=0.0;
			//packing forces should be reset at the beginning of the step, but not during NR iterations
			if (resetPacking){
				//PackingForcesTwoStepsAgoStep[j][k] = PackingForcesPreviousStep[j][k];
				//PackingForcesPreviousStep[j][k] = PackingForces[j][k];
				PackingForces[j][k]=0.0;
			}
			if (recordForcesOnFixedNodes){
				FixedNodeForces[j][k]=0.0;
			}
		}
	}
}


void Simulation::calculateApicalSize(){
	double sizeLim[2][2] = {{0.0,0.0},{0.0,0.0}};
	bool found[2][2] = {{false,false},{false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 0 && (*itNode)->tissuePlacement == 1){ //checking only columnar apical layer nodes
			for (int i=0; i<2; ++i){
				if ( (*itNode)->Position[i] < sizeLim[0][i] ){
					sizeLim[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>sizeLim[1][i]){
					sizeLim[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	if (!found[0][0] && !found[0][1] && !found[1][0] && !found[1][1]){
		cerr<<" error in apical bounding  box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[1][0]<<" "<<found[1][1]<<endl;
	}
	double DV = sizeLim[1][0] - sizeLim[0][0];
	double AP = sizeLim[1][1] - sizeLim[0][1];
	outputFile<<"At time: "<<currSimTimeSec<<" apical bounding box size: "<<DV<<" "<<AP<<endl;
}

void Simulation::calculateBoundingBox(){
	boundingBox[0][0] =  100000.0;	//lower left x
	boundingBox[0][1] =  100000.0;	//lower left y
	boundingBox[0][2] =  100000.0;	//lower z
	boundingBox[1][0] = -100000.0;	//upper right x
	boundingBox[1][1] = -100000.0;	//upper right y
	boundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		//if (!(*itNode)->allOwnersAblated){
			//There is at least one element owning this node that is not ablated
		//if ((*itNode)->tissueType == 0 || (*itNode)->tissueType==1) {	//only consider peripodial or columnar nodes, not the linkers
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < boundingBox[0][i] ){
					boundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>boundingBox[1][i]){
					boundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		//}
		//}
	}
	//cout<<"calculating bounding box, symmetricity x & y: "<<symmetricX<<" "<<symmetricY<<endl;
	//cout<<"bounding box before update: "<<boundingBox[0][0]<<" "<<boundingBox[0][1]<<" "<<boundingBox[1][0]<<" "<<boundingBox[1][1]<<endl;
	if (symmetricY){
		boundingBox[0][1] = (-1.0)*boundingBox[1][1]; //if there is Y symmetricity, then the bounding box is extended to double the size in y
	}
	cout<<"bounding box after update: "<<boundingBox[0][0]<<" "<<boundingBox[0][1]<<" "<<boundingBox[1][0]<<" "<<boundingBox[1][1]<<endl;
	for (int i=0; i<3; ++i){
		boundingBoxSize[i] = boundingBox[1][i] - boundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}


void Simulation::bringMyosinStimuliUpToDate(){
	//need to go throug all myosin functions and apply them to here:
	//for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	//	(*itElement)->bringMyosinStimuliUpToDate();
	//}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->adjustCMyosinFromSave();
	}
}

/*
void Simulation::calculateColumnarLayerBoundingBox(){
	columnarBoundingBox[0][0] =  100000.0;	//lower left x
	columnarBoundingBox[0][1] =  100000.0;	//lower left y
	columnarBoundingBox[0][2] =  100000.0;	//lower z
	columnarBoundingBox[1][0] = -100000.0;	//upper right x
	columnarBoundingBox[1][1] = -100000.0;	//upper right y
	columnarBoundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 0){ //checking only columnar layer nodes
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < columnarBoundingBox[0][i] ){
					columnarBoundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>columnarBoundingBox[1][i]){
					columnarBoundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	for (int i=0; i<3; ++i){
		columnarBoundingBoxSize[i] = columnarBoundingBox[1][i] - columnarBoundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}

void Simulation::calculatePeripodialBoundingBox(){
	peripodialBoundingBox[0][0] =  100000.0;	//lower left x
	peripodialBoundingBox[0][1] =  100000.0;	//lower left y
	peripodialBoundingBox[0][2] =  100000.0;	//lower z
	peripodialBoundingBox[1][0] = -100000.0;	//upper right x
	peripodialBoundingBox[1][1] = -100000.0;	//upper right y
	peripodialBoundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType != 0){ //checking any node that is not columnar
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < peripodialBoundingBox[0][i] ){
					peripodialBoundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>peripodialBoundingBox[1][i]){
					peripodialBoundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	for (int i=0; i<3; ++i){
		peripodialBoundingBoxSize[i] = peripodialBoundingBox[1][i] - peripodialBoundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}
*/
void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<currSimTimeSec<<" sec ("<<currSimTimeSec/3600<<" hr )"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
    writeGrowth();
	writeForces();
	writePacking();
	writeProteins();
	writePhysicalProp();
}

void Simulation::writeSaveFileStepHeader(){
    saveFileMesh<<"=============== TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
	saveFileMesh<<"==================================================="<<endl;
}

void Simulation::writeSaveFileStepFooter(){
	saveFileMesh<<"=============== END OF TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
	saveFileMesh<<"============================================"<<endl;
}

void Simulation::writeNodes(){
	saveFileMesh<<nNodes<<endl;
	for (int i = 0; i<nNodes; ++i){
		if(Nodes[i]->nDim==2){
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[0];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[1];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<0.0;

		}
		else{
			int ndim = Nodes[i]->nDim;
			for (int j=0;j<ndim;++j){
				saveFileMesh.precision(10);saveFileMesh.width(20);
				saveFileMesh<<Nodes[i]->Position[j];
			}

		}
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissuePlacement;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissueType;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->atCircumference;
		saveFileMesh<<endl;
	}
}

void Simulation::writeElements(){
	saveFileMesh<<nElements<<endl;
	for (int i = 0; i<nElements; ++i){
		int shapetype = Elements[i]->getShapeType();
		saveFileMesh.width(4);
		saveFileMesh<<shapetype;
		if(shapetype == 4 ){
			double height = Elements[i]->getElementHeight();
			saveFileMesh.precision(5);saveFileMesh.width(12);
			saveFileMesh<<height;
		}
		saveFileMesh.width(4);
		saveFileMesh<<Elements[i]->IsAblated;
		int nodeNumber = Elements[i]->getNodeNumber();
		int*  NodeIds = Elements[i]->getNodeIds();
		for (int j = 0; j<nodeNumber; ++j ){
			saveFileMesh.width(9);saveFileMesh<<NodeIds[j];
		}
		int dim  = Elements[i]->getDim();
		double** refPos = Elements[i]->getReferencePos();
		for (int j = 0; j<nodeNumber; ++j ){
			for (int k = 0; k<dim; ++k ){
				saveFileMesh.precision(4);saveFileMesh.width(15);
				saveFileMesh<<refPos[j][k];
			}
		}
		saveFileMesh<<endl;
	}
}

void Simulation::writeTensionCompression(){
	//for (int i=0;i<6;++i){
	//	cout<<" at timestep :"<< timestep<<" the plastic strains of element 0:	"<<Elements[0]->PlasticStrain(i)<<"	normal strain: 	"<<Elements[i]->Strain(0)<<endl;
	//}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        for (int j=0; j<6; ++j){
            double S = gsl_matrix_get((*itElement)->Strain,j,0);
            saveFileTensionCompression.write((char*) &S, sizeof S);
        }
    }
	saveFileTensionCompression.flush();
}

void Simulation::writeGrowth(){
    //for (int i=0;i<6;++i){
    //	cout<<" at timestep :"<< timestep<<" the plastic strains of element 0:	"<<Elements[0]->PlasticStrain(i)<<"	normal strain: 	"<<Elements[i]->Strain(0)<<endl;
    //}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        gsl_matrix* currFg = (*itElement)->getFg();
        double* growthRate = (*itElement)->getGrowthRate();
        for (int j=0; j<3; ++j){
            for (int k=0; k<3; ++k){
                double Fgjk = gsl_matrix_get(currFg,j,k);
                saveFileGrowth.write((char*) &Fgjk, sizeof Fgjk);
            }
            saveFileGrowthRate.write((char*) &growthRate[j], sizeof growthRate[j]);
        }
        gsl_matrix_free(currFg);
    }
    saveFileGrowth.flush();
    saveFileGrowthRate.flush();
}



void Simulation::writeForces(){
	//Write system forces first, on a nodal basis
	for (int i=0;i<nNodes;++i){
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileForces.write((char*) &SystemForces[i][1], sizeof SystemForces[i][1]);
		saveFileForces.write((char*) &SystemForces[i][2], sizeof SystemForces[i][2]);
	}
	//Then write myosin forces
	for (int i=0;i<nElements;++i){
		int n=Elements[i]->getNodeNumber();
		for (int j=0; j<n; ++j){
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][0], sizeof Elements[i]->MyoForce[j][0]);
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][1], sizeof Elements[i]->MyoForce[j][1]);
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][2], sizeof Elements[i]->MyoForce[j][2]);
		}
	}
	saveFileForces.flush();
}

void Simulation::writePacking(){
	int n = pacingNodeCouples0.size();
	vector <int>::iterator itId0;
	vector <int>::iterator itId1;
	itId1=pacingNodeCouples1.begin();
	saveFilePacking.write((char*) &n, sizeof n);
	for (itId0 = pacingNodeCouples0.begin(); itId0 < pacingNodeCouples0.end(); ++itId0){
		saveFilePacking.write((char*) &(*itId0), sizeof (*itId0));
		saveFilePacking.write((char*) &(*itId1), sizeof (*itId1));
		itId1++;
	}
	for (int i=0;i<n;++i){
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][0], sizeof PackingForces[pacingNodeCouples0[i]][0]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][1], sizeof PackingForces[pacingNodeCouples0[i]][1]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][2], sizeof PackingForces[pacingNodeCouples0[i]][2]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][0], sizeof PackingForces[pacingNodeCouples1[i]][0]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][1], sizeof PackingForces[pacingNodeCouples1[i]][1]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][2], sizeof PackingForces[pacingNodeCouples1[i]][2]);
	}
	saveFilePacking.flush();
}

void Simulation::writeProteins(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double* cMyo = new double[4];
		double* cMyoEq = new double[4];
		(*itElement)->getMyosinLevels(cMyo);
		(*itElement)->getEquilibriumMyosinLevels(cMyoEq);
		saveFileProteins.write((char*) &cMyo[0], sizeof cMyo[0]);
		saveFileProteins.write((char*) &cMyo[1], sizeof cMyo[1]);
		saveFileProteins.write((char*) &cMyo[2], sizeof cMyo[2]);
		saveFileProteins.write((char*) &cMyo[3], sizeof cMyo[3]);
		saveFileProteins.write((char*) &cMyoEq[0], sizeof cMyoEq[0]);
		saveFileProteins.write((char*) &cMyoEq[1], sizeof cMyoEq[1]);
		saveFileProteins.write((char*) &cMyoEq[2], sizeof cMyoEq[2]);
		saveFileProteins.write((char*) &cMyoEq[3], sizeof cMyoEq[3]);
		delete[] cMyo;
		delete[] cMyoEq;
	}
	saveFileProteins.flush();
}

void Simulation::writePhysicalProp(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double E = (*itElement)->getYoungModulus();
		double interalVisc	= (*itElement)->getInternalViscosity();
		double zRemodellingSoFar = (*itElement)->getZRemodellingSoFar();
		saveFilePhysicalProp.write((char*) &E, sizeof E);
		saveFilePhysicalProp.write((char*) &interalVisc, sizeof interalVisc);
		saveFilePhysicalProp.write((char*) &zRemodellingSoFar, sizeof zRemodellingSoFar);
	}
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[0], sizeof (*itNode)->externalViscosity[0]);
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[1], sizeof (*itNode)->externalViscosity[1]);
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[2], sizeof (*itNode)->externalViscosity[2]);
	}
	saveFilePhysicalProp.flush();
}

void Simulation::calculateMyosinForces(){
	//cout<<"Entered calculateMyosinForces"<<endl;
	cleanUpMyosinForces();
	bool basedOnArea = false;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updateMyosinConcentration(dt, kMyo, thereIsMyosinFeedback, MyosinFeedbackCap);
		if (basedOnArea){
			//The forces are based on the total apical/basal area of the element, it may not fluctuate with rapid shape changes driven by myosin itself.
			(*itElement)->calculateMyosinForcesAreaBased(forcePerMyoMolecule);
		}
		else{
			//The forces are based on the total size of the element, it does not fluctuate with rapid shape changes, based on the grown volume
			(*itElement)->calculateMyosinForcesTotalSizeBased(forcePerMyoMolecule);
		}
	}
}

void Simulation::cleanUpMyosinForces(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->cleanMyosinForce();
	}
}

void Simulation::checkForMyosinUpdates(){
	for (int i=0; i<nMyosinFunctions; ++i){
		if(currSimTimeSec == myosinFunctions[i]->initTime  ){ //the application time of the signal is given in seconds.
			updateEquilibriumMyosinsFromInputSignal(myosinFunctions[i]);
		}
	}
}


void Simulation::updateEquilibriumMyosinsFromInputSignal(MyosinFunction* currMF){
	//cout<<"inside updateEquilibriumMyosinsFromInputSignal "<<endl;
	int nGridX = currMF->getGridX();
	int nGridY = currMF->getGridY();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((currMF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currMF->applyToPeripodialMembrane && (*itElement)->tissueType == 1) ){//|| (*itElement)->tissueType == 2){
			if ((*itElement)->spansWholeTissue || (currMF->isApical && (*itElement)->tissuePlacement == 1) || (!currMF->isApical && (*itElement)->tissuePlacement == 0)){
				// 1) The element spans the whole tissue therefore both apical and basal responses should be applied
				// 2) The myosin response is applicable to apical surface, and the tissue placement of the element is apical,
				// 3) The myosin response is applicable to basal surface, and the tissue placement of the lement is basal.
				if (currMF->manualStripes){
					//Here, I need to check if all the nodes of the element fall into the stripe of myosin up-regulations.
					//double stripeSize1, stripeSize2;
					//double initialPoint, endPoint;
					bool inActiveZone = (*itElement)->calculateIfInsideActiveStripe(currMF->initialPoint,currMF->endPoint,currMF->stripeSize1,currMF->stripeSize2);
					if (inActiveZone){
						if (currMF->isPolarised){
							(*itElement)->updateUnipolarEquilibriumMyosinConcentration(currMF->isApical,currMF->manualCMyoEq, currMF->manualOrientation[0], currMF->manualOrientation[1]);
						}
						else{
							(*itElement)->updateUniformEquilibriumMyosinConcentration(currMF->isApical, currMF->manualCMyoEq);
						}
					}
				}
				else{
					//calculating the grid indices:
					double* ReletivePos = new double[2];
					//normalising the element centre position with bounding box
					(*itElement)->getRelativePosInBoundingBox(ReletivePos);
					/*if ((*itElement)->tissueType == 0){
						(*itElement)->getRelativePosInColumnarBoundingBox(ReletivePos);
					}
					else{
						(*itElement)->getRelativePosInPeripodialBoundingBox(ReletivePos);
					}*/
					int indexX, indexY;
					double fracX, fracY;
					(*itElement)->convertRelativePosToGridIndex(ReletivePos, indexX, indexY, fracX, fracY, nGridX, nGridY);
					//reading the equilibrium myosin value
					double cEqYmid[2]= {0.0,0.0};
					double cEq;
					cEqYmid[0] = currMF->getEquilibriumMyoMatrixElement(indexX,indexY)*(1.0-fracX) + currMF->getEquilibriumMyoMatrixElement(indexX+1,indexY)*fracX;
					cEqYmid[1] = currMF->getEquilibriumMyoMatrixElement(indexX,indexY+1)*(1.0-fracX) + currMF->getEquilibriumMyoMatrixElement(indexX+1,indexY+1)*fracX;
					cEq = cEqYmid[0]*(1.0-fracY) + cEqYmid[1]*fracY;
					if (currMF->isPolarised){
						//If the function is polarised, reading the orientation
						double orientation[2];
						for (int axis=0; axis<2; ++axis){
							cEqYmid[0] = currMF->getOrientationMatrixElement(indexX,indexY,axis)*(1.0-fracX) + currMF->getOrientationMatrixElement(indexX+1,indexY,axis)*fracX;
							cEqYmid[1] = currMF->getOrientationMatrixElement(indexX,indexY+1,axis)*(1.0-fracX) + currMF->getOrientationMatrixElement(indexX+1,indexY+1,axis)*fracX;
							orientation[axis] = cEqYmid[0]*(1.0-fracY) + cEqYmid[1]*fracY;
					}
					//updating the values of the shape for polarised myosin:
					(*itElement)->updateUnipolarEquilibriumMyosinConcentration(currMF->isApical, cEq, orientation[0], orientation[1]);
					}
					else{
						//updating the values of the shape for uniform contractile:
						(*itElement)->updateUniformEquilibriumMyosinConcentration(currMF->isApical, cEq);
					}
					delete[] ReletivePos;
				}
			}
		}
	}
	//cout<<"finalised updateEquilibriumMyosinsFromInputSignal "<<endl;
}


void Simulation::calculateGrowth(){
	//cout<<"Calculating Growth"<<endl;
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctions[i]->Type == 1){
			//cout<<"Calculating Uniform Growth, function: "<<i<<endl;
			calculateGrowthUniform(GrowthFunctions[i]);
		}
		else if(GrowthFunctions[i]->Type == 2){
			calculateGrowthRing(GrowthFunctions[i]);
		}
		else if(GrowthFunctions[i]->Type == 3){
			calculateGrowthGridBased(GrowthFunctions[i]);
		}
	}
	//calculateGrowthGridBased();
}

void Simulation::calculateShapeChange(){
	cleanUpShapeChangeRates();
	for (int i=0; i<nShapeChangeFunctions; ++i){
		if (GrowthFunctions[i]->Type == 1){
			//cout<<"Calculating Uniform Growth"<<endl;
			calculateShapeChangeUniform(ShapeChangeFunctions[i]);
		}
	}
}


void Simulation::cleanUpGrowthRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setGrowthRate(dt,0.0,0.0,0.0);
		(*itElement)->updateGrowthIncrementFromRate();
	}
}

void Simulation::cleanUpShapeChangeRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setShapeChangeRate(0.0,0.0,0.0,0.0,0.0,0.0);
	}
}

void Simulation::calculateShapeChangeUniform (GrowthFunctionBase* currSCF){
	if(currSimTimeSec >= currSCF->initTime && currSimTimeSec < currSCF->endTime ){
		//cout<<"calculating shape change uniform"<<endl;
		double *maxValues;
        maxValues = new double[3];
        currSCF->getGrowthRate(maxValues);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			//tissue type == 0 is columnar layer, ==1 is peripodial membrane, ==2 id linker zone
			if ( currSCF->applyToColumnarLayer){
				if ((*itElement)->tissueType == 0){ //columnar layer, grow directly
					(*itElement)->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getColumnarness();
					(*itElement)->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
			if ( currSCF->applyToPeripodialMembrane){
				if ((*itElement)->tissueType == 1){ //peripodial membrane, grow directly
					(*itElement)->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getPeripodialness();
					(*itElement)->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
		}
		delete[] maxValues;
	}
}

void Simulation::calculateGrowthUniform(GrowthFunctionBase* currGF){
	//cout<<"inside uniform growth function, initTime: "<<currGF->initTime <<" endtime: "<<currGF->endTime<<" simTime"<<simTime<<endl;
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
		double *growthRates = new double[3];
		currGF->getGrowthRate(growthRates);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    		if (!thereIsExplicitECM || !(*itElement)->isECMMimicing){
    			//tissue type == 0 is columnar layer, ==1 is peripodial membrane, ==2 id linker zone
				gsl_matrix_set_identity(columnarFgIncrement);
				gsl_matrix_set_identity(peripodialFgIncrement);
				if (currGF->applyToColumnarLayer){
					(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
				}
				if (currGF->applyToPeripodialMembrane){
					(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
				}
				(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
    		}
    	}
		delete[] growthRates;
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::calculateGrowthRing(GrowthFunctionBase* currGF){
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		//The growth function is active at current time, now I will grow the elements.
		//First get the remaining data from the growth function parameters
		float centre[2];
		currGF->getCentre(centre[0], centre[1]);
		float innerRadius = currGF->getInnerRadius();
		float outerRadius = currGF->getOuterRadius();
		double* maxValues;
        maxValues = new double[3];
		currGF->getGrowthRate(maxValues);
		float innerRadius2 = innerRadius*innerRadius;
		float outerRadius2 = outerRadius*outerRadius;
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    		if (!thereIsExplicitECM || !(*itElement)->isECMMimicing){
				double* Elementcentre = new double[3];
				Elementcentre = (*itElement)->getCentre();
				//the distance is calculated in the x-y projection
				double d[2] = {centre[0] - Elementcentre[0], centre[1] - Elementcentre[1]};
				double dmag2 = d[0]*d[0] + d[1]*d[1];
				if (dmag2 > innerRadius2 && dmag2 < outerRadius2){
					//the element is within the growth zone.
					float distance = pow(dmag2,0.5);
					//calculating the growth rate: as a fraction increase within this time point
					double sf = (1.0 - (distance - innerRadius) / (outerRadius - innerRadius) );
					double growthscale[3] = {maxValues[0]*sf,maxValues[1]*sf,maxValues[2]*sf};
					gsl_matrix_set_identity(columnarFgIncrement);
					gsl_matrix_set_identity(peripodialFgIncrement);
					if (currGF->applyToColumnarLayer){
						(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
					}
					if (currGF->applyToPeripodialMembrane){
						(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
					}
					(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
				}
				delete[] Elementcentre;
    		}
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
		delete[] maxValues;
	}
}

void Simulation::setUpECMMimicingElements(){
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->tissuePlacement == 0 ){
			//basal elements are mimicing the ECM:
			(*itElement)->setECMMimicing(true);
		}
		if ( (*itElement)->tissuePlacement == 2 && (*itElement)->spansWholeTissue ){
			//elemetns that  are called mid -line, as they
			//span the whole tissue, are treated as ECM mimicing
			(*itElement)->setECMMimicing(true);
		}
	}
}

void Simulation::setUpActinMimicingElements(){
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->tissuePlacement == 1 ){
			//apical elements are mimicing actin:
			(*itElement)->setActinMimicing(true);
		}
		if ( (*itElement)->tissuePlacement == 2 && (*itElement)->spansWholeTissue ){
			//elemetns that  are called mid -line, as they
			//span the whole tissue, are treated as ECM mimicing
			(*itElement)->setActinMimicing(true);
		}
	}
}

void Simulation::calculateGrowthGridBased(GrowthFunctionBase* currGF){
	int nGridX = currGF->getGridX();
	int nGridY = currGF->getGridY();
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
		for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			gsl_matrix_set_identity(columnarFgIncrement);
			gsl_matrix_set_identity(peripodialFgIncrement);

			int IndexX = 0.0, IndexY = 0.0;
			double FracX = 1.0,  FracY = 1.0;
			//(*itElement)->getTissuePositionWeigths(columnarnessWeight, peripodialnessWeight);
			if (GridGrowthsPinnedOnInitialMesh){
				(*itElement)->getInitialRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
			else{
				(*itElement)->getRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
			if (!thereIsExplicitECM || !(*itElement)->isECMMimicing){
				//There is either no explicit ECM definition, or the element is not ECM mimicing.
				//If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
				//If there is no explicit ecm, then all should proceed as usual.
				if (currGF->applyToColumnarLayer){
					(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
				}
				if (currGF->applyToPeripodialMembrane){
					(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, peripodialFgIncrement, 1, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 1 for peripodial membrane
				}
			}
			(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::TissueAxisPositionDisplay(){
	cerr<<"DV border: "<<endl;
	for (int i=0;i<nNodes;++i){
		double x= Nodes[i]->Position[0];
		if (x < 0.2 && x > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
	cerr<<"AP border: "<<endl;
	for (int i=0;i<nNodes;++i){
		double y= Nodes[i]->Position[1];
		if (y < 0.2 && y > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
}

void Simulation::coordinateDisplay(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		int type =(*itElement)-> getShapeType();
		if(type ==1){
			for(int j=0;j<3;++j){
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		int type =(*itElement)-> getShapeType();
		if(type ==1){
			for(int j=3;j<6;++j){
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
}

void Simulation::setStretch(){
	cout<<"setting the stretcher"<<endl;
	recordForcesOnFixedNodes = true;
	for (int i=0;i<3;i++){
		leftClampForces[i] = 0.0;
		rightClampForces[i] = 0.0;
	}
	vector <int> clampedNodeIds;
	double distance = 0;


	if (DVClamp){
		distance = fabs(Nodes[ventralTipIndex]->Position[0] - Nodes[dorsalTipIndex]->Position[0]);
		cerr<<"Total DV distance: "<<distance<<" ";
		distanceIndex = 0; //if the clamp is on DV axis, then the direction of interest is x, index is 0;
	}
	else{
		distance = fabs(Nodes[anteriorTipIndex]->Position[1] - Nodes[posteriorTipIndex]->Position[1]);
		cerr<<"Total AP distance: "<<distance<<" ";
		distanceIndex = 1; //if the clamp is on AP axis, then the direction of interest is y, index is 1.
	}
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[distanceIndex]> StretchMax || Nodes[i]->Position[distanceIndex] < StretchMin){
			Nodes[i]->FixedPos[0]=1;
			Nodes[i]->FixedPos[1]=1;
			Nodes[i]->FixedPos[2]=1;
			clampedNodeIds.push_back(Nodes[i]->Id);
		}
	}
	setUpClampBorders(clampedNodeIds);
	//the distance that is to be moved:
	distance *= StretchStrain;
	cerr<<"the distance that is to be moved: "<<distance<<" ";
	//the time steps that the stretch operation should take place in:
	double StretchTimeSteps = (StretchEndTime - StretchInitialTime)/dt;
	cerr<<"stretchTimeSteps: "<<StretchTimeSteps<<" ";
	StretchDistanceStep = 0.5* (distance / StretchTimeSteps);
	cerr<<"StretchDistanceStep: "<<StretchDistanceStep<<endl;
}

void Simulation::setUpClampBorders(vector<int>& clampedNodeIds){
	int* nodeIds;
	int n = clampedNodeIds.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		nodeIds = (*itElement)->getNodeIds();
		bool hasClampedNode = false;
		bool hasNonClampedNode = false;
		bool leftHandSide = false;
		vector <int> clampedBorderNodes;
		int nNodes = (*itElement)->getNodeNumber();
		for (int j=0; j< nNodes; j++){
			bool lastNodeWasClamped = false;
			for (int k=0; k<n; k++){
				if (nodeIds[j] == clampedNodeIds[k]){
					hasClampedNode = true;
					lastNodeWasClamped = true;
					//check if node is recorded:
					bool alreadyRecorded = false;
					int nLeftClamp = leftClampBorder.size();
					for (int m=0; m<nLeftClamp; m++){
						if (leftClampBorder[m] == nodeIds[j]){
							alreadyRecorded = true;
							break;
						}
					}
					if (!alreadyRecorded){
						int nRightClamp = rightClampBorder.size();
						for (int m=0; m<nRightClamp; m++){
							if (rightClampBorder[m] == nodeIds[j]){
								alreadyRecorded = true;
								break;
							}
						}
					}
					if (!alreadyRecorded){
						clampedBorderNodes.push_back(nodeIds[j]);
					}
					if (Nodes[nodeIds[j]]->Position[distanceIndex] < StretchMin){
						leftHandSide = true;
					}
					break;
				}
			}
			if (lastNodeWasClamped == false){
				hasNonClampedNode = true;
			}
		}
		if (hasClampedNode && hasNonClampedNode){
			cout<<" Element "<<(*itElement)->Id<<" is at the border"<<endl;
			int nClamp = clampedBorderNodes.size();
			if(leftHandSide){
				cout<<"Element is on left hand side"<<endl;
				for (int k=0; k<nClamp; k++){
					leftClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
			else {
				cout<<"Element is on right hand side"<<endl;
				for (int k=0; k<nClamp; k++){
					rightClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
		}
	}int nLeftClamp = leftClampBorder.size();
	for (int k=0; k<nLeftClamp; k++){
		cout<<"left clamp border nodes: "<<leftClampBorder[k]<<endl;
	}
	int nRightClamp = rightClampBorder.size();
	for (int k=0; k<nRightClamp; k++){
		cout<<"right clamp border nodes: "<<rightClampBorder[k]<<endl;
	}
}

void Simulation::recordForcesOnClampBorders(){
	if (recordForcesOnFixedNodes){
		for (int i=0;i<3;i++){
			leftClampForces[i] = 0.0;
			rightClampForces[i] = 0.0;
		}
		int nLeftClamp = leftClampBorder.size();
		for (int k=0; k<nLeftClamp; k++){
			for (int j=0; j<3; ++j){
				leftClampForces[j] += FixedNodeForces[leftClampBorder[k]][j];
			}
		}
		int nRightClamp = rightClampBorder.size();
		for (int k=0; k<nRightClamp; k++){
			for (int j=0; j<3; ++j){
				rightClampForces[j] += FixedNodeForces[rightClampBorder[k]][j];
			}
		}
	}
	outputFile<<"Forces on clamps lhs: "<<leftClampForces[0]<<" "<<leftClampForces[1]<<" "<<leftClampForces[2]<<" rhs: "<<rightClampForces[0]<<" "<<rightClampForces[1]<<" "<<rightClampForces[2]<<endl;
}

void Simulation::moveClampedNodesForStretcher(){
	if (currSimTimeSec>=StretchInitialTime && currSimTimeSec<StretchEndTime){
		for (int i=0; i<nNodes; ++i){
			if (Nodes[i]->Position[distanceIndex]> StretchMax){
				Nodes[i]->Position[distanceIndex] += StretchDistanceStep;
			}
			else if( Nodes[i]->Position[distanceIndex] < StretchMin ){
				Nodes[i]->Position[distanceIndex] -= StretchDistanceStep;
			}
		}
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			(*itElement)->updatePositions(Nodes);
        }
	}
}

void Simulation::setupPipetteExperiment(){
	pipetteRadiusSq = pipetteRadius*pipetteRadius;
	effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
	effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;
	//Now I set up the upper and lower boundaries to cover both apical and basal suction, but I need to extend the boundary
	//inside the tube, such that the nodes will not stop being suced in after they go inside a certain length into the tube
	//This means increasing the upper boundary for apical suction, and reducing the lower boundary for basal suction:
	if(ApicalSuction){
		effectLimitsInZ[1] += 1000;
	}
	else{
		effectLimitsInZ[0] -= 1000;
	}
	//cout<<"set the system, fixing nodes:"<<endl;
    //Now I am sticking the other side of the tissue to a surface
	if (ApicalSuction){
		for (int i=0; i<nNodes; ++i){
			//fix basal nodes of columnar layer:
			if(Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0) {
				fixAllD(i, false); //this is fixing with adhesives, should be a hard fix at all times
			}
		}
	}
	else{
		for (int i=0; i<nNodes; ++i){
			//fix apical nodes and all peripodial membrane nodes:
			if(Nodes[i]->tissuePlacement == 1 || Nodes[i]->tissueType == 1) {
				fixAllD(i, false); //this is fixing with adhesives, should be a hard fix at all times
			}
		}
	}
}

void Simulation::addMyosinForces(gsl_matrix* gExt){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    int* nodeIds = (*itElement)->getNodeIds();
	    int nNodes= (*itElement)->getNodeNumber();
	    for (int j=0; j<nNodes; ++j){
			double Fx = (*itElement)->MyoForce[j][0];
			double Fy = (*itElement)->MyoForce[j][1];
			double Fz = (*itElement)->MyoForce[j][2];
			int indice = nodeIds[j]*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
			Fz += gsl_matrix_get(gExt,indice+2,0);
			gsl_matrix_set(gExt,indice+2,0,Fz);
	    }
	}
}

void Simulation::addPipetteForces(gsl_matrix* gExt){
    if      (timestep == 0) {SuctionPressure[2] = 0;}
    else if (timestep == 1) {SuctionPressure[2] = 100;}
    else if (timestep == 2) {SuctionPressure[2] = 200;}
    else if (timestep == 3) {SuctionPressure[2] = 300;}
    else if (timestep == 4) {SuctionPressure[2] = 400;}
    else if (timestep == 5) {SuctionPressure[2] = 500;}
    else if (timestep == 6) {SuctionPressure[2] = 600;}

    //cout<<"in add pipette forces, pipette pos: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<endl;
    int dim = 3;
	for (int i=0; i<nNodes; ++i){
        //cout<<"Node "<<i<<" z pos: "<<Nodes[i]->Position[2]<<" effectLimitsInZ: "<<effectLimitsInZ[0]<<" "<<effectLimitsInZ[1]<<endl;
		if (Nodes[i]->Position[2]> effectLimitsInZ[0] &&  Nodes[i]->Position[2]< effectLimitsInZ[1]){
            //cout<<"Node "<<i<<" is within z range"<<endl;
			double dx = pipetteCentre[0] - Nodes[i]->Position[0];
			double dy = pipetteCentre[1] - Nodes[i]->Position[1];
			//cout<<"dx: "<<dx<<" dy: "<<dy<<" dx*dx+dy*dy: "<<dx*dx+dy*dy<<endl;
			double d2 = dx*dx+dy*dy ;
			if (d2 < pipetteRadiusSq){
                //cout<<"Node "<<i<<" is within planar range"<<endl;
                double multiplier = 1.0;
                bool scalePressure = false;
				if(scalePressure){ multiplier = (1 - d2/pipetteRadiusSq);}
				double LocalPressure[3] = { multiplier*SuctionPressure[0], multiplier*SuctionPressure[1],multiplier*SuctionPressure[2]};
				//cout<<"Node: "<<i<<" being sucked into pipette, force: "<<SuctionPressure[0]*Nodes[i]->surface<<" "<<SuctionPressure[1]*Nodes[i]->surface<<" "<<SuctionPressure[2]*Nodes[i]->surface<<endl;
				//node is within suciton range
				if(!Nodes[i]->FixedPos[0]){
                    double value = gsl_matrix_get(gExt,i*dim,0);
                    value +=LocalPressure[0]*Nodes[i]->zProjectedArea;
                    gsl_matrix_set(gExt,i*dim,0,value);
                }
				if(!Nodes[i]->FixedPos[1]){
                    double value = gsl_matrix_get(gExt,i*dim+1,0);
                    value +=LocalPressure[1]*Nodes[i]->zProjectedArea;
                    gsl_matrix_set(gExt,i*dim+1,0,value);
                }
				if(!Nodes[i]->FixedPos[2]){
                    double value = gsl_matrix_get(gExt,i*dim+2,0);
                    value +=LocalPressure[2]*Nodes[i]->zProjectedArea;
                    gsl_matrix_set(gExt,i*dim+2,0,value);
                }
                //Elements[0]->displayMatrix(gExt,"gExtInsidePipetteFunc");
			}
		}
	}
}

void Simulation::packToPipetteWall(){
	//Fixing the z-height of nodes near the pipette wall:
	double pipLim2[2] = {0.9*pipetteRadius, 1.4*(pipetteRadius+pipetteThickness)};
	pipLim2[0] *= pipLim2[0];
	pipLim2[1] *= pipLim2[1];
	int nZfix = TransientZFixListForPipette.size();
	for (int zfixiter =0; zfixiter <nZfix; zfixiter++){
		Nodes[TransientZFixListForPipette[zfixiter]]->FixedPos[2] = 0;
	}
	TransientZFixListForPipette.clear();
	for (int currId=0; currId<nNodes;currId++){
		if (Nodes[currId]->FixedPos[2] == 0){
			//the node is not already fixed in z
			//I will check if it is in the range of the tip of pipette:
			double zLimits[2] = {-2.0,2.0}; //the range in z to freeze nodes
			double v0[3] = {Nodes[currId]->Position[0]-pipetteCentre[0], Nodes[currId]->Position[1]-pipetteCentre[1], Nodes[currId]->Position[2]-pipetteCentre[2]};
			if (v0[2] < zLimits[1] && v0[2]> zLimits[0]){
				if (currId == 61) {cout<<"Node[61] is in z range"<<endl;}
				//node is in Z range, is it in x-y range?
				double xydist2 = v0[0]*v0[0]+v0[1]*v0[1];
				if (xydist2>pipLim2[0] && xydist2<pipLim2[1]){
					 //node is within x-y range, freeze!
					Nodes[currId]->FixedPos[2] = 1;
					TransientZFixListForPipette.push_back(currId);
				}
			}
		}
	}
}

void Simulation::laserAblateTissueType(int ablationType){
	vector <int> AblatedElements;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->tissueType == ablationType && !Nodes[i]->hasLateralElementOwner){
			//I want to ablate a whole tisse type, such as ablating the disc proper at the beginning of simulation.
			//BUT, I want the nodes that are owned by the linker nodes to stay intact. Otherwise, I will loose part of the unintended tissue.
			AblatedNodes.push_back(i);
		}
		else if (Nodes[i]->Position[2] > 14.25){//I do not want any node above 14.25;
			AblatedNodes.push_back(i);
		}
	}

	int nAN = AblatedNodes.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					//cerr<<"Ablating element:" <<Elements[i]->Id<<endl;
					break;
				}
			}
		}
	}
	//some nodes are left with zero mass, which will cause problems in later calculations:
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::laserAblate(double OriginX, double OriginY, double Radius){
	vector <int> AblatedNodes;
	vector <int> AblatedElements;
	double thres2 = Radius*Radius;
	for (int i=0; i<nNodes; ++i){
		double dx = Nodes[i]->Position[0]- OriginX;
		double dy = Nodes[i]->Position[1]- OriginY;
		double d2 = dx *dx + dy*dy;
		if (d2 < thres2){
			AblatedNodes.push_back(i);
		}
		//if (i != 0 && i != 1 && i != 22 && i != 88 && i != 66 && i != 67  ){
			//AblatedNodes.push_back(i);
		//}
	}

	int nAN = AblatedNodes.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					cerr<<"Ablating element:" <<(*itElement)->Id<<endl;
					break;
				}
			}
		}
	}
	//some nodes are ledt with zero mass, which will cause problems in later calculations:
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::updateElementVolumesAndTissuePlacements(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<"updating element: "<<(*itElement)->Id<<endl;
		(*itElement)->updateElementVolumesAndTissuePlacementsForSave(Nodes);
	}
}

void Simulation::clearNodeMassLists(){
	for (int i=0 ;i<nNodes;++i){
		Nodes[i]->connectedElementIds.size();
		Nodes[i]->connectedElementIds.clear();
		Nodes[i]->connectedElementWeights.clear();
		Nodes[i]->mass=0.0;
		//Nodes[i]->surface=0.0;
	}
}

void Simulation::clearLaserAblatedSites(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->IsAblated){
			(*itElement)->removeMassFromNodes(Nodes);
		}
	}
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::setupYsymmetricity(){
	double yLimPos = 0.1;
	double yLimNeg = (-1.0)*yLimPos;
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[1]< yLimPos){
			if (Nodes[i]->Position[1] > yLimNeg){
				//symmetricYBoundaryNodes.push_back(Nodes[i]);
				Nodes[i]->atSymmetricityBorder = true;
				fixY(Nodes[i],false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
			else{
				AblatedNodes.push_back(i);
				fixAllD(Nodes[i], false); //this is fixing for ablated nodes, no need for calculations
				//setSymmetricNode(Nodes[i],yLimPos);
			}
		}
	}
	int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					break;
				}
			}
		}
	}
}

void Simulation::setupXsymmetricity(){
	double xLimPos = 0.1;
	double xLimNeg = (-1.0)*xLimPos;
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[0]< xLimPos){
			if (Nodes[i]->Position[0] > xLimNeg){
				//symmetricXBoundaryNodes.push_back(Nodes[i]);
				Nodes[i]->atSymmetricityBorder = true;
				fixX(Nodes[i],false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
			else{
				AblatedNodes.push_back(i);
				fixAllD(Nodes[i], false); //this is fixing for ablated nodes, no need for calculations
				//setSymmetricNode(Nodes[i],yLimPos);
			}
		}
	}
	int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					break;
				}
			}
		}
	}
}

void Simulation::ablateSpcific(){
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		fixAllD(Nodes[i],  false); //this is fixing for ablated nodes, no need for calculations);
	}
	//int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double* c = new double[3];
		c = (*itElement)->getCentre();
		//if ((*itElement)->Id != 38 && (*itElement)->Id != 39 /*&& (*itElement)->Id != 30 && (*itElement)->Id != 41*/){
		if (c[0]<20 || c[1] >10 || c[1]<-10){
			(*itElement)->removeMassFromNodes(Nodes);
			(*itElement)->IsAblated = true;
		}
		else{
			for (int i=1; i<nNodes; ++i){
				if ((*itElement)->DoesPointBelogToMe(i)){
					Nodes[i]->FixedPos[0]= false;
					Nodes[i]->FixedPos[1]= false;
					Nodes[i]->FixedPos[2]= false;
				}
			}
		}
		delete[] c;
	}
	Nodes[1]->FixedPos[2]= true;
	Nodes[18]->FixedPos[1]= true;
	Nodes[18]->FixedPos[2]= true;
}

void Simulation::pokeElement(int elementId, double dx, double dy, double dz){
	cout<<" poking element: "<<elementId<<endl;
	int* nodeIds = Elements[elementId]->getNodeIds();
	int nNodes= Elements[elementId]->getNodeNumber();
	for (int j=0; j<nNodes; ++j){
		Nodes[nodeIds[j]]->Position[0] += dx;
		Nodes[nodeIds[j]]->Position[1] += dy;
		Nodes[nodeIds[j]]->Position[2] += dz;
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
    }
}

void Simulation::writeMeshRemovingAblatedRegions(){
	cout<<"writing non-ablated mesh"<<endl;
	int nonAblatedNodeMap[( const int ) nNodes];
	int nNonAblatedNode = 0;
	for (int i=0; i<nNodes; ++i){
		nonAblatedNodeMap[i] = -10;
		if (Nodes[i]->mass > 0){
			nonAblatedNodeMap[i] = nNonAblatedNode;
			nNonAblatedNode++;
		}
	}
	cout<<"got non-ablated node number and map "<<endl;
	int nNonAblatedElements = 0;
	for (int i=0; i<nElements; ++i){
		if (!Elements[i]->IsAblated){
			nNonAblatedElements++;
		}
	}
	cout<<"opening file"<<endl;
	string meshSaveString = saveDirectory +"/MeshFromNonAblated.mesh";
	const char* name_meshSaveString = meshSaveString.c_str();;
	ofstream file;
	file.open(name_meshSaveString, ofstream::out);
	file<<nNonAblatedNode;
	file<<endl;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass > 0){
			file << Nodes[i]->Position[0];
			file<<" \t";
			file << Nodes[i]->Position[1];
			file<<" \t";
			file << Nodes[i]->Position[2];
			file<<" \t";
			file << Nodes[i]->tissuePlacement;
			file<<" \t";
			file << Nodes[i]->tissueType;
			file<<" \t";
			file << Nodes[i]->atCircumference;
			file << endl;
		}
	}
	file<<nNonAblatedElements;
	file<<endl;
	for (int i=0; i<nElements; ++i){
		if (!Elements[i]->IsAblated){
			file<< Elements[i]->getShapeType();
			file<<" \t";
			const int n = Elements[i]->getNodeNumber();
			int* NodeIds;
			NodeIds = new int[n];
			NodeIds = Elements[i]->getNodeIds();
			for (int j=0; j<n; ++j){
				file<< nonAblatedNodeMap[NodeIds[j]];
				file<<" \t";
			}
			int dim  = Elements[i]->getDim();
			double** refPos = Elements[i]->getReferencePos();
			for (int j = 0; j<6; ++j ){
				for (int k = 0; k<dim; ++k ){
					file.precision(5);file.width(12);
					file<<refPos[j][k];
				}
			}
			file << endl;
		}
	}
	//recording tissue weights
	file <<1<<endl;
	for (int i=0; i<nElements; ++i){
			if (!Elements[i]->IsAblated){
				file <<Elements[i]->getPeripodialness() << endl;
			}
	}
	file.close();
}
