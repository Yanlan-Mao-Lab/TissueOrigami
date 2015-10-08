
#include "Simulation.h"
#include "Prism.h"
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
	reachedEndOfSaveFile = false;
	AddPeripodialMembrane = false;
	lumenHeight = -20;
	columnarBoundingBoxSize[0]=1000.0; columnarBoundingBoxSize[1]=1000.0; columnarBoundingBoxSize[2]=1000.0;
	peripodialBoundingBoxSize[0]=1000.0; peripodialBoundingBoxSize[1]=1000.0; peripodialBoundingBoxSize[2]=1000.0;
	ContinueFromSave = false;
    growthRotationUpdateFrequency = 60.0/dt;
    nElements = 0;
    nNodes = 0;
    if (growthRotationUpdateFrequency<1) {growthRotationUpdateFrequency =1;}
	setDefaultParameters();

	//double GrowthMatrix[3][3][3] = {
	//							{{1.00, 0.50, 0.00}, {1.00, 0.50, 0.00}, {1.00, 0.50, 0.00}},
	//							{{1.00, 0.75, 0.00}, {1.00, 0.625,0.00}, {1.00, 0.50, 0.00}},
	//							{{1.00, 1.00, 0.00}, {1.00, 0.75, 0.00}, {1.00, 0.50, 0.00}}
	//							};
}

Simulation::~Simulation(){
    cerr<<"destructor for simulation called"<<endl;
	delete ModInp;
	//n nodes
	for (int j=0;j<nNodes;++j){
		delete[] SystemForces[j];
		delete[] PackingForces[j];
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
    delete[] PackingForces;
    delete[] FixedNodeForces;
    cout<<"deleting elements"<<endl;
	while(!Elements.empty()){
		ShapeBase* tmp_pt;
		tmp_pt = Elements.back();
		Elements.pop_back();
		delete tmp_pt;
        //cerr<<"Element list size: "<<Elements.size()<<endl;
	}
	while(!Nodes.empty()){
		Node* tmp_pt;
		tmp_pt = Nodes.back();
		Nodes.pop_back();
		delete tmp_pt;
        //cerr<<"Node list size: "<<Nodes.size()<<endl;
	}
	while(!GrowthFunctions.empty()){
		GrowthFunctionBase* tmp_GF;
		tmp_GF = GrowthFunctions.back();
		GrowthFunctions.pop_back();
		delete tmp_GF;
        //cerr<<"GrowtFuncitons list size: "<<GrowthFunctions.size()<<endl;
	}
}

void Simulation::setDefaultParameters(){
	dt = 0.01;				//sec
	SimLength = 10.0/dt; 	//timesteps wothr 10 sec of simulation
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
	ApicalVisc = 0.0;
	BasalVisc = 0.0;
	zeroViscosity = true;
	memset(noiseOnPysProp,0,4*sizeof(int));
	// The default input is a calculated mesh of width 4 elements, each element being 2.0 unit high
	// and having 1.0 unit sides of the triangles.
	MeshType = 2;
	Row = 4;
	Column = Row-2;
	SideLength=1.0;
	zHeight = 2.0;
	for (int i=0; i<3; ++i){
		ApicalNodeFix[i]= false;
		BasalNodeFix[i]= false;
		for (int j=0; j<3; j++){
			CircumferentialNodeFix[i][j]= false;
		}
	}
	nGrowthFunctions = 0;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
    GrowthSaved = true;
	ForcesSaved = true;
	VelocitiesSaved = true;
	ProteinsSaved = true;
	PeripodialElasticity = 0.0;
	PeripodialViscosity = ApicalVisc;
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
	StretchInitialStep = -100;
	StretchEndStep = -100;
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

	columnarBoundingBox[0][0] =  1000.0;	//left x
	columnarBoundingBox[0][1] =  1000.0;	//low y
	columnarBoundingBox[0][2] =  1000.0;	//bottom z
	columnarBoundingBox[1][0] = -1000.0;	//right x
	columnarBoundingBox[1][1] = -1000.0;	//high y
	columnarBoundingBox[1][2] = -1000.0;	//top z
	peripodialBoundingBox[0][0] =  1000.0;	//left x
	peripodialBoundingBox[0][1] =  1000.0;	//low y
	peripodialBoundingBox[0][2] =  1000.0;	//bottom z
	peripodialBoundingBox[1][0] = -1000.0;	//right x
	peripodialBoundingBox[1][1] = -1000.0;	//high y
	peripodialBoundingBox[1][2] = -1000.0;	//top z
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
	bool ZerothFrame = true;
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
		if (!success){
			return false;
		}
		readElementDataToContinueFromSave();

		//assignPhysicalParameters();

		if (TensionCompressionSaved){
			readTensionCompressionToContinueFromSave();
		}
        if (GrowthSaved){
            readGrowthToContinueFromSave();
        }
        if (ProteinsSaved){
        	readProteinsToContinueFromSave();
		}
		if (ZerothFrame){
			ZerothFrame = false;
		}
		else{
			timestep = timestep + dataSaveInterval;
		}
		cout<<"read time point: "<<timestep<<" this is "<<timestep*dt<<" sec"<<endl;
		//skipping the footer:
		getline(saveFileToDisplayMesh,currline);
		while (currline.empty() && !saveFileToDisplayMesh.eof()){
			//skipping empty line
			getline(saveFileToDisplayMesh,currline);
		}
	}
	//Outside the loop for reading all time frames:
	updateElementVolumesAndTissuePlacements();
    //cleanMatrixUpdateData();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	calculateColumnarLayerBoundingBox();
	if (AddPeripodialMembrane){
		calculatePeripodialBoundingBox();
	}
	//bring the time step and data save stime steps to the main modelinput:
	dataSaveInterval = dataSaveIntervalCurrentSim;
	dt = timeStepCurrentSim;
	//During a simulation, the data is saved after the step is run, and the time step is incremented after the save. Now I have read in the final step,
	//I need to increment my time step to continue the next time step from here.
	timestep++;
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
				cerr<<"There is no peripodial membrane, while growth function "<<i<<" is applicable to peropodial membrane"<<endl;
				return false;
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
		//cout<<"mesh type : 4"<<endl;
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	Success = calculateTissueHeight(); //Calculating how many layers the columnar layer has, and what the actual height is.
	if (!Success){
		return Success;
	}
	if (AddPeripodialMembrane){
		Success = addPeripodialMembraneToTissue();
	}
	checkForNodeFixing();
	assignTips();
	if (!Success){
		return Success;
	}
	fillInNodeNeighbourhood();
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
	checkForZeroViscosity();
    //calculateStiffnessMatrices();
    calculateShapeFunctionDerivatives();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	//for (int i=0; i<nNodes;++i){
	//	cout<<"Node: "<<i<<endl;
	//	Nodes[i]->displayConnectedElementIds();
	//	Nodes[i]->displayConnectedElementWeights();
	//}
    alignTissueDVToXPositive();
    calculateColumnarLayerBoundingBox();
	if (AddPeripodialMembrane){
		calculatePeripodialBoundingBox();
	}
    calculateDVDistance();
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateRelativePosInBoundingBox(columnarBoundingBox[0][0],columnarBoundingBox[0][1],columnarBoundingBoxSize[0],columnarBoundingBoxSize[1], peripodialBoundingBox[0][0],peripodialBoundingBox[0][1],peripodialBoundingBoxSize[0],peripodialBoundingBoxSize[1]);
    }
	if (stretcherAttached){
		setStretch();
	}
    cout<<"setting the pipette"<<endl;
	if(PipetteSuction){
		setupPipetteExperiment();
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
	return Success;
}

void Simulation::assignNodeMasses(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    (*itElement)->assignVolumesToNodes(Nodes);
		(*itElement)->assignSurfaceAreaToNodes(Nodes);
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

        //Growth information at each step:
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

		//Velocity information at each step:
		saveFileString = saveDirectory +"/Save_Velocity";
		const char* name_saveFileVelocities = saveFileString.c_str();
		saveFileVelocities.open(name_saveFileVelocities, ofstream::binary);
		if (saveFileVelocities.good() && saveFileVelocities.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileVelocities<<endl;
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
		//
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
		cerr<<"at step: "<<timestep<<" could not open file: "<<name_outputFile<<endl;
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
	saveFileSimulationSummary<<"Mesh Type:  ";
	saveFileSimulationSummary<<MeshType<<endl;
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
	int currIndex = 0;
	for (int i=0; i<nGrowthFunctions; ++i){
		GrowthFunctions[i]->writeSummary(saveFileSimulationSummary,dt);
	}
}

void Simulation::writeMyosinSummary(){
	int currIndex = 0;
	for (int i=0; i<nMyosinFunctions; ++i){
		myosinFunctions[i]->writeSummary(saveFileSimulationSummary);
	}
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
	cout<<"the file opened for tension and compression: "<<name_saveFileToDisplayTenComp<<endl;
	saveFileToDisplayTenComp.open(name_saveFileToDisplayTenComp, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayTenComp<<endl;
		TensionCompressionSaved = false;
	}

    saveFileString = saveDirectoryToDisplayString +"/Save_Growth";
    const char* name_saveFileToDisplayGrowth = saveFileString.c_str();
    cout<<"the file opened for tension and compression: "<<name_saveFileToDisplayGrowth<<endl;
    saveFileToDisplayGrowth.open(name_saveFileToDisplayGrowth, ifstream::in);
    if (!(saveFileToDisplayGrowth.good() && saveFileToDisplayGrowth.is_open())){
        cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowth<<endl;
        GrowthSaved = false;
    }

	saveFileString = saveDirectoryToDisplayString +"/Save_Velocity";
	const char* name_saveFileToDisplayVel = saveFileString.c_str();;
	saveFileToDisplayVel.open(name_saveFileToDisplayVel, ifstream::in);
	if (!(saveFileToDisplayVel.good() && saveFileToDisplayVel.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayVel<<endl;
		VelocitiesSaved = false;
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
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if (VelocitiesSaved){
		updateVelocitiesFromSave();
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
	return true;
}

bool Simulation::readNodeDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
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
	//for (int i=0; i<n; ++i){
	//	cout<<"node "<<i<<"pos: "<<
	//}
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
		return false;
	}
	for (int i=0; i<nElements; ++i){
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
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
	delete[] PackingForces;
	delete[] FixedNodeForces;
	//reinitiating with the new size:
	const int n = nNodes;
	SystemForces = new double*[n];
	PackingForces = new double*[n];
	FixedNodeForces = new double*[n];
	for (int j=0;j<n;++j){
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
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
}

void Simulation::updateVelocitiesFromSave(){
	for (int i=0;i<nNodes;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			saveFileToDisplayVel.read((char*) &Nodes[i]->Velocity[j], sizeof Nodes[i]->Velocity[j]);
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
		gsl_matrix * currFg = gsl_matrix_calloc(3,3);
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
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        gsl_matrix * currFg = gsl_matrix_calloc(3,3);
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
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	vector<ShapeBase*>::iterator it = Elements.begin();
	it += k;
	Elements.insert(it,PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
}


void Simulation::updateOneStepFromSave(){
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
	//cout<<"skipped header: "<<currline<<endl;
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
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if(VelocitiesSaved){
		updateVelocitiesFromSave();
	}
	if(ProteinsSaved){
		updateProteinsFromSave();
	}
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
	calculateColumnarLayerBoundingBox();
	if (AddPeripodialMembrane){
		calculatePeripodialBoundingBox();
	}
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
	if ( Nodes.size() != nNodes){
		//the node number is change, I yupdated the node list, now I need to fix system forces:
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

bool Simulation::initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight ){
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
	else if ( MeshType ==4 ){
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
}

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
		initiateNodesFromMeshInput();
		initiateElementsFromMeshInput();
		saveFileToDisplayMesh.close();
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
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
		return false;
	}
	return true;
}

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

	/*cout<<"ColumnarCircumferencialNodeList before sort"<<endl;
	for (int j =0 ; j<n; ++j){
		cout<<j<<"of "<<n<<": "<<  ColumnarCircumferencialNodeList[j]<<" "<<angles[j]<<endl;
	}*/
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

	/*cout<<"ColumnarCircumferencialNodeList after sort"<<endl;
	for (int j =0 ; j<n; ++j){
		cout<<ColumnarCircumferencialNodeList[j]<<" "<<angles[j]<<endl;
	}*/
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
			if((*itNode)->tissuePlacement == 0){ //Node is basal
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
	return true;
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

void Simulation::fixAllD(Node* currNode){
	for (int j =0 ; j<currNode->nDim; ++j){
		currNode->FixedPos[j]=true;
	}
}
void Simulation::fixAllD(int i){
	for (int j =0 ; j<Nodes[i]->nDim; ++j){
		Nodes[i]->FixedPos[j]=true;
	}
}

void Simulation::fixX(Node* currNode){
	if(currNode->nDim>0){
		currNode->FixedPos[0]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have x-dimension"<<endl;
	}
}
void Simulation::fixX(int i){
	if(Nodes[i]->nDim>0){
		Nodes[i]->FixedPos[0]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have x-dimension"<<endl;
	}
}

void Simulation::fixY(Node* currNode){
	if(currNode->nDim>1){
		currNode->FixedPos[1]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have y-dimension"<<endl;
	}
}
void Simulation::fixY(int i){
	if(Nodes[i]->nDim>1){
		Nodes[i]->FixedPos[1]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have y-dimension"<<endl;
	}
}

void Simulation::fixZ(Node* currNode){
	if(currNode->nDim>2){
		currNode->FixedPos[2]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have z-dimension"<<endl;
	}
}

void Simulation::fixZ(int i){
	if(Nodes[i]->nDim>2){
		Nodes[i]->FixedPos[2]=true;
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
	FixedNodeForces = new double*[n];

	for (int j=0;j<n;++j){
		//3 dimensions
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
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
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	fixZ(0);
	fixZ(1);
	fixZ(2);
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

void Simulation::checkForNodeFixing(){
	//Are there any circumferential node fixing options enabled:
	bool thereIsCircumFix = false;
	for (int i=0;i<3; ++i){
		for (int j=0;j<3; ++j){
			if (CircumferentialNodeFix[i][j] == true){
				thereIsCircumFix = true;
				break;
			}
		}
	}
	//If there is any circumferential Node fixing:
	if (thereIsCircumFix){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->atCircumference){
				//Node is at the circumference, now checking for all possibilities:
				// if i == 0 , I am checking for apical circumference
				// if i == 1 , I am checking for basal  circumference
				// if i == 2 , I am checking for all    circumference
				for (int i=0;i<3; ++i){
					if ( (i == 0 && (*itNode)->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
						 (i == 1 && (*itNode)->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
						 (i == 2)){										  //tissuePlacement is irrelevant, fixing all

						//The node is at circumference; if
						if (CircumferentialNodeFix[i][0]){
							fixX((*itNode));
						}
						if (CircumferentialNodeFix[i][1]){
							fixY((*itNode));
						}
						if (CircumferentialNodeFix[i][2]){
							fixZ((*itNode));
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
					fixX((*itNode));
				}
				if (BasalNodeFix[1]){
					fixY((*itNode));
				}
				if (BasalNodeFix[2]){
					fixZ((*itNode));
				}
			}
		}
	}
	if (ApicalNodeFix[0] || ApicalNodeFix[1]  || ApicalNodeFix[2]){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 1){
				if (ApicalNodeFix[0]){
					fixX((*itNode));
				}
				if (ApicalNodeFix[1]){
					fixY((*itNode));
				}
				if (ApicalNodeFix[2]){
					fixZ((*itNode));
				}
			}
		}
	}
}

/*
void Simulation::GenerateCircumferencialNodeList(vector<int> &NodesToFix, int nLastRow){
	//Now I have a list of nodes to fix, written in a specific order.
	//I would like reorder these, to generate a continuous strip around the tissue
	vector <int> list1, list2, list3, list4;
    int n = NodesToFix.size();
    int i=0;
    list1.push_back(NodesToFix[i]);
    i++;
    list2.push_back(NodesToFix[i]);
    i++;
    cout<<"generating Circumferencial node list = before while loop"<<endl;
    while(i<n-2*nLastRow){
    	list1.push_back(NodesToFix[i]);
    	list2.push_back(NodesToFix[i+1]);
    	list3.push_back(NodesToFix[i+2]);
    	list4.push_back(NodesToFix[i+3]);
    	i += 4;
    }
    cout<<"generating Circumferencial node list - after while loop"<<endl;
    for ( int k=0; k<nLastRow; ++k){
    	list1.push_back(NodesToFix[i]);
    	i++;
    }
    for ( int k=0; k<nLastRow; ++k){
    	list3.push_back(NodesToFix[i]);
        i++;
    }
     //now I will need to merge these 4 lists in the order:
    //list1 + inverse of list2 + list4 + inverse of list3
    //list1:
    for(vector <int>::iterator it = list1.begin(); it<list1.end(); ++it){
    	CircumferencialNodeList.push_back(*it);
    }
    //inverse of list2:
    for(vector <int>::iterator it = list2.end()-1; it>=list2.begin(); --it){
    	CircumferencialNodeList.push_back(*it);
    }
    //list4:
    for(vector <int>::iterator it = list4.begin(); it<list4.end(); ++it){
    	CircumferencialNodeList.push_back(*it);
    }
    //inverse of list3:
    for(vector <int>::iterator it = list3.end()-1; it>=list3.begin(); --it){
    	CircumferencialNodeList.push_back(*it);
    }
    cout<<"merged lists"<<endl;
}
*/
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
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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

void Simulation::assignPhysicalParameters(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double r = (rand() % 200) / 100.0;	//random number between 0.00 and 2.00
		r = r - 1.0; 						//random number between -1.00 and 1.00
		float noise1 = r*noiseOnPysProp[0];	//percent noise on current element
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise2 = r*noiseOnPysProp[1];
		if ((*itElement)->tissueType == 0){ //Element is on the columnar layer
			(*itElement)->setElasticProperties(EApical*(1 + noise1/100.0),EBasal*(1 + noise1/100.0),EMid*(1 + noise1/100.0),poisson*(1 + noise2/100));
		}
		else if ((*itElement)->tissueType == 1){ //Element is on the peripodial membrane
			double currE = PeripodialElasticity*(1 + noise1/100.0);
			(*itElement)->setElasticProperties(currE,currE,currE,poisson*(1 + noise2/100));
		}
		else if ((*itElement)->tissueType ==2 ){ //Element is on the linker Zone, I will weight the values:
			double currPeripodialE = PeripodialElasticity*(1 + noise1/100.0);
			double currEApical = EApical*(1 + noise1/100.0);
			double currEBasal = EBasal*(1 + noise1/100.0);
			double currEMid = EMid*(1 + noise1/100.0);
			double periWeight = (*itElement)->getPeripodialness();
			double colWeight = (*itElement)->getColumnarness();
			currEApical = colWeight * currEApical + periWeight * currPeripodialE;
			currEBasal  = colWeight * currEBasal  + periWeight * currPeripodialE;
			currEMid    = colWeight * currEMid    + periWeight * currPeripodialE;
			(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,poisson*(1 + noise2/100));
		}

	}
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		(*itNode)->setViscosity(ApicalVisc*noise3, BasalVisc*noise3, PeripodialViscosity*noise3);
	}
}

void Simulation::checkForZeroViscosity(){
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->Viscosity > 0){
			zeroViscosity = false;
			break;
		}
	}
    if (zeroViscosity){
    	if(nNodes < 3) {
    		cerr<<"Cannot run zero viscosity simulation with less than 3 nodes, run simulation at your own risk!"<<endl;
    	}
        //fix x,y,z for 0
        Nodes[0]->FixedPos[0] = true;
        Nodes[0]->FixedPos[1] = true;
        Nodes[0]->FixedPos[2] = true;
        //fix x and z for 1
        Nodes[1]->FixedPos[0] = true;
        Nodes[1]->FixedPos[2] = true;
        //fix z for 2
        Nodes[2]->FixedPos[2] = true;
    }
}

void Simulation::manualPerturbationToInitialSetup(bool deform, bool rotate){
    if(timestep==0){
        double scaleX = 1.0;
        double scaleY = 1.0;
        double scaleZ = 1.0;

        double PI = 3.14159265359;
        double tetX = 0 *PI/180.0;
        double tetY = 0 *PI/180.0;
        double tetZ = 45 *PI/180.0;
        if(deform){
        	vector<Node*>::iterator itNode;
        	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
        		(*itNode)->Position[0] *=scaleX;
        		(*itNode)->Position[1] *=scaleY;
                (*itNode)->Position[2] *=scaleZ;
            }
        }

        if (rotate){
        	//cerr<<"Rotating system"<<endl;
            double R[3][3];
            double c = cos(tetX), s = sin(tetX);
            double Rx[3][3] = {{1,0,0},{0,c,-s},{0,s,c}};
            c = cos(tetY); s = sin(tetY);
            double Ry[3][3] = {{c,0,s},{0,1,0},{-s,0,c}};
            c = cos(tetZ), s = sin(tetZ);
            double Rz[3][3] = {{c,-s,0},{s,c,0},{0,0,1}};
            for (int j =0; j<3;++j){
                for(int k=0; k<3;++k){
                    R[j][k] = Rz[j][k];
                }
            }
        	vector<Node*>::iterator itNode;
        	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
                double x = (*itNode)->Position[0]*R[0][0] + (*itNode)->Position[1]*R[0][1] + (*itNode)->Position[2]*R[0][2];
                double y = (*itNode)->Position[0]*R[1][0] + (*itNode)->Position[1]*R[1][1] + (*itNode)->Position[2]*R[1][2];
                double z = (*itNode)->Position[0]*R[2][0] + (*itNode)->Position[1]*R[2][1] + (*itNode)->Position[2]*R[2][2];
                (*itNode)->Position[0]=x;
                (*itNode)->Position[1]=y;
                (*itNode)->Position[2]=z;
            }
        }
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
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

void Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double* pos, double sideThickness){
	double* vec0;
	vec0 = new double[3];
	double midpoint[3] = {0,0,0};
	for (int j=0; j<Nodes[nodeId0]->nDim; j++){
		vec0[j] = (Nodes[nodeId1]->Position[j] - Nodes[nodeId0]->Position[j])/2.0;
		midpoint[j] = Nodes[nodeId0]->Position[j] + vec0[j];
	}
	Elements[0]->normaliseVector3D(vec0);
	//The list is sorted counter-cock-wise, to point out, I will roate normalised vector v0 -90 degrees on z axis:
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
		cerr<<"WARNING, the lateral connection thickness between the peripodial membrane and the columnar layer is too different than aerage alament size (more than 5 fold diference)"<<endl;
	}
	//Now I need the average side of an element, to add new nodes accordingly:
	int nCircumference = ColumnarBasedNodeArray.size();
	for (int i=0; i<nCircumference; ++i){
		//cout<<"at node in the list: "<<i<<endl;
		int nodeId0 = ColumnarBasedNodeArray[i][0];
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
	for (int i=0;i<nCircumference; ++i){
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
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
						//found the current node on hte existing list
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
		int indiceTri0Corner0, indiceTri0Corner1,indiceTri0Corner2;
		int n = PeripodialCapNodeArray.size();
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
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
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
	//I have a list of the Ids fot the OuterNodes I have added, changeing their circumferential node flags to true:
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

bool Simulation::addPeripodialMembraneToTissue(){
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

void Simulation::runOneStep(){
    cout<<"entered run one step"<<endl;
    manualPerturbationToInitialSetup(false,false); //bool deform, bool rotate
    cleanGrowthData();
    resetForces();
    int freq = 10.0/dt ;
    if (freq <1 ){freq =1;}
    if ((timestep - 1)% freq  == 0){
        cout<<"At time -- "<<dt*timestep<<" sec ("<<dt*timestep/3600<<" hours - "<<timestep<<" timesteps)"<<endl;
        alignTissueDVToXPositive();
        //alignTissueAPToXYPlane();
        calculateColumnarLayerBoundingBox();
        calculateApicalSize();
    	if (AddPeripodialMembrane){
    		calculatePeripodialBoundingBox();
    	}
        calculateDVDistance();
        vector<ShapeBase*>::iterator itElement;
        for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        	(*itElement)->calculateRelativePosInBoundingBox(columnarBoundingBox[0][0],columnarBoundingBox[0][1],columnarBoundingBoxSize[0],columnarBoundingBoxSize[1], peripodialBoundingBox[0][0],peripodialBoundingBox[0][1],peripodialBoundingBoxSize[0],peripodialBoundingBoxSize[1]);
        }
    }
    checkForExperimentalSetupsBeforeIteration();
    if(nMyosinFunctions > 0){
    	checkForMyosinUpdates();
    }
    if(nGrowthFunctions>0 || nShapeChangeFunctions >0){
        if ((timestep - 1)% growthRotationUpdateFrequency  == 0){
            updateGrowthRotationMatrices();
        }
        //outputFile<<"calculating growth"<<endl;
        if(nGrowthFunctions>0){
        	calculateGrowth();
        }
        if(nShapeChangeFunctions>0){
        	calculateShapeChange();
        }
    }
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->growShapeByFg(dt);
        	(*itElement)->changeShapeByFsc(dt);
        }
    }
    updateNodeMasses();
    updateElementToConnectedNodes(Nodes);
    calculateMyosinForces();
    updateStepNR();
    processDisplayDataAndSave();
    cout<<"finalised run one step"<<endl;
    calculateColumnarLayerBoundingBox();
	if (AddPeripodialMembrane){
		calculatePeripodialBoundingBox();
	}
    //cout<<" step: "<<timestep<<" Pressure: "<<SuctionPressure[2]<<" Pa, maximum z-height: "<<boundingBox[1][2]<<" L/a: "<<(boundingBox[1][2]-50)/(2.0*pipetteRadius)<<endl;
    timestep++;
    //for (int i=0; i<nNodes; ++i){
    //	cout<<" Nodes["<<i<<"]->Position[0]="<<Nodes[i]->Position[0]<<"; Nodes["<<i<<"]->Position[1]="<<Nodes[i]->Position[1]<<";  Nodes["<<i<<"]->Position[2]="<<Nodes[i]->Position[2]<<"; "<<endl;
    //}
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

void Simulation::updateStepNR(){
    int dim  = 3;
    int iteratorK = 0;
    bool converged = false;
    gsl_matrix* un = gsl_matrix_calloc(dim*nNodes,1);
    constructUnMatrix(un);
    gsl_matrix* mviscdt = gsl_matrix_calloc(dim*nNodes,dim*nNodes);
    constructLumpedMassViscosityDtMatrix(mviscdt);
    gsl_matrix* ge = gsl_matrix_calloc(dim*nNodes,1);
    gsl_matrix* gv = gsl_matrix_calloc(dim*nNodes,1);
    gsl_matrix* gExt = gsl_matrix_calloc(dim*nNodes,1);
    gsl_vector* gSum = gsl_vector_calloc(dim*nNodes);
    gsl_matrix* uk = gsl_matrix_calloc(dim*nNodes,1);
    gsl_vector* deltaU = gsl_vector_calloc(dim*nNodes);
    Elements[0]->createMatrixCopy(uk,un);
    gsl_matrix* K = gsl_matrix_calloc(dim*nNodes,dim*nNodes);
    //Elements[0]->displayMatrix(mviscdt,"mviscdt");
    //cout<<" checking for Pipette forces: "<<PipetteSuction<<" Pipette time: "<<PipetteInitialStep<<" "<<PipetteEndStep<<" timestep "<<timestep<<endl;
    while (!converged){
        cout<<"iteration: "<<iteratorK<<endl;
        //cout<<"Element 0 BasalArea: "<<Elements[0]->ReferenceShape->BasalArea<<endl;
        resetForces();
        gsl_matrix_set_zero(ge);
        gsl_matrix_set_zero(gv);
        gsl_matrix_set_zero(K);
        gsl_vector_set_zero(gSum);

        //Elements[0]->displayMatrix(gExt,"gExt");
        //cout<<"calculating elastic forces"<<endl;
        calculateElasticForcesForNR(ge);
        //cout<<"calculating viscous forces"<<endl;
        calculateViscousForcesForNR(gv,mviscdt,uk,un);
        calculateImplicitKElastic(K);
        //calculateImplucitKElasticNumerical(K, ge);
        calculateImplucitKViscous(K,mviscdt);
        //calculateImplucitKViscousNumerical(mviscdt,un,uk);
        for (int i=0; i<dim*nNodes; ++i){
            gsl_vector_set(gSum,i,gsl_matrix_get(ge,i,0)+gsl_matrix_get(gv,i,0));
        }
        gsl_matrix_set_zero(gExt);
		if (PipetteSuction && timestep >= PipetteInitialStep && timestep<PipetteEndStep){
			packToPipetteWall();
			calculateZProjectedAreas();
			addPipetteForces(gExt);
		}
		addMyosinForces(gExt);
        for (int i=0; i<dim*nNodes; ++i){
            gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
        }
        checkForExperimentalSetupsWithinIteration();
        //Adding external forces:
        //double F = 5.509e+04 ;
        double F = 0 ;
        const int nNodesToPush = 1;
        int NodesToPush[nNodesToPush] = {5};
        int AxesToPush[nNodesToPush] = {0};
        double ForceDir[nNodesToPush] = {1};
        for(int jj = 0; jj<nNodesToPush; ++jj){
            double value = gsl_vector_get(gSum,NodesToPush[jj]*dim+AxesToPush[jj]);
            value += F*ForceDir[jj];
            gsl_vector_set(gSum,NodesToPush[jj]*dim+AxesToPush[jj],value);
        }
        //calculating external force with stress on node 0:
        /*gsl_matrix* Externalstress = gsl_matrix_calloc(6,1);
        gsl_matrix* ExternalNodalForces = gsl_matrix_calloc(3,1);
        gsl_matrix_set(Externalstress,2,0,50.0);
        Elements[0]->calculateForceFromStress(0, Externalstress, ExternalNodalForces);
        for(int jj = 0; jj<nNodesToPush; ++jj){
            for (int kk=0; kk<dim; ++kk){
                double value = gsl_vector_get(gSum,NodesToPush[jj]*dim+kk);
                value -= gsl_matrix_get(ExternalNodalForces,kk,0);
                gsl_vector_set(gSum,NodesToPush[jj]*dim+kk,value);
            }
        }*/
        //Finalised adding external forces
        calcutateFixedK(K,gSum);
        //converged = checkConvergenceViaForce(gSum);
        if (converged){
            break;
        }
        //Elements[0]->displayMatrix(K,"K");
        //Elements[0]->displayMatrix(ge,"ElasticForces");
        //Elements[0]->displayMatrix(gv,"ViscousForces");
        //Elements[0]->displayMatrix(gSum,"TotalForces");
        //Elements[0]->displayMatrix(mviscdt,"mviscdt");
        //Elements[0]->displayMatrix(un,"un");
        //Elements[0]->displayMatrix(uk,"uk");
        solveForDeltaU(K,gSum,deltaU);
        converged = checkConvergenceViaDeltaU(deltaU);
        //Elements[0]->displayMatrix(deltaU,"deltaU");
        updateUkInNR(uk,deltaU);
        updateElementPositionsinNR(uk);
        updateNodePositionsNR(uk);
        iteratorK ++;
        //Elements[0]->displayMatrix(deltaU,"MovementInIteration");
        //Elements[0]->displayMatrix(uk,"newPosiitons");
        if (iteratorK > 1000){
            cerr<<"Error: did not converge!!!"<<endl;
            converged = true;
        }
        /*if (iteratorK  ){
            cout<<"System Forces - Displacements"<<endl;
            for (int i=0; i<nNodes; i++){
                cout<<SystemForces[i][0]<<" "<<gsl_vector_get(deltaU,3*i)<<endl;
                cout<<SystemForces[i][1]<<" "<<gsl_vector_get(deltaU,3*i+1)<<endl;
                cout<<SystemForces[i][2]<<" "<<gsl_vector_get(deltaU,3*i+2)<<endl;
            }
        }*/
        //cout<<"Node 0 position: "<<Nodes[0]->Position[0]<<" "<<Nodes[0]->Position[1]<<" "<<Nodes[0]->Position[2]<<endl;
    }
    checkForExperimentalSetupsAfterIteration();
    //Now the calculation is converged, I update the node positions with the latest positions uk:
    updateNodePositionsNR(uk);
    //Element positions are already up to date.
    gsl_matrix_free(un);
    gsl_matrix_free(mviscdt);
    gsl_matrix_free(ge);
    gsl_matrix_free(gv);
    gsl_matrix_free(gExt);
    gsl_vector_free(gSum);
    gsl_matrix_free(uk);
    gsl_vector_free(deltaU);
    gsl_matrix_free(K);
}

void Simulation::calculateImplucitKElasticNumerical(gsl_matrix* K, gsl_matrix* geNoPerturbation){
    int dim = 3;
    double perturbation = 1E-6;
    gsl_matrix* KPerturbed = gsl_matrix_calloc(nNodes*dim,nNodes*dim);
    gsl_matrix* geWithPerturbation = gsl_matrix_calloc(nNodes*dim,1);
    gsl_matrix* geDiff = gsl_matrix_calloc(nNodes*dim,1);

    for (int i =0 ; i<nNodes; i++){
        for (int j=0; j<dim; j++){
            gsl_matrix_set_zero(geWithPerturbation);
            gsl_matrix_set_zero(geDiff);
            resetForces();            
            Elements[0]->Positions[i][j] += perturbation;
            calculateElasticForcesForNR(geWithPerturbation);
            for(int k=0; k<nNodes*dim; ++k){
                gsl_matrix_set(geDiff,k,0,gsl_matrix_get(geWithPerturbation,k,0));
            }
            gsl_matrix_sub(geDiff,geNoPerturbation);
            for(int k=0; k<nNodes*dim; ++k){
                double numericalValue = gsl_matrix_get(geNoPerturbation,k,0) - gsl_matrix_get(geWithPerturbation,k,0);
                numericalValue /= perturbation;
                gsl_matrix_set(KPerturbed,k,i*dim+j,numericalValue);
            }
            //cout<<"Perturbed node: "<<i<<" dimention: "<<j<<endl;
            //Elements[0]->displayMatrix(geWithPerturbation,"geWithPerturbation");
            //Elements[0]->displayMatrix(geNoPerturbation,"geNoPerturbation");
            //Elements[0]->displayMatrix(geDiff,"geWith-NoPerturb");
            Elements[0]->Positions[i][j] -= perturbation;
        }
    }
    //Elements[0]->displayMatrix(KPerturbed,"KeNumerical");
    //Elements[0]->displayMatrix(K,"KeAnalytical");

    bool useNumericalK = true;
    if (useNumericalK){
        for (int i =0 ; i<nNodes*dim; i++){
            for (int j =0 ; j<nNodes*dim; j++){
                gsl_matrix_set(K,i, j, gsl_matrix_get(KPerturbed,i,j));
            }
        }
    }
    gsl_matrix_free(KPerturbed);
    gsl_matrix_free(geWithPerturbation);
    gsl_matrix_free(geDiff);
}

void Simulation::calcutateFixedK(gsl_matrix* K, gsl_vector* g){
    int dim = 3;
    int Ksize = K->size1;
    //int nFixed = FixedNodes.size();
    //Elements[0]->displayMatrix(K,"Kinitial");
    //for (int iter = 0; iter<nFixed; ++iter){
        //int i = FixedNodes[iter];
        //bool Xfixed = true;
        //bool Yfixed = true;
        //bool Zfixed = true;
        //if (iter == 1) { Yfixed = false;}
        //if (iter == 2) { Xfixed = false;Yfixed = false;}
    for(int i=0; i<nNodes; i++){
        for (int j=0; j<dim; ++j){
            //if ( (Xfixed && j == 0) || (Yfixed && j == 1) || (Zfixed && j == 2) ){
            if (Nodes[i]->FixedPos[j]){
                int index1 = i*dim+j;
                gsl_vector_set(g,index1,0.0); // making the forces zero
                for (int k =0; k<Ksize; ++k){
                    double value =0.0;
                    if (index1 == k ){value =1.0;}
                    gsl_matrix_set(K, index1, k, value);
                    gsl_matrix_set(K, k, index1, value); //K is symmetric;
                }
            }
        }
    }
    //Elements[0]->displayMatrix(g,"gfixed");
    //Elements[0]->displayMatrix(K,"Kfixed");
}

void Simulation::updateNodePositionsNR(gsl_matrix* uk){
    int dim = 3;
    for (int i = 0; i<nNodes; ++i){
        for (int j=0; j<dim; ++j){
            Nodes[i]->Position[j]=gsl_matrix_get(uk,dim*i+j,0);
        }
    }
    //cout<<"finised node pos update"<<endl;
}

void Simulation::calculateImplucitKViscous(gsl_matrix* K, gsl_matrix*  mviscdt){
    gsl_matrix_add(K,mviscdt);
}

void Simulation::calculateImplucitKViscousNumerical(gsl_matrix*  mviscdt, gsl_matrix*  un, gsl_matrix* uk){
    int dim = 3;
    gsl_matrix* ViscousForces0 = gsl_matrix_calloc(nNodes*dim,1);
    gsl_matrix* ViscousForces1 = gsl_matrix_calloc(nNodes*dim,1);
    gsl_matrix* ukepsilon = gsl_matrix_calloc(nNodes*dim,1);
    double epsilon = 0.5;
    for (int i=0; i<nNodes*dim; i++){
        double value = gsl_matrix_get(uk,i,0)+epsilon;
        gsl_matrix_set(ukepsilon,i,0,value);
    }
    for (int i=0; i<nNodes;++i){
        for (int j=0; j<dim; j++){
            int k = dim*i+j;
            double d0 = gsl_matrix_get(un,k,0) - gsl_matrix_get(uk,k,0);
            double d1 = gsl_matrix_get(un,k,0) - gsl_matrix_get(ukepsilon,k,0);
            double v0 = d0 / dt;
            double v1 = d1 / dt;
            double F0 = Nodes[i]->mass*Nodes[i]->Viscosity*v0;
            double F1 = Nodes[i]->mass*Nodes[i]->Viscosity*v1;
            gsl_matrix_set(ViscousForces0,k,0,F0);
            gsl_matrix_set(ViscousForces1,k,0,F1);
        }
    }
    gsl_matrix* NumericKv = gsl_matrix_calloc(nNodes*dim,nNodes*dim);
    for (int i=0; i<nNodes*dim; i++){
        double value = gsl_matrix_get(ViscousForces1,i,0) - gsl_matrix_get(ViscousForces0,i,0) ;
        value /= epsilon;
        gsl_matrix_set(NumericKv,i,i,value);
    }
    Elements[0]->displayMatrix(NumericKv,"NumericKv");
    Elements[0]->displayMatrix(mviscdt,"mviscdt");
    Elements[0]->displayMatrix(ViscousForces0,"ViscousForces0");
    Elements[0]->displayMatrix(ViscousForces1,"ViscousForces1");
    gsl_matrix_free(ViscousForces0);
    gsl_matrix_free(ViscousForces1);
    gsl_matrix_free(ukepsilon);

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

void Simulation::updateUkInNR(gsl_matrix* uk, gsl_vector* deltaU){
    int n = uk->size1;
    for (int i=0; i<n;++i){
        double newValue = gsl_matrix_get(uk,i,0)+gsl_vector_get(deltaU,i);
        gsl_matrix_set(uk,i,0,newValue);
    }
}

void Simulation::solveForDeltaU(gsl_matrix* K, gsl_vector* g, gsl_vector* deltaU){
    int dim = 3;
    const int nmult  = dim*nNodes;

    int *ia = new int[nmult+1];
    double *b = new double[nmult];
    vector <int> ja_vec;
    vector <double> a_vec;

    constructiaForPardiso(K, ia, nmult, ja_vec, a_vec);
    const int nNonzero = ja_vec.size();
    int* ja = new int[nNonzero];
    double* a = new double [nNonzero];
    writeKinPardisoFormat(nNonzero, ja_vec, a_vec, ja, a);
    writeginPardisoFormat(g,b,nmult);
    int error = solveWithPardiso(a, b, ia, ja, deltaU , nmult);
    if (error != 0){cerr<<"Pardiso solver did not return success!!"<<endl;}
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
}

#include <math.h>
/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, double *, int    *,    int *, int *,   int *, int *,   int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


int Simulation::solveWithPardiso(double* a, double*b, int* ia, int* ja, gsl_vector* deltaU ,const int n_variables){

    // I am copying my libraries to a different location for this to work:
    // On MAC:
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgfortran.3.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgomp.1.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libquadmath.0.dylib /usr/local/lib/
    // cp libpardiso500-MACOS-X86-64.dylib usr/local/lib
    //
    // compilation:
    // g++ pardiso_sym.cpp -o pardiso_sym  -L./ -L/usr/local/lib -L/usr/lib/  -lpardiso500-MACOS-X86-64 -llapack


    // On ubuntu,
    // cp libpardiso500-GNU461-X86-64.so /usr/lib/
    //
    // sometimes linux cannot recognise liblapack.so.3gf or liblapack.so.3.0.1 or others like this, are essentially liblapack.so
    // on ubuntu you can get this solved by installing liblapack-dev:
    // sudo apt-get install liblapack-dev
    //
    // compilation:
    // gcc test.cpp -o testexe  -L/usr/lib/  -lpardiso500-GNU461-X86-64  -fopenmp  -llapack

    //
    // also for each terminal run:
    // export OMP_NUM_THREADS=1
    // For mkl this is :
    // export MKL_PARDISO_OOC_MAX_CORE_SIZE=10000
    // export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
    //
    // MSGLVL: the level of verbal output, 0 is no output.

    int    n = n_variables;
    int    nnz = ia[n];
    int    mtype = 11;        /* Real symmetric matrix */

    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    double   x[n_variables], diag[n_variables];
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i, k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */

    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        //printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		    /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    iparm[10] = 0; /* no scaling  */
    iparm[12] = 0; /* no matching */

    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    bool carryOutDebuggingChecks = false;
    if (carryOutDebuggingChecks){
        pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
        if (error != 0) {
            printf("\nERROR in consistency of matrix: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    if (carryOutDebuggingChecks){
        pardiso_chkvec (&n, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR  in right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */
    if (carryOutDebuggingChecks){
        pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11;
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//cout<<"symbolic factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */
    phase = 22;
    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//cout<<"numerical factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    //printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */
   /* phase = 33;

    iparm[7] = 1;       // Max numbers of iterative refinement steps.

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    bool displayResult = false;
    if (displayResult){
        printf("\nSolve completed ... ");
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n");
    }
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }
    */
/* -------------------------------------------------------------------- */
/* ..  Back substitution with tranposed matrix A^t x=b                  */
/* -------------------------------------------------------------------- */

	phase = 33;

	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
	iparm[11] = 1;       /* Solving with transpose matrix. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			 &n, a, ia, ja, &idum, &nrhs,
			 iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	bool displayResult = false;
	if (displayResult){
		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		printf ("\n");
	}
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }

/* -------------------------------------------------------------------- */
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    return 0;
}
/*
void Simulation::constructiaForPardiso(gsl_matrix* K, int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
    double negThreshold = -1E-14, posThreshold = 1E-14;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=i; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //cout<<"writing is for row "<<i<<" column is: "<<j<<endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}
*/
void Simulation::constructiaForPardiso(gsl_matrix* K, int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
    double negThreshold = -1E-14, posThreshold = 1E-14;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=0; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //cout<<"writing is for row "<<i<<" column is: "<<j<<endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}
void Simulation::writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}

void Simulation::writeginPardisoFormat(gsl_vector* g, double* b, const int n){
    for (int i=0; i<n; ++i){
        b[i] = gsl_vector_get(g,i);
    }
}

void Simulation::calculateImplicitKElastic(gsl_matrix* K){
    //cout<<"calculateImplicitKElastic for thw whole system"<<endl;
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateImplicitKElastic();
    }
    //writing all elements K values into big K matrix:
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->writeKelasticToMainKatrix(K);
    }
}

void Simulation::calculateElasticForcesForNR(gsl_matrix* ge){
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
            //cout<<"calculating forces for element:  "<<i<<endl;
        	(*itElement)->calculateForces(SystemForces, Nodes, recordForcesOnFixedNodes, FixedNodeForces, outputFile);
            /*
            cout<<"Element "<<i<<"system forces: "<<endl;
            for (int nn=0; nn< nNodes; ++nn){
                 cout<<"Node: "<<nn<<" "<<SystemForces[nn][0]<<" "<<SystemForces[nn][1]<<" "<<SystemForces[nn][2]<<endl;
            }*/
        }
    }
    //now all the forces are written on SysyemForces
    //I will add them into ge, this step can be made faster by separating calculate forces function into two,
    //and filling up either ge or System forces depending on hte solution method:
    for (int i = 0; i< nNodes; ++i){
        for ( int j=0; j<3; j++){
            gsl_matrix_set(ge, 3*i+j,0,SystemForces[i][j]);
        }
    }
    //cout<<"finalised elastic force calculation of the system"<<endl;
}

void Simulation::calculateViscousForcesForNR(gsl_matrix* gv, gsl_matrix* mviscdt, gsl_matrix* uk, gsl_matrix* un){
    int dim  = 3;
    gsl_matrix* displacement = gsl_matrix_calloc(dim*nNodes,1);
    Elements[0]->createMatrixCopy(displacement,un);
    gsl_matrix_sub(displacement,uk);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mviscdt, displacement,0.0, gv);
    gsl_matrix_free(displacement);
}

void Simulation::constructUnMatrix(gsl_matrix* un){
    for (int i = 0; i<nNodes; ++i ){
        for (int j=0; j<3; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void Simulation::constructLumpedMassViscosityDtMatrix(gsl_matrix* mviscdt){
    for (int i = 0; i<nNodes; ++i ){
        double matrixValue = Nodes[i]->mass*Nodes[i]->Viscosity / dt;
        for (int j=0; j<3; ++j){
            gsl_matrix_set(mviscdt,3*i+j,3*i+j,matrixValue);
        }
    }
    //cout<<" Node 0 - mass: "<<Nodes[0]->mass<<" visc: "<<Nodes[0]->Viscosity<<" matrixvalue: "<<gsl_matrix_get(mviscdt,0,0)<<endl;
}

bool Simulation::checkConvergenceViaDeltaU(gsl_vector* deltaU){
    bool converged = true;
    double Threshold = 1E-10;
    double d = gsl_blas_dnrm2 (deltaU);

    if (d>Threshold){
        converged = false;
        cout<<" not  yet converged via du: norm "<<d<<endl;
    }
    else{
        cout<<"converged with displacement: norm"<<d<<endl;
    }
    return converged;
}

bool Simulation::checkConvergenceViaForce(gsl_vector* gSum){
    bool converged = true;
    double Threshold = 1E-10;
    double d = gsl_blas_dnrm2 (gSum);
    if (d>Threshold){
        converged = false;
        cout<<" not  yet converged via forces: norm "<<d<<endl;
    }
    else{
        cout<<"converged with forces: norm"<<d<<endl;
    }
    return converged;
}


void Simulation::processDisplayDataAndSave(){
    //if (displayIsOn && !DisplaySave){
        //The simulation is not displaying a saved setup, it is running and displaying
        //I need to correct the values to be displayed, and store averages, otherwise
        //the displayed values will be from artificial setups of different RK steps. (RK1 or RK4 depending on parameter)
        //updateDisplaySaveValuesFromRK();
    //}
    if (saveData && timestep % dataSaveInterval == 0){
        //updateDisplaySaveValuesFromRK();
        saveStep();
    }
}

void Simulation::updateNodeMasses(){
    vector<Node*>::iterator itNode;
    for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
    	(*itNode)->mass = 0;
    }
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->assignVolumesToNodes(Nodes);
        }
    }
}

void 	Simulation:: updateElementToConnectedNodes(vector <Node*>& Nodes){
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
	    int n = (*itNode)->connectedElementIds.size();
        for (int i=0; i<n; ++i){
            //if(!Elements[Nodes[j]->connectedElementIds[i]] -> IsAblated){
        	(*itNode)->connectedElementWeights[i] = Elements[(*itNode)->connectedElementIds[i]]->VolumePerNode/(*itNode)->mass;
            //}
        }
    }
}


void Simulation::fillInNodeNeighbourhood(){
	vector<ShapeBase*>::iterator itEle;
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->fillNodeNeighbourhood(Nodes);
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

void Simulation::getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner, bool& bothperipodial){
	if (NodePointer->tissueType == 1 && NodePointer->tissueType == 1){
		//both element and the node are on the peripodial membrane
		//I will update the normal of the element here, and flag the interaction to be between two peripodial membrane entities.
		//This will mean, they will push away from each other, regardless of direction.
		//The flag will be checked to reverse the normal direction in case the dot product is negative. I am not doing it here
		//as it will be a duplicate calculation.
		getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
		bothperipodial = true;
	}
	else if (NodePointer->tissueType == 0 ){
		//the node is on the columnar layer. apical nodes should pack apically and basal nodes should pack basally, the element will give the correct vectors
		//even if it is on peripodial membrane
		if (NodePointer->tissuePlacement == 1){  //node is apical, should pack to apical nodes only
			getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
		}
		else if (NodePointer->tissuePlacement == 0){ //node is basal, should pack to apical nodes only
			getBasalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
		}
	}
	else if (ElementPointer->tissueType == 0 && NodePointer->tissueType == 1){
		//Element is on the columnar layer and the node is on the peripodial membrane
		//if the tissue has multiple discretisation layers, the peripodial node should be pushed out of the tissue side the elemnt belongs to
		//(if the element is apical, push it out towards the apical surface and vice versa.
		if (TissueHeightDiscretisationLayers>1){
			if (ElementPointer->tissuePlacement == 0 ){ //element is basal:
				getBasalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
			}
			else if (ElementPointer->tissuePlacement == 1 ){//element is apical:
				getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
			}
			else{cerr<<"Element has error in tissue placement, packing, Element id: "<<ElementPointer->Id<<endl;}
		}
		else{
			//tissue does not have multiple layers, I will push out from the apical surface
			//if your tissue is complicated enough to have wrapping around peripodial membrane,
			//it will have flips with a single layer of columnar elements anyway.
			//USE multi-layer discretisation for columnar layer if you have packing issues.
			getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
		}
	}
}

inline void Simulation::CapPackingForce(double& Fmag){
	double Fcap = 1e4;
	if (Fmag > Fcap){
		Fmag = Fcap;
	}
}

void Simulation::checkPackingToPipette(bool& packsToPip, double* pos, double* pipF, double mass, int id){
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

void Simulation::calculatePacking(double PeriThreshold, double ColThreshold){
	//cout<<"inside calculate packing"<<endl;
	double multiplier = 0.05; //just a term to scale the forces down
	double threshold = 5;	 //packing forces start at 3 microns
	double t2 = threshold*threshold;

	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	//resetting all the normal update flags
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->ApicalNormalForPackingUpToDate = false;
		(*itEle)->BasalNormalForPackingUpToDate = false;
	}
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
		if (NodeHasPacking){
			//if ((*itNode)->Id == 358) {
			//	cout<<"calculating packing for node: "<<(*itNode)->Id<<endl;
			//}
			double* pos;
			pos = new double[3];
			(*itNode)->getCurrentPosition(pos);
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this element
				bool PackingToThisElement =  (*itEle)->checkPackingToThisNodeViaState(TissueHeightDiscretisationLayers, (*itNode));
				if (PackingToThisElement){
					//the node does not belong to the element, and placed correctly. Lets have a preliminary distance check:
					PackingToThisElement = (*itEle)->IsPointCloseEnoughForPacking((*itNode)->Position, PeriThreshold, ColThreshold,(*itNode)->tissuePlacement, (*itNode)->tissueType);
				}
				if (PackingToThisElement){
					//All position and state conditions are satisfied, and the point is close enough to the element for packing.
					double* normalForPacking;
					normalForPacking = new double[3];
					normalForPacking[0]=0.0;normalForPacking[1]=0.0;normalForPacking[2]=0.0;
					double* posCorner;
					posCorner = new double[3];
					bool bothperipodial = false;
					getNormalAndCornerPosForPacking((*itNode),(*itEle),normalForPacking,posCorner,bothperipodial);
					double dProjectionOffset = 0.0;
					for (int i=0; i<3; ++i){
						dProjectionOffset += (pos[i] - posCorner[i])*normalForPacking[i];
					}
					//cout<<"dProjectionOffset: "<<dProjectionOffset<<endl;
					if (bothperipodial && dProjectionOffset<0){
						//both node and element lie on the peripodial membrane, the directionality does not come from tissue placement, they should always push away.
						dProjectionOffset *= -1.0;
						normalForPacking[0] *= -1.0;
						normalForPacking[1] *= -1.0;
						normalForPacking[2] *= -1.0;
					}
					//cout<<" calculated distace between node and element: "<<dProjectionOffset<<" threshold : "<<threshold<<endl;
					//In the case that there is packing (the point will lie on the triangle, the distance
					//I will be interested in is dProjectionOffet. I will check if this is below threshold,
					//then I will bother with the following calculation if necessary
					if ((dProjectionOffset > 0 && dProjectionOffset < threshold) || (dProjectionOffset < 0 && dProjectionOffset > (-1.0)*threshold)){
						//cout<<"checking element "<<(*itEle)->Id<<" for node: "<<(*itNode)->Id<<endl;
						double VecprojectedPoint[3];
						VecprojectedPoint[0] = (-1.0)*dProjectionOffset*(normalForPacking[0]);
						VecprojectedPoint[1] = (-1.0)*dProjectionOffset*(normalForPacking[1]);
						VecprojectedPoint[2] = (-1.0)*dProjectionOffset*(normalForPacking[2]);
						double projectedPoint[3] = {pos[0] + VecprojectedPoint[0], pos[1] + VecprojectedPoint[1], pos[2] + VecprojectedPoint[2]};
						//cout<<"projected point: "<<projectedPoint[0]<<" "<<projectedPoint[1]<<projectedPoint[2]<<endl;
						//check if the projected point lies on the triangle, if so the distance is (-1.0)*dProjectionOffet
						bool pointInsideTriangle = (*itEle)->IspointInsideTriangle((*itNode)->tissuePlacement,projectedPoint[0],projectedPoint[1],projectedPoint[2]);
						if (pointInsideTriangle){
							//there is packing:
							//cout<<" Node : " <<(*itNode)->Id <<" element: "<<(*itEle)->Id<<" projected point: "<<projectedPoint[0]<<" "<<projectedPoint[1]<<projectedPoint[2]<<endl;
							float d2 = dProjectionOffset*dProjectionOffset;
							//cout<<"	point is inside triangle: d2: "<<d2 <<endl;
							//normalising force direction
							//cout<<" Node : " <<(*itNode)->Id <<" element: "<<(*itEle)->Id<<"	- before normalisation VecprojectedPoint: "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							if (dProjectionOffset<0){
								VecprojectedPoint[0] /= dProjectionOffset;
								VecprojectedPoint[1] /= dProjectionOffset;
								VecprojectedPoint[2] /= dProjectionOffset;
							}
							else{
								VecprojectedPoint[0] /= (-1.0)*dProjectionOffset;
								VecprojectedPoint[1] /= (-1.0)*dProjectionOffset;
								VecprojectedPoint[2] /= (-1.0)*dProjectionOffset;
							}
							//cleaning up noise to allow accurate calculation of strong pushing forces
							//cout<<" Node : " <<(*itNode)->Id <<" element: "<<(*itEle)->Id<<"	- after normalisation VecprojectedPoint: "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							for (int i=0;i<3; i++){
								if(VecprojectedPoint[i]<1E-6 && VecprojectedPoint[i]>-1E-6){
									VecprojectedPoint[i]=0.0;
								}
							}
							//cout<<"Normal for packing for element: "<<(*itEle)->Id<<" vec: "<<normalForPacking[0]<<" "<<normalForPacking[1]<<" "<<normalForPacking[2]<<endl;
							//cout<<"Node "<<(*itNode)->Id<<" position: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<" corner position: "<<posCorner[0]<<" "<<posCorner[1]<<" "<<posCorner[2]<<endl;
							//cout<<"	normalised VecprojectedPoint "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							double averageMass = 0.5*((*itEle)->VolumePerNode + (*itNode)->mass);
							double Fmag = averageMass* multiplier * (1.0/d2 - 1.0/t2);
							CapPackingForce(Fmag);
							//cout<<"	d2  "<<d2 <<endl;
							//cout<<"	Fmag  "<<Fmag <<endl;
							double F[3];
							F[0] = Fmag * VecprojectedPoint[0];
							F[1] = Fmag * VecprojectedPoint[1];
							F[2] = Fmag * VecprojectedPoint[2];
							//cout<<"Node "<<(*itNode)->Id<<" is packing to Element: "<<(*itEle)->Id<<" mag: "<<Fmag<<" Force: "<<F[0] <<" "<<F[1]<<" " <<F[2]<<endl;
							//cout<<" Force: "<<F[0] <<" "<<F[1]<<" " <<F[2]<<endl;
							//direction correction with (-) sign
							for(int j=0; j<3; ++j){
								if (!(*itNode)->FixedPos[j]){
									SystemForces[(*itNode)->Id][j] += F[j];
									PackingForces[(*itNode)->Id][j] += F[j];
								}
								else if (recordForcesOnFixedNodes) {
									FixedNodeForces[(*itNode)->Id][j] += F[j];
								}
							}
							//cout<<"updated the system forces: "<<SystemForces[(*itNode)->Id][0]<<" "<<SystemForces[(*itNode)->Id][1]<<" "<<SystemForces[(*itNode)->Id][2]<<endl;
							//add opposite force to the element nodes:
							(*itEle)->AddPackingToSurface((*itNode)->tissuePlacement, F[0],F[1],F[2], SystemForces,PackingForces, Nodes);
							//cout<<"added packing to surface"<<endl;

							//if(/*(*itNode)->Id == 13 && (*itEle)->Id == 44*/ memberOf44 && (*itEle)->Id>592){
							//if ((*itNode)->Id == 358) {
							//	cout<<" threshold: "<<threshold<<" Node: "<<(*itNode)->Id<<" Element: "<<(*itEle)->Id<<" forcemag: "<<Fmag;
							//}
							//cout<<"threshold: "<<threshold<<" distance: "<<dProjectionOffet<<" Node: "<<(*itNode)->Id<<" Element: "<<(*itEle)->Id<<" forcemag: "<<Fmag;

							//	cout<<" Fdir: "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							//	cout<<"	pos:"<<(*itNode)->Position[0]<<" "<<(*itNode)->Position[1]<<" "<<(*itNode)->Position[2]<<endl;
							//}
						}
					}
					delete[] posCorner;
				}
			}
			if (PipetteSuction && timestep> PipetteInitialStep && timestep<PipetteEndStep){
				if ( (ApicalSuction && (*itNode)->tissuePlacement == 1) || (!ApicalSuction && (*itNode)->tissuePlacement == 0)){
					//checking packing of node to the pipette:
					double* pipF;
					pipF = new double[3];
					bool packsToPip = false;
					checkPackingToPipette(packsToPip, pos, pipF,(*itNode)->mass,(*itNode)->Id);
					if (packsToPip){
						for (int j=0;j<3;++j){
							if (!(*itNode)->FixedPos[j]){
								SystemForces[(*itNode)->Id][j] += pipF[j];
								PackingForces[(*itNode)->Id][j] += pipF[j];
							}
							else if(recordForcesOnFixedNodes){
								FixedNodeForces[(*itNode)->Id][j] += pipF[j];
							}
						}
					}
					delete[] pipF;
				}
			}
			delete[] pos;
		}
	}
	//cout<<"finalised calculate packing"<<endl;
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
	//cout<<"DV Tip node indexes: "<<dorsalTipIndex<<" "<<ventralTipIndex<<endl;
}

void Simulation::alignTissueDVToXPositive(){
	//ventralTipIndex = 0;
	//dorsalTipIndex = 88;
	double* u = new double[3];
	double* v = new double[3];
	//For simulations with no viscosity, the position of node 0 is fixed. Node 0 is the dorsal tip of the tissue.
	// To keep the fixed position of hte
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
	/*for(int i=0;i<n;++i){
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			u[j] = Nodes[i]->Position[j];
		}
		Elements[0]->rotateVectorByRotationMatrix(u,rotMat);
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			Nodes[i]->Position[j] = u[j];
		}
	}*/
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
	for (int i=0;i<3;++i){
		d[i] = Nodes[dorsalTipIndex]->Position[i] - Nodes[ventralTipIndex]->Position[i];
		d[i] *=d[i];

	}
	double dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	//outputFile<<"time: "<<timestep * dt<<" DV distance is: "<<dmag<<endl;
}

void Simulation::cleanGrowthData(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->resetCurrStepGrowthData();
		//Elements[i]->resetCurrStepShapeChangeData();
	}
}

void Simulation::resetForces(){
	int dim = 3;
	//n nodes
	for (int j=0;j<nNodes;++j){
		//3 dimensions
		for (int k=0;k<dim;++k){
			SystemForces[j][k]=0.0;
			if (recordForcesOnFixedNodes){
				FixedNodeForces[j][k]=0.0;
			}
			//PackingForces[j][k]=0.0;
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
	outputFile<<"At time: "<<dt*timestep<<" apical bounding box size: "<<DV<<" "<<AP<<endl;
}

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

void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<timestep*dt<<" sec ("<<timestep*dt/3600<<" hr )"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
    writeGrowth();
	writeForces();
	writeVelocities();
	writeProteins();
}

void Simulation::writeSaveFileStepHeader(){
    saveFileMesh<<"=============== TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<timestep*dt;
	saveFileMesh<<"==================================================="<<endl;
}

void Simulation::writeSaveFileStepFooter(){
	saveFileMesh<<"=============== END OF TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<timestep*dt;
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
			saveFileMesh.width(5);saveFileMesh<<NodeIds[j];
		}
		int dim  = Elements[i]->getDim();
		double** refPos = Elements[i]->getReferencePos();
		for (int j = 0; j<nodeNumber; ++j ){
			for (int k = 0; k<dim; ++k ){
				saveFileMesh.precision(5);saveFileMesh.width(12);
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
        for (int j=0; j<3; ++j){
            for (int k=0; k<3; ++k){
                double Fgjk = gsl_matrix_get(currFg,j,k);
                saveFileGrowth.write((char*) &Fgjk, sizeof Fgjk);
            }
        }
        gsl_matrix_free(currFg);

    }
    saveFileGrowth.flush();
}



void Simulation::writeForces(){
	for (int i=0;i<nNodes;++i){
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileForces.write((char*) &SystemForces[i][1], sizeof SystemForces[i][1]);
		saveFileForces.write((char*) &SystemForces[i][2], sizeof SystemForces[i][2]);
	}
	saveFileForces.flush();
}

void Simulation::writeVelocities(){
	for (int i=0;i<nNodes;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			saveFileVelocities.write((char*) &Nodes[i]->Velocity[j], sizeof Nodes[i]->Velocity[j]);
		}
	}
	saveFileVelocities.flush();
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
	}
	saveFileProteins.flush();
}


void Simulation::calculateMyosinForces(){
	//cout<<"Entered calculateMyosinForces"<<endl;
	cleanUpMyosinForces();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updateMyosinConcentration(dt, kMyo);
		(*itElement)->calculateMyosinForces(forcePerMyoMolecule);
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
		if(timestep == myosinFunctions[i]->initTime  ){ //the application time of the signal is given in simulation time steps.
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
		if ((currMF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currMF->applyToPeripodialMembrane && (*itElement)->tissueType == 1) || (*itElement)->tissueType == 2){
			//calculating the grid indices:
			double* ReletivePos = new double[2];
			//normalising the element centre position with bounding box
			if ((*itElement)->tissueType == 0){
				(*itElement)->getRelativePosInColumnarBoundingBox(ReletivePos);
			}
			else{
				(*itElement)->getRelativePosInPeripodialBoundingBox(ReletivePos);
			}
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
	//cout<<"finalised updateEquilibriumMyosinsFromInputSignal "<<endl;
}


void Simulation::calculateGrowth(){
	//cout<<"Calculating Growth"<<endl;
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctions[i]->Type == 1){
			//cout<<"Calculating Uniform Growth"<<endl;
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

/*
void Simulation::changeCellShapesInSystem(){
	int currIndexForParameters = 0;
	cleanUpShapeChangeRates();
	for (int i=0; i<nShapeChangeFunctions; ++i){
		if(ShapeChangeFunctionTypes[i] == 1){
			changeCellShapeRing(currIndexForParameters);
			currIndexForParameters += 8;
		}
	}
}
*/

void Simulation::cleanUpGrowthRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setGrowthRate(0.0,0.0,0.0);
	}
}

void Simulation::cleanUpShapeChangeRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setShapeChangeRate(0.0,0.0,0.0,0.0,0.0,0.0);
	}
}

void Simulation::calculateShapeChangeUniform (GrowthFunctionBase* currSCF){
	float simTime = dt*timestep;
	if(simTime >= currSCF->initTime && simTime < currSCF->endTime ){
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
	float simTime = dt*timestep;
	//cout<<"inside uniform growth function, initTime: "<<currGF->initTime <<" endtime: "<<currGF->endTime<<" simTime"<<simTime<<endl;
	if(simTime >= currGF->initTime && simTime < currGF->endTime ){
		//cout<<"calculating growth"<<endl;
		double *maxValues;
        maxValues = new double[3];
		currGF->getGrowthRate(maxValues);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			//tissue type == 0 is columnar layer, ==1 is peripodial membrane, ==2 id linker zone
			if ( currGF->applyToColumnarLayer){
				if ((*itElement)->tissueType == 0){ //columnar layer, grow directly
					(*itElement)->updateGrowthRate(maxValues[0],maxValues[1],maxValues[2]);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getColumnarness();
					(*itElement)->updateGrowthRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2]);
				}
			}
			if ( currGF->applyToPeripodialMembrane){
				if ((*itElement)->tissueType == 1){ //peripodial membrane, grow directly
					(*itElement)->updateGrowthRate(maxValues[0],maxValues[1],maxValues[2]);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getPeripodialness();
					(*itElement)->updateGrowthRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2]);
				}
			}
			//if ((currGF->applyToColumnarLayer && Elements[i]->tissueType == 0) || (currGF->applyToPeripodialMembrane && Elements[i]->tissueType == 1)){
			//	 //cout<<"updating growth for element: "<<Elements[i]->Id<<endl;
			//	Elements[i]->updateGrowthRate(maxValues[0],maxValues[1],maxValues[2]);
			//}
		}
		delete[] maxValues;
	}
}

void Simulation::calculateGrowthRing(GrowthFunctionBase* currGF){
	float simTime = dt*timestep;
	if(simTime >= currGF->initTime && simTime < currGF->endTime ){
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
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			if ((currGF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currGF->applyToPeripodialMembrane && (*itElement)->tissueType == 1) || (*itElement)->tissueType == 2){
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
					//growing the shape
                    if ( currGF->applyToColumnarLayer){
						if ((*itElement)->tissueType == 0){ //columnar layer, grow directly
							(*itElement)->updateGrowthRate(growthscale[0],growthscale[1],growthscale[2]);
						}
						else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
							double weight = (*itElement)->getColumnarness();
							(*itElement)->updateGrowthRate(weight*growthscale[0],weight*growthscale[1],weight*growthscale[2]);
						}
                    }
                    if ( currGF->applyToPeripodialMembrane){
						if ((*itElement)->tissueType == 1){ //peripodial membrane, grow directly
							(*itElement)->updateGrowthRate(growthscale[0],growthscale[1],growthscale[2]);
						}
						else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
							double weight = (*itElement)->getPeripodialness();
							(*itElement)->updateGrowthRate(weight*growthscale[0],weight*growthscale[1],weight*growthscale[2]);
						}
					}

				}
				delete[] Elementcentre;
			}
		}
		delete[] maxValues;
	}
}


void Simulation::calculateGrowthGridBased(GrowthFunctionBase* currGF){
	int nGridX = currGF->getGridX();
	int nGridY = currGF->getGridY();
	float simTime = dt*timestep;
	//cout<<"calculating growth grid based, initTime: "<<currGF->initTime<<" endTime: "<< currGF->endTime<<" simTime: "<<simTime<<endl;
	if(simTime >= currGF->initTime && simTime < currGF->endTime ){
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			//cout<<"calculating the element: "<<(*itElement)->Id<<endl;
			if ((currGF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currGF->applyToPeripodialMembrane && (*itElement)->tissueType == 1) || (*itElement)->tissueType == 2){
				double* columnarReletivePos = new double[2];
				double* peripodialRelativePos = new double[2];
				int columnarIndexX = 0.0, peripodialIndexX = 0.0, columnarIndexY = 0.0, peripodialIndexY = 0.0;
				double columnarFracX = 1.0,  peripodialFracX = 1.0, columnarFracY = 1.0, peripodialFracY = 1.0;
				//normalising the element centre position with bounding box
				if ((*itElement)->tissueType == 0){//columnar layer element
					(*itElement)->getRelativePosInColumnarBoundingBox(columnarReletivePos);
					(*itElement)->convertRelativePosToGridIndex(columnarReletivePos, columnarIndexX, columnarIndexY, columnarFracX, columnarFracY, nGridX, nGridY);
					//cout<<"columnarReletivePos: "<<columnarReletivePos[0]<<" "<<columnarReletivePos[1]<<endl;
				}
				else if((*itElement)->tissueType == 1){// peripodial layer element
					(*itElement)->getRelativePosInPeripodialBoundingBox(peripodialRelativePos);
					(*itElement)->convertRelativePosToGridIndex(peripodialRelativePos, peripodialIndexX,peripodialIndexY, peripodialFracX, peripodialFracY, nGridX, nGridY);
					//cout<<"peripodialRelativePos: "<<peripodialRelativePos[0]<<" "<<peripodialRelativePos[1]<<endl;
				}
				else {
					(*itElement)->getRelativePosInColumnarBoundingBox(columnarReletivePos);
					(*itElement)->getRelativePosInPeripodialBoundingBox(peripodialRelativePos);
					(*itElement)->convertRelativePosToGridIndex(columnarReletivePos, columnarIndexX, columnarIndexY, columnarFracX, columnarFracY, nGridX, nGridY);
					(*itElement)->convertRelativePosToGridIndex(peripodialRelativePos, peripodialIndexX,peripodialIndexY, peripodialFracX, peripodialFracY, nGridX, nGridY);
					//cout<<"columnarReletivePos: "<<columnarReletivePos[0]<<" "<<columnarReletivePos[1]<<endl;
					//cout<<"peripodialRelativePos: "<<peripodialRelativePos[0]<<" "<<peripodialRelativePos[1]<<endl;
				}
				//cout<<"The grid size: "<<nGridX<<" "<<nGridY<<endl;
				//cout<<"The indices  : cX: "<<columnarIndexX<<" pX: "<<peripodialIndexX<<" cY: "<<columnarIndexY<<" pY: "<<peripodialIndexY<<endl;
				//cout<<"The fractions: cX "<<columnarFracX<<" pX: "<<peripodialFracX<<" cY: "<<columnarFracY<<" pY: "<<peripodialFracY<<endl;

                double growthYmid[2][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0}};
                double columnarGrowthscale[3] = {0.0,0.0,0.0}, peripodialGrowthscale[3] = {0.0,0.0,0.0};
				for (int axis = 0; axis<3; ++axis){
					if (currGF->applyToColumnarLayer){
						growthYmid[0][axis] = currGF->getGrowthMatrixElement(columnarIndexX,columnarIndexY,axis)*(1.0-columnarFracX) + currGF->getGrowthMatrixElement(columnarIndexX+1,columnarIndexY,axis)*columnarFracX;
						//cout<<"growthYmid[0]["<<axis<<"]: "<<growthYmid[0][axis]<<endl;
						growthYmid[1][axis] = currGF->getGrowthMatrixElement(columnarIndexX,columnarIndexY+1,axis)*(1.0-columnarFracX) + currGF->getGrowthMatrixElement(columnarIndexX+1,columnarIndexY+1,axis)*columnarFracX;
						//cout<<"growthYmid[1]["<<axis<<"]: "<<growthYmid[1][axis]<<endl;
						columnarGrowthscale[axis] = growthYmid[0][axis]*(1.0-columnarFracY) + growthYmid[1][axis]*columnarFracY;
						//cout<<"columnarGrowthscale["<<axis<<"]: "<<columnarGrowthscale[axis]<<endl;
					}
					if (currGF->applyToPeripodialMembrane){
						growthYmid[0][axis] = currGF->getGrowthMatrixElement(peripodialIndexX,peripodialIndexY,axis)*(1.0-peripodialFracX) + currGF->getGrowthMatrixElement(peripodialIndexX+1,peripodialIndexY,axis)*peripodialFracX;
						//cout<<"growthYmid[0]["<<axis<<"]: "<<growthYmid[0][axis]<<endl;
						growthYmid[1][axis] = currGF->getGrowthMatrixElement(peripodialIndexX,peripodialIndexY+1,axis)*(1.0-peripodialFracX) + currGF->getGrowthMatrixElement(peripodialIndexX+1,peripodialIndexY+1,axis)*peripodialFracX;
						//cout<<"growthYmid[1]["<<axis<<"]: "<<growthYmid[1][axis]<<endl;
						peripodialGrowthscale[axis] = growthYmid[0][axis]*(1.0-peripodialFracY) + growthYmid[1][axis]*peripodialFracY;
						//cout<<"peripodialGrowthscale["<<axis<<"]: "<<peripodialGrowthscale[axis]<<endl;
					}
				}
				//growing the shape
				//cout<<"peripodialGrowthscale: "<<peripodialGrowthscale[0]<<" "<<peripodialGrowthscale[1]<<" "<<peripodialGrowthscale[2]<<endl;
				//cout<<"columnarGrowthscale: "<<columnarGrowthscale[0]<<" "<<columnarGrowthscale[1]<<" "<<columnarGrowthscale[2]<<endl;

				if ( currGF->applyToColumnarLayer){
					if ((*itElement)->tissueType == 0){ //columnar layer, grow directly
						(*itElement)->updateGrowthRate(columnarGrowthscale[0],columnarGrowthscale[1],columnarGrowthscale[2]);
					}
					else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
						double weight = (*itElement)->getColumnarness();
						(*itElement)->updateGrowthRate(weight*columnarGrowthscale[0],weight*columnarGrowthscale[1],weight*columnarGrowthscale[2]);
					}
				}
				if ( currGF->applyToPeripodialMembrane){
					if ((*itElement)->tissueType == 1){ //peripodial membrane, grow directly
						(*itElement)->updateGrowthRate(peripodialGrowthscale[0],peripodialGrowthscale[1],peripodialGrowthscale[2]);
					}
					else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
						double weight = (*itElement)->getPeripodialness();
						(*itElement)->updateGrowthRate(weight*peripodialGrowthscale[0],weight*peripodialGrowthscale[1],weight*peripodialGrowthscale[2]);
					}
				}
				delete[] columnarReletivePos;
				delete[] peripodialRelativePos;
			}
		}
	}
	//cout<<"Finished growth rate calculation"<<endl;
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
	double StretchTimeSteps = StretchEndStep - StretchInitialStep;
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
					for (int m=0; m<leftClampBorder.size(); m++){
						if (leftClampBorder[m] == nodeIds[j]){
							alreadyRecorded = true;
							break;
						}
					}
					if (!alreadyRecorded){
						for (int m=0; m<rightClampBorder.size(); m++){
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
			if(leftHandSide){
				cout<<"Element is on left hand side"<<endl;
				for (int k=0; k<clampedBorderNodes.size(); k++){
					leftClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
			else {
				cout<<"Element is on right hand side"<<endl;
				for (int k=0; k<clampedBorderNodes.size(); k++){
					rightClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
		}
	}
	for (int k=0; k<leftClampBorder.size(); k++){
		cout<<"left clamp border nodes: "<<leftClampBorder[k]<<endl;
	}
	for (int k=0; k<rightClampBorder.size(); k++){
		cout<<"right clamp border nodes: "<<rightClampBorder[k]<<endl;
	}
}

void Simulation::recordForcesOnClampBorders(){
	if (recordForcesOnFixedNodes){
		for (int i=0;i<3;i++){
			leftClampForces[i] = 0.0;
			rightClampForces[i] = 0.0;
		}
		for (int k=0; k<leftClampBorder.size(); k++){
			for (int j=0; j<3; ++j){
				leftClampForces[j] += FixedNodeForces[leftClampBorder[k]][j];
			}
		}
		for (int k=0; k<rightClampBorder.size(); k++){
			for (int j=0; j<3; ++j){
				rightClampForces[j] += FixedNodeForces[rightClampBorder[k]][j];
			}
		}
	}
	outputFile<<"Forces on clamps lhs: "<<leftClampForces[0]<<" "<<leftClampForces[1]<<" "<<leftClampForces[2]<<" rhs: "<<rightClampForces[0]<<" "<<rightClampForces[1]<<" "<<rightClampForces[2]<<endl;
}

void Simulation::moveClampedNodesForStretcher(){
	if (timestep>=StretchInitialStep && timestep<StretchEndStep){
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
				fixAllD(i);
			}
		}
	}
	else{
		for (int i=0; i<nNodes; ++i){
			//fix apical nodes and all peripodial membrane nodes:
			if(Nodes[i]->tissuePlacement == 1 || Nodes[i]->tissueType == 1) {
				fixAllD(i);
			}
		}
	}
}

void Simulation::addMyosinForces(gsl_matrix* gExt){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    int* nodeIds = (*itElement)->getNodeIds();
	    int nNodes= (*itElement)->getNodeNumber();
	    double sumMyo[3] = {0.0,0.0,0.0};
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
    else if (timestep == 1) {SuctionPressure[2] = 0;}
    else if (timestep == 2) {SuctionPressure[2] = 100;}
    else if (timestep == 3) {SuctionPressure[2] = 200;}
    else if (timestep == 4) {SuctionPressure[2] = 300;}
    else if (timestep == 5) {SuctionPressure[2] = 400;}
    else if (timestep == 6) {SuctionPressure[2] = 500;}

    //cout<<"in add pipette forces, pipette pos: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<endl;
    int dim = 3;
	for (int i=0; i<nNodes; ++i){
        cout<<"Node "<<i<<" z pos: "<<Nodes[i]->Position[2]<<" effectLimitsInZ: "<<effectLimitsInZ[0]<<" "<<effectLimitsInZ[1]<<endl;
		if (Nodes[i]->Position[2]> effectLimitsInZ[0] &&  Nodes[i]->Position[2]< effectLimitsInZ[1]){
            cout<<"Node "<<i<<" is within z range"<<endl;
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
		(*itElement)->updateElementVolumesAndTissuePlacementsForSave(Nodes);
	}
}

void Simulation::clearNodeMassLists(){
	for (int i=0 ;i<nNodes;++i){
		Nodes[i]->connectedElementIds.size();
		Nodes[i]->connectedElementIds.clear();
		Nodes[i]->connectedElementWeights.clear();
		Nodes[i]->mass=0.0;
		Nodes[i]->surface=0.0;
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
