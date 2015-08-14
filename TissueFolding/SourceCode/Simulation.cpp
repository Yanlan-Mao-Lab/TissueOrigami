
#include "Simulation.h"
#include "Prism.h"
#include "Triangle.h"
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
	PeripodialMembraneType = -1;
	lumenHeight = -20;
	BoundingBoxSize[0]=1000.0; BoundingBoxSize[1]=1000.0; BoundingBoxSize[2]=1000.0;
	ContinueFromSave = false;
    growthRotationUpdateFrequency = 60.0/dt;
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
	int n = Nodes.size();
	//4 RK steps
	for (int i = 0; i<4; ++i){
		//n nodes
		for (int j=0;j<n;++j){
			delete[] SystemForces[i][j];
			delete[] PackingForces[i][j];
		}
	}
	//4 RK steps
	for (int i = 0; i<4; ++i){
		delete[] SystemForces[i];
		delete[] PackingForces[i];
	}
	delete[] SystemForces;
    delete[] PackingForces;
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
	ApicalVisc = 10.0;
	BasalVisc = 100.0;
	memset(noiseOnPysProp,0,4*sizeof(int));
	// The default input is a calculated mesh of width 4 elements, each element being 2.0 unit high
	// and having 1.0 unit sides of the triangles.
	MeshType = 2;
	Row = 4;
	Column = Row-2;
	SideLength=1.0;
	zHeight = 2.0;
	ApicalNodeFix[0]= false;
	ApicalNodeFix[1]= false;
	BasalNodeFix[0]= false;
	BasalNodeFix[1]= false;
	nGrowthFunctions = 0;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
    GrowthSaved = true;
	ForcesSaved = true;
	VelocitiesSaved = true;
	PeripodialElasticity = 0.0;
	PeripodialViscosity = ApicalVisc;
	PeripodialThicnessScale = 1.0;
	lumenHeightScale = 0.3;
	DVRight = 0;
	DVLeft = 1;
	stretcherAttached = false;
	StretchInitialStep = -100;
	StretchEndStep = -100;
	PipetteSuction = false;
	StretchVelocity = 0.0;
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
	//This is updating positions from save, I am only interested in normal positions, no Runge-Kutta steps. This corresponds to RK step 4, RKId = 3
	updateElementPositions(3);
	calculateColumnarLayerBoundingBox();
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
	if(!ApicalNodeFix[0] && ApicalNodeFix[1]){
		cerr <<"Apical nodes: There is no option to fix x and y movement while keeping z free for pinned nodes, correct mesh options"<<endl;
		ApicalNodeFix[0] = true;
		return false;
	}
	if(!BasalNodeFix[0] && BasalNodeFix[1]){
		cerr <<"Basal nodes: There is no option to fix x and y movement while keeping z free for pinned nodes, correct mesh options"<<endl;
		BasalNodeFix[0] = true;
		return false;
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
		cout<<"mesh type : 4"<<endl;
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	Success = CalculateTissueHeight(); //Calculating how many layers the columnar layer has, and what the actual height is.
	if (!Success){
		return Success;
	}
	if (AddPeripodialMembrane){
		Success = addPeripodialMembraneToTissue();
	}
	if (!Success){
		return Success;
	}
	fillInNodeNeighbourhood();
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
    //calculateStiffnessMatrices();
    calculateShapeFunctionDerivatives();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	//for (int i=0; i<Nodes.size();++i){
	//	cout<<"Node: "<<i<<endl;
	//	Nodes[i]->displayConnectedElementIds();
	//	Nodes[i]->displayConnectedElementWeights();
	//}
	if (AddPeripodialMembrane){
		assignMassWeightsDueToPeripodialMembrane();
	}
	if (stretcherAttached){
		setStretch();
	}
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
	return Success;
}

void Simulation::assignNodeMasses(){
	int n= Elements.size();
	for (int i=0; i<n; ++i){
		Elements[i]->assignVolumesToNodes(Nodes);
		Elements[i]->assignSurfaceAreaToNodes(Nodes);
	}
}

void Simulation::assignConnectedElementsAndWeightsToNodes(){
	int n= Elements.size();
	for (int i=0; i<n; ++i){
		Elements[i]->assignElementToConnectedNodes(Nodes);
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
	//This is updating positions from save, I am only interested in normal positions, no Runge-Kutta steps. This corresponds to RK step 4, RKId = 3
	updateElementPositions(3);
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
	if(Nodes.size() != n){
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
	int n;
	saveFileToDisplayMesh >> n;
	cout<<"number of nodes: "<<n<<endl;
	Node* tmp_nd;
	for (int i=0; i<n; ++i){
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
	cout<<"number of nodes: "<<Nodes.size()<<endl;
}

void Simulation::initiateNodesFromMeshInput(){
	int n;
	vector <int> NodesToFix;
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
		if (atCircumference ==1) {
			tmp_nd->atCircumference = true;
			NodesToFix.push_back(Nodes.size()); //If there is any fixing of the circumference, then this is a node that should be fixed
		}
		Nodes.push_back(tmp_nd);
		delete[] pos;
	}
	if (BasalNodeFix[0] || BasalNodeFix[1] || ApicalNodeFix[0] || ApicalNodeFix[1] ){
		fixApicalBasalNodes(NodesToFix);
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
		else if (shapeType == 4){
			initiateTriangleFromMeshInput();
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
		else if (shapeType == 4){
			double height = 0.0;
			saveFileToDisplayMesh >> height;
			initiateTriangleFromSave(height);
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}
	}
	cout<<"number of elements: "<<Elements.size()<<endl;
}

bool Simulation::readElementDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	if (Elements.size() != n){
		cerr<<"The element number from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	for (int i=0; i<n; ++i){
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

void Simulation::initiateTriangleFromSave(double height){
	//inserts a new triangle at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[3];
	for (int i =0 ;i<3; ++i){
		NodeIds[i] = 0;
	}
	Triangle* TrianglePnt01;
	TrianglePnt01 = new Triangle(NodeIds, Nodes, currElementId, height);
	TrianglePnt01->updateShapeFromSave(saveFileToDisplayMesh);
	calculateSystemCentre();
	TrianglePnt01->AlignReferenceApicalNormalToZ(SystemCentre);  //correcting the alignment of the triangular element such that the apical side will be aligned with (+)ve z
	calculateSystemCentre();
	Elements.push_back(TrianglePnt01);
	currElementId++;
	delete[] NodeIds;
	//Elements[Elements.size()-1]->displayReferencePositions();
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
	currElementId++;
	delete[] NodeIds;
	//Elements[Elements.size()-1]->displayReferencePositions();
}

void Simulation::initiateTriangleFromMeshInput(){
	int* NodeIds;
	NodeIds = new int[3];
	for (int i =0 ;i<3; ++i){
		int savedId;
		saveFileToDisplayMesh >> savedId;
		NodeIds[i] = savedId;
	}
	double height;
	saveFileToDisplayMesh >> height;
	Triangle* TrianglePnt01;
	TrianglePnt01 = new Triangle(NodeIds, Nodes, currElementId, height);
	TrianglePnt01->updateReferencePositionMatrixFromMeshInput(saveFileToDisplayMesh);
	calculateSystemCentre();
	//CORRECT THIS LATER ON!!!
	//SystemCentre[2] -= 10; // The apical surface is assumed to look towards (-)ve z
	TrianglePnt01->AlignReferenceApicalNormalToZ(SystemCentre);  //correcting the alignment of the triangular element such that the apical side will be aligned with (+)ve z
	calculateSystemCentre();
	Elements.push_back(TrianglePnt01);
	currElementId++;
	delete[] NodeIds;
	//Elements[Elements.size()-1]->displayReferencePositions();
}

void Simulation::reInitiateSystemForces(int oldSize){
	//deleting the old system forces:
	for (int i = 0; i<4; ++i){
		for (int j=0;j<oldSize;++j){
			delete[] SystemForces[i][j];
		}
	}
	for (int i = 0; i<4; ++i){
		delete[] SystemForces[i];
	}
	delete[] SystemForces;

	//reinitiating with the new size:
	const int n = Nodes.size();
	//4 RK steps
	SystemForces =  new double**[4];
	PackingForces =  new double**[4];
	for (int i=0;i<4;++i){
		SystemForces[i] = new double*[n];
		PackingForces[i] = new double*[n];
		for (int j=0;j<n;++j){
			SystemForces[i][j] = new double[3];
			PackingForces[i][j] = new double[3];
			SystemForces[i][j][0]=0.0;
			SystemForces[i][j][1]=0.0;
			SystemForces[i][j][2]=0.0;
			PackingForces[i][j][0]=0.0;
			PackingForces[i][j][1]=0.0;
			PackingForces[i][j][2]=0.0;
		}
	}
}

void Simulation::updateForcesFromSave(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		saveFileToDisplayForce.read((char*) &SystemForces[0][i][0], sizeof SystemForces[0][i][0]);
		saveFileToDisplayForce.read((char*) &SystemForces[0][i][1], sizeof SystemForces[0][i][1]);
		saveFileToDisplayForce.read((char*) &SystemForces[0][i][2], sizeof SystemForces[0][i][2]);
	}
}

void Simulation::updateVelocitiesFromSave(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			saveFileToDisplayVel.read((char*) &Nodes[i]->Velocity[j], sizeof Nodes[i]->Velocity[j]);
		}
	}
}

void Simulation::updateTensionCompressionFromSave(){
	int n = Elements.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<6; ++j){
            double S = gsl_matrix_get(Elements[i]->Strain,j,0);
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set(Elements[i]->Strain,j,0,S);
        }
	}
}

void Simulation::updateGrowthFromSave(){
    int n = Elements.size();
    for (int i=0;i<n;++i){
        gsl_matrix * currFg = gsl_matrix_calloc(3,3);
        for (int j=0; j<3; ++j){
            for(int k=0; k<3; ++k){
                double Fgjk;
                saveFileToDisplayGrowth.read((char*) &Fgjk, sizeof Fgjk);
                gsl_matrix_set(currFg,j,k,Fgjk);
            }
        }
        Elements[i]->setFg(currFg);
        gsl_matrix_free(currFg);
    }
}

void Simulation::readTensionCompressionToContinueFromSave(){
	int n = Elements.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<6; ++j){
            double S;
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set(Elements[i]->Strain,j,0,S);
		}
	}
}

void Simulation::readGrowthToContinueFromSave(){
    int n = Elements.size();
    for (int i=0;i<n;++i){
        gsl_matrix * currFg = gsl_matrix_calloc(3,3);
        for (int j=0; j<3; ++j){
            for(int k=0; k<3; ++k){
                double Fgjk;
                saveFileToDisplayGrowth.read((char*) &Fgjk, sizeof Fgjk);
                gsl_matrix_set(currFg,j,k,Fgjk);
            }
        }
        Elements[i]->setFg(currFg);
        gsl_matrix_free(currFg);
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
	currElementId++;
	delete[] NodeIds;
}

void Simulation::initiateTriangleFromSaveForUpdate(int k, double height){
	//inserts a new triangle at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[3];
	for (int i =0 ;i<3; ++i){
		//the node ids will be updated in function: updateShapeFromSave
		NodeIds[i] = 0;
	}
	Triangle* TrianglePnt01;
	TrianglePnt01 = new Triangle(NodeIds, Nodes, currElementId, height);
	TrianglePnt01->updateShapeFromSave(saveFileToDisplayMesh);
	calculateSystemCentre();
	//CORRECT THIS LATER ON!!!
	SystemCentre[2] -= 10; // The apical surface is assumed to took towards (-)ve z
	TrianglePnt01->AlignReferenceApicalNormalToZ(SystemCentre);  //correcting the alignment of the triangular element such that the apical side will be aligned with (+)ve z
	calculateSystemCentre();
	vector<ShapeBase*>::iterator it = Elements.begin();
	it += k;
	Elements.insert(it,TrianglePnt01);
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
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		//This is updating positions from save, I am only interested in normal positions, no Runge-Kutta steps. This corresponds to RK step 4, RKId = 3
		Elements[i]->updatePositions(3, Nodes);
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
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
	calculateColumnarLayerBoundingBox();
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
	int n;
	//cout<<"Is save file open: "<<saveFileToDisplay.is_open()<<" is file good? "<<saveFileToDisplay.good()<<endl;
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
			delete[] pos;
		}
	}
	else{
		for (int i = 0; i<(currNodeNumber-n); ++i){
			Node* tmp_nd;
			tmp_nd = Nodes.back();
			Nodes.pop_back();
			delete tmp_nd;
		}
	}
	if ( currNodeNumber != n ){
		//the node number is change, I yupdated the node list, now I need to fix system forces:
		reInitiateSystemForces(currNodeNumber);
	}
	//cout<<"end of funciton, number of nodes from save file: "<<n <<" number of nodes on the vector: "<<Nodes.size()<<endl;
}

void  Simulation::updateElementStatesFromSave(){
	//cout<<"Updating element states from save"<<endl;
	int n;
	saveFileToDisplayMesh >> n;
	int currElementNumber = Elements.size();
	//The elements list is bigger than necessary, I am deleting elements form the end of the list
	//cout<<"number of elements from save file: "<<n <<" number of element on the Elements vector: "<<currElementNumber<<endl;
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
			if (shapeType == 4){
				//the new shape is a triangle, I will add it now;
				double height;
				saveFileToDisplayMesh >> height;
				initiateTriangleFromSaveForUpdate(i,height);
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
		if (shapeType == 4){
			int i = Elements.size()-1;
			//the new shape is a triangle, I will add it now;
			double height;
			saveFileToDisplayMesh >> height;
			initiateTriangleFromSaveForUpdate(i,height);
			currElementNumber = Elements.size();
		}
	}
}

void Simulation::removeElementFromEndOfList(){
	ShapeBase* tmp_element;
	tmp_element = Elements.back();
	//int shapetype = tmp_element->getShapeType();
	Elements.pop_back();
	delete tmp_element;
}

void Simulation::updateNodePositionsFromSave(){
	//cout<<"Updating node positions from save"<<endl;
	//I have already read the current number of nodes in function "updateNodeNumberFromSave",
	//and I updated the number of nodes where necessary
	//the "cursor" in the file progressed, and is at the beginning of positions now
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
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



bool Simulation::addPeripodialMembraneToTissue(){
	bool Success = true;
	Success = generateColumnarCircumferenceNodeList();
	if (!Success){
		return Success;
	}
	calculateSystemCentre();
	sortColumnarCircumferenceNodeList();
	//if (PeripodialMembraneType == 1){
		//2D triangular dome of peripodial membrane, attached to the midline of the tissue
		vector <int*> trianglecornerlist;
		double d=0.0, dummy =0.0;
		getAverageSideLength(dummy,d);	//first term will get you the average side length of the peripodial membrane elements, second is the columnar elements
		if (!Success){
			return Success;
		}
		lumenHeight = TissueHeight*lumenHeightScale;
		addPeripodialMembraneNodes(trianglecornerlist, TissueHeight, d);
		FillNodeAssociationDueToPeripodialMembrane();
		addPeripodialMembraneElements(trianglecornerlist, PeripodialThicnessScale*TissueHeight);
	//}
	/*else if (PeripodialMembraneType == 2){
		//2D triangular dome of peripodial membrane, attached to the apical side of the tissue
		vector <int*> trianglecornerlist;
		double d=0.0, dummy =0.0;
		getAverageSideLength(dummy,d);	//first term will get you the average side length of the peripodial membrane elements, second is the columnar elements
		if (!Success){
			return Success;
		}
		addPeripodialMembraneNodes(trianglecornerlist, TissueHeight, d);
		FillNodeAssociationDueToPeripodialMembrane();
		addPeripodialMembraneElements(trianglecornerlist, TissueHeight);
	}*/
	return Success;
}

bool Simulation::generateColumnarCircumferenceNodeList(){
	//generating a list of nodes that are at the circumference and at the basal surface
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
		if (Nodes[i]->atCircumference && Nodes[i]->tissuePlacement == 0){ // tissuePlacement = 0 -> basal node
			ColumnarCircumferencialNodeList.push_back(i);

		}
		if (Nodes[i]->atCircumference && Nodes[i]->tissuePlacement == 1){ // tissuePlacement = 0 -> apical node
			ApicalColumnarCircumferencialNodeList.push_back(i);
		}
	}
	n = ColumnarCircumferencialNodeList.size();
	if (n<=0){
		cerr<<"No circumferncial nodes indicated! Cannot generate PeripodialMembrane"<<endl;
		AddPeripodialMembrane = false;
		return false;
	}
	return true;
}

void Simulation::sortColumnarCircumferenceNodeList(){
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
				//also change the apical list!
				temp=ApicalColumnarCircumferencialNodeList[i-1];
				ApicalColumnarCircumferencialNodeList[i-1]=ApicalColumnarCircumferencialNodeList[i];
				ApicalColumnarCircumferencialNodeList[i]=temp;
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


void Simulation::calculateCentreOfNodes(double* centre){
	centre[0] = 0.0;
	centre[1] = 0.0;
	centre[2] = 0.0;
	int n = ColumnarCircumferencialNodeList.size();
	double x =0.0, y = 0.0, z= 0.0;
	for (int i=0; i<n; ++i){
		for (int j=0; j<Nodes[ColumnarCircumferencialNodeList[i]]->nDim; ++j){
			centre[j] += Nodes[ColumnarCircumferencialNodeList[i]]->Position[j];
		}
	}
	for (int j=0; j<3; ++j){
		centre[j] /= n;
	}
}

void Simulation::AddPeripodialMembraneCircumference(double height, int& index_begin, int &index_end){
	double zoffset = 0.0;
	if (PeripodialMembraneType == 0) {
		zoffset = height + lumenHeight; //the offset of the circumferential peripodial nodes from the basal layer, this is a hoovering membrane, it should be higher than the actual tissue
	}
	if (PeripodialMembraneType == 1) {
		//adding the nodes to midzone, the offset should be half the height
		zoffset = height/2.0;
	}
	if (PeripodialMembraneType == 2) {
		//adding the nodes to apical layer, (the forces still apply to whole column) the offset should be equal to height
		zoffset = height;
	}
	int n = ColumnarCircumferencialNodeList.size();
	for (int i=0; i<n; ++i){
		int index = ColumnarCircumferencialNodeList[i];
		double* pos;
		pos = new double[3];
		pos[0] = Nodes[index]->Position[0];
		pos[1] = Nodes[index]->Position[1];
		pos[2] = Nodes[index]->Position[2] + zoffset;
		Node* tmp_nd;
		tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is PeripodialMembrane, node type is basal
		Nodes.push_back(tmp_nd);
		if (i==0){index_begin = tmp_nd->Id;}
		else if (i==n-1){index_end = tmp_nd->Id;}
		PeripodialMembraneCircumferencialNodeList.push_back(tmp_nd->Id);
		tmp_nd->atPeripodialCircumference = true;
		//cerr<<"NodeId: "<<tmp_nd->Id<<" pos: "<<tmp_nd->Position[0]<<" "<<tmp_nd->Position[1]<<" "<<tmp_nd->Position[2]<<endl;
		delete[] pos;
	}
}

void Simulation::AddHorizontalRowOfPeripodialMembraneNodes(vector <int*> &trianglecornerlist, double d, int &index_begin, int &index_end){
	int tmp_begin=0, tmp_end=0;
	for (int i = index_begin; i<=index_end; ++i){
		int index1 = i;
		int index2 = i+1;
		if (i==index_end) {
			//the connected node is node zero for the last node on the list, for all else, it is (i+1)th node on the list
			index2 = index_begin;
		}
		double MidPoint[3];
		MidPoint[0] = 0.5 * (Nodes[index1]->Position[0] + Nodes[index2]->Position[0]);
		MidPoint[1] = 0.5 * (Nodes[index1]->Position[1] + Nodes[index2]->Position[1]);
		MidPoint[2] = 0.5 * (Nodes[index1]->Position[2] + Nodes[index2]->Position[2]);
		double vec[2] = { Nodes[index2]->Position[0] - Nodes[index1]->Position[0], Nodes[index2]->Position[1] - Nodes[index1]->Position[1]};
		//rotate the vector from node 0 to node 1 cw 90 degrees
		double* norm;
		norm = new double[3];
		norm[0] = vec[1];
		norm[1] = -1.0*vec[0];
		norm[2] = 0.0;
		//cout<<"	1 norm: "<<norm[0]<<" "<<norm[1]<<" "<<norm[2]<<endl;
		//normalise the vector, and scale with average side length:
		double mag = pow((norm[0]*norm[0] + norm[1]*norm[1]),0.5);
		mag /=d;
		norm[0] /= mag; norm[1] /=mag;
		//calculate x,y position of the node:
		norm[0] += MidPoint[0];
		norm[1] += MidPoint[1];
		norm[2] += MidPoint[2];
		Node* tmp_nd;
		//Adding Peroipodium node:
		tmp_nd = new Node(Nodes.size(), 3, norm,0,1);  //tissue type is PeripodialMembrane, node type is basal
		Nodes.push_back(tmp_nd);
		if (i==index_begin){tmp_begin = tmp_nd->Id;}
		else if (i==index_end){tmp_end = tmp_nd->Id;}
		int* TriNodeIds0;
		TriNodeIds0 = new int[3];
		TriNodeIds0[0]=index1;
		TriNodeIds0[1]=index2;
		TriNodeIds0[2]=tmp_nd->Id;
		trianglecornerlist.push_back(TriNodeIds0);
		//Adding node to list prior to creation of the node itself, this is the second layer of triangles
		int* TriNodeIds1;
		TriNodeIds1 = new int[3];
		TriNodeIds1[0]=tmp_nd->Id;
		TriNodeIds1[1]=index2;
		if (index2 ==index_begin){TriNodeIds1[2]=index_end+1;}
		else{TriNodeIds1[2]=tmp_nd->Id+1;}
		trianglecornerlist.push_back(TriNodeIds1);

		//cerr<<"Triangle Corner Ids: "<<trianglecornerlist[trianglecornerlist.size()-1][0]<<" "<<trianglecornerlist[trianglecornerlist.size()-1][1]<<" "<<trianglecornerlist[trianglecornerlist.size()-1][2]<<endl;
		delete[] norm;
	}
	index_begin = tmp_begin;
	index_end = tmp_end;
}

void Simulation::AddVerticalRowOfPeripodialMembraneNodes(int& layerCount, int nLayers, vector <int*> &trianglecornerlist, double height, int &index_begin, int &index_end){
	int tmp_begin=0, tmp_end=0;
	int counter = 0;
	for (int i = index_begin; i<=index_end; ++i){
		int index1 = i;
		int index2 = i+1;
		if (i==index_end) {
			//the connected node is node zero for the last node on the list, for all else, it is (i+1)th node on the list
			index2 = index_begin;
		}
		double MidPoint[3];
		MidPoint[0] = 0.5 * (Nodes[index1]->Position[0] + Nodes[index2]->Position[0]);
		MidPoint[1] = 0.5 * (Nodes[index1]->Position[1] + Nodes[index2]->Position[1]);
		MidPoint[2] = 0.5 * (Nodes[index1]->Position[2] + Nodes[index2]->Position[2]);
		double targetpoint[3] = {0.0,0.0,0.0};
		int n = ColumnarCircumferencialNodeList.size();
		if(layerCount %2 ==0 ){ //even layer count, starting from zero, the target is the node on circumference
			int idx = counter + layerCount*0.5 +1;
			//cout<<"layerCount: "<<layerCount<<" index1 "<<index1<<" counter: "<<counter<<endl;
			if (idx >= n){
				idx -= n;
			}
			int targetindex = ColumnarCircumferencialNodeList[idx];
			targetpoint[0] = Nodes[targetindex]->Position[0];
			targetpoint[1] = Nodes[targetindex]->Position[1];
			targetpoint[2] = Nodes[targetindex]->Position[2] + height + lumenHeight*0.5 ;
		}
		else{
			//the target point should be the midpoint of the two connected nodes:
			int idx1 = 0;
			int idx2 = 0;
			idx1 = counter + (layerCount-1)*0.5 +1;
			idx2 = idx1+1;
			if (idx1 >= n){
				idx1 -= n;
			}
			if (idx2 >= n){
				idx2 -= n;
			}
			int targetindex1 = ColumnarCircumferencialNodeList[idx1];
			int targetindex2 = ColumnarCircumferencialNodeList[idx2];
			targetpoint[0] = 0.5 * (Nodes[targetindex1]->Position[0] + Nodes[targetindex2]->Position[0] );
			targetpoint[1] = 0.5 * (Nodes[targetindex1]->Position[1] + Nodes[targetindex2]->Position[1] );
			targetpoint[2] = 0.5 * (Nodes[targetindex1]->Position[2] + Nodes[targetindex2]->Position[2] ) + height + lumenHeight*0.5 ;
		}
		double vec[3] = {targetpoint[0] - MidPoint[0], targetpoint[1] - MidPoint[1], targetpoint[2] - MidPoint[2]};
		double mag = pow((vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]),0.5);
		vec[0] /= mag;
		vec[1] /= mag;
		vec[2] /= mag;
		//The height spanned in this direction should be (half height of the tissue + half lumen height) / number of layers:
		double z = (height*0.5 + lumenHeight*0.5) / (float) nLayers;
		double multiplier = z / vec[2];
		double* norm;
		norm = new double[3];
		norm[0] = MidPoint[0]+vec[0]*multiplier;
		norm[1] = MidPoint[1]+vec[1]*multiplier;
		norm[2] = MidPoint[2]+z;
		//cout<<"layerCount: "<<layerCount<<"  MidPoint          : "<<MidPoint[0]<<" "<<MidPoint[1]<<" "<<MidPoint[2]<<endl;
		//cout<<"layerCount: "<<layerCount<<"  z: "<<z<<" vec  : "<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<endl;
		//cout<<"layerCount: "<<layerCount<<" norm for added node: "<<norm[0]<<" "<<norm[1]<<" "<<norm[2]<<endl;
		Node* tmp_nd;
		//Adding Peroipodium node:
		tmp_nd = new Node(Nodes.size(), 3, norm,0,1);  //tissue type is PeripodialMembrane, node type is basal
		Nodes.push_back(tmp_nd);
		if (i==index_begin){tmp_begin = tmp_nd->Id;}
		else if (i==index_end){tmp_end = tmp_nd->Id;}
		int* TriNodeIds0;
		TriNodeIds0 = new int[3];
		TriNodeIds0[0]=index1;
		TriNodeIds0[1]=index2;
		TriNodeIds0[2]=tmp_nd->Id;
		trianglecornerlist.push_back(TriNodeIds0);
		//Adding node to list prior to creation of the node itself, this is the second layer of triangles
		int* TriNodeIds1;
		TriNodeIds1 = new int[3];
		TriNodeIds1[0]=tmp_nd->Id;
		TriNodeIds1[1]=index2;
		if (index2 ==index_begin){TriNodeIds1[2]=index_end+1;}
		else{TriNodeIds1[2]=tmp_nd->Id+1;}
		trianglecornerlist.push_back(TriNodeIds1);

		//cerr<<"Triangle Corner Ids: "<<trianglecornerlist[trianglecornerlist.size()-1][0]<<" "<<trianglecornerlist[trianglecornerlist.size()-1][1]<<" "<<trianglecornerlist[trianglecornerlist.size()-1][2]<<endl;
		counter ++;
		delete[] norm;
	}
	index_begin = tmp_begin;
	index_end = tmp_end;
	layerCount++;
}

void Simulation::AddPeripodialMembraneCapToMidAttached(int layerCount,  vector <int*> &trianglecornerlist, double height, int index_begin, int index_end){
	//Now I have the indices of the nodes specifying the last row.
	//I want to cap the tissue, with the topology of the apical surfaces of the columnar layer
	vector <int> PeripodialMembraneNodeId;
	vector <int> CorrespondingApicalNodeId;
	int n = ApicalColumnarCircumferencialNodeList.size();
	//map the circumference to PeripodialMembrane nodes:
	int counter =0;
	for (int i = index_begin; i <= index_end; ++i){
		int idx = counter + layerCount*0.5 +1;
		if (idx >= n){
			idx -= n;
		}
		counter++;
		PeripodialMembraneNodeId.push_back(Nodes[i]->Id);
		CorrespondingApicalNodeId.push_back(Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id);
		//cout<<"PeripodialMembrane Node: "<<Nodes[i]->Id<<" pos: "<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<" corr. node id: "<<Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id<<endl;
	}
	// The end point is 0.5*lumenHeight above the apical surface.
	//I want to add the remaining 50% height of the lumen as a curvature coered by the first layer of PeripodialMembrane nodes
	// the layer count is equal to (nLayers -1 ) at the moment, as I am at the topmost layer
	// I can obtain the increment I need from adding the sum of these as a z offset:
	//double zOffset = height*0.5 + height / (float) (layerCount+1);
	double zOffset = lumenHeight;
	for (int i = 0; i<Nodes.size(); ++i){
		if (Nodes[i]->tissuePlacement == 1){ //Node is apical
			int id = Nodes[i]->Id;
			bool AtCircumference = false;
			for (int j =0 ; j< n; ++j){
				if (id == Nodes[ApicalColumnarCircumferencialNodeList[j]]->Id ){
					//if it is not on the circumference of PeripodialMembrane, skip
					AtCircumference = true;
					break;
				}
			}
			if (!AtCircumference){
				double* pos;
				pos = new double[3];
				pos[0] = Nodes[i]->Position[0];
				pos[1] = Nodes[i]->Position[1];
				pos[2] = Nodes[i]->Position[2] + zOffset;
				Node* tmp_nd;
				//Adding Peroipodium node:
				tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is PeripodialMembrane, node type is basal
				Nodes.push_back(tmp_nd);
				PeripodialMembraneNodeId.push_back(tmp_nd->Id);
				CorrespondingApicalNodeId.push_back(Nodes[i]->Id);
				//cout<<"temp Node: "<<tmp_nd->Id<<" pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" corr. node id: "<<Nodes[i]->Id<<endl;
			}
		}
	}
	/*cout<<"apical node - PeripodialMembrane node couples: "<<endl;
	for (int a = 0; a< PeripodialMembraneNodeId.size(); a++){
		cout<<PeripodialMembraneNodeId[a]<<" "<<CorrespondingApicalNodeId[a]<<endl;
	}*/
	//Now I have generated the nodes and have the corresponding node mapping, I can add the triangles to list:
	int initialpointfordisplay = trianglecornerlist.size();
	for (int i = 0; i<Elements.size(); ++i){
		vector <int> ApicalTriangles;
		Elements[i]->getApicalTriangles(ApicalTriangles);
		int nList = ApicalTriangles.size();
		for (int k=0; k<nList-2; ++k){
			if (Nodes[ApicalTriangles[k]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+1]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+2]]->tissuePlacement == 1){
				int* TriNodeIds;
				TriNodeIds = new int[3];
				int nDict = CorrespondingApicalNodeId.size();
				for (int dictionaryIndex =0; dictionaryIndex<nDict; ++dictionaryIndex){
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k]){
						TriNodeIds[0]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+1]){
						TriNodeIds[1]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+2]){
						TriNodeIds[2]=PeripodialMembraneNodeId[dictionaryIndex];
					}
				}
				trianglecornerlist.push_back(TriNodeIds);
			}
		}
	}
	/*cout<<"triangle corners for peripodial membrane cap: "<<endl;
	for (int a =initialpointfordisplay; a< trianglecornerlist.size(); a++){
		cout<<trianglecornerlist[a][0]<<" "<<trianglecornerlist[a][1]<<" "<<trianglecornerlist[a][2]<<endl;
	}*/
}



void Simulation::AddPeripodialMembraneCapToApicalAttached(int layerCount,  vector <int*> &trianglecornerlist, double height, int index_begin, int index_end){
	//Now I have the indices of the nodes specifying the last row.
	//I want to cap the tissue, with the topology of the apical surfaces of the columnar layer
	vector <int> PeripodialMembraneNodeId;
	vector <int> CorrespondingApicalNodeId;
	int n = ApicalColumnarCircumferencialNodeList.size();
	//map the circumference to peripodial membrane nodes:
	int counter =0;
	for (int i = index_begin; i <= index_end; ++i){
		PeripodialMembraneNodeId.push_back(Nodes[i]->Id);
		CorrespondingApicalNodeId.push_back(Nodes[ApicalColumnarCircumferencialNodeList[i-index_begin]]->Id);
		//cout<<"Peripodial membrane Node: "<<Nodes[i]->Id<<" pos: "<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<" corr. node id: "<<Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id<<endl;
	}
	// The end point is on hte apical surface
	//I want to add the height of the lumen as a curvature covered by the first layer of peripodial membrane nodes
	double zOffset = lumenHeight;
	for (int i = 0; i<Nodes.size(); ++i){
		if (Nodes[i]->tissuePlacement == 1){ //Node is apical
			int id = Nodes[i]->Id;
			bool AtCircumference = false;
			for (int j =0 ; j< n; ++j){
				if (id == Nodes[ApicalColumnarCircumferencialNodeList[j]]->Id ){
					//if it is not on the circumference of Peripodial membrane, skip
					AtCircumference = true;
					break;
				}
			}
			if (!AtCircumference){
				double* pos;
				pos = new double[3];
				pos[0] = Nodes[i]->Position[0];
				pos[1] = Nodes[i]->Position[1];
				pos[2] = Nodes[i]->Position[2] + zOffset;
				Node* tmp_nd;
				//Adding Peroipodium node:
				tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is Peripodial Membrane, node type is basal
				Nodes.push_back(tmp_nd);
				PeripodialMembraneNodeId.push_back(tmp_nd->Id);
				CorrespondingApicalNodeId.push_back(Nodes[i]->Id);
				//cout<<"temp Node: "<<tmp_nd->Id<<" pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" corr. node id: "<<Nodes[i]->Id<<endl;
			}
		}
	}
	/*cout<<"apical node - peripodial membrane node couples: "<<endl;
	for (int a = 0; a< PeripodialMembraneNodeId.size(); a++){
		cout<<PeripodialMembraneNodeId[a]<<" "<<CorrespondingApicalNodeId[a]<<endl;
	}*/
	//Now I have generated the nodes and have the corresponding node mapping, I can add the triangles to list:
	int initialpointfordisplay = trianglecornerlist.size();
	for (int i = 0; i<Elements.size(); ++i){
		vector <int> ApicalTriangles;
		Elements[i]->getApicalTriangles(ApicalTriangles);
		int nList = ApicalTriangles.size();
		for (int k=0; k<nList-2; ++k){
			if (Nodes[ApicalTriangles[k]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+1]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+2]]->tissuePlacement == 1){
				//This is an element that is adding a peripodial membrane triangle,
				//This will be the base node of the said triangle, if the triangle is tilted
				//I will record all the information now, and apply the changes to tilted triangles i
				int* TriNodeIds;
				TriNodeIds = new int[3];
				int nDict = CorrespondingApicalNodeId.size();
				for (int dictionaryIndex =0; dictionaryIndex<nDict; ++dictionaryIndex){
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k]){
						TriNodeIds[0]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+1]){
						TriNodeIds[1]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+2]){
						TriNodeIds[2]=PeripodialMembraneNodeId[dictionaryIndex];
					}
				}
				trianglecornerlist.push_back(TriNodeIds);
			}
		}
	}
	/*cout<<"triangle corners for Peripodial membrane cap: "<<endl;
	for (int a =initialpointfordisplay; a< trianglecornerlist.size(); a++){
		cout<<trianglecornerlist[a][0]<<" "<<trianglecornerlist[a][1]<<" "<<trianglecornerlist[a][2]<<endl;
	}*/
}

void Simulation::AddPeripodialMembraneCapToHoovering(int layerCount,  vector <int*> &trianglecornerlist, double height, int index_begin, int index_end){
	//Now I have the indices of the nodes specifying the last row.
	//I want to cap the tissue, with the topology of the apical surfaces of the columnar layer
	vector <int> PeripodialMembraneNodeId;
	vector <int> CorrespondingApicalNodeId;
	int n = ApicalColumnarCircumferencialNodeList.size();
	//map the circumference to peripodial membrane nodes:
	int counter =0;
	for (int i = index_begin; i <= index_end; ++i){
		PeripodialMembraneNodeId.push_back(Nodes[i]->Id);
		CorrespondingApicalNodeId.push_back(Nodes[ApicalColumnarCircumferencialNodeList[i-index_begin]]->Id);
		//cout<<"Peripodial membrane Node: "<<Nodes[i]->Id<<" pos: "<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<" corr. node id: "<<Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id<<endl;
	}
	// The end point is on hte apical surface
	//I want to add the height of the lumen as a curvature covered by the first layer of peripodial membrane nodes
	double zOffset = lumenHeight;
	for (int i = 0; i<Nodes.size(); ++i){
		if (Nodes[i]->tissuePlacement == 1){ //Node is apical
			int id = Nodes[i]->Id;
			bool AtCircumference = false;
			for (int j =0 ; j< n; ++j){
				if (id == Nodes[ApicalColumnarCircumferencialNodeList[j]]->Id ){
					//if it is not on the circumference of peripodial membrane, skip
					AtCircumference = true;
					break;
				}
			}
			if (!AtCircumference){
				double* pos;
				pos = new double[3];
				pos[0] = Nodes[i]->Position[0];
				pos[1] = Nodes[i]->Position[1];
				pos[2] = Nodes[i]->Position[2] + zOffset;
				Node* tmp_nd;
				//Adding Peroipodium node:
				tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is peripodial membrane, node type is basal
				Nodes.push_back(tmp_nd);
				PeripodialMembraneNodeId.push_back(tmp_nd->Id);
				CorrespondingApicalNodeId.push_back(Nodes[i]->Id);
				//cout<<"temp Node: "<<tmp_nd->Id<<" pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" corr. node id: "<<Nodes[i]->Id<<endl;
			}
		}
	}
	/*cout<<"apical node - peripodial membrane node couples: "<<endl;
	for (int a = 0; a< PeripodialMembranemNodeId.size(); a++){
		cout<<PeripodialMembraneNodeId[a]<<" "<<CorrespondingApicalNodeId[a]<<endl;
	}*/
	//Now I have generated the nodes and have the corresponding node mapping, I can add the triangles to list:
	int initialpointfordisplay = trianglecornerlist.size();
	for (int i = 0; i<Elements.size(); ++i){
		vector <int> ApicalTriangles;
		Elements[i]->getApicalTriangles(ApicalTriangles);
		int nList = ApicalTriangles.size();
		for (int k=0; k<nList-2; ++k){
			if (Nodes[ApicalTriangles[k]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+1]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+2]]->tissuePlacement == 1){
				//This is an element that is adding a peripodial membrane triangle,
				//This will be the base node of the said triangle, if the triangle is tilted
				//I will record all the information now, and apply the changes to tilted triangles i
				int* TriNodeIds;
				TriNodeIds = new int[3];
				int nDict = CorrespondingApicalNodeId.size();
				for (int dictionaryIndex =0; dictionaryIndex<nDict; ++dictionaryIndex){
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k]){
						TriNodeIds[0]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+1]){
						TriNodeIds[1]=PeripodialMembraneNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+2]){
						TriNodeIds[2]=PeripodialMembraneNodeId[dictionaryIndex];
					}
				}
				trianglecornerlist.push_back(TriNodeIds);
			}
		}
	}
	/*cout<<"triangle corners for peripodial membrane cap: "<<endl;
	for (int a =initialpointfordisplay; a< trianglecornerlist.size(); a++){
		cout<<trianglecornerlist[a][0]<<" "<<trianglecornerlist[a][1]<<" "<<trianglecornerlist[a][2]<<endl;
	}*/
}


void Simulation::FillNodeAssociationDueToPeripodialMembrane(){
	//Take the peripodial membrane circumferential list
	//Start from the associated Basal columnar circumferential node
	//Move apically through the elements until you reach the apical surface - an apical node
	int n = PeripodialMembraneCircumferencialNodeList.size();
	int nE = Elements.size();
	for (int i=0; i< n; ++i){
		int index = PeripodialMembraneCircumferencialNodeList[i];
		int currNodeId = ColumnarCircumferencialNodeList[i];
		Nodes[index]->AssociatedNodesDueToPeripodialMembrane.push_back(ColumnarCircumferencialNodeList[i]);
		Nodes[ColumnarCircumferencialNodeList[i]]->LinkedPeripodialNodeId = Nodes[index]->Id;
		while(Nodes[currNodeId]->tissuePlacement != 1){ //While I have not reached the apical node
			for (int j= 0; j<nE; ++j){
				bool IsBasalOwner = Elements[j]->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					currNodeId = Elements[j]->getCorrecpondingApical(currNodeId);
					break;
				}
			}
			Nodes[index]->AssociatedNodesDueToPeripodialMembrane.push_back(currNodeId);
			if (currNodeId != Nodes[currNodeId]->Id ){cerr<<"Error in node association index"<<endl;}
			Nodes[currNodeId]->LinkedPeripodialNodeId = Nodes[index]->Id;
		}
		//cout<<"Node: "<<index<<" pos : "<<Nodes[index]->Position[0]<<" "<<Nodes[index]->Position[1]<<" "<<Nodes[index]->Position[2]<<endl;
		//for(int j = 0; j< Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size(); ++j){
		//	int a = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
		//	cout<<"	Associated Node id: "<<Nodes[a]->Id<<" Pos: "<<Nodes[a]->Position[0]<<" "<<Nodes[a]->Position[1]<<" "<<Nodes[a]->Position[2]<<endl;
		//}
	}
}


void Simulation::assignMassWeightsDueToPeripodialMembrane(){
	//Take the peripodial membrane circumferential list
	//calculate the sum of associated node masses
	//Add the weigthing fractions to AssociatedNodeWeightsDueToPeripodialMembrane of each Node
	int n = PeripodialMembraneCircumferencialNodeList.size();
	for (int i=0; i< n; ++i){
		int index = PeripodialMembraneCircumferencialNodeList[i];
		int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
		double weightSum = 0.0;
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
			double w = Nodes[index2]->mass;
			weightSum += w;
			Nodes[index]->AssociatedNodeWeightsDueToPeripodialMembrane.push_back(w);
		}
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
			Nodes[index]->AssociatedNodeWeightsDueToPeripodialMembrane[j] /= weightSum;
			//Distributing the weight of this node onto the associated nodes
			Nodes[index2]->mass += Nodes[index]->mass*Nodes[index]->AssociatedNodeWeightsDueToPeripodialMembrane[j];
		}
	}
}


void Simulation::addPeripodialMembraneNodes(vector <int*> &trianglecornerlist, double height, double d){
	//cerr<<"Adding peripodial membrane nodes"<<endl;
	//int n = ColumnarCircumferencialNodeList.size();
	int index_begin = 0, index_end =0;
	//Adding a midline range of nodes
	AddPeripodialMembraneCircumference(height, index_begin, index_end);
	int layerCount = 0;
	if (PeripodialMembraneType == 0){
		AddPeripodialMembraneCapToHoovering(layerCount, trianglecornerlist, height, index_begin, index_end);
	}
	else if (PeripodialMembraneType == 1){
		double triangleHeight = 0.866*d; //0.866 is square-root(3)/2, this is the height of the triangle I am adding,
		AddHorizontalRowOfPeripodialMembraneNodes(trianglecornerlist, triangleHeight, index_begin, index_end);
		//calculating how many layers of triangles I need for spanning the necessary height, I want the side layer to go up
		// 50% of tissue height (to reach the same level as the top of columnar layer)
		// + 50% of the lumen size. The remaining 50% of lumen size will be spanned by the first row of cap elements.
		int nLayers = ceil((height*0.5 + lumenHeight*0.5) / triangleHeight );
		//I need an odd number of layers, so that the final shape will be the same as the circumference of columnar layer:
		if (nLayers %2 == 0){
			nLayers --;
		}
		//now I want nLayers many layers to cover exactly the height of the tissue, the distance should be the hipotenus divided by nLayers
		//I will adjust the triangle height
		triangleHeight = pow((triangleHeight*triangleHeight + height*height),0.5) / (float)nLayers;
		for (int i=0; i<nLayers; ++i){
			AddVerticalRowOfPeripodialMembraneNodes(layerCount, nLayers, trianglecornerlist, height, index_begin, index_end);
		}
		AddPeripodialMembraneCapToMidAttached(layerCount, trianglecornerlist, height, index_begin, index_end);
	}
	else if (PeripodialMembraneType == 2){
		AddPeripodialMembraneCapToApicalAttached(layerCount, trianglecornerlist, height, index_begin, index_end);
	}
	//AssignAssociatedNodesToPeripodialMembraneCircumference();
}

void Simulation::getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength){
	double dsumPeri =0.0, dsumCol = 0.0;
	int colCounter =0, periCounter=0;
	int n = Elements.size();
	for(int i=0; i<n; ++i){
		if (!Elements[i]->IsAblated){
			//do not count ablated elements
			if (Elements[i]->tissueType==0){ //element belongs to columnar layer
				dsumCol += Elements[i]->getApicalSideLengthAverage();
				colCounter++;
			}
			else{
				dsumPeri += Elements[i]->getApicalSideLengthAverage();
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

bool Simulation::CalculateTissueHeight(){
	//check if the columnar layer is made of 3D elements:
	bool ColumnarLayer3D = isColumnarLayer3D();
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

void Simulation::addPeripodialMembraneElements(vector <int*> &trianglecornerlist, double height){
	int n = trianglecornerlist.size();
	calculateSystemCentre();
	//SystemCentre[2] -= 10; // The apical surface is assumed to took towards (-)ve z
	for (int i=0; i<n; ++i){
		int* NodeIds;
		NodeIds = new int[3];
		NodeIds[0] = trianglecornerlist[i][0];
		NodeIds[1] = trianglecornerlist[i][1];
		NodeIds[2] = trianglecornerlist[i][2];
		Triangle* TrianglePnt01;
		//cout<<"NodeIds: "<<NodeIds[0]<<" "<<NodeIds[1]<<" "<<NodeIds[2]<<endl;
		TrianglePnt01 = new Triangle(NodeIds, Nodes, currElementId, height);
		TrianglePnt01->AlignReferenceApicalNormalToZ(SystemCentre);  //correcting the alignment of the triangular element such that the apical side will be aligned with (+)ve z
		Elements.push_back(TrianglePnt01);
		currElementId++;
		delete[] NodeIds;
	}
	calculateSystemCentre();
}

void Simulation::calculateStiffnessMatrices(){
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		cout<<" setting up element :  "<<i<<" of "<<n<<endl;
		Elements[i]->calculateReferenceStiffnessMatrix();
	}
}

void Simulation::calculateShapeFunctionDerivatives(){
    int n = Elements.size();
    for (int i=0; i<n; ++i){
        cout<<" setting up element :  "<<i<<" of "<<n<<endl;
        Elements[i]->calculateElementShapeFunctionDerivatives();
    }
}

void Simulation::fixAllD(int i){
	for (int j =0 ; j<Nodes[i]->nDim; ++j){
		Nodes[i]->FixedPos[j]=true;
	}
}

void Simulation::fixZ(int i){
	if(Nodes[i]->nDim>2){
		Nodes[i]->FixedPos[2]=true;
	}
	else{
		cerr<<"ERROR: Node : "<<i<<" does not have z-dimension"<<endl;
	}
}

void Simulation::zeroForcesOnNode(int RKId, int i){
	double ForceBalance[3];
	ForceBalance[0] = SystemForces[RKId][i][0];
	ForceBalance[1] = SystemForces[RKId][i][1];
	ForceBalance[2] = SystemForces[RKId][i][2];
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		SystemForces[RKId][i][0]-=ForceBalance[0];
		SystemForces[RKId][i][1]-=ForceBalance[1];
		SystemForces[RKId][i][2]-=ForceBalance[2];
	}
}

void Simulation::initiateSystemForces(){
	const int n = Nodes.size();
	//4 RK steps
	SystemForces = new double**[4];
	PackingForces = new double**[4];
	for (int i=0;i<4;++i){
		//n nodes
		SystemForces[i] = new double*[n];
		PackingForces[i] = new double*[n];
		for (int j=0;j<n;++j){
			//3 dimensions
			SystemForces[i][j] = new double[3];
			PackingForces[i][j] = new double[3];
			SystemForces[i][j][0]=0.0;
			SystemForces[i][j][1]=0.0;
			SystemForces[i][j][2]=0.0;
			PackingForces[i][j][0]=0.0;
			PackingForces[i][j][1]=0.0;
			PackingForces[i][j][2]=0.0;
			//cout<<"systemforces[i][j]: "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<endl;
		}
	}
}

void Simulation::initiateSinglePrismNodes(float zHeight){
	double *pos = new double[3];
	Node* tmp_nd;
	pos[0]=0;pos[1]=1;pos[2]=0;
	tmp_nd = new Node(0, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	pos[0]=1;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(1, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(2, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=1;pos[2]=zHeight;
	tmp_nd = new Node(3, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	pos[0]=1;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(4, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(5, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
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
	currElementId++;
	fixZ(0);
	fixZ(1);
	fixZ(2);
}


void Simulation::initiateNodesByRowAndColumn(int Row, int Column, float SideLength, float zHeight){
	DVRight = Row;
	DVLeft = 0;
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
	}
	//Adding the apical level, all will form columnar elements:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = zHeight;
		tmp_nd = new Node(n+i, 3, pos,1,0);
		Nodes.push_back(tmp_nd);
	}

	if (BasalNodeFix[0] || BasalNodeFix[1] || ApicalNodeFix[0] || ApicalNodeFix[1] ){
		fixApicalBasalNodes(NodesToFix);
	}

	delete[] pos;
}

void Simulation::fixApicalBasalNodes(vector<int> &NodesToFix){
	int nNodesToFix = NodesToFix.size();
	int n = Nodes.size()/2;
	for (int i=0; i<nNodesToFix; ++i){
		if(BasalNodeFix[0]){
			//The nodes on the circumference on the basal side (bottom)
			//have their z position fixed.
			fixZ(NodesToFix[i]);
		}
		if(BasalNodeFix[1]){
			//The nodes on the circumference on the basal side (bottom)
			//have their x & y positions fixed.
			fixAllD(NodesToFix[i]);
		}
		if(ApicalNodeFix[0]){
			//The nodes on the circumference on the apical side (top)
			//have their z position fixed.
			fixZ(NodesToFix[i]+n);
		}
		if(ApicalNodeFix[1]){
			//The nodes on the circumference on the apical side (top)
			//have their x & y positions fixed.
			fixAllD(NodesToFix[i]+n);
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
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
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
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
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
	int n = Nodes.size();
	for (int i = 0; i< n; ++i){
		for (int j =0; j<Nodes[i]->nDim; ++j){
			SystemCentre[j] += Nodes[i]->Position[j];
		}
	}
	SystemCentre[0]= SystemCentre[0]/n;
	SystemCentre[1]= SystemCentre[1]/n;
	SystemCentre[2]= SystemCentre[2]/n;
}

void Simulation::assignPhysicalParameters(){
	for (int i=0; i<Elements.size();++i){
		double r = (rand() % 200) / 100.0;	//random number between 0.00 and 2.00
		r = r - 1.0; 						//random number between -1.00 and 1.00
		float noise1 = r*noiseOnPysProp[0];	//percent noise on current element
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise2 = r*noiseOnPysProp[1];
		if (Elements[i]->tissueType == 0){ //Element is on the columnar layer
			Elements[i]->setElasticProperties(EApical*(1 + noise1/100.0),EBasal*(1 + noise1/100.0),EMid*(1 + noise1/100.0),poisson*(1 + noise2/100));
		}
		if (Elements[i]->tissueType == 1){ //Element is on the peripodial membrane
			double currE = PeripodialElasticity*(1 + noise1/100.0);
			Elements[i]->setElasticProperties(currE,currE,currE,poisson*(1 + noise2/100));
		}
	}
	for (int i=0; i<Nodes.size(); ++i){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		Nodes[i]->setViscosity(ApicalVisc*noise3, BasalVisc*noise3, PeripodialViscosity*noise3);
	}
}

void Simulation::manualPerturbationToInitialSetup(bool deform, bool rotate){
    if(timestep==0){
        double scaleX = 1.0;
        double scaleY = 1.0;
        double scaleZ = 2.0;

        double PI = 3.14159265359;
        double tetX = -45 *PI/180.0;
        double tetY = -30 *PI/180.0;
        double tetZ = 90 *PI/180.0;
        if(deform){
            for(int i=0;i<Nodes.size();++i){
                Nodes[i]->Position[0] *=scaleX;
                Nodes[i]->Position[1] *=scaleY;
                Nodes[i]->Position[2] *=scaleZ;
            }
            Nodes[3]->Position[1] -= 15.0;
            Nodes[5]->Position[0] -= 15.0;
        }

        if (rotate){
            double R[3][3];
            double c = cos(tetX), s = sin(tetX);
            double Rx[3][3] = {{1,0,0},{0,c,-s},{0,s,c}};
            c = cos(tetY); s = sin(tetY);
            double Ry[3][3] = {{c,0,s},{0,1,0},{-s,0,c}};
            c = cos(tetZ), s = sin(tetZ);
            double Rz[3][3] = {{c,-s,0},{s,c,0},{0,0,1}};
            for (int j =0; j<3;++j){
                for(int k=0; k<3;++k){
                    R[j][k] = Ry[j][k];
                }
            }
            for(int i=0;i<Nodes.size();++i){
                double x = Nodes[i]->Position[0]*R[0][0] + Nodes[i]->Position[1]*R[0][1] + Nodes[i]->Position[2]*R[0][2];
                double y = Nodes[i]->Position[0]*R[1][0] + Nodes[i]->Position[1]*R[1][1] + Nodes[i]->Position[2]*R[1][2];
                double z = Nodes[i]->Position[0]*R[2][0] + Nodes[i]->Position[1]*R[2][1] + Nodes[i]->Position[2]*R[2][2];
                Nodes[i]->Position[0]=x;
                Nodes[i]->Position[1]=y;
                Nodes[i]->Position[2]=z;
            }
        }

        for(int i=0;i<Elements.size();++i){
            Elements[i]->updatePositions(3,Nodes);
        }
    }
}

void Simulation::updateGrowthRotationMatrices(){
    int nElement = Elements.size();
    for (int i=0; i<nElement; ++i){
        if (!Elements[i]->IsAblated){
            Elements[i]->CalculateGrowthRotationByF();
        }
    }
}

void Simulation::AddPeripodialMembraneNewScratch(){
    bool Success = true;
    Success = generateColumnarCircumferenceNodeList();
    if (!Success){
        cerr<<"Error!! circumferential nodes not sorted"<<endl;
    }
    calculateSystemCentre();
    sortColumnarCircumferenceNodeList();
    double avrSide=0.0, dummy =0.0;
    getAverageSideLength(dummy,avrSide);	//first term will get you the average side length of the peripodial membrane elements, second is the columnar elements

    cout<<"Sorted list for columnar basal circumferential nodes"<<endl;
    /*for (int i=0; i<ColumnarCircumferencialNodeList.size(); ++i){
        cout<<ColumnarCircumferencialNodeList[i]<<endl;
    }*/
    vector <int> NewNodeIds;
    vector <int> CorrespondingOldNodeIds;
    for (int i=0; i<ColumnarCircumferencialNodeList.size(); ++i){
        cout<<i<<endl;
        double* vec0;
        double* vec1;
        vec0 = new double[3];
        vec1 = new double[3];
        int indice0 = i-1;
        int indice1 = i+1;
        if (i == 0){
            indice0 = ColumnarCircumferencialNodeList.size() -1;
        }
        else if( i == ColumnarCircumferencialNodeList.size() -1){
            indice1 = 0;
        }
        int nodeId0 = ColumnarCircumferencialNodeList[indice0];
        int nodeId1 = ColumnarCircumferencialNodeList[indice1];
        int nodeIdi = ColumnarCircumferencialNodeList[i];
        for (int j=0; j<Nodes[nodeIdi]->nDim; j++){
            vec0[j] = Nodes[nodeId0]->Position[j] - Nodes[nodeIdi]->Position[j];
            vec1[j] = Nodes[nodeId1]->Position[j] - Nodes[nodeIdi]->Position[j];
        }
        Elements[0]->normaliseVector3D(vec0);
        Elements[0]->normaliseVector3D(vec1);
        for (int j=0; j<Nodes[nodeIdi]->nDim; j++){
            vec0[j] += vec1[j];
        }
        delete[] vec1;
        Elements[0]->normaliseVector3D(vec0);
        double* pos;
        pos = new double[3];
        for (int j=0; j<Nodes[nodeIdi]->nDim; j++){
            pos[j] = Nodes[nodeIdi]->Position[j]+avrSide*vec0[j];
        }
        Node* tmp_nd;
        //tmp_nd = new Node(Nodes.size(), 3, pos,Nodes[nodeIdi]->tissuePlacement, 1);
        NewNodeIds.push_back(Nodes.size());
        CorrespondingOldNodeIds.push_back(nodeIdi);

        int currNodeId = nodeIdi;
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
            currNodeId = (*itElement)->getCorrecpondingApical(currNodeId); //have the next node
            for (int j=0; j<Nodes[currNodeId]->nDim; j++){
                pos[j] = Nodes[currNodeId]->Position[j]+avrSide*vec0[j];
            }
            //tmp_nd = new Node(Nodes.size(), 3, pos,Nodes[nodeIdi]->tissuePlacement, 1);
            NewNodeIds.push_back(Nodes.size());
            CorrespondingOldNodeIds.push_back(currNodeId);
            if (Nodes[currNodeId]->tissuePlacement == 1) {
                //Added the apical node now, I can stop:
                finishedTissueThicness = true;
            }
        }
        //Nodes.push_back(tmp_nd);
        delete[] vec0;
        delete[] pos;
    }
    //lumenHeight = TissueHeight*lumenHeightScale;bibap
    //Now I have all the nodes added, I need to add the elements:
}

void Simulation::runOneStep(){
    //cout<<"entered run one step"<<endl;
    manualPerturbationToInitialSetup(false,false); //bool deform, bool rotate
    cleanGrowthData();
    resetForces();
    int freq = 10.0/dt ;
    if (freq <1 ){freq =1;}
    if ((timestep - 1)% freq  == 0){
        for (int i= 0; i<Nodes.size(); ++i){
            if (Nodes[i]->atCircumference){
                for (int j=0; j<3; j++){
                    Nodes[i]->FixedPos[0] = true;
                }
            }
        }
        cout<<"At time -- "<<dt*timestep<<" sec ("<<dt*timestep/3600<<" hours - "<<timestep<<" timesteps)"<<endl;
        calculateColumnarLayerBoundingBox();
        calculateDVDistance();
        int nElement = Elements.size();
        for (int i=0; i<nElement; ++i){
            Elements[i]->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],BoundingBoxSize[0],BoundingBoxSize[1]);
        }
    }
    if(nGrowthFunctions>0){
        if ((timestep - 1)% growthRotationUpdateFrequency  == 0){
            updateGrowthRotationMatrices();
        }
        //outputFile<<"calculating growth"<<endl;
        calculateGrowth();
    }
    int nElement = Elements.size();
    for (int i=0; i<nElement; ++i){
        if (!Elements[i]->IsAblated){
            Elements[i]->growShapeByFg(dt);
        }
    }
    updateNodeMasses();
    updateElementToConnectedNodes(Nodes);
    bool ExplicitCalculation = false;
    if (ExplicitCalculation){
        updateStep4RK();
    }
    else{
        updateStepNR();
    }
    processDisplayDataAndSave();
    //cout<<"finalised run one step"<<endl;
    calculateColumnarLayerBoundingBox();
    cout<<" step: "<<timestep<<" Pressure: "<<SuctionPressure[2]<<" Pa, maximum z-height: "<<boundingBox[1][2]<<" L/a: "<<(boundingBox[1][2]-50)/(2.0*pipetteRadius)<<endl;
    timestep++;
}

void Simulation::updateStep4RK(){
    for (int RKId = 0; RKId<4; ++RKId){
        //cout<<"entered RK: "<<RKId<<endl;
        int nElement = Elements.size();
        for (int i=0; i<nElement; ++i){
            if (!Elements[i]->IsAblated){
                //cout<<"calculating forces for element:  "<<i<<endl;
                Elements[i]->calculateForces(RKId, SystemForces, Nodes, outputFile);
                if(RKId == 0 ){
                    Elements[i]->createMatrixCopy(Elements[i]->RK1Strain,Elements[i]->Strain);
                }
            }
        }
        /*if (RKId ==0 ){
            cout<<"System Forces - Displacements"<<endl;
            for (int i=0; i<Nodes.size(); i++){
                cout<<SystemForces[0][i][0]<<" "<<SystemForces[0][i][0]/(Nodes[i]->Viscosity*Nodes[i]->mass)*dt<<endl;
                cout<<SystemForces[0][i][1]<<" "<<SystemForces[0][i][1]/(Nodes[i]->Viscosity*Nodes[i]->mass)*dt<<endl;
                cout<<SystemForces[0][i][2]<<" "<<SystemForces[0][i][2]/(Nodes[i]->Viscosity*Nodes[i]->mass)*dt<<endl;
            }
        }*/
        //cout<<"finalised elements & RK: "<<RKId<<endl;
        updateNodePositions(RKId);
        updateElementPositions(RKId);
    }
}

void Simulation::clearProjectedAreas(){
    int n = Nodes.size();
    for (int i=0;i<n; ++i){
        Nodes[i]->zProjectedArea = 0.0;
    }
}

void Simulation::correctzProjectedAreaForMidNodes(){
    int n = Nodes.size();
    for (int i=0;i<n; ++i){
        if (Nodes[i]->tissuePlacement == 2 ){ // the node is on midlayer
            //For the midline nodes, the area is added from apical and basal surfaces of elemetns on both sides.
            Nodes[i]->zProjectedArea /= 2.0;
        }
    }
}

void Simulation::calculateZProjectedAreas(){
    int n = Elements.size();
    clearProjectedAreas();
    for (int i=0; i<n; ++i){
        Elements[i]->calculateZProjectedAreas();
        Elements[i]->assignZProjectedAreas(Nodes);
    }
    correctzProjectedAreaForMidNodes();
}

void Simulation::updateStepNR(){
    int nNodes = Nodes.size();
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
    /*if (PipetteSuction){
        gsl_matrix_set_zero(gExt);
        packToPipetteWall(gExt);
        if (timestep>= PipetteInitialStep && timestep<PipetteEndStep){
            calculateZProjectedAreas();
            addPipetteForces(gExt);
        }
    }*/
    bool zeroViscosity = true;
    if (zeroViscosity){
        //fix x,y,z for 0
        Nodes[0]->FixedPos[0] = true;
        Nodes[0]->FixedPos[1] = true;
        Nodes[0]->FixedPos[2] = true;
        //fix x and z for 1
        Nodes[1]->FixedPos[0] = true;
        Nodes[1]->FixedPos[2] = true;
        //fix z for 2
        Nodes[2]->FixedPos[2] = true;

        /*Nodes[1]->FixedPos[1] = true;
        Nodes[2]->FixedPos[0] = true;
        Nodes[2]->FixedPos[1] = true;
        Nodes[3]->FixedPos[0] = true;
        Nodes[3]->FixedPos[1] = true;
        Nodes[3]->FixedPos[2] = true;
        Nodes[4]->FixedPos[0] = true;
        Nodes[4]->FixedPos[1] = true;
        Nodes[4]->FixedPos[2] = true;*/
    }
    vector <int> TransientZFixList;
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
        if (PipetteSuction){
            gsl_matrix_set_zero(gExt);
            //Z-fix trial"
            double pipLim2[2] = {0.9*pipetteRadius, 1.4*(pipetteRadius+pipetteThickness)};
            pipLim2[0] *= pipLim2[0];
            pipLim2[1] *= pipLim2[1];
            int nZfix = TransientZFixList.size();
            for (int zfixiter =0; zfixiter <nZfix; zfixiter++){
                Nodes[TransientZFixList[zfixiter]]->FixedPos[2] = 0;
            }
            TransientZFixList.clear();
            for (int currId=0; currId<Nodes.size();currId++){
                if (Nodes[currId]->FixedPos[2] == 0){
                    //the node is not already fixed in z
                    //I will check if it is in the range of the tip of pipette:
                    double zLimits[2] = {-2.0,2.0}; //the range in z to freeze nodes
                    double v0[3] = {Nodes[currId]->Position[0]-pipetteCentre[0], Nodes[currId]->Position[1]-pipetteCentre[1], Nodes[currId]->Position[2]-pipetteCentre[2]};
                    //if (currId == 61) {cout<<"Node[61] z: "<<Nodes[currId]->Position[2]<<"v0[2]: "<<v0[2]<<endl;}
                    if (v0[2] < zLimits[1] && v0[2]> zLimits[0]){
                        if (currId == 61) {cout<<"Node[61] is in z range"<<endl;}
                        //node is in Z range, is it in x-y range?
                        double xydist2 = v0[0]*v0[0]+v0[1]*v0[1];
                        //if (currId == 61) {cout<<"Node[61] xy2: "<<xydist2<<" lims: "<<pipLim2[0]<<" "<<pipLim2[1]<<endl;}
                        if (xydist2>pipLim2[0] && xydist2<pipLim2[1]){
                            //if (currId == 61) {cout<<"Node[61] is in xy range"<<endl;}
                             //node is within x-y range, freeze!
                            Nodes[currId]->FixedPos[2] = 1;
                            TransientZFixList.push_back(currId);
                        }
                    }
                }
            }
            // end of Zfix trial
            //packToPipetteWall(gExt,gSum);
            //Elements[0]->displayMatrix(gExt,"gExt-packed");
            if (timestep>= PipetteInitialStep && timestep<PipetteEndStep){
                calculateZProjectedAreas();
                addPipetteForces(gExt);
            }
            //Elements[0]->displayMatrix(gExt,"gExt-sucked");
        }
        for (int i=0; i<dim*nNodes; ++i){
            gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
        }

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
        //EliminateNode0Displacement(deltaU);
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
            for (int i=0; i<Nodes.size(); i++){
                cout<<SystemForces[0][i][0]<<" "<<gsl_vector_get(deltaU,3*i)<<endl;
                cout<<SystemForces[0][i][1]<<" "<<gsl_vector_get(deltaU,3*i+1)<<endl;
                cout<<SystemForces[0][i][2]<<" "<<gsl_vector_get(deltaU,3*i+2)<<endl;
            }
        }*/
    }
    //Now the calculation is converged, I update the node positions with the latest positions uk:
    updateNodePositionsNR(uk);
    //Element positions are already up to date.
}

void Simulation::calculateImplucitKElasticNumerical(gsl_matrix* K, gsl_matrix* geNoPerturbation){
    int n = Nodes.size();
    int dim = 3;
    double perturbation = 1E-6;
    gsl_matrix* KPerturbed = gsl_matrix_calloc(n*dim,n*dim);
    gsl_matrix* geWithPerturbation = gsl_matrix_calloc(n*dim,1);
    gsl_matrix* geDiff = gsl_matrix_calloc(n*dim,1);

    for (int i =0 ; i<n; i++){
        for (int j=0; j<dim; j++){
            gsl_matrix_set_zero(geWithPerturbation);
            gsl_matrix_set_zero(geDiff);
            resetForces();            
            Elements[0]->Positions[i][j] += perturbation;
            calculateElasticForcesForNR(geWithPerturbation);
            for(int k=0; k<n*dim; ++k){
                gsl_matrix_set(geDiff,k,0,gsl_matrix_get(geWithPerturbation,k,0));
            }
            gsl_matrix_sub(geDiff,geNoPerturbation);
            for(int k=0; k<n*dim; ++k){
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
        for (int i =0 ; i<n*dim; i++){
            for (int j =0 ; j<n*dim; j++){
                gsl_matrix_set(K,i, j, gsl_matrix_get(KPerturbed,i,j));
            }
        }
    }
}

//void Simulation::calcutateFixedK(vector<int> FixedNodes, gsl_matrix* K, gsl_vector* g){
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
    int n = Nodes.size();
    for(int i=0; i<n; i++){
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
    int n = Nodes.size();
    int dim = 3;
    for (int i = 0; i<n; ++i){
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
    int n = Nodes.size();
    int dim = 3;
    gsl_matrix* ViscousForces0 = gsl_matrix_calloc(n*dim,1);
    gsl_matrix* ViscousForces1 = gsl_matrix_calloc(n*dim,1);
    gsl_matrix* ukepsilon = gsl_matrix_calloc(n*dim,1);
    double epsilon = 0.5;
    for (int i=0; i<n*dim; i++){
        double value = gsl_matrix_get(uk,i,0)+epsilon;
        gsl_matrix_set(ukepsilon,i,0,value);
    }
    for (int i=0; i<n;++i){
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
    gsl_matrix* NumericKv = gsl_matrix_calloc(n*dim,n*dim);
    for (int i=0; i<n*dim; i++){
        double value = gsl_matrix_get(ViscousForces1,i,0) - gsl_matrix_get(ViscousForces0,i,0) ;
        value /= epsilon;
        gsl_matrix_set(NumericKv,i,i,value);
    }
    Elements[0]->displayMatrix(NumericKv,"NumericKv");
    Elements[0]->displayMatrix(mviscdt,"mviscdt");
    Elements[0]->displayMatrix(ViscousForces0,"ViscousForces0");
    Elements[0]->displayMatrix(ViscousForces1,"ViscousForces1");


}

void Simulation::updateElementPositionsinNR(gsl_matrix* uk){
    int dim = 3;
    int n = Elements.size();
    for (int i=0;i<n;i++){
        int* nodeIds = Elements[i]->getNodeIds();
        int nNodes= Elements[i]->getNodeNumber();
        for (int j=0; j<nNodes; ++j){
            double x = gsl_matrix_get(uk,dim*nodeIds[j],0);
            double y = gsl_matrix_get(uk,dim*nodeIds[j]+1,0);
            double z = gsl_matrix_get(uk,dim*nodeIds[j]+2,0);
            Elements[i]->Positions[j][0] = x;
            Elements[i]->Positions[j][1] = y;
            Elements[i]->Positions[j][2] = z;
        }
    }
}

void Simulation::EliminateNode0Displacement(gsl_vector* deltaU){
    double dU0[3] = {gsl_vector_get(deltaU,0),gsl_vector_get(deltaU,1),gsl_vector_get(deltaU,2)};
    int n = deltaU->size;
    for (int i=0; i<n; i=i+3){
        gsl_vector_set(deltaU,i,  gsl_vector_get(deltaU,i)-dU0[0]);
        gsl_vector_set(deltaU,i+1,gsl_vector_get(deltaU,i+1)-dU0[1]);
        gsl_vector_set(deltaU,i+2,gsl_vector_get(deltaU,i+2)-dU0[2]);
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
    int nNodes = Nodes.size();
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
    int    mtype = -2;        /* Real symmetric matrix */

    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    double   x[n_variables];
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
    int      i;

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
    phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

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
    int nElement = Elements.size();
    for (int i=0; i<nElement; ++i){
        Elements[i]->calculateImplicitKElastic();
    }
    //writing all elements K values into big K matrix:
    for (int i=0; i<nElement; ++i){
        Elements[i]->writeKelasticToMainKatrix(K);
    }
}

void Simulation::calculateElasticForcesForNR(gsl_matrix* ge){
    int nElement = Elements.size();
    for (int i=0; i<nElement; ++i){
        if (!Elements[i]->IsAblated){
            //cout<<"calculating forces for element:  "<<i<<endl;
            Elements[i]->calculateForces(0, SystemForces, Nodes, outputFile);
            /*int nNodes = Nodes.size();
            cout<<"Element "<<i<<"system forces: "<<endl;
            for (int nn=0; nn< nNodes; ++nn){
                 cout<<"Node: "<<nn<<" "<<SystemForces[0][nn][0]<<" "<<SystemForces[0][nn][1]<<" "<<SystemForces[0][nn][2]<<endl;
            }*/
        }
    }
    //now all the forces are written on SysyemForces RK step = 0
    //I will add them into ge, this step can be made faster by separating calculate forces function into two,
    //and filling up either ge or System forces depending on hte solution method:
    int nNodes = Nodes.size();
    for (int i = 0; i< nNodes; ++i){
        for ( int j=0; j<3; j++){
            gsl_matrix_set(ge, 3*i+j,0,SystemForces[0][i][j]);
        }
    }
    //cout<<"finalised elastic force calculation of the system"<<endl;
}

void Simulation::calculateViscousForcesForNR(gsl_matrix* gv, gsl_matrix* mviscdt, gsl_matrix* uk, gsl_matrix* un){
    int nNodes = Nodes.size();
    int dim  = 3;
    gsl_matrix* displacement = gsl_matrix_calloc(dim*nNodes,1);
    Elements[0]->createMatrixCopy(displacement,un);
    gsl_matrix_sub(displacement,uk);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mviscdt, displacement,0.0, gv);

/*
    gsl_matrix* dispPerturb = gsl_matrix_calloc(dim*nNodes,1);
    gsl_matrix* ukPerturb = gsl_matrix_calloc(dim*nNodes,1);
    gsl_matrix* gvPerturb = gsl_matrix_calloc(dim*nNodes,1);
    Elements[0]->createMatrixCopy(dispPerturb,un);
    Elements[0]->createMatrixCopy(ukPerturb,uk);
    double perturb = 5.0;
    double value = gsl_matrix_get(ukPerturb,9,0) + perturb;
    gsl_matrix_set(ukPerturb,9,0,value);
    gsl_matrix_sub(dispPerturb,ukPerturb);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mviscdt, dispPerturb,0.0, gvPerturb);
    Elements[0]->displayMatrix(displacement,"displacement");
    Elements[0]->displayMatrix(gv,"gv");
    Elements[0]->displayMatrix(gvPerturb,"gvPerturb");
    gsl_matrix_sub(gvPerturb,gv);
    Elements[0]->displayMatrix(gvPerturb,"PerturbDiff");
    cout<<gsl_matrix_get(gvPerturb,9,0)/perturb<<" mass value:  "<<gsl_matrix_get(mviscdt,9,9)<<endl;

*/
}

void Simulation::constructUnMatrix(gsl_matrix* un){
    int nNodes = Nodes.size();
    for (int i = 0; i<nNodes; ++i ){
        for (int j=0; j<3; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void Simulation::constructLumpedMassViscosityDtMatrix(gsl_matrix* mviscdt){
    int nNodes = Nodes.size();
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
    if (displayIsOn && !DisplaySave){
        //The simulation is not displaying a saved setup, it is running and displaying
        //I need to correct the values to be displayed, and store averages, otherwise
        //the displayed values will be from artificial setups of different RK steps. (RK1 or RK4 depending on parameter)
        updateDisplaySaveValuesFromRK();
    }
    if (saveData && timestep % dataSaveInterval == 0){
        updateDisplaySaveValuesFromRK();
        saveStep();
    }
}

void Simulation::updateNodeMasses(){
    int nNodes = Nodes.size();
    for (int i=0; i<nNodes; ++i){
        Nodes[i]->mass = 0;
    }
    int nElements = Elements.size();
    for (int i=0; i<nElements; ++i){
        if (!Elements[i]->IsAblated){
            Elements[i]->assignVolumesToNodes(Nodes);
        }
    }
}

void 	Simulation:: updateElementToConnectedNodes(vector <Node*>& Nodes){
    for (int j=0; j<Nodes.size(); ++j){
        int n = Nodes[j]->connectedElementIds.size();
        for (int i=0; i<n; ++i){
            //if(!Elements[Nodes[j]->connectedElementIds[i]] -> IsAblated){
                Nodes[j]->connectedElementWeights[i] = Elements[Nodes[j]->connectedElementIds[i]]->VolumePerNode/Nodes[j]->mass;
            //}
        }
    }
}

void Simulation::smallStrainrunOneStep(){
	//cout<<"time step: "<<timestep<<endl;
	if(timestep==0){
		calculateColumnarLayerBoundingBox();
		calculateDVDistance();
		//outputFile<<"calculating element health"<<endl;
		int nElement = Elements.size();
		for (int i=0; i<nElement; ++i){
			Elements[i]->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],BoundingBoxSize[0],BoundingBoxSize[1]);
		}
		//for (int i=0; i<nElement; ++i){
		//	Elements[i]->calculatGrowthScalingMatrices();
		//}
		//cout<<"In loop for t=0: "<<endl;
		//Elements[0]->displayRelativePosInBoundingBox();
		/*
		for(int i=0;i<Nodes.size();++i){
			//if (Nodes[i]->atCircumference){
			//	Nodes[i]->FixedPos[0] = true;
			//	Nodes[i]->FixedPos[1] = true;
			//	Nodes[i]->FixedPos[2] = true;
			//}
			Nodes[i]->Position[0] *=2.0;
			Nodes[i]->Position[1] *=1.0;
			Nodes[i]->Position[2] *=1.0;
		}
		double R[3][3];
		double Rx[3][3] = {{1,0,0},{0,0,-1},{0,1,0}};
		double Ry[3][3] = {{0,0,1},{0,1,0},{-1,0,0}};
		double Rz[3][3] = {{0,-1,0},{1,0,0},{0,0,1}};
		for (int j =0; j<3;++j){
			for(int k=0; k<3;++k){
				R[j][k] = Ry[j][k];
			}
		}
		for(int i=0;i<Nodes.size();++i){
			double x = Nodes[i]->Position[0]*R[0][0] + Nodes[i]->Position[1]*R[0][1] + Nodes[i]->Position[2]*R[0][2];
			double y = Nodes[i]->Position[0]*R[1][0] + Nodes[i]->Position[1]*R[1][1] + Nodes[i]->Position[2]*R[1][2];
			double z = Nodes[i]->Position[0]*R[2][0] + Nodes[i]->Position[1]*R[2][1] + Nodes[i]->Position[2]*R[2][2];
			Nodes[i]->Position[0]=x;
			Nodes[i]->Position[1]=y;
			Nodes[i]->Position[2]=z;
		}
		for(int i=0;i<Elements.size();++i){
			Elements[i]->updatePositions(3,Nodes);
		}*/
		//for (int i=0; i<Nodes.size();++i){
		//	cout<<"Node: "<<i<<" pos: "<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		//}

	}
	//cout<<"outisde loop for t=0: "<<endl;
	//Elements[0]->displayRelativePosInBoundingBox();
	if(timestep==-10/dt){
		LaserAblate(0.0,0.0,1.8);
	}
	int displayfreq = 60/dt;
	if (timestep%displayfreq == 0){
		//Success = reOpenOutputFile();
		//if(!Success){cerr<<" outputfile not open, time : "<<timestep * dt<<endl;}
		outputFile<<"time : "<<timestep * dt<<endl;
	}
	int perturbstep = -1200/dt;
	//cout<<"perturbstep: "<<perturbstep<<endl;
	if (timestep == perturbstep ){
		cerr<<"Perturbing the system, pushing middle nodes down by 0.5 units"<<endl;
		for(int i=0;i<Nodes.size();++i){
			if (Nodes[i]->Position[0]<20 && Nodes[i]->Position[0]>-20 && Nodes[i]->Position[1]>-10 && Nodes[i]->Position[1]<-5 ){
				Nodes[i]->Position[2] += -1.0;
			}
		}
	}
	//if(timestep==0){Elements[0]->LocalGrowthStrainsMat(0,0) = 1.0;}
	//cleanreferenceupdates();
    //cleanMatrixUpdateData();
	cleanGrowthData();
	resetForces();
	//alignTissueDVToXPositive();
	int freq = 10.0/dt ;
	if ((timestep - 1)% freq  == 0){
		calculateColumnarLayerBoundingBox();
		calculateDVDistance();
		int nElement = Elements.size();
		for (int i=0; i<nElement; ++i){
			Elements[i]->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],BoundingBoxSize[0],BoundingBoxSize[1]);
		}
	}
	int nElement = Elements.size();
	if(nGrowthFunctions>0){
		//outputFile<<"calculating growth"<<endl;
		calculateGrowth();
	}
	//cout<<"after calculate growth: "<<endl;
	//Elements[0]->displayRelativePosInBoundingBox();
	//if(nShapeChangeFunctions>0){
	//	changeCellShapesInSystem();
	//}
	//outputFile<<"calculating alignment of reference"<<endl;
	for (int i=0; i<nElement; ++i){
		//deletion here:
		//if (Elements[i]->IsGrowing && Elements[i]->tissueType == 0){ //only columnar layer is grown this way, peripodial membrane is already grown without alignment
		if (Elements[i]->IsGrowing){
            //Elements[i]->growShape();
		}
	}

	double periPackingThreshold = 1000, colPackingThreshold = 1000;
	getAverageSideLength(periPackingThreshold,colPackingThreshold);
	//packingThreshold *=0.7; //I am adding a 40% safety to average side length, and then taking half of it as threshold for packing
	for (int RKId = 0; RKId<4; ++RKId){
		//outputFile<<"started RK: "<<RKId<<endl;
		calculatePacking(RKId,periPackingThreshold,colPackingThreshold);
		for (int i=0; i<nElement; ++i){
			if (!Elements[i]->IsAblated){
				Elements[i]->calculateForces(RKId, SystemForces, Nodes, outputFile);
			}
		}
		//outputFile<<"     calculated forces"<<endl;
		if (stretcherAttached && timestep>StretchInitialStep && timestep<StretchEndStep){
			addStretchForces(RKId);
		}

        //if (PipetteSuction && timestep> PipetteInitialStep && timestep<PipetteEndStep){
        //	addPipetteForces(RKId);
        //}
		redistributePeripodialMembraneForces(RKId);
		updateNodePositions(RKId);
		//outputFile<<"     updated node pos"<<endl;
		updateElementPositions(RKId);
		//for (int i=0; i<Nodes.size(); ++i){
		//	cout<<"RK: "<<RKId<<" Node: "<<i<<"	PackingForces: "<<PackingForces[RKId][i][0]<<"	"<<PackingForces[RKId][i][1]<<"	"<<PackingForces[RKId][i][2]<<endl;
		//}
		//outputFile<<"     updated element pos"<<endl;
	}
	if (displayIsOn && !DisplaySave){
		//The simulation is not displaying a saved setup, it is running and displaying
		//I need to correct the values to be displayed, and store averages, otherwise
		//the displayed values will be from artificial setups of different RK steps. (RK1 or RK4 depending on parameter)
		updateDisplaySaveValuesFromRK();
	}
	if (saveData && timestep % dataSaveInterval == 0){
		updateDisplaySaveValuesFromRK();
		saveStep();
	}
	timestep++;
	//outputFile<<"finished runonestep"<<endl;
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

void Simulation::calculatePacking(int RKId, double PeriThreshold, double ColThreshold){
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
			(*itNode)->getCurrentRKPosition(RKId,pos);
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
									SystemForces[RKId][(*itNode)->Id][j] += F[j];
									PackingForces[RKId][(*itNode)->Id][j] += F[j];
								}
							}
							//cout<<"updated the system forces: "<<SystemForces[RKId][(*itNode)->Id][0]<<" "<<SystemForces[RKId][(*itNode)->Id][1]<<" "<<SystemForces[RKId][(*itNode)->Id][2]<<endl;
							//add opposite force to the element nodes:
							(*itEle)->AddPackingToSurface((*itNode)->tissuePlacement, F[0],F[1],F[2],RKId, SystemForces,PackingForces, Nodes);
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
								SystemForces[RKId][(*itNode)->Id][j] += pipF[j];
								PackingForces[RKId][(*itNode)->Id][j] += pipF[j];
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

/*	float multiplier = 100.0;
	vector<Node*>::iterator itNode0;
	vector<Node*>::iterator itNode1;
	float dmin = threshold;
	float dminNeg =(-1.0)*threshold;
	float t2 = threshold*threshold;
	for (itNode0=Nodes.begin(); itNode0<Nodes.end(); ++itNode0){
		if (!(*itNode0)->atPeripodialMembraneCircumference){
			for (itNode1=itNode0+1; itNode1<Nodes.end(); ++itNode1){
				if (!(*itNode1)->atPeripodialMembraneCircumference){
					float dx =100.0, dy = 100.0, dz = 100.0;
					if (RKId ==0){
						dx = (*itNode1)->Position[0]-(*itNode0)->Position[0];
						dy = (*itNode1)->Position[1]-(*itNode0)->Position[1];
						dz = (*itNode1)->Position[2]-(*itNode0)->Position[2];
					}
					else{
						dx = (*itNode1)->RKPosition[0]-(*itNode0)->RKPosition[0];
						dy = (*itNode1)->RKPosition[1]-(*itNode0)->RKPosition[1];
						dz = (*itNode1)->RKPosition[2]-(*itNode0)->RKPosition[2];
					}
					if ((dx >0 && dx < dmin) || (dx <0 && dx >dminNeg)){
						if ((dy >0 && dy < dmin) || (dy <0 && dy >dminNeg)){
							if ((dz >0 && dz < dmin) || (dz <0 && dz >dminNeg)){
								//distance is close enough for further investigation, but are these nodes neighbours?
								bool AreNeigs = (*itNode0)->checkIfNeighbour((*itNode1)->Id);
								if (!AreNeigs){
									float d2 = dx*dx+dy*dy+dz*dz;
									if (d2 < t2){
										float Fmag = multiplier * (1.0/d2 - 1.0/t2);
										float d = pow(d2,0.5);
										//normalising
										dx /= d;
										dy /= d;
										dz /= d;
										//vector id from node 0 to node 1, force should be in the opposite direction
										SystemForces[RKId][(*itNode0)->Id][0] -= dx*Fmag;
										SystemForces[RKId][(*itNode0)->Id][1] -= dy*Fmag;
										SystemForces[RKId][(*itNode0)->Id][2] -= dz*Fmag;
										//add opposite force to the second node:
										SystemForces[RKId][(*itNode1)->Id][0] += dx*Fmag;
										SystemForces[RKId][(*itNode1)->Id][1] += dy*Fmag;
										SystemForces[RKId][(*itNode1)->Id][2] += dz*Fmag;
										if ((*itNode0)->Id == 497 || (*itNode1)->Id == 497){
											cout<<"threshold: "<<threshold<<" distance: "<<d<<" Nodes: "<<(*itNode0)->Id<<" "<<(*itNode1)->Id<<" force: "<<Fmag<<endl;
											cout<<"	pos:"<<(*itNode0)->Position[0]<<" "<<(*itNode0)->Position[1]<<" "<<(*itNode0)->Position[2]<<endl;
											cout<<"	pos:"<<(*itNode1)->Position[0]<<" "<<(*itNode1)->Position[1]<<" "<<(*itNode1)->Position[2]<<endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
*/


void Simulation::redistributePeripodialMembraneForces(int RKId){
	int n = PeripodialMembraneCircumferencialNodeList.size();
	for (int i=0; i<n; ++i){
		int index = PeripodialMembraneCircumferencialNodeList[i];
		double F[3] = {SystemForces[RKId][index][0], SystemForces[RKId][index][1], SystemForces[RKId][index][2]};
		int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
			double w = Nodes[index]->AssociatedNodeWeightsDueToPeripodialMembrane[j];
			if(!Nodes[index2]->FixedPos[0]){
				SystemForces[RKId][index2][0] += F[0]*w;
			}
			if(!Nodes[index2]->FixedPos[1]){
				SystemForces[RKId][index2][1] += F[1]*w;
			}
			if(!Nodes[index2]->FixedPos[2]){
				SystemForces[RKId][index2][2] += F[2]*w;
			}
		}
		SystemForces[RKId][index][0] = 0.0;
		SystemForces[RKId][index][1] = 0.0;
		SystemForces[RKId][index][2] = 0.0;
	}
}



void Simulation::updateNodePositions(int RKId){
	//Update Node positions:
	int n = Nodes.size();
	if (RKId < 3){
		//the first 3 RK steps, the velocity will be calculated, and the positions will be updated from normal positions to RKPositions, with half dt:
		double multiplier=0.0;
		if (RKId<2){
			multiplier =0.5;
		}
		else{
			multiplier =1.0;
		}
        double v0[3] = {0.0,0.0,0.0};
        for (int j=0; j<Nodes[0]->nDim; ++j){
            Nodes[0]->Velocity[RKId][j] = SystemForces[RKId][0][j]/ (Nodes[0]->Viscosity*Nodes[0]->mass) ;
            v0[j] = Nodes[0]->Velocity[RKId][j];
            Nodes[0]->Velocity[RKId][j] = 0.0;
            Nodes[0]->RKPosition[j] = Nodes[0]->Position[j];
        }
        for (int i=1;i<n;++i){
			for (int j=0; j<Nodes[i]->nDim; ++j){
					Nodes[i]->Velocity[RKId][j] = SystemForces[RKId][i][j]/ (Nodes[i]->Viscosity*Nodes[i]->mass) ;
                    Nodes[i]->Velocity[RKId][j] -= v0[j];
                    Nodes[i]->RKPosition[j] = Nodes[i]->Position[j] + Nodes[i]->Velocity[RKId][j]*multiplier*dt;
			}
            //cout<<"RK: "<<RKId<<" node: "<<i<<"mass: "<<Nodes[i]->mass<<" visc: "<<Nodes[i]->Viscosity<<endl;
            //for (int j=0; j<Nodes[i]->nDim; ++j){
            //    cout<<"	"<<Nodes[i]->Velocity[RKId][j]<<" ";
            //}
            //cout<<endl;
			//cout<<"	Old pos: ";
			//for (int j=0; j<Nodes[i]->nDim; ++j){
			//	cout<<Nodes[i]->Position[j]<<" ";
			//}
			//cout<<endl<<"	New pos: ";
			//for (int j=0; j<Nodes[i]->nDim; ++j){
			//	cout<<Nodes[i]->RKPosition[j]<<" ";
			//}
			//cout<<endl;
		}
	}
	else{
		//this is the last RK step, I need to update the velocity only with RK, then I need to calculate the final positions
		//from 4 RK velocities:
        for (int i=1;i<n;++i){
            //cout<<"Nodes "<<i<<" velocity: ";
            double v0[3] = {0.0,0.0,0.0};
            for (int j=0; j<Nodes[0]->nDim; ++j){
                Nodes[0]->Velocity[RKId][j] = SystemForces[RKId][0][j]/ (Nodes[0]->Viscosity*Nodes[0]->mass) ;
                v0[j] = Nodes[0]->Velocity[RKId][j];
                Nodes[0]->Velocity[RKId][j] = 0.0;
            }
            for (int j=0; j<Nodes[i]->nDim; ++j){
				Nodes[i]->Velocity[RKId][j] = SystemForces[RKId][i][j]/(Nodes[i]->Viscosity*Nodes[i]->mass);
                Nodes[i]->Velocity[RKId][j] -= v0[j];
                //now I have 4 velocity data (corresponding to Runge-Kutta  k1, k2, k3, and k4)
				//writing  the velocity into v[0]
				//cout<<Nodes[i]->Velocity[0][j]<<" "<<Nodes[i]->Velocity[1][j]<<" "<<Nodes[i]->Velocity[2][j]<<" "<<Nodes[i]->Velocity[3][j]<<" ";
				Nodes[i]->Velocity[0][j] = 1.0/6.0 * (Nodes[i]->Velocity[0][j] + 2.0 * (Nodes[i]->Velocity[1][j] + Nodes[i]->Velocity[2][j]) + Nodes[i]->Velocity[3][j]);
				Nodes[i]->Position[j] += Nodes[i]->Velocity[0][j]*dt;
			}
		}
	}
	updateNodePositionsForPeripodialMembraneCircumference(RKId);
}

void Simulation::realignPositionsForMidAttachedPeripodialMembrane(int RKId){
	int n = PeripodialMembraneCircumferencialNodeList.size();
	if (RKId < 3){
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			double RKsumV[3] = {0.0,0.0,0.0};
			double RKsumPos[3] = {0.0,0.0,0.0};
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				for (int k=0; k<Nodes[index2]->nDim; ++k){
					RKsumV[k] += Nodes[index2]->Velocity[RKId][k];
					RKsumPos[k] += Nodes[index2]->RKPosition[k];
				}
			}
			for (int k=0; k<Nodes[index]->nDim; ++k){
				Nodes[index]->Velocity[RKId][k] = RKsumV[k]/nA;
				Nodes[index]->RKPosition[k] = RKsumPos[k]/nA;
			}
		}
	}
	else{
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			double sumV[3] = {0.0,0.0,0.0};
			double RKsumV[3] = {0.0,0.0,0.0};
			double sumPos[3] = {0.0,0.0,0.0};
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				for (int k=0; k<Nodes[index2]->nDim; ++k){
					RKsumV[k] += Nodes[index2]->Velocity[RKId][k];
					sumV[k] += Nodes[index2]->Velocity[0][k];
					sumPos[k] += Nodes[index2]->Position[k];
				}
			}
			for (int k=0; k<Nodes[index]->nDim; ++k){
				Nodes[index]->Velocity[RKId][k] = RKsumV[k]/nA;
				Nodes[index]->Velocity[0][k] = sumV[k]/nA;
				Nodes[index]->Position[k] = sumPos[k]/nA;
			}
		}
	}
}

void Simulation::realignPositionsForApicalAttachedPeripodialMembrane(int RKId){
	int n = PeripodialMembraneCircumferencialNodeList.size();
	if (RKId < 3){
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				if(Nodes[index2] -> tissuePlacement ==1 ){
					//found the apical node!
					for (int k=0; k<Nodes[index2]->nDim; ++k){
						Nodes[index]->Velocity[RKId][k] = Nodes[index2]->Velocity[RKId][k];
						Nodes[index]->RKPosition[k] =  Nodes[index2]->RKPosition[k];
					}
					break;
				}
			}
		}
	}
	else{
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				if(Nodes[index2] -> tissuePlacement ==1 ){
					//found Apical Node!
					for (int k=0; k<Nodes[index2]->nDim; ++k){
						Nodes[index]->Velocity[RKId][k] = Nodes[index2]->Velocity[RKId][k];
						Nodes[index]->Velocity[0][k] = Nodes[index2]->Velocity[0][k];
						Nodes[index]->Position[k] =  Nodes[index2]->Position[k];
					}
					break;
				}
			}
		}
	}
}


void Simulation::realignPositionsForHooveringPeripodialMembrane(int RKId){
	int n = PeripodialMembraneCircumferencialNodeList.size();
	if (RKId < 3){
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				if(Nodes[index2] -> tissuePlacement ==1 ){
					//found the apical node!
					for (int k=0; k<Nodes[index2]->nDim; ++k){
						Nodes[index]->Velocity[RKId][k] = Nodes[index2]->Velocity[RKId][k];
						Nodes[index]->RKPosition[k] =  Nodes[index2]->RKPosition[k];
					}
					Nodes[index]->RKPosition[2] += lumenHeight;
					break;
				}
			}
		}
	}
	else{
		for (int i=0; i<n; ++i){
			int index = PeripodialMembraneCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodialMembrane.size();
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodialMembrane[j];
				if(Nodes[index2] -> tissuePlacement ==1 ){
					//found Apical Node!
					for (int k=0; k<Nodes[index2]->nDim; ++k){
						Nodes[index]->Velocity[RKId][k] = Nodes[index2]->Velocity[RKId][k];
						Nodes[index]->Velocity[0][k] = Nodes[index2]->Velocity[0][k];
						Nodes[index]->Position[k] =  Nodes[index2]->Position[k];
					}
					Nodes[index]->Position[2] += lumenHeight;
					break;
				}
			}
		}
	}
}


void Simulation::updateNodePositionsForPeripodialMembraneCircumference(int RKId){
	if (PeripodialMembraneType==0){
		realignPositionsForHooveringPeripodialMembrane(RKId);
	}
	if(PeripodialMembraneType == 2) {
		realignPositionsForApicalAttachedPeripodialMembrane(RKId);
	}
	if (PeripodialMembraneType == 1) {
		realignPositionsForMidAttachedPeripodialMembrane(RKId);
	}
}


void Simulation::updateElementPositions(int RKId){
	for (int i=0;i<Elements.size(); ++i ){
		Elements[i]->updatePositions(RKId, Nodes);
	}
}

void Simulation::updateElementPositionsSingle(int RKId, int i ){
	Elements[i]->updatePositions(RKId, Nodes);
}

void Simulation::alignTissueDVToXPositive(){
	double* u = new double[3];
	double* v = new double[3];
	for (int i=0;i<3;++i){
		u[i] = Nodes[DVRight]->Position[i] - Nodes[DVLeft]->Position[i];
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
	int n = Nodes.size();
	for(int i=0;i<n;++i){
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			u[j] = Nodes[i]->Position[j];
		}
		Elements[0]->rotateVectorByRotationMatrix(u,rotMat);
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			Nodes[i]->Position[j] = u[j];
		}
	}
	for(int i=0;i<Elements.size();++i){
		//this will need to update the positions of the elements using the actual positions, which corresponds to RK ste 4, RKId = 3;
		Elements[i]->updatePositions(3, Nodes);
	}
	delete[] rotAx;
	delete[] rotMat;
	delete[] u;
	delete[] v;
}

void Simulation::calculateDVDistance(){
	double d[3];
	for (int i=0;i<3;++i){
		d[i] = Nodes[DVRight]->Position[i] - Nodes[DVLeft]->Position[i];
		d[i] *=d[i];

	}
	double dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	//outputFile<<"time: "<<timestep * dt<<" DV distance is: "<<dmag<<endl;
}

void Simulation::cleanGrowthData(){
	int nElement = Elements.size();
	for (int i=0; i<nElement; ++i){
		Elements[i]->resetCurrStepGrowthData();
		//Elements[i]->resetCurrStepShapeChangeData();
	}
}

void Simulation::resetForces(){
	int n = Nodes.size();
	int dim = 3;
	//4 RK steps
	for (int i=0;i<4;++i){
		//n nodes
		for (int j=0;j<n;++j){
			//3 dimensions
			for (int k=0;k<dim;++k){
				SystemForces[i][j][k]=0.0;
                //PackingForces[i][j][k]=0.0;
			}
		}
	}
}

void Simulation::calculateColumnarLayerBoundingBox(){
	boundingBox[0][0] =  100000.0;	//lower left x
	boundingBox[0][1] =  100000.0;	//lower left y
	boundingBox[0][2] =  100000.0;	//lower z
	boundingBox[1][0] = -100000.0;	//upper right x
	boundingBox[1][1] = -100000.0;	//upper right y
	boundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
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
	}
	for (int i=0; i<3; ++i){
		BoundingBoxSize[i] = boundingBox[1][i] - boundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}

void Simulation::updateDisplaySaveValuesFromRK(){
	//in display and save, I am using the storage for RK step 1. I need to make this the average values of 4 steps, otherwise
	//it will show the misleading RK 1 values only.
	int n = Nodes.size();
	for (int j=0;j<n;++j){
		SystemForces[0][j][0] = 1.0/6.0 * (SystemForces[0][j][0] + 2 * (SystemForces[1][j][0] + SystemForces[2][j][0]) + SystemForces[3][j][0]);
		SystemForces[0][j][1] = 1.0/6.0 * (SystemForces[0][j][1] + 2 * (SystemForces[1][j][1] + SystemForces[2][j][1]) + SystemForces[3][j][1]);
		SystemForces[0][j][2] = 1.0/6.0 * (SystemForces[0][j][2] + 2 * (SystemForces[1][j][2] + SystemForces[2][j][2]) + SystemForces[3][j][2]);

        //PackingForces[0][j][0] = 1.0/6.0 * (PackingForces[0][j][0] + 2 * (PackingForces[1][j][0] + PackingForces[2][j][0]) + PackingForces[3][j][0]);
        //PackingForces[0][j][1] = 1.0/6.0 * (PackingForces[0][j][1] + 2 * (PackingForces[1][j][1] + PackingForces[2][j][1]) + PackingForces[3][j][1]);
        //PackingForces[0][j][2] = 1.0/6.0 * (PackingForces[0][j][2] + 2 * (PackingForces[1][j][2] + PackingForces[2][j][2]) + PackingForces[3][j][2]);
		//I do not need to do velocities, as I already calculate the average velocity on storage space of RK step 1, for posiiton update
	}
	n = Elements.size();
    for (int j=0;j<n;++j){
        Elements[j]->createMatrixCopy(Elements[j]->Strain, Elements[j]->RK1Strain);
    }
}

void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<timestep*dt<<" sec"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
    writeGrowth();
	writeForces();
	writeVelocities();
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
	int n = Nodes.size();
	saveFileMesh<<n<<endl;
	for (int i = 0; i<n; ++i){
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
	int n= Elements.size();
	saveFileMesh<<n<<endl;
	for (int i = 0; i<n; ++i){
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
	int n = Elements.size();
	for (int i=0;i<n;++i){
        for (int j=0; j<6; ++j){
            double S = gsl_matrix_get(Elements[i]->Strain,j,0);
            saveFileTensionCompression.write((char*) &S, sizeof S);
        }
    }
	saveFileTensionCompression.flush();
}

void Simulation::writeGrowth(){
    //for (int i=0;i<6;++i){
    //	cout<<" at timestep :"<< timestep<<" the plastic strains of element 0:	"<<Elements[0]->PlasticStrain(i)<<"	normal strain: 	"<<Elements[i]->Strain(0)<<endl;
    //}
    int n = Elements.size();
    for (int i=0;i<n;++i){
        gsl_matrix* currFg = Elements[i]->getFg();
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
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		saveFileForces.write((char*) &SystemForces[0][i][0], sizeof SystemForces[0][i][0]);
		saveFileForces.write((char*) &SystemForces[0][i][1], sizeof SystemForces[0][i][1]);
		saveFileForces.write((char*) &SystemForces[0][i][2], sizeof SystemForces[0][i][2]);
	}
	saveFileForces.flush();
}

void Simulation::writeVelocities(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			saveFileVelocities.write((char*) &Nodes[i]->Velocity[j], sizeof Nodes[i]->Velocity[j]);
		}
	}
	saveFileVelocities.flush();
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
	int  n = Elements.size();
	for ( int i = 0; i < n; ++i ){
		Elements[i]->setGrowthRate(0.0,0.0,0.0);
	}
}

void Simulation::cleanUpShapeChangeRates(){
	int  n = Elements.size();
	for ( int i = 0; i < n; ++i ){
		Elements[i]->setShapeChangeRate(0.0,0.0,0.0);
	}
}

void Simulation::calculateGrowthUniform(GrowthFunctionBase* currGF){
	float simTime = dt*timestep;
	//cout<<"inside uniform growth function, initTime: "<<currGF->initTime <<" endtime: "<<currGF->endTime<<" simTime"<<simTime<<endl;
	if(simTime > currGF->initTime && simTime < currGF->endTime ){
		//cout<<"calculating growth"<<endl;
		double *maxValues;
        maxValues = new double[3];
		currGF->getGrowthRate(maxValues);
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			//tissue type == 0 is columnar layer, ==1 is peripodial membrane
			if ((currGF->applyToColumnarLayer && Elements[i]->tissueType == 0) || (currGF->applyToPeripodialMembrane && Elements[i]->tissueType == 1)){
				 //cout<<"updating growth for element: "<<Elements[i]->Id<<endl;
                Elements[i]->updateGrowthRate(maxValues[0],maxValues[1],maxValues[2]);
           }
		}
		delete[] maxValues;
	}
}

void Simulation::calculateGrowthRing(GrowthFunctionBase* currGF){
	float simTime = dt*timestep;
	if(simTime > currGF->initTime && simTime < currGF->endTime ){
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
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			if ((currGF->applyToColumnarLayer && Elements[i]->tissueType == 0) || (currGF->applyToPeripodialMembrane && Elements[i]->tissueType == 1)){
				double* Elementcentre = new double[3];
				Elementcentre = Elements[i]->getCentre();
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
                    Elements[i]->updateGrowthRate(growthscale[0],growthscale[1],growthscale[2]);
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
	if(simTime > currGF->initTime && simTime < currGF->endTime ){
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			if ((currGF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currGF->applyToPeripodialMembrane && (*itElement)->tissueType == 1)){
				double* ReletivePos = new double[2];
				//normalising the element centre position with bounding box
				(*itElement)->getRelativePosInBoundingBox(ReletivePos);
				ReletivePos[0] *= (float) (nGridX-1);
				ReletivePos[1] *= (float) (nGridY-1);
				int indexX = floor(ReletivePos[0]);
				double fracX  = ReletivePos[0] - indexX;
				if (indexX == nGridX) { //this is for the point that is exactly the point determining the bounding box high end in X
					indexX--;
					fracX = 1.0;
				}
				int indexY = floor(ReletivePos[1]);
				double fracY  = ReletivePos[1] - indexY;
				if (indexY == nGridY) { //this is for the point that is exactly the point determining the bounding box high end in Y
					indexY--;
					fracY = 1.0;
				}
                double growthYmid[2][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0}};
                double growthscale[3];
				for (int axis = 0; axis<3; ++axis){
                    growthYmid[0][axis] = currGF->getGrowthMatrixElement(indexX,indexY,axis)*(1.0-fracX) + currGF->getGrowthMatrixElement(indexX+1,indexY,axis)*fracX;
                    growthYmid[1][axis] = currGF->getGrowthMatrixElement(indexX,indexY+1,axis)*(1.0-fracX) + currGF->getGrowthMatrixElement(indexX+1,indexY+1,axis)*fracX;
					growthscale[axis] = growthYmid[0][axis]*(1.0-fracY) + growthYmid[1][axis]*fracY;
				}
				//growing the shape
                (*itElement)->updateGrowthRate(growthscale[0],growthscale[1],growthscale[2]);
				delete[] ReletivePos;
			}
		}
	}
}

/*
void Simulation::changeCellShapeRing(int currIndex){
	//cout<<"entered shape change with ring function"<<endl;
	float initTime = ShapeChangeParameters[currIndex];
	float endTime = ShapeChangeParameters[currIndex+1];
	float simTime = dt*timestep;
	if(simTime > initTime && simTime < endTime ){
		float centre[2] = {ShapeChangeParameters[currIndex+2], ShapeChangeParameters[currIndex+3]};
		float innerRadius = {ShapeChangeParameters[currIndex+4]};
		float outerRadius = {ShapeChangeParameters[currIndex+5]};
		int axis = {(int) ShapeChangeParameters[currIndex+6]};
		float maxValue = {ShapeChangeParameters[currIndex+7]};
		float innerRadius2 = innerRadius*innerRadius;
		float outerRadius2 = outerRadius*outerRadius;
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			double* Elementcentre = new double[3];
			Elementcentre = Elements[i]->getCentre();
			//the distance is calculated in the x-y projection
			double d[2] = {centre[0] - Elementcentre[0], centre[1] - Elementcentre[1]};
			double dmag2 = d[0]*d[0] + d[1]*d[1];
			if (dmag2 > innerRadius2 && dmag2 < outerRadius2){
				float distance = pow(dmag2,0.5);
				double sf = (1.0 - (distance - innerRadius) / (outerRadius - innerRadius) );
				double shapeChangeScale = maxValue* sf;
				Elements[i]->changeShape(shapeChangeScale,axis);
				float timescale = 60*60/dt;
				Elements[i]->updateShapeChangeRate(shapeChangeScale*timescale,axis);
			}
		}
	}
}
*/
void Simulation::TissueAxisPositionDisplay(){
	cerr<<"DV border: "<<endl;
	for (int i=0;i<Nodes.size();++i){
		double x= Nodes[i]->Position[0];
		if (x < 0.2 && x > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
	cerr<<"AP border: "<<endl;
	for (int i=0;i<Nodes.size();++i){
		double y= Nodes[i]->Position[1];
		if (y < 0.2 && y > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
}

void Simulation::CoordinateDisplay(){
	for (int i=0;i<Elements.size();++i){
		int type =Elements[i]-> getShapeType();
		if(type ==1){
			for(int j=0;j<3;++j){
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
	for (int i=0;i<Elements.size();++i){
		int type =Elements[i]-> getShapeType();
		if(type ==1){
			for(int j=3;j<6;++j){
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[Elements[i]->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
}

void Simulation::setStretch(){
	double DVmin = 1000.0;
	double DVmax = -1000.0;
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
		if (Nodes[i]->Position[0]> StretchMax || Nodes[i]->Position[0] < StretchMin){
			Nodes[i]->FixedPos[0]=1;
			Nodes[i]->FixedPos[1]=1;
			Nodes[i]->FixedPos[2]=1;
			if (Nodes[i]->Position[0]<DVmin){
				DVmin = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[0]>DVmax){
				DVmax = Nodes[i]->Position[0];
			}
		}
	}
	double distance = DVmax - DVmin;
	cerr<<"Total DV distance: "<<distance<<" ";
	//the distance that is to be moved:
	distance *= StretchStrain;
	cerr<<"the distance that is to be moved: "<<distance<<" ";
	//the time steps that the stretch operation should take place in:
	double StretchTimeSteps = StretchEndStep - StretchInitialStep;
	cerr<<"stretchTimeSteps: "<<StretchTimeSteps<<" ";
	StretchVelocity = distance / StretchTimeSteps;
	cerr<<"stretchVelocity: "<<StretchVelocity<<endl;

}

void Simulation::addStretchForces(int RKId){
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
		if (Nodes[i]->Position[0]> StretchMax){
			SystemForces[RKId][i][0]=StretchVelocity*Nodes[i]->Viscosity*Nodes[i]->mass/dt;
			SystemForces[RKId][i][1]=0.0;
			SystemForces[RKId][i][2]=0.0;
		}
		else if( Nodes[i]->Position[0] < StretchMin ){
			SystemForces[RKId][i][0]=(-1.0)*StretchVelocity*Nodes[i]->Viscosity*Nodes[i]->mass/dt;
			SystemForces[RKId][i][1]=0.0;
			SystemForces[RKId][i][2]=0.0;
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
	int n = Nodes.size();
	if (ApicalSuction){
		for (int i=0; i<n; ++i){
			//fix basal nodes of columnar layer:
			if(Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0) {
				fixAllD(i);
			}
		}
	}
	else{
		for (int i=0; i<n; ++i){
			//fix apical nodes and all peripodial membrane nodes:
			if(Nodes[i]->tissuePlacement == 1 || Nodes[i]->tissueType == 1) {
				fixAllD(i);
			}
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

    //cout<<"in add pipette forces, pipette pos: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<endl;
    int dim = 3;
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
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

void Simulation::packToPipetteWall(gsl_matrix* gExt, gsl_vector* gSum){
    double zeroThres = 1E-5;
    double threshold  = 1.2*pipetteThickness/2.0;	 //pipette packing forces start at 2 microns
    double t2 = threshold*threshold;
    double multiplier = 10.0;
    double cap = 1E3;
    double pipRange2[2] = {pipetteRadius, pipetteRadius+2.0*threshold};
    double alignmentZ = 0.0;
    pipRange2[0] *= pipRange2[0];
    pipRange2[1] *= pipRange2[1];
    int n = Nodes.size();
    int dim = 3;
    for (int i=0; i<n; ++i){
        bool thereIsPacking = false;
        //v0 is the vector frpm the pipete centre to the point
        double v0[3] = {Nodes[i]->Position[0]-pipetteCentre[0], Nodes[i]->Position[1]-pipetteCentre[1], Nodes[i]->Position[2]-pipetteCentre[2]};       
        double v02[3] = {v0[0]*v0[0], v0[1]*v0[1],v0[2]*v0[2]};
        double xydist2 = v02[0]+v02[1];
        if (xydist2 > pipRange2[0] && xydist2<pipRange2[1]){
            //cout<<"Node "<<id<<" is within range - pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
            //the node is within the packing ring range
            // Now check if the node is inside the pipette or not:
            // if the suction is apical, to be outside, the z position should be lower than pipette tip, v0[2]<0,
            // for basal suction, the z-position should be higher than pipette tip, v0[2]>0
            if((ApicalSuction && v0[2]>0) || (!ApicalSuction && v0[2]<0)){
                //node falls towards the inside the pipette, the alignment z should be the z height of the poitn:
                alignmentZ = v0[2]+0.5;
                if(ApicalSuction){
                  alignmentZ += 0.5;   //always pushing it slightly down as well as to the side
                }
                else{
                    alignmentZ -=0.5;
                }

                thereIsPacking = true;
                //cout<<" Pipette packing to node: "<<i<<" inside the pipette"<<endl;
            }
            else{
                //node falls towards the outside the pipette
                //there will be packing only if the node is close enough in z
                //if there is packing, z-alignment should be to pipette tip, i.e the difference being zero:
                if (v02[2]<t2){
                    thereIsPacking = true;
                    alignmentZ = 0;
                    //cout<<" Pipette packing to node: "<<i<<" outside the pipette"<<endl;
                }
            }
            if (thereIsPacking){
                //if suction strated, I will increase pushing forces accordingly:
                //if (timestep>= PipetteInitialStep && timestep<PipetteEndStep){
                //    multiplier = 0.01*SuctionPressure[2]*Nodes[i]->zProjectedArea;
                //}
                //Now I need to calculate the point on the circle (the pipette surface cross-section)
                double v1[3] = {v0[0],v0[1],alignmentZ};
                double v1norm = v1[0]*v1[0] + v1[1]*v1[1]; //magnitude in xy plane.
                v1norm = pow(v1norm,0.5);
                v1norm /= (pipetteRadius+threshold); // packing from the centre of the pipette wall
                v1[0] /= v1norm;
                v1[1] /= v1norm; //do not scale the z. the point should be at the height of zAlignment.
                //now v1 is the point on pipette surface circle, that is closest to the point.
                //calculating the vector d, from the closest poitn on pipette surface to the point:
                //the force vector should point away from the pipette wall, therefore calculated as v0-v1:
                double d[3] = {v0[0]-v1[0], v0[1]-v1[1], v0[2]-v1[2]};
                //cout<<"Node: "<<i<<" v0: "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<endl;
                //cout<<"Node: "<<i<<" v1: "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<endl;
                //cout<<"Node: "<<i<<" d:  "<<d[0] <<" "<<d[1] <<" "<<d[2] <<endl;

                double dmag2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                double dmag = pow(dmag2,0.5);
                if (dmag < zeroThres){
                    dmag = zeroThres;
                    dmag2 = zeroThres*zeroThres;
                    d[0]= 0; d[1]= 0; d[2] = zeroThres;
                    if (ApicalSuction) {d[2] *= -1.0;}
                }
                //cout<<"vec0: "<<v0[0]<<" "<<v0[1]<<" "<<v0[2]<<endl;
                //cout<<"vec1: "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<endl;

                //cout<<"distance vec: "<<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;
                //double Fmag = multiplier * (1.0/dmag2 - 1.0/t2);

                /* coverintg for new method trial:
                double Fmag = multiplier * (1.0/dmag - 1.0/threshold);
                //cout<<"dmag: "<<dmag<<" threshold: "<<threshold<<" 1/dmag: "<<1.0/dmag<<" 1/th "<<1.0/threshold<<" Fmag: "<<Fmag<<endl;
                if (Fmag > cap){
                    Fmag = cap;
                }
                //cout<<"Fmag: "<<Fmag<<" Force: "<<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;
                //cout<<" (1.0/dmag - 1.0/threshold) "<<(1.0/dmag - 1.0/threshold)<<" multiplier "<<multiplier<<endl;
                Fmag *= Nodes[i]->mass;
                //now converting d to the force vector, divide by dmag and multiply by Fmag ( or divide by dmag/Fmag)

                dmag /=Fmag;
                d[0] /= dmag;
                d[1] /= dmag;
                d[2] /= dmag;
                PackingForces[0][i][0] = d[0];
                PackingForces[0][i][1] = d[1];
                PackingForces[0][i][2] = d[2];
                */
                //new method trial:
                //making d the unit vector from pipette to node
                d[0] /= dmag;
                d[1] /= dmag;
                d[2] /= dmag;
                //get the total current force on node:
                double FNode[3] ={gsl_vector_get(gSum,i*3), gsl_vector_get(gSum,i*3+1), gsl_vector_get(gSum,i*3+2)};
                //project the force on the unit vector d:
                double Fcostet = FNode[0]*d[0] + FNode[1]*d[1]+FNode[2]*d[2];
                //now generate the force vector in the opposite direction of the projection:
                if (Fcostet < 0){  //only if the node is trying to move towards the inside of the pipette:
                    Fcostet *= -1.0*(1.2*threshold-dmag)/dmag; //increasing the force %20 for safety
                }
                else{
                    Fcostet = 0.0;
                }
                d[0] *= Fcostet;
                d[1] *= Fcostet;
                d[2] *= Fcostet;
                PackingForces[0][i][0] = d[0];
                PackingForces[0][i][1] = d[1];
                PackingForces[0][i][2] = d[2];
                //now write the force to be added (the negative force) to external forces:
                //cout<<"Force of node["<<i<<"] : "<<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;
                for (int j=0; j<dim; j++){
                    double value = gsl_matrix_get(gExt,i*dim+j,0);
                    value +=d[j];
                    gsl_matrix_set(gExt,i*dim+j,0,value);
                }
            }
        }
    }
}

void Simulation::LaserAblate(double OriginX, double OriginY, double Radius){
	vector <int> AblatedNodes;
	vector <int> AblatedElements;
	int n = Nodes.size();
	double thres2 = Radius*Radius;
	for (int i=0; i<n; ++i){
		double dx = Nodes[i]->Position[0]- OriginX;
		double dy = Nodes[i]->Position[1]- OriginY;
		double d2 = dx *dx + dy*dy;
		if (d2 < thres2){
			AblatedNodes.push_back(i);
		}
	}
	n = Elements.size();
	int nAN = AblatedNodes.size();
	for (int i=0; i<n; ++i){
		if(!Elements[i]->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = Elements[i]->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					Elements[i]->removeMassFromNodes(Nodes);
					Elements[i]->IsAblated = true;
					//cerr<<"Ablating element:" <<Elements[i]->Id<<endl;
					break;
				}
			}
		}
	}
	//some nodes are ledt with zero mass, which will cause problems in later calculations:
	n = Nodes.size();
	for (int i=0; i<n; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}

}



void Simulation::updateElementVolumesAndTissuePlacements(){
	int n = Elements.size();
	for (int i=0; i<n;++i){
		Elements[i]->updateElementVolumesAndTissuePlacementsForSave(Nodes);
	}
}

void Simulation::clearNodeMassLists(){
	int n = Nodes.size();
	for (int i=0 ;i<n;++i){
		Nodes[i]->connectedElementIds.size();
		Nodes[i]->connectedElementIds.clear();
		Nodes[i]->connectedElementWeights.clear();
		Nodes[i]->mass=0.0;
		Nodes[i]->surface=0.0;
	}
}

void Simulation::clearLaserAblatedSites(){
	int nE = Elements.size();
	for (int i=0; i<nE; ++i){
		if (Elements[i]->IsAblated){
			Elements[i]->removeMassFromNodes(Nodes);
		}
	}
	int n = Nodes.size();
	for (int i=0; i<n; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}
