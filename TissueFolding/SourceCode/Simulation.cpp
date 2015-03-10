
#include "Simulation.h"
#include "Prism.h"
#include "Triangle.h"
#include <string.h>
#include <algorithm>    // std::sort
#include <vector>

using namespace std;

Simulation::Simulation(){
	currElementId = 0;
	ModInp = new ModelInputObject();
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	TissueHeight = 0.0;
	TissueHeightDiscretisationLayers= 1;
	timestep = 0;
	reachedEndOfSaveFile = false;
	AddPeripodium = false;
	setDefaultParameters();
};

Simulation::~Simulation(){
	//cerr<<"destructor for simulation called"<<endl;
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

	while(!Elements.empty()){
		ShapeBase* tmp_pt;
		tmp_pt = Elements.back();
		Elements.pop_back();
		delete tmp_pt;
	}
	while(!Nodes.empty()){
		Node* tmp_pt;
		tmp_pt = Nodes.back();
		Nodes.pop_back();
		delete tmp_pt;
	}
	//cerr<<"destructor for simulation finalised"<<endl;
};

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
	ForcesSaved = true;
	VelocitiesSaved = true;
	PeripodiumElasticity = 0.0;
	DVRight = 0;
	DVLeft = 1;
	stretcherAttached = false;
	StretchInitialStep = -100;
	StretchEndStep = -100;
	StretchVelocity = 0.0;
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
		cerr<<" input the mode of simulation: {DisplaySave, SimulationOnTheGo, Default}"<<endl;
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
	else{
		cerr<<"Please provide input mode: -mode {DisplaySave, SimulationOnTheGo, Default}";
		return false;
	}
}

bool Simulation::readParameters(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the model input file"<<endl;
		return false;
	}
	cerr<<"Reading parameter input file"<<endl;
	const char* inpstring = argv[i];
	ModInp->Sim=this;
	ModInp->parameterFileName =  inpstring;
	cerr<<" Reading parameters from file: "<<ModInp->parameterFileName<<endl;
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

bool Simulation::checkInputConsistency(){
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
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	if (AddPeripodium){
		Success = addPeripodiumToTissue();
	}
	if (!Success){
		return Success;
	}
	fillInNodeNeighbourhood();
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
	calculateStiffnessMatrices();
	assignNodeMasses();
	assignConnectedElementsAndWightsToNodes();
	//for (int i=0; i<Nodes.size();++i){
	//	cout<<"Node: "<<i<<endl;
	//	Nodes[i]->displayConnectedElementIds();
	//	Nodes[i]->displayConnectedElementWeights();
	//}
	if (AddPeripodium){
		assignMassWeightsDueToPeripodium();
	}
	if (stretcherAttached){
		setStretch();
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
	}
}

void Simulation::assignConnectedElementsAndWightsToNodes(){
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
		saveFileTensionCompression.open(name_saveFileTenComp, ofstream::binary);
		if (saveFileTensionCompression.good() && saveFileTensionCompression.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileTenComp<<endl;
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
	}
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
		if (GrowthFunctionTypes[i] == 1){
			saveFileSimulationSummary<<"Growth Type:  Uniform (1)"<<endl;
			saveFileSimulationSummary<<"	Initial time(sec): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex];
			saveFileSimulationSummary<<"	FinalTime time(sec): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+1];
			saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+2]/dt*3600.0;
			saveFileSimulationSummary<<"  ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+3]/dt*3600.0;
			saveFileSimulationSummary<<"  ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+4]/dt*3600.0;
			saveFileSimulationSummary<<endl;
			currIndex += 5;
		}
		else if(GrowthFunctionTypes[i] == 2){
			saveFileSimulationSummary<<"Growth Type:  Ring (2)"<<endl;
			saveFileSimulationSummary<<"	Initial time(sec): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex];
			saveFileSimulationSummary<<"	FinalTime time(sec): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+1];
			saveFileSimulationSummary<<"	Centre(micron): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+2];
			saveFileSimulationSummary<<"  ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+3];
			saveFileSimulationSummary<<"	Inner radius(micron): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+4];
			saveFileSimulationSummary<<"	Outer radius(micron): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+5];
			saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+6]/dt*3600.0;
			saveFileSimulationSummary<<"  ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+7]/dt*3600.0;
			saveFileSimulationSummary<<"  ";
			saveFileSimulationSummary<<GrowthParameters[currIndex+8]/dt*3600.0;
			saveFileSimulationSummary<<endl;
		}
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
	const char* name_saveFileToDisplayTenComp = saveFileString.c_str();;
	saveFileToDisplayTenComp.open(name_saveFileToDisplayTenComp, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayTenComp<<endl;
		TensionCompressionSaved = false;
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
	initiateSystemForces();
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if (VelocitiesSaved){
		updateVelocitiesFromSave();
	}
	updateElementVolumesAndTissuePlacements();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWightsToNodes();
	clearLaserAblatedSites();
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
	dataSaveInterval =  dummydouble/dt;

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
	//CORRECT THIS LATER ON!!!
	SystemCentre[2] -= 1; // The apical surface is assumed to took towards (-)ve z
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
	SystemCentre[2] -= 1; // The apical surface is assumed to took towards (-)ve z
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
		for (int j=0; j<3; ++j){
			saveFileToDisplayTenComp.read((char*) &Elements[i]->StrainTissueMat(j,j), sizeof Elements[i]->StrainTissueMat(j,j));
		}
		for (int j=0; j<3; ++j){
			saveFileToDisplayTenComp.read((char*) &Elements[i]->CurrPlasticStrainsInTissueCoordsMat(j,j), sizeof Elements[i]->CurrPlasticStrainsInTissueCoordsMat(j,j));
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
	SystemCentre[2] -= 1; // The apical surface is assumed to took towards (-)ve z
	TrianglePnt01->AlignReferenceApicalNormalToZ(SystemCentre);  //correcting the alignment of the triangular element such that the apical side will be aligned with (+)ve z
	calculateSystemCentre();
	vector<ShapeBase*>::iterator it = Elements.begin();
	it += k;
	Elements.insert(it,TrianglePnt01);
	currElementId++;
	delete[] NodeIds;
}


void Simulation::updateOneStepFromSave(){
	//cout<<"Updating step from save:"<<endl;
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
	//cout<<"skipped header: "<<currline<<endl;
	updateNodeNumberFromSave();
	cout<<"updated node number: "<<Nodes.size();
	updateNodePositionsFromSave();
	updateElementStatesFromSave();
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		//This is updating positions from save, I am only interested in normal positions, no Runge-Kutta steps. This corresponds to RK step 4, RKId = 3
		Elements[i]->updatePositions(3, Nodes);
	}
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if(VelocitiesSaved){
		updateVelocitiesFromSave();
	}
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWightsToNodes();
	clearLaserAblatedSites();
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
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
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
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}



bool Simulation::addPeripodiumToTissue(){
	bool Success = true;
	Success = generateColumnarCircumferenceNodeList();
	if (!Success){
		return Success;
	}
	calculateSystemCentre();
	sortColumnarCircumferenceNodeList();
	vector <int*> trianglecornerlist;
	double d = getAverageSideLength();
	Success = CalculateTissueHeight();
	if (!Success){
		return Success;
	}
	addPeripodiumNodes(trianglecornerlist, TissueHeight, d);
	FillNodeAssociationDueToPeripodium();
	addPeripodiumElements(trianglecornerlist, TissueHeight);
	return Success;
	//addMassToPeripodiumNodes();
	//distributeCircumferenceMass();
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
		cerr<<"No circumferncial nodes indicated! Cannot generate peripodium"<<endl;
		AddPeripodium = false;
		return false;
	}
}

void Simulation::sortColumnarCircumferenceNodeList(){
	//ordering the circumferencial nodes of the basal surface in clockwise rotation
	//sort (ColumnarCircumferencialNodeList.begin(), ColumnarCircumferencialNodeList.end(), bind(&Simualtion::sortcomparison, this, _1, _2));
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
	/*
	cout<<"ColumnarCircumferencialNodeList after sort"<<endl;
	for (int j =0 ; j<n; ++j){
		cout<<ColumnarCircumferencialNodeList[j]<<" "<<angles[j]<<endl;
	}*/
}

bool Simulation::sortcomparison(int i, int j){
/*	if ((Nodes[i]->Position[0] - SystemCentre[0] ) >= 0 && ( Nodes[j]->Position[0] - SystemCentre[0] )< 0)
		return true;
	if ((Nodes[i]->Position[0] - SystemCentre[0] ) < 0 && ( Nodes[j]->Position[0] - SystemCentre[0] ) >= 0)
		return false;
	if ((Nodes[i]->Position[0] - SystemCentre[0] ) == 0 && (Nodes[j]->Position[0] - SystemCentre[0]) == 0) {
		if ( (Nodes[i]->Position[1] - SystemCentre[1]) >= 0 || (Nodes[j]->Position[1] - SystemCentre[1]) >= 0)
			return Nodes[i]->Position[1] > Nodes[j]->Position[1];
		return Nodes[j]->Position[1] > Nodes[i]->Position[1];
	}

	// compute the cross product of vectors (center -> a) x (center -> b)
	int det = (Nodes[i]->Position[0] - SystemCentre[0]) * (Nodes[j]->Position[1] - SystemCentre[1]) - (Nodes[j]->Position[0] - SystemCentre[0]) * (Nodes[i]->Position[1] - SystemCentre[1]);
	if (det < 0)
		return true;
	if (det > 0)
		return false;

	// points a and b are on the same line from the center
	// check which point is closer to the center
	int d1 = (Nodes[i]->Position[0] - SystemCentre[0]) * (Nodes[i]->Position[0] - SystemCentre[0]) + (Nodes[i]->Position[1] - SystemCentre[1]) * (Nodes[i]->Position[1] - SystemCentre[1]);
	int d2 = (Nodes[j]->Position[0] - SystemCentre[0]) * (Nodes[j]->Position[0] - SystemCentre[0]) + (Nodes[j]->Position[1] - SystemCentre[1]) * (Nodes[j]->Position[1] - SystemCentre[1]);
	return d1 > d2;
*/
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

void Simulation::AddPeripodiumCircumference(double height, int& index_begin, int &index_end){
	int n = ColumnarCircumferencialNodeList.size();
	for (int i=0; i<n; ++i){
		int index = ColumnarCircumferencialNodeList[i];
		double* pos;
		pos = new double[3];
		pos[0] = Nodes[index]->Position[0];
		pos[1] = Nodes[index]->Position[1];
		pos[2] = Nodes[index]->Position[2] + height/2.0;
		Node* tmp_nd;
		tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is peripodium, node type is basal
		Nodes.push_back(tmp_nd);
		if (i==0){index_begin = tmp_nd->Id;}
		else if (i==n-1){index_end = tmp_nd->Id;}
		PeripodiumCircumferencialNodeList.push_back(tmp_nd->Id);
		tmp_nd->atPeripodiumCircumference = true;
		//cerr<<"NodeId: "<<tmp_nd->Id<<" pos: "<<tmp_nd->Position[0]<<" "<<tmp_nd->Position[1]<<" "<<tmp_nd->Position[2]<<endl;
		delete[] pos;
	}
}

void Simulation::AddHorizontalRowOfPeripodiumNodes(vector <int*> &trianglecornerlist, double d, int &index_begin, int &index_end){
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
		tmp_nd = new Node(Nodes.size(), 3, norm,0,1);  //tissue type is peripodium, node type is basal
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

void Simulation::AddVerticalRowOfPeripodiumNodes(int& layerCount, int nLayers, vector <int*> &trianglecornerlist, double height,  int &index_begin, int &index_end){
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
			targetpoint[2] = Nodes[targetindex]->Position[2] + 1.5*height ;
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
			targetpoint[2] = 0.5 * (Nodes[targetindex1]->Position[2] + Nodes[targetindex2]->Position[2] ) + 1.5*height ;
		}
		double vec[3] = {targetpoint[0] - MidPoint[0], targetpoint[1] - MidPoint[1], targetpoint[2] - MidPoint[2]};
		double mag = pow((vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]),0.5);
		vec[0] /= mag;
		vec[1] /= mag;
		vec[2] /= mag;
		//The height spanned in this direction should be height of the tissue / number of layers:
		double z = height / (float) nLayers;
		double multiplier = z / vec[2];
		double* norm;
		norm = new double[3];
		norm[0] = MidPoint[0]+vec[0]*multiplier;
		norm[1] = MidPoint[1]+vec[1]*multiplier;
		norm[2] = MidPoint[2]+z;
		Node* tmp_nd;
		//Adding Peroipodium node:
		tmp_nd = new Node(Nodes.size(), 3, norm,0,1);  //tissue type is peripodium, node type is basal
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

void Simulation::AddPeripodiumCap(int layerCount,  vector <int*> &trianglecornerlist, double height, int index_begin, int index_end){
	//Now I have the indices of the nodes specifying the last row.
	//I want to cap the tissue, with the topology of the apical surfaces of the columnar layer
	vector <int> PeripodiumNodeId;
	vector <int> CorrespondingApicalNodeId;
	int n = ApicalColumnarCircumferencialNodeList.size();
	//map the circumference to peripodium nodes:
	int counter =0;
	for (int i = index_begin; i <= index_end; ++i){
		int idx = counter + layerCount*0.5 +1;
		if (idx >= n){
			idx -= n;
		}
		counter++;
		PeripodiumNodeId.push_back(Nodes[i]->Id);
		CorrespondingApicalNodeId.push_back(Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id);
		cout<<"peripodium Node: "<<Nodes[i]->Id<<" pos: "<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<" corr. node id: "<<Nodes[ApicalColumnarCircumferencialNodeList[idx]]->Id<<endl;
	}
	// The end point is 0.5*height above the apical surface.
	// I also want to add one more layer of curvature, same height as each triangle in the side peripodium spans
	// the layer count is equal to (nLayers -1 ) at the moment, as I am at the topmost layer
	// I can obtain the increment I need from adding the sum of these as a z offset:
	double zOffset = height*0.5 + height / (float) (layerCount+1);
	for (int i = 0; i<Nodes.size(); ++i){
		if (Nodes[i]->tissuePlacement == 1){ //Node is apical
			int id = Nodes[i]->Id;
			bool AtCircumference = false;
			for (int j =0 ; j< n; ++j){
				if (id == Nodes[ApicalColumnarCircumferencialNodeList[j]]->Id ){
					//if it is not on the circumference of peripodium, skip
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
				tmp_nd = new Node(Nodes.size(), 3, pos,0,1);  //tissue type is peripodium, node type is basal
				Nodes.push_back(tmp_nd);
				PeripodiumNodeId.push_back(tmp_nd->Id);
				CorrespondingApicalNodeId.push_back(Nodes[i]->Id);
				cout<<"temp Node: "<<tmp_nd->Id<<" pos: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" corr. node id: "<<Nodes[i]->Id<<endl;
			}
		}
	}
	/*cout<<"apical node - peripodium node couples: "<<endl;
	for (int a = 0; a< PeripodiumNodeId.size(); a++){
		cout<<PeripodiumNodeId[a]<<" "<<CorrespondingApicalNodeId[a]<<endl;
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
						TriNodeIds[0]=PeripodiumNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+1]){
						TriNodeIds[1]=PeripodiumNodeId[dictionaryIndex];
					}
					if (CorrespondingApicalNodeId[dictionaryIndex] == ApicalTriangles[k+2]){
						TriNodeIds[2]=PeripodiumNodeId[dictionaryIndex];
					}
				}
				trianglecornerlist.push_back(TriNodeIds);
			}
		}
	}
	/*cout<<"triangle corners for peripodium cap: "<<endl;
	for (int a =initialpointfordisplay; a< trianglecornerlist.size(); a++){
		cout<<trianglecornerlist[a][0]<<" "<<trianglecornerlist[a][1]<<" "<<trianglecornerlist[a][2]<<endl;
	}*/
}

void Simulation::FillNodeAssociationDueToPeripodium(){
	//Take the peripodium circumferential list
	//Start from the associated Basal columnar circumferential node
	//Move apically through the elements until you reach the apical surface - an apical node
	int n = PeripodiumCircumferencialNodeList.size();
	int nE = Elements.size();
	for (int i=0; i< n; ++i){
		int index = PeripodiumCircumferencialNodeList[i];
		int currNodeId = ColumnarCircumferencialNodeList[i];
		Nodes[index]->AssociatedNodesDueToPeripodium.push_back(ColumnarCircumferencialNodeList[i]);
		while(Nodes[currNodeId]->tissuePlacement != 1){ //While I have not reached the apical node
			for (int j= 0; j<nE; ++j){
				bool IsBasalOwner = Elements[j]->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					currNodeId = Elements[j]->getCorrecpondingApical(currNodeId);
					break;
				}
			}
			Nodes[index]->AssociatedNodesDueToPeripodium.push_back(currNodeId);
		}
		/*cout<<"Node: "<<index<<" pos : "<<Nodes[index]->Position[0]<<" "<<Nodes[index]->Position[1]<<" "<<Nodes[index]->Position[2]<<endl;
		for(int j = 0; j< Nodes[index]->AssociatedNodesDueToPeripodium.size(); ++j){
			int a = Nodes[index]->AssociatedNodesDueToPeripodium[j];
			cout<<"	Associated Node Pos: "<<Nodes[a]->Position[0]<<" "<<Nodes[a]->Position[1]<<" "<<Nodes[a]->Position[2]<<endl;
		}*/
	}
}


void Simulation::assignMassWeightsDueToPeripodium(){
	//Take the peripodium circumferential list
	//calculate the sum of associated node masses
	//Add the weigthing fractions to AssociatedNodeWeightsDueToPeripodium of each Node
	int n = PeripodiumCircumferencialNodeList.size();
	for (int i=0; i< n; ++i){
		int index = PeripodiumCircumferencialNodeList[i];
		int nA = Nodes[index]->AssociatedNodesDueToPeripodium.size();
		double weightSum = 0.0;
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodium[j];
			double w = Nodes[index2]->mass;
			weightSum += w;
			Nodes[index]->AssociatedNodeWeightsDueToPeripodium.push_back(w);
		}
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodium[j];
			Nodes[index]->AssociatedNodeWeightsDueToPeripodium[j] /= weightSum;
			//Distributing the weight of this node onto the associated nodes
			Nodes[index2]->mass += Nodes[index]->mass*Nodes[index]->AssociatedNodeWeightsDueToPeripodium[j];
		}
	}
}


void Simulation::addPeripodiumNodes(vector <int*> &trianglecornerlist, double height, double d){
	cerr<<"Adding peripodium nodes"<<endl;
	//int n = ColumnarCircumferencialNodeList.size();
	int index_begin = 0, index_end =0;
	//Adding a midline range of nodes
	AddPeripodiumCircumference(height, index_begin, index_end);
	double triangleHeight = 0.866*d; //0.866 is square-root(3)/2, this is the height of the triangle I am adding,
	AddHorizontalRowOfPeripodiumNodes(trianglecornerlist, triangleHeight, index_begin, index_end);
	//calculating how many layers of triangles I need for spanning the height of the tissue:
	int nLayers = ceil(height / triangleHeight );
	//I need an odd number of layers, so that the final shape will be the same as the circumference of aolumnar layer:
	if (nLayers %2 == 0){
		nLayers --;
	}
	//now I want nLayers many layers to cover exactly the height of the tissue, the distance should be the hipotenus divided by nLayers
	//I will adjust the triangle height
	triangleHeight = pow((triangleHeight*triangleHeight + height*height),0.5) / (float)nLayers;
	int layerCount = 0;
	for (int i=0; i<nLayers; ++i){
		AddVerticalRowOfPeripodiumNodes(layerCount, nLayers, trianglecornerlist, height, index_begin, index_end);
	}
	AddPeripodiumCap(layerCount, trianglecornerlist, height, index_begin, index_end);
	//AssignAssociatedNodesToPeripodiumCircumference();
}

double Simulation::getAverageSideLength(){
	double dsum =0.0;
	int n = Elements.size();
	for(int i=0; i<n; ++i){
		dsum += Elements[i]->getApicalSideLengthAverage();
	}
	dsum /= (double)n;
	return dsum;
}

bool Simulation::CalculateTissueHeight(){
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
	while(Nodes[currNodeId]->tissuePlacement != 1 && foundElement){ //while the node is not apical, and I could fing the next element
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
	return true;
}

void Simulation::addPeripodiumElements(vector <int*> &trianglecornerlist, double height){
	int n = trianglecornerlist.size();
	calculateSystemCentre();
	SystemCentre[2] -= 1; // The apical surface is assumed to took towards (-)ve z
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
		if (Elements[i]->tissueType == 1){ //Element is on the peripodium
			double currE = PeripodiumElasticity*(1 + noise1/100.0);
			Elements[i]->setElasticProperties(currE,currE,currE,poisson*(1 + noise2/100));
		}
	}
	for (int i=0; i<Nodes.size(); ++i){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		Nodes[i]->setViscosity(ApicalVisc*noise3, BasalVisc*noise3);
	}
}

void Simulation::runOneStep(){
	/*if(timestep==0){
		for(int i=0;i<Nodes.size();++i){
			if (Nodes[i]->atCircumference){
				Nodes[i]->FixedPos[0] = true;
				Nodes[i]->FixedPos[1] = true;
				Nodes[i]->FixedPos[2] = true;
			}
			Nodes[i]->Position[0] *=1.5;
			Nodes[i]->Position[1] *=1.5;
			Nodes[i]->Position[2] *=1.0;
		}
		double R[3][3];
		double Rx[3][3] = {{1,0,0},{0,0,-1},{0,1,0}};
		double Ry[3][3] = {{0,0,1},{0,1,0},{-1,0,0}};
		double Rz[3][3] = {{0,-1,0},{1,0,0},{0,0,1}};
		for (int j =0; j<3;++j){
			for(int k=0; k<3;++k){
				R[j][k] = Rx[j][k];
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
		}
	}*/
	if(timestep==-10/dt){
		LaserAblate(0.0,0.0,1.8);
	}
	int displayfreq = 60/dt;
	if (timestep%displayfreq == 0){
		//reOpenOutputFile();
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
	cleanMatrixUpdateData();
	cleanGrowthData();
	resetForces();
	//alignTissueDVToXPositive();
	int freq = 10.0/dt + 1 ;
	if (timestep% freq  == 0){
		calculateDVDistance();
		//outputFile<<"calculating element health"<<endl;
		int nElement = Elements.size();
		for (int i=0; i<nElement; ++i){
			Elements[i]->checkHealth();
		}
	}
	int nElement = Elements.size();
	if(nGrowthFunctions>0){
		//outputFile<<"calculating growth"<<endl;
		calculateGrowth();
	}
	//if(nShapeChangeFunctions>0){
	//	changeCellShapesInSystem();
	//}
	//outputFile<<"calculating alignment of reference"<<endl;
	for (int i=0; i<nElement; ++i){
		Elements[i]->alignElementOnReference();
		if (Elements[i]->IsGrowing){
			Elements[i]->growShape();
		}
	}

	double packingThreshold = getAverageSideLength();
	packingThreshold *=0.7; //I am adding a 40% safety to average side length, and then taking half of it as threshold for packing
	for (int RKId = 0; RKId<4; ++RKId){
		//outputFile<<"started RK: "<<RKId<<endl;
		//calculatePacking(RKId,packingThreshold);
		//cout<<"Packing Forces: "<<endl;
		//for (int p = 648; p<651; ++p){
		//	cout<<"Node: "<<p<<" RK: "<<RKId<<" pos: "<<Nodes[p]->Position[0]<<" "<<Nodes[p]->Position[1]<<" "<<Nodes[p]->Position[2]<<" f: "<<PackingForces[RKId][p][0]<<" "<<PackingForces[RKId][p][1]<<" "<<PackingForces[RKId][p][2]<<endl;
		//}
		//cout<<"After Packing, systemForces:"<<endl;
		//cout<<"node 497 - at peri: "<<Nodes[497]->atPeripodiumCircumference<<" pos: "<<Nodes[497]->Position[0]<<" "<<Nodes[497]->Position[1]<<" "<<Nodes[497]->Position[2]<<endl;
		//cout<<"	forces: "<<SystemForces[RKId][497][0]<<" "<<SystemForces[RKId][497][1]<<" "<<SystemForces[RKId][497][2]<<endl;
		//cout<<"node 394 - at peri: "<<Nodes[394]->atPeripodiumCircumference<<" pos: "<<Nodes[394]->Position[0]<<" "<<Nodes[394]->Position[1]<<" "<<Nodes[394]->Position[2]<<endl;
		//cout<<"	forces: "<<SystemForces[RKId][394][0]<<" "<<SystemForces[RKId][394][1]<<" "<<SystemForces[RKId][394][2]<<endl;

		for (int i=0; i<nElement; ++i){
			if (!Elements[i]->IsAblated){
				Elements[i]->calculateForces(RKId, SystemForces, Nodes, outputFile);
			}
		}
		//outputFile<<"     calculated forces"<<endl;
		if (stretcherAttached && timestep>StretchInitialStep && timestep<StretchEndStep){
			addStretchForces(RKId);
		}
		redistributePeripodiumForces(RKId);
		updateNodePositions(RKId);
		//outputFile<<"     updated node pos"<<endl;
		updateElementPositions(RKId);
		//outputFile<<"     updated element pos"<<endl;
		/*cout<<"After elastic calc, systemForces:"<<endl;
		cout<<"node 497 - at peri: "<<Nodes[497]->atPeripodiumCircumference<<" pos: "<<Nodes[497]->Position[0]<<" "<<Nodes[497]->Position[1]<<" "<<Nodes[497]->Position[2]<<endl;
		cout<<"	forces: "<<SystemForces[RKId][497][0]<<" "<<SystemForces[RKId][497][1]<<" "<<SystemForces[RKId][497][2]<<endl;
		cout<<"node 197 - pos: "<<Nodes[197]->Position[0]<<" "<<Nodes[197]->Position[1]<<" "<<Nodes[197]->Position[2]<<endl;
		cout<<"	forces: "<<SystemForces[RKId][197][0]<<" "<<SystemForces[RKId][197][1]<<" "<<SystemForces[RKId][197][2]<<endl;
		 */
	}
	/*for (int RKId = 0; RKId<4; ++RKId){
		for (int i=0; i<nElement; ++i){
			Elements[i]->calculateForces(RKId, SystemForces, Nodes);
			Elements[i]->updateElementsNodePositions(RKId, SystemForces, Nodes, dt);
			updateElementPositions(RKId);

		}
	}*/
	/*cout<<"Forces on node 0: "<<endl;
	for (int i=0;i<4;++i){
		cout<<" RK: "<<i<<" ";
		for (int j=0;j<3;++j){
			cout<<SystemForces[i][0][j]<<" ";
		}
		cout<<endl;
	}
	*/

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

void Simulation::calculatePacking(int RKId, double threshold){
	cout<<"inside calculate packing"<<endl;
	//threshold = 2.1;
	double multiplier = 0.0001;
	double t2 = threshold*threshold;
	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	//resetting all the normal update flags
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->NormalForPackingUpToDate = false;
	}
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if (!(*itNode)->atPeripodiumCircumference){
			bool memberOf44 = Elements[44]->DoesPointBelogToMe((*itNode)->Id);
			double* pos;
			pos = new double[3];
			if (RKId == 0){
				pos[0] = (*itNode)->Position[0];
				pos[1] = (*itNode)->Position[1];
				pos[2] = (*itNode)->Position[2];
			}
			else{
				pos[0]= (*itNode)->RKPosition[0];
				pos[1]= (*itNode)->RKPosition[1];
				pos[2]= (*itNode)->RKPosition[2];
			}
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this element
				bool IsNodeMemberOfThisElement = (*itEle)->DoesPointBelogToMe((*itNode)->Id);
				bool checkPackingToThisElement = false;
				//the node does not belong to the element, lets have a preliminary distance check:
				if (!IsNodeMemberOfThisElement){
					checkPackingToThisElement = (*itEle)->IsPointCloseEnoughForPacking((*itNode)->Position, threshold);
				}
				if (checkPackingToThisElement){
					if (!(*itEle)->NormalForPackingUpToDate){
						(*itEle)->calculateNormalForPacking();
					}
					double* posCorner;
					posCorner = new double[3];
					(*itEle)->getApicalNodePos(posCorner);
					double dProjectionOffet = 0.0;
					for (int i=0; i<3; ++i){
						dProjectionOffet += (pos[i] - posCorner[i])*(*itEle)->normalForPacking[i];
					}
					//In the case that there is packing (the point will lie on the triangle, the distance
					//I will be interested in is dProjectionOffet. I will check if this is below threshold,
					//then I will bother with the following calculation if necessary
					if ((dProjectionOffet > 0 && dProjectionOffet < threshold) || (dProjectionOffet < 0 && dProjectionOffet > (-1.0)*threshold)){
						//cout<<"checking element "<<(*itEle)->Id<<" for node: "<<(*itNode)->Id<<endl;
						double VecprojectedPoint[3];
						VecprojectedPoint[0] = (-1.0)*dProjectionOffet*((*itEle)->normalForPacking[0]);
						VecprojectedPoint[1] = (-1.0)*dProjectionOffet*((*itEle)->normalForPacking[1]);
						VecprojectedPoint[2] = (-1.0)*dProjectionOffet*((*itEle)->normalForPacking[2]);
						double projectedPoint[3] = {pos[0] + VecprojectedPoint[0], pos[1] + VecprojectedPoint[1], pos[2] + VecprojectedPoint[2]};
						//check if the projected point lies on the triangle, if so the distance is (-1.0)*dProjectionOffet
						bool pointInsideTriangle = (*itEle)->IspointInsideApicalTriangle(projectedPoint[0],projectedPoint[1],projectedPoint[2]);
						/*if(memberOf44 &&  (*itEle)->Id==592){
							cout<<"checking packing to element: "<<(*itEle)->Id<<endl;
							cout<<"	normalForPacking: "<<(*itEle)->normalForPacking[0]<<" "<<(*itEle)->normalForPacking[1]<<" "<<(*itEle)->normalForPacking[2]<<endl;
							cout<<"	posCorner: "<<posCorner[0]<<" "<<posCorner[1]<<" "<<posCorner[2]<<endl;
							cout<<"	dProjectionOffet: "<<dProjectionOffet<<endl;
							cout<<"	normalForPacking: "<<(*itEle)->normalForPacking[0]<<" "<<(*itEle)->normalForPacking[1]<<" "<<(*itEle)->normalForPacking[2]<<endl;
							cout<<"	VecprojectedPoint comp "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							cout<<"	VecprojectedPoint manu "<<(-1.0)*dProjectionOffet*((*itEle)->normalForPacking[0])<<" "<<(-1.0)*dProjectionOffet*((*itEle)->normalForPacking[1])<<" "<<(-1.0)*dProjectionOffet*((*itEle)->normalForPacking[2])<<endl;
							cout<<"	projectedPoint "<<projectedPoint[0]<<" "<<projectedPoint[1]<<" "<<projectedPoint[2]<<endl;
							cout<<"	point inside triangle: "<<pointInsideTriangle<<endl;
						}*/
						if (pointInsideTriangle){
							//there is packing:
							float d2 = dProjectionOffet*dProjectionOffet;
							//normalising force direction
							if (dProjectionOffet<0){
								VecprojectedPoint[0] /= (-1.0)*dProjectionOffet;
								VecprojectedPoint[1] /= (-1.0)*dProjectionOffet;
								VecprojectedPoint[2] /= (-1.0)*dProjectionOffet;
							}
							else{
								VecprojectedPoint[0] /= dProjectionOffet;
								VecprojectedPoint[1] /= dProjectionOffet;
								VecprojectedPoint[2] /= dProjectionOffet;
							}
							double Fmag = multiplier * (1.0/d2 - 1.0/t2);
							double F[3];
							F[0] = Fmag * VecprojectedPoint[0];
							F[1] = Fmag * VecprojectedPoint[1];
							F[2] = Fmag * VecprojectedPoint[2];
							//direction correction with (-) sign
							for(int j=0; j<3; ++j){
								if (!(*itNode)->FixedPos[j]){
									SystemForces[RKId][(*itNode)->Id][j] -= F[j];
									PackingForces[RKId][(*itNode)->Id][j] -= F[j];
								}
							}
							//add opposite force to the element nodes:
							(*itEle)->AddPackingToApicalSurface(F[0],F[1],F[2],RKId, SystemForces,PackingForces, Nodes);
							//if(/*(*itNode)->Id == 13 && (*itEle)->Id == 44*/ memberOf44 && (*itEle)->Id>592){
							//	cout<<"threshold: "<<threshold<<" distance: "<<dProjectionOffet<<" Node: "<<(*itNode)->Id<<" Element: "<<(*itEle)->Id<<" forcemag: "<<Fmag;
							//	cout<<" Fdir: "<<VecprojectedPoint[0]<<" "<<VecprojectedPoint[1]<<" "<<VecprojectedPoint[2]<<endl;
							//	cout<<"	pos:"<<(*itNode)->Position[0]<<" "<<(*itNode)->Position[1]<<" "<<(*itNode)->Position[2]<<endl;
							//}
						}
					}
					delete[] posCorner;
				}
			}
			delete[] pos;
		}
	}
	cout<<"finalised calculate packing"<<endl;
}

/*	float multiplier = 100.0;
	vector<Node*>::iterator itNode0;
	vector<Node*>::iterator itNode1;
	float dmin = threshold;
	float dminNeg =(-1.0)*threshold;
	float t2 = threshold*threshold;
	for (itNode0=Nodes.begin(); itNode0<Nodes.end(); ++itNode0){
		if (!(*itNode0)->atPeripodiumCircumference){
			for (itNode1=itNode0+1; itNode1<Nodes.end(); ++itNode1){
				if (!(*itNode1)->atPeripodiumCircumference){
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


void Simulation::redistributePeripodiumForces(int RKId){
	int n = PeripodiumCircumferencialNodeList.size();
	for (int i=0; i<n; ++i){
		int index = PeripodiumCircumferencialNodeList[i];
		double F[3] = {SystemForces[RKId][index][0], SystemForces[RKId][index][1], SystemForces[RKId][index][2]};
		int nA = Nodes[index]->AssociatedNodesDueToPeripodium.size();
		for(int j = 0; j < nA; ++j){
			int index2 = Nodes[index]->AssociatedNodesDueToPeripodium[j];
			double w = Nodes[index]->AssociatedNodeWeightsDueToPeripodium[j];
			SystemForces[RKId][index2][0] += F[0]*w;
			SystemForces[RKId][index2][1] += F[1]*w;
			SystemForces[RKId][index2][2] += F[2]*w;
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
		for (int i=0;i<n;++i){
			for (int j=0; j<Nodes[i]->nDim; ++j){
					Nodes[i]->Velocity[RKId][j] = SystemForces[RKId][i][j]/ (Nodes[i]->Viscosity*Nodes[i]->mass) ;
					Nodes[i]->RKPosition[j] = Nodes[i]->Position[j] + Nodes[i]->Velocity[RKId][j]*multiplier*dt;
			}
			//cout<<"RK: "<<RKId<<" node: "<<i<<"mass: "<<Nodes[i]->mass<<" visc: "<<Nodes[i]->Viscosity<<endl;
			//for (int j=0; j<Nodes[i]->nDim; ++j){
			//	cout<<"	"<<Nodes[i]->Velocity[RKId][j]<<" ";
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
		for (int i=0;i<n;++i){
		//	cout<<"Nodes "<<i<<" velocity: ";
			for (int j=0; j<Nodes[i]->nDim; ++j){
				Nodes[i]->Velocity[RKId][j] = SystemForces[RKId][i][j]/(Nodes[i]->Viscosity*Nodes[i]->mass);
				//now I have 4 velocity data (corresponding to Runge-Kutta  k1, k2, k3, and k4)
				//writing  the velocity into v[0]
				//cout<<Nodes[i]->Velocity[0][j]<<" "<<Nodes[i]->Velocity[1][j]<<" "<<Nodes[i]->Velocity[2][j]<<" "<<Nodes[i]->Velocity[3][j]<<" ";
				Nodes[i]->Velocity[0][j] = 1.0/6.0 * (Nodes[i]->Velocity[0][j] + 2.0 * (Nodes[i]->Velocity[1][j] + Nodes[i]->Velocity[2][j]) + Nodes[i]->Velocity[3][j]);
				Nodes[i]->Position[j] += Nodes[i]->Velocity[0][j]*dt;
			}
		//	cout<<endl;
		}
	}
	updateNodePositionsForPeripodiumCircumference(RKId);
};

void Simulation::updateNodePositionsForPeripodiumCircumference(int RKId){
	int n = PeripodiumCircumferencialNodeList.size();
	if (RKId < 3){
		for (int i=0; i<n; ++i){
			int index = PeripodiumCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodium.size();
			double RKsumV[3] = {0.0,0.0,0.0};
			double RKsumPos[3] = {0.0,0.0,0.0};
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodium[j];
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
			int index = PeripodiumCircumferencialNodeList[i];
			int nA = Nodes[index]->AssociatedNodesDueToPeripodium.size();
			double sumV[3] = {0.0,0.0,0.0};
			double RKsumV[3] = {0.0,0.0,0.0};
			double sumPos[3] = {0.0,0.0,0.0};
			for(int j = 0; j < nA; ++j){
				int index2 = Nodes[index]->AssociatedNodesDueToPeripodium[j];
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
	Elements[0]->calculateRotationAxis(u,v,rotAx);	//calculating the rotation axis that is perpendicular to both u and v
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
//cleanreferenceupdates();
void Simulation::cleanMatrixUpdateData(){
	int nElement = Elements.size();
	for (int i=0; i<nElement; ++i){
		Elements[i]->WorldToTissueRotMatUpToDate=false;
		Elements[i]->GrowthStrainsRotMatUpToDate= false;
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
				PackingForces[i][j][k]=0.0;
			}
		}
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

		PackingForces[0][j][0] = 1.0/6.0 * (PackingForces[0][j][0] + 2 * (PackingForces[1][j][0] + PackingForces[2][j][0]) + PackingForces[3][j][0]);
		PackingForces[0][j][1] = 1.0/6.0 * (PackingForces[0][j][1] + 2 * (PackingForces[1][j][1] + PackingForces[2][j][1]) + PackingForces[3][j][1]);
		PackingForces[0][j][2] = 1.0/6.0 * (PackingForces[0][j][2] + 2 * (PackingForces[1][j][2] + PackingForces[2][j][2]) + PackingForces[3][j][2]);
		//I do not need to do velocities, as I already calculate the average velocity on storage space of RK step 1, for posiiton update
	}
	n = Elements.size();
	for (int j=0;j<n;++j){
		Elements[j]->Strain = Elements[j]->RK1Strain;
	}

}
void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<timestep*dt<<" sec"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
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
			saveFileMesh.precision(3);saveFileMesh.width(10);
			saveFileMesh<<Nodes[i]->Position[0];
			saveFileMesh.precision(3);saveFileMesh.width(10);
			saveFileMesh<<Nodes[i]->Position[1];
			saveFileMesh.precision(3);saveFileMesh.width(10);
			saveFileMesh<<0.0;

		}
		else{
			int ndim = Nodes[i]->nDim;
			for (int j=0;j<ndim;++j){
				saveFileMesh.precision(3);saveFileMesh.width(10);
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
			saveFileMesh.precision(3);saveFileMesh.width(10);
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
				saveFileMesh.precision(3);saveFileMesh.width(10);
				saveFileMesh<<refPos[j][k];
			}
		}
		saveFileMesh<<endl;
	}
}

void Simulation::writeTensionCompression(){
	int n = Elements.size();
	for (int i=0;i<n;++i){
		saveFileTensionCompression.write((char*) &Elements[i]->Strain(0), sizeof Elements[i]->Strain(0));
		saveFileTensionCompression.write((char*) &Elements[i]->Strain(1), sizeof Elements[i]->Strain(1));
		saveFileTensionCompression.write((char*) &Elements[i]->Strain(2), sizeof Elements[i]->Strain(2));
		saveFileTensionCompression.write((char*) &Elements[i]->PlasticStrain(0), sizeof Elements[i]->PlasticStrain(0));
		saveFileTensionCompression.write((char*) &Elements[i]->PlasticStrain(1), sizeof Elements[i]->PlasticStrain(1));
		saveFileTensionCompression.write((char*) &Elements[i]->PlasticStrain(2), sizeof Elements[i]->PlasticStrain(2));
	}
}


void Simulation::writeForces(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		saveFileForces.write((char*) &SystemForces[0][i][0], sizeof SystemForces[0][i][0]);
		saveFileForces.write((char*) &SystemForces[0][i][1], sizeof SystemForces[0][i][1]);
		saveFileForces.write((char*) &SystemForces[0][i][2], sizeof SystemForces[0][i][2]);
	}
}

void Simulation::writeVelocities(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			saveFileVelocities.write((char*) &Nodes[i]->Velocity[j], sizeof Nodes[i]->Velocity[j]);
		}
	}
}

void Simulation::calculateGrowth(){
	int currIndexForParameters = 0;
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctionTypes[i] == 1){
			calculateGrowthUniform(currIndexForParameters);
			currIndexForParameters += 5;
		}
		else if(GrowthFunctionTypes[i] == 2){
			calculateGrowthRing(currIndexForParameters);
			currIndexForParameters += 9;
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

void Simulation::calculateGrowthUniform(int currIndex){
	float initTime = GrowthParameters[currIndex];
	float endTime = GrowthParameters[currIndex+1];
	float simTime = dt*timestep;
	if(simTime > initTime && simTime < endTime ){
		double MaxValue[3] = {GrowthParameters[currIndex+2],GrowthParameters[currIndex+3],GrowthParameters[currIndex+4]};
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			if (Elements[i]->tissueType == 0){ //grow columnar layer
				Elements[i]->updateGrowthToAdd(MaxValue);
				//This value is stored as fraction per hour, conversion is done by a time scale variable:
				float timescale = 60*60/dt;
				Elements[i]->updateGrowthRate(MaxValue[0]*timescale,MaxValue[1]*timescale,MaxValue[2]*timescale);
			}
		}
	}
}

void Simulation::calculateGrowthRing(int currIndex){
	float initTime = GrowthParameters[currIndex];
	float endTime = GrowthParameters[currIndex+1];
	float simTime = dt*timestep;
	if(simTime > initTime && simTime < endTime ){
		//The growth function is active at current time, now I will grow the elements.
		//First get the remaining data from the growth function parameters
		float centre[2] = {GrowthParameters[currIndex+2], GrowthParameters[currIndex+3]};
		float innerRadius = {GrowthParameters[currIndex+4]};
		float outerRadius = {GrowthParameters[currIndex+5]};
		float maxValue[3] = {GrowthParameters[currIndex+6],GrowthParameters[currIndex+7],GrowthParameters[currIndex+8]};
		float innerRadius2 = innerRadius*innerRadius;
		float outerRadius2 = outerRadius*outerRadius;
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			if (Elements[i]->tissueType == 0){ //grow columnar layer
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
					double growthscale[3] = {maxValue[0] * sf, maxValue[1] * sf, maxValue[2] * sf};
					//growing the shape
					Elements[i]->updateGrowthToAdd(growthscale);
					//This value is stored as fraction per hour, conversion is done by a time scale variable:
					float timescale = 60*60/dt;
					Elements[i]->updateGrowthRate(maxValue[0]*sf*timescale,maxValue[1]*sf*timescale,maxValue[2]*sf*timescale);
				}
				delete[] Elementcentre;
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

