
#include "Simulation.h"
#include "Prism.h"
#include "PrismLateral.h"
#include <string.h>

using namespace std;

Simulation::Simulation(){
	currElementId = 0;
	ModInp = new ModelInputObject();
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	timestep = 0;
	ReferencePeripodiumArea = 0.0;
	PeripodiumStrain = 0.0;
	reachedEndOfSaveFile = false;
	AddLateralNodes = false;
	AddPeripodialArea = false;
	setDefaultParameters();
};

Simulation::~Simulation(){
	//cerr<<"destructor for simulation called"<<endl;
	delete ModInp;
	int n = Nodes.size();

	for (int i=0;i<n;++i){
		delete[] SystemForces[i];
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
	E = 10.0;
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
	ApicalNodeFix[0]= true;
	ApicalNodeFix[1]= true;
	BasalNodeFix[0]= true;
	BasalNodeFix[1]= true;
	LateralNodeFix[0]= false;
	LateralNodeFix[1]= false;
	nGrowthFunctions = 0;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
	ForcesSaved = true;
	VelocitiesSaved = true;
	nLateralNodes = 0;
	AddPeripodialArea = false;
	PeripodiumElasticity = 0.0;
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
	if(AddLateralNodes == false && AddPeripodialArea){
		cerr <<"You need to add lateral nodes to have peripodium effects, correct mesh options"<<endl;
		AddPeripodialArea=false;
		return false;
	}
	if(AddLateralNodes == false && (LateralNodeFix[0] || LateralNodeFix[1] )){
		cerr <<"You need to add lateral nodes to be able to fix them in any coordinate, correct mesh options"<<endl;
		LateralNodeFix[0]=false;
		LateralNodeFix[1]=false;
		return false;
	}
	if(!LateralNodeFix[0] && LateralNodeFix[1]){
		cerr <<"Lateral nodes: There is no option to fix x and y movement while keeping z free for pinned nodes, correct mesh options"<<endl;
		LateralNodeFix[0]=true;
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
		initiateMesh(MeshType);
	}
	else if (MeshType == 2){
		initiateMesh(MeshType, Row, Column,  SideLength,  zHeight);
	}
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
	CalculateStiffnessMatrices();
	if (saveData){
		cout<<"writing the summary current simulation parameters"<<endl;
		writeSimulationSummary();
	}

	return Success;
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
	string outputFileString;
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

void Simulation::writeSimulationSummary(){
	saveFileSimulationSummary<<"TimeStep(sec):  ";
	saveFileSimulationSummary<<dt<<endl;
	saveFileSimulationSummary<<"DataSaveInterval(sec):  ";
	saveFileSimulationSummary<<dataSaveInterval*dt<<endl;
	saveFileSimulationSummary<<"ModelinputName  ";
	saveFileSimulationSummary<<ModInp->parameterFileName<<endl;

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
	//skipping the footer:
	getline(saveFileToDisplayMesh,currline);
	getline(saveFileToDisplayMesh,currline);
	cout<<"skipped footer: "<<currline<<endl;
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
	Node* tmp_nd;
	for (int i=0; i<n; ++i){
		double* pos = new double[3];
		saveFileToDisplayMesh >> pos[0];
		saveFileToDisplayMesh >> pos[1];
		saveFileToDisplayMesh >> pos[2];
		tmp_nd = new Node(i, 3, pos);
		Nodes.push_back(tmp_nd);
		delete[] pos;
	}
}

void Simulation::initiateElementsFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	for (int i=0; i<n; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			initiatePrismFromSave();
		}
		else if (shapeType == 2){
			initiateLateralPrismFromSave();
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}
	}
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
}

void Simulation::initiateLateralPrismFromSave(){
	//inserts a new prism at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave

	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		NodeIds[i] = 0;
	}
	PrismLateral* PrismPnt01;
	PrismPnt01 = new PrismLateral(NodeIds, Nodes, currElementId);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	Elements.push_back(PrismPnt01);
	currElementId++;
	delete[] NodeIds;
}

void Simulation::reInitiateSystemForces(int oldSize){
	//deleting the old system forces:
	for (int i=0;i<oldSize;++i){
		delete[] SystemForces[i];
	}
	delete[] SystemForces;

	//reinitiating with the new size:
	const int n = Nodes.size();
	SystemForces = new double*[n];
	for (int i=0;i<n;++i){
		SystemForces[i] = new double[3];
		SystemForces[i][0]=0.0;
		SystemForces[i][1]=0.0;
		SystemForces[i][2]=0.0;
	}
}

void Simulation::updateForcesFromSave(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		saveFileToDisplayForce.read((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][0], sizeof SystemForces[i][1]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][0], sizeof SystemForces[i][2]);
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
	saveFileToDisplayTenComp.read((char*) &PeripodiumStrain, sizeof PeripodiumStrain);
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

void Simulation::updateOneStepFromSave(){
	cout<<"Updating step from save:"<<endl;
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
	cout<<"skipped header: "<<currline<<endl;
	updateNodeNumberFromSave();
	updateNodePositionsFromSave();
	updateElementStatesFromSave();
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		Elements[i]->updatePositions(Nodes);
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

	//skipping the footer:
	getline(saveFileToDisplayMesh,currline);
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
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
			tmp_nd = new Node(i, 3, pos);
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
		if (shapeType == currShapeType){
			//cout<<"shape type is correct, moving on to update"<<endl;
			//The current shape on the list, and the shape I am reading are of same type, I can read it:
			Elements[i]->updateShapeFromSave(saveFileToDisplayMesh);
		}
		else{
			//cout<<"shape type is wrong"<<endl;
			//the current element is a different type, I need to insert generate a new element, and insert it here:
			if (shapeType == 1){

				//the new shape is a prism, I will add it now;
				initiatePrismFromSaveForUpdate(i);
				removeElementFromEndOfList();
				currElementNumber = Elements.size();
			}
		}
	}
	while(n>currElementNumber){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			int i = Elements.size()-1;
			//this will initiate a prism at the current point in the Elemetns vector
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
	}
}

void Simulation::initiateMesh(int MeshType){
	if (MeshType == 1 ){
		initiateSinglePrismNodes();
		initiateSinglePrismElement();
	}
	else if (MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
	}
	else if ( MeshType == 3 || MeshType ==4 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
	}
}

void Simulation::initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
	}
	if ( MeshType == 2){
		initiateNodesByRowAndColumn(Row,Column,SideLength,zHeight);
		initiateElementsByRowAndColumn(Row,Column);
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
	}
	else if ( MeshType ==4 ){
			cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
	}
}

void Simulation::initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight ){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
	}
	else if ( MeshType == 3 ){
		//this will be inputting circumference of the tissue, and the sidelength and z-height
		//generate mesh by triangulation
	}
	else if ( MeshType ==4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		//this will be reading full mesh data
		//read mesh data from file
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
	}
}

void Simulation::initiateMesh(int MeshType, string inputtype ){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
	}
	else if ( MeshType ==4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		//this will be reading full mesh data
		//read mesh data from file
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
	}
}

void Simulation::CalculateStiffnessMatrices(){
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		cout<<" setting up element :  "<<i<<" of "<<n<<endl;
		Elements[i]->calculateReferenceStiffnessMatrix();
		Elements[i]->setStiffnessMatrixBuffers();
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

void Simulation::zeroForcesOnNode(int i){
	double ForceBalance[3];
	ForceBalance[0] = SystemForces[i][0];
	ForceBalance[1] = SystemForces[i][1];
	ForceBalance[2] = SystemForces[i][2];
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		SystemForces[i][0]-=ForceBalance[0];
		SystemForces[i][1]-=ForceBalance[1];
		SystemForces[i][2]-=ForceBalance[2];
	}
}

void Simulation::initiateSystemForces(){
	const int n = Nodes.size();
	SystemForces = new double*[n];
	for (int i=0;i<n;++i){
		SystemForces[i] = new double[3];
		SystemForces[i][0]=0.0;
		SystemForces[i][1]=0.0;
		SystemForces[i][2]=0.0;
		//cout<<"systemforces[i][j]: "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<endl;
	}
}

void Simulation::initiateSinglePrismNodes(){
	double *pos = new double[3];
	Node* tmp_nd;
	pos[0]=0;pos[1]=1;pos[2]=0;
	tmp_nd = new Node(0, 3, pos);
	Nodes.push_back(tmp_nd);
	pos[0]=1;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(1, 3, pos);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(2, 3, pos);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=1;pos[2]=1;
	tmp_nd = new Node(3, 3, pos);
	Nodes.push_back(tmp_nd);
	pos[0]=1;pos[1]=0;pos[2]=1;
	tmp_nd = new Node(4, 3, pos);
	Nodes.push_back(tmp_nd);
	pos[0]=0;pos[1]=0;pos[2]=1;
	tmp_nd = new Node(5, 3, pos);
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
	//Adding the lower level of nodes:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = 0.0;
		tmp_nd = new Node(i, 3, pos);
		Nodes.push_back(tmp_nd);
	}
	//Adding the upper level:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = zHeight;
		tmp_nd = new Node(n+i, 3, pos);
		Nodes.push_back(tmp_nd);
	}

	if (BasalNodeFix[0] || BasalNodeFix[1] || ApicalNodeFix[0] || ApicalNodeFix[1] ){
		fixApicalBasalNodes(NodesToFix);
	}
	if (AddLateralNodes){
		cout<<"generating lateral node list"<<endl;
        GenerateLateralNodeList(NodesToFix,toprowcounter);
        GenerateLateralNodes();
	}
	if (LateralNodeFix[0] || LateralNodeFix[1]){
		fixLateralNodes();
	}
	if (AddPeripodialArea){
		ReferencePeripodiumArea = calculatePeripodiumArea();
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

void Simulation::fixLateralNodes(){
	for (int i=0; i<nLateralNodes; ++i){
		if(LateralNodeFix[0]){
			//The nodes on the circumference on the basal side (bottom)
			//have their z position fixed.
			fixZ(PeripodiumAnchorNodeList[i]);
		}
		if(LateralNodeFix[1]){
			//The nodes on the circumference on the basal side (bottom)
			//have their x & y positions fixed.
			fixAllD(PeripodiumAnchorNodeList[i]);
		}
	}
}
void Simulation::GenerateLateralNodeList(vector<int> &NodesToFix, int nLastRow){
	//Now I have a list of nodes to fix, written in a specific order.
	//I would like reorder these, to generate a continuous strip around the tissue
	cout<<"generating lateral node list - inside"<<endl;
	vector <int> list1, list2, list3, list4;
    int n = NodesToFix.size();
    int i=0;
    list1.push_back(NodesToFix[i]);
    i++;
    list2.push_back(NodesToFix[i]);
    i++;
    cout<<"generating lateral node list = before while loop"<<endl;
    while(i<n-2*nLastRow){
    	list1.push_back(NodesToFix[i]);
    	list2.push_back(NodesToFix[i+1]);
    	list3.push_back(NodesToFix[i+2]);
    	list4.push_back(NodesToFix[i+3]);
    	i += 4;
    }
    cout<<"generating lateral node list - after while loop"<<endl;
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
    	LateralNodeList.push_back(*it);
    }
    //inverse of list2:
    for(vector <int>::iterator it = list2.end()-1; it>=list2.begin(); --it){
    	LateralNodeList.push_back(*it);
    }
    //list4:
    for(vector <int>::iterator it = list4.begin(); it<list4.end(); ++it){
    	LateralNodeList.push_back(*it);
    }
    //inverse of list3:
    for(vector <int>::iterator it = list3.end()-1; it>=list3.begin(); --it){
    	LateralNodeList.push_back(*it);
    }
    cout<<"merged lists"<<endl;
    //for (i=0; i<LateralNodeList.size(); ++i){
    //	cout<<"i: "<<i<<" LateralNodeList[i]: "<<LateralNodeList[i]<<endl;
    //}
}

void Simulation::GenerateLateralNodes(){
	//I have the list of nodes, that are the apical base for the nodes I am going to add
	//as lateral tips. The list of nodes starts form the first node, rotating clockwise
	//the position of each added node will be pointing out from the tissue
	int n = LateralNodeList.size();
	double vec1[3], vec2[3];
	double* dir = new double[3];
	//calculating the vector pointing out, the first connection is between the beginning
	// and the end of the list:
    for(int i=0; i < n; ++i){
		int index0 = LateralNodeList[i];
		int index1, index2;
		if (i==0){
			index1 = LateralNodeList[i+1];
			index2 = LateralNodeList[n-1];
		}
		else if(i==n-1){
			index1 = LateralNodeList[0];
			index2 = LateralNodeList[i-1];
		}
		else{
			index1 = LateralNodeList[i+1];
			index2 = LateralNodeList[i-1];
		}
		double norm2 = 0.0;
		for (int j=0;j<3; ++j){
			vec1[j]=Nodes[index1]->Position[j] - Nodes[index0]->Position[j];
			vec2[j]=Nodes[index2]->Position[j] - Nodes[index0]->Position[j];
			dir[j] = -vec1[j] - vec2[j];
			norm2 +=dir[j]*dir[j];
		}
		if (norm2 <1E-8){
			//the vectors form a straight line, I will need to rotate one to get the perpendicular:
			//rotating -90 degrees around z:
			dir[0]=(1.0)*vec2[1];
			dir[1]=(-1.0)*vec2[0];
			dir[2]=vec2[2];
			norm2 = 0;
			for (int j=0;j<3; ++j){
				norm2 += vec1[j]*vec1[j];
			}
		}
		//normalise vec1:
		double norm = pow(norm2,0.5);
		for (int j=0; j<3; ++j){
			dir[j] /= norm;
		}
		//now I will get the point that is moving in the dir direction form my current node, and also migrating apically,
		//half the length of the tissue:
		double sqrt3 = 1.7321;
		float h = sqrt3/2*SideLength;
		dir[0] = dir[0]*h + Nodes[index0]->Position[0];
		dir[1] = dir[1]*h + Nodes[index0]->Position[1];
		dir[2] = dir[2]*h + Nodes[index0]->Position[2] + zHeight/2.0;
		Node* tmp_nd;
		int nNodes = Nodes.size();
		tmp_nd = new Node(nNodes, 3, dir);
		Nodes.push_back(tmp_nd);
		PeripodiumAnchorNodeList.push_back(nNodes);
		//fixZ(nNodes);
		//fixAllD(nNodes);
		nLateralNodes++;
    }
    delete[] dir;
    cout<<"finished generating lateral nodes"<<endl;
}

void Simulation::initiateElementsByRowAndColumn(int Row, int Column){
	int xinit1 = 0;
	int xinit2 = xinit1+Row+1;
	int xinit3 = 0;
	int xinit4 = xinit1+2*(Row+1)-1;
	int n = (Nodes.size() - nLateralNodes) /2.0;
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
	//initialising the lateral elements of the circumference:
	cout<<"initiating the lateral elements"<<endl;
	for (int i =nLateralNodes-1 ;i>=0; --i){
		int secondIndex = i-1;
		if(i==0){
			secondIndex = nLateralNodes-1;
		}
		int* NodeIds;
		NodeIds = new int[6];
		NodeIds[0] = LateralNodeList[i];
		NodeIds[1] = LateralNodeList[i]+n;
		NodeIds[2] = 2*n+i;
		NodeIds[3] = LateralNodeList[secondIndex];
		NodeIds[4] = LateralNodeList[secondIndex]+n;
		NodeIds[5] = 2*n+secondIndex;
		PrismLateral* PrismPnt01;
		PrismPnt01 = new PrismLateral(NodeIds, Nodes, currElementId);
		//Prism* PrismPnt01;
		//PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
		Elements.push_back(PrismPnt01);
		currElementId++;
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
		Elements[i]->setElasticProperties(E*(1 + noise1/100.0),poisson*(1 + noise2/100));
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		Elements[i]->setViscosity(ApicalVisc*noise3, BasalVisc*noise3, Nodes);
	}
}

void Simulation::runOneStep(){
	int displayfreq = 60/dt;
	/*if(timestep==0){
		for(int i=0;i<Nodes.size();++i){
			Nodes[i]->Position[0] *=2.0;
			Nodes[i]->Position[1] *=2.0;
			Nodes[i]->Position[2] *=2.0;
		}
		double R[3][3];
		double Rx[3][3] = {{1,0,0},{0,0,-1},{0,1,0}};
		double Ry[3][3] = {{0,0,1},{0,1,0},{-1,0,0}};
		double Rz[3][3] = {{0,-1,0},{1,0,0},{0,0,1}};
		for (int j =0; j<3;++j){
			for(int k=0; k<3;++k){
				R[j][k] = Rz[j][k];
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
			Elements[i]->updatePositions(Nodes);
		}
	}*/
	if (timestep%displayfreq == 0){
		cout<<"time : "<<timestep * dt<<endl;
		float dx = Nodes[Elements[0]->NodeIds[0]]->Position[0] - Nodes[Elements[4]->NodeIds[1]]->Position[0];
		float dy = Nodes[Elements[0]->NodeIds[0]]->Position[1] - Nodes[Elements[4]->NodeIds[1]]->Position[1];
		float dz = Nodes[Elements[0]->NodeIds[0]]->Position[2] - Nodes[Elements[4]->NodeIds[1]]->Position[2];
		double d2 = dx*dx + dy*dy + dz*dz;
		float d = pow(d2,0.5);
		cout<<"distance: "<<d<<endl;
	}
	cleanreferenceupdates();
	cleanGrowthData();
	resetForces();

	int nElement = Elements.size();
	for (int i=0; i<nElement; ++i){
		//aligning the centres:
		Elements[i]->alignReference();
		if (Elements[i]->updatedReference){
			Elements[i]->calculateReferenceStiffnessMatrix();
		}
	}
	if(nGrowthFunctions>0){
		growSystem();
	}
	if(nShapeChangeFunctions>0){
		changeCellShapesInSystem();
	}
	for (int i=0; i<nElement; ++i){
		//cout<<"calculating element: "<<i<<endl;
		//updating the plastic strains with the calculations on growth and shape change:
		Elements[i]->calculatePlasticStrain();
		//calculating the forces:
		Elements[i]->calculateForces(1,SystemForces, Nodes);
	}

	if(AddPeripodialArea){
		addPeripodiumResistance();
	}
	//Fixing Node 0 in space:
	//zeroForcesOnNode(0);

	//Update Node positions:
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			Nodes[i]->Velocity[j] = SystemForces[i][j]/Nodes[i]->Viscosity;
			Nodes[i]->Position[j] += Nodes[i]->Velocity[j]*dt;
		}
	}
	for (int i=0;i<Elements.size();++i){
		Elements[i]->updatePositions(Nodes);
	}
	if (saveData && timestep % dataSaveInterval == 0){
		saveStep();
	}
	timestep++;

}

void Simulation::cleanGrowthData(){
	int nElement = Elements.size();
	for (int i=0; i<nElement; ++i){
		Elements[i]->resetCurrStepGrowthData();
		Elements[i]->resetCurrStepShapeChangeData();
	}
}

void Simulation::resetForces(){
	int n = Nodes.size();
	int dim = 3;
	for (int i=0;i<n;++i){
		for (int j=0;j<dim;++j){
			SystemForces[i][j]=0.0;
		}
	}
}

void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<timestep*dt<<" sec"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
	writePeripodiumTensionCompression();
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
			saveFileMesh<<0.0<<endl;

		}
		else{
			int ndim = Nodes[i]->nDim;
			for (int j=0;j<ndim;++j){
				saveFileMesh.precision(3);saveFileMesh.width(10);
				saveFileMesh<<Nodes[i]->Position[j];
			}
			saveFileMesh<<endl;
		}
	}
}

void Simulation::writeElements(){
	int n= Elements.size();
	saveFileMesh<<n<<endl;
	for (int i = 0; i<n; ++i){
		int shapetype = Elements[i]->getShapeType();
		saveFileMesh.width(4);
		saveFileMesh<<shapetype;
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

void Simulation::writePeripodiumTensionCompression(){
	saveFileTensionCompression.write((char*) &PeripodiumStrain, sizeof PeripodiumStrain);
}

void Simulation::writeForces(){
	int n = Nodes.size();
	for (int i=0;i<n;++i){
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][1]);
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][2]);
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

void Simulation::growSystem(){
	int currIndexForParameters = 0;
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctionTypes[i] == 1){
			growSystemUniform(currIndexForParameters);
			currIndexForParameters += 5;
		}
		else if(GrowthFunctionTypes[i] == 2){
			growSystemRing(currIndexForParameters);
			currIndexForParameters += 9;
		}
	}
}

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


void Simulation::cleanreferenceupdates(){
	int nElement = Elements.size();
	for (int i=0; i<nElement; ++i){
		Elements[i]->updatedReference = false;
	}
}

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

void Simulation::growSystemUniform(int currIndex){
	float initTime = GrowthParameters[currIndex];
	float endTime = GrowthParameters[currIndex+1];
	float simTime = dt*timestep;
	if(simTime > initTime && simTime < endTime ){
		double MaxValue[3] = {GrowthParameters[currIndex+2],GrowthParameters[currIndex+3],GrowthParameters[currIndex+4]};
		int  n = Elements.size();
		for ( int i = 0; i < n; ++i ){
			Elements[i]->updateGrowthToAdd(MaxValue);
			//This value is stored as fraction per hour, conversion is done by a time scale variable:
			float timescale = 60*60/dt;
			Elements[i]->updateGrowthRate(MaxValue[0]*timescale,MaxValue[1]*timescale,MaxValue[2]*timescale);
		}
	}
}

void Simulation::growSystemRing(int currIndex){
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

double Simulation::calculatePeripodiumArea(){
	//approximating the area of peripodial membrane from lateral nodes
	double sum =0;
	for (int i=0; i< nLateralNodes; ++i){
		int index1 = PeripodiumAnchorNodeList[i];
		int index0,index2;
		if (index1 ==0 ){
			index0 = PeripodiumAnchorNodeList[nLateralNodes-1];
			index2 = PeripodiumAnchorNodeList[i+1];

		}
		else if(index1 == nLateralNodes-1){
			index0 = PeripodiumAnchorNodeList[i-1];
			index2 = PeripodiumAnchorNodeList[0];
		}
		else{
			index0 = PeripodiumAnchorNodeList[i-1];
			index2 = PeripodiumAnchorNodeList[i+1];
		}
		double x = Nodes[index1]->Position[0];
		double y0 = Nodes[index0]->Position[1];
		double y2 = Nodes[index2]->Position[1];
		sum += x * (y2 - y0);
	}
	sum *=0.5;
	return sum;
};

double Simulation::calculatePeripodiumResistanceForce(){
	double CurrentPeripodiumArea = calculatePeripodiumArea();
	PeripodiumStrain = (CurrentPeripodiumArea - ReferencePeripodiumArea) / ReferencePeripodiumArea;
	double ForceMagnitude = PeripodiumStrain * PeripodiumElasticity;
	//cout<<"referenceperipodium: "<<ReferencePeripodiumArea<<" Current Area: "<< CurrentPeripodiumArea<<" PeripodiumStrain: "<< PeripodiumStrain<<" ForceMagnitude: "<<ForceMagnitude<<endl;
	return ForceMagnitude;
}

void Simulation::addPeripodiumResistance(){
	double ForceMagnitude = calculatePeripodiumResistanceForce();
	if (ForceMagnitude <1E-6 &&  ForceMagnitude>-1E-6){
		return;
	}
	ForceMagnitude /= nLateralNodes;
	calculateSystemCentre();
	for (int i=0; i< nLateralNodes; ++i){
		double x = Nodes[PeripodiumAnchorNodeList[i]]->Position[0];
		double y = Nodes[PeripodiumAnchorNodeList[i]]->Position[1];
		double z = Nodes[PeripodiumAnchorNodeList[i]]->Position[2];
		x = SystemCentre[0]-x;
		y = SystemCentre[1]-y;
		z = SystemCentre[2]-z;
		double mag = x*x +y*y+z*z;
		if (mag > 1E-6){
			mag = pow(mag,0.5);
			x /= mag;
			y /= mag;
			z /= mag;
		}
		SystemForces[i][0] += x * ForceMagnitude;
		SystemForces[i][1] += y * ForceMagnitude;
		SystemForces[i][2] += z * ForceMagnitude;
		//cout<<"i: "<<i<<" (x,y,z): "<<x<<" , "<<y<<" , "<<z<<", ForceMag(x,y,z): "<<x * ForceMagnitude<<" , "<<y * ForceMagnitude<<" , "<<z * ForceMagnitude<<endl;
	}
}
