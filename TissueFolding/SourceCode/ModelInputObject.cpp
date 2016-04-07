/*
 * ModelInputObject.cpp
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#include "ModelInputObject.h"
#include "Simulation.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"

#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

ModelInputObject::ModelInputObject(){
	meshFileName = "NA";
};

ModelInputObject::~ModelInputObject(){
	//cout<<"called the destructor for ModelInputObject class"<<endl;
	//cout<<"finalised the destructor for ModelInputObject class"<<endl;
};

bool ModelInputObject::readParameters(){
	/**
	 *  This function will read all available model inputs from the file ModelInputObject#parameterFileName. \n
	 *  It will start by opening the model input file, after each attempt to open a file, there will be a health check to ensure the file
	 *  could be opened. In case there are issues with the file (most common one being the file is not opened due to a path error),
	 *  the function will throw an error with corresponding explanatory error message, and quit the simulation.
	 */
	bool Success = true;
	ifstream parametersFile;
	parametersFile.open(parameterFileName, ifstream::in);
	Success = checkFileStatus(parametersFile,parameterFileName);
	if (!Success){
		return Success;
	}
	/**
	 *  After successfully opening the input file, the function will read it until it reaches to the end of the file.
	 *  This will involve a series of private functions, which are thoroughly documented in source code, while
	 *  the processed documentation may or may not be available in this user interface documentation structure.
	 */
	while (!parametersFile.eof()){
		string currline;
		getline(parametersFile,currline);
		if (!parametersFile.eof()){
			if (currline.empty()){
				//skipping empty line
				continue;
			}
			istringstream currSStrem(currline);
			string currParameterHeader;
			currSStrem >> currParameterHeader;
			if (currParameterHeader == ""){
				//just in case an empty line stores two consecutive line ends (/n/n),
				//and thus, currline is not technically empty, while the line is still empty
				continue;
			}
			/**
			 * Depending on the header the function it encounters it will read:
			 */
			if(currParameterHeader == "InputMeshParameters:"){
				/**
				 * Mesh geometry related parameters through the private function ModelInputObject#readMeshParameters
				 */
				Success  = readMeshParameters(parametersFile);
			}
			else if(currParameterHeader == "PeripodialMembraneParameters:"){
				/**
				 * Peripodial membrane structure related parameters through the private function ModelInputObject#readPeripodialMembraneParameters
				 */
				Success  = readPeripodialMembraneParameters(parametersFile);
			}
			else if(currParameterHeader == "LinkerZoneParameters:"){
				/**
				 * Linker zone physical properties through the private function ModelInputObject#readLinkerZoneParameters
				 */
				Success  = readLinkerZoneParameters(parametersFile);
			}
			else if(currParameterHeader == "NodeFixingOptions:"){
				/**
				 * Inputs relating to fixing the nodes of the tissue through the private function ModelInputObject#readNodeFixingParameters
				 */
				Success  = readNodeFixingParameters(parametersFile);
			}
			else if(currParameterHeader == "Manipulations:"){
				/**
				 * Inputs relating to manual manipulations to tissue after the mesh is read in
				 */
				Success  = readManupulationParamters(parametersFile);
			}
			else if(currParameterHeader == "TimeParameters:"){
				/**
				 * Simulation time and time step related parameters through the private function ModelInputObject#readTimeParameters
				 */
				Success  = readTimeParameters(parametersFile);
			}
			else if(currParameterHeader == "PysicalProperties:"){
				/**
				 * Physical parameters of the tissue  through the private function ModelInputObject#readPysicalProperties
				 */
				Success  = readPysicalProperties(parametersFile);
			}
			else if (currParameterHeader == "SaveOptions:"){
				/**
				 * Save options of the simulation through the private function ModelInputObject#readSaveOptions
				 */
				Success  = readSaveOptions(parametersFile);
			}
			else if (currParameterHeader == "GrowthOptions:"){
				/**
				 * Growth functions and related parameters of the simulation through the private function ModelInputObject#readGrowthOptions
				 */
				Success  = readGrowthOptions(parametersFile);
			}
			else if(currParameterHeader == "ShapeChangeOptions:"){
				/**
				 * Shape change functions and related parameters of the simulation through the private function ModelInputObject#readShapeChangeOptions
				 */
				Success  = readShapeChangeOptions(parametersFile);
			}
			else if (currParameterHeader == "PlasticDeformationOptions:"){
				/**
				* Plastic deformation options through the private function ModelInputObject#readPlasticDeformationOptions
				*/
				Success  = readPlasticDeformationOptions(parametersFile);
			}
			else if (currParameterHeader == "MyosinOptions:"){
				/**
				 * Myosin concentrations and related parameters of the simulation through the private function ModelInputObject#readMyosinOptions
				 */
				Success  = readMyosinOptions(parametersFile);
			}
			else if(currParameterHeader == "Stretcher:"){
				/**
				 * Stretcher experimental setup parameters of the simulation through the private function ModelInputObject#readStretcherSetup
				 */
				Success  = readStretcherSetup(parametersFile);
			}
			else if(currParameterHeader == "Pipette_Aspiration:"){
				/**
				 * Pipette aspiration experimental setup parameters of the simulation through the private function ModelInputObject#readPipetteSetup
				 */
				Success  = readPipetteSetup(parametersFile);
			}
			else {
				/**
				 * In the case that the function, or any of the above listed parameter reading functions, encounters an
				 * unexpected line in the model input file, it will throw an error with a corresponding explanatory message,
				 * and quit the simulation.
				 */
				cerr<<"Unidentified parameter input line: "<<endl;
				cerr<<"		"<<currParameterHeader<<endl;
				return false;
			}
			if (!Success){
				return false;
			}
		}
	}
	return Success;
}

bool ModelInputObject::checkFileStatus(ifstream& file, string fileName){
	if (!file.is_open()) {
		cout<<"Cannot open parameter input file, "<<fileName<<endl;
		return false;
	}
	if (!file.good()) {
		cout<<"File does not exist, "<<fileName<<endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readGrowthOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	int n;
	cout<<"entered read growth options, current header: "<<currHeader<<endl;
	if(currHeader == "NumberofGrowthFunctions(int):"){
		file >> n;
		Sim->nGrowthFunctions = n;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: NumberofGrowthFunctions(int):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "GridGrowthsPinnedOnInitialMesh(bool):"){
		bool tmpBool;
		file >> tmpBool;
		Sim->GridGrowthsPinnedOnInitialMesh = tmpBool;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: GridGrowthsPinnedOnInitialMesh(bool):" <<endl;
		return false;
	}
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		cout<<"inside the loop, read growth options, current header: "<<currHeader<<endl;
		if(currHeader == "GrowthFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
			cerr<<"Error in reading growth type, curr string: "<<currHeader<<", should have been: GrowthFunctionType(int-seeDocumentation):" <<endl;
			return false;
		}
		bool Success = true;
		if (type == 1){
			Success = readGrowthType1(file);
		}
		else if (type == 2){
			Success = readGrowthType2(file);
		}
		else if (type == 3){
			Success = readGrowthType3(file);
		}
		else{
			cerr<<"Error in reading growth type, please enter a valid type: {1, 2, 3}, current type: "<<type<<endl;
			return false;
		}
		if (!Success){
			return false;
		}
	}
	cout<<"finalised reading growths"<<endl;
	return true;
}

bool ModelInputObject::readGrowthType1(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 1, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	float DVRate;
	float APRate;
	float ABRate;
	float angle;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
		//Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> DVRate;
		DVRate *= timeMultiplier;
		file >> APRate;
		APRate *= timeMultiplier;
		file >> ABRate;
		ABRate *= timeMultiplier;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Angle(degrees):"){
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degreesPerHour):" <<endl;
		return false;
	}

	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 1
	GSBp = new UniformGrowthFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, DVRate, APRate,  ABRate, angle);
	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readGrowthType2(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 2, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	float CentreX,CentreY;
	float innerR, outerR;
	float DVRate;
	float APRate;
	float ABRate;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre:"){
		file >> CentreX;
		file >> CentreY;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Centre:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InnerRadius:"){
		file >> innerR;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "OuterRadius:"){
		file >> outerR;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> DVRate;
		DVRate *= timeMultiplier;
		file >> APRate;
		APRate *= timeMultiplier;
		file >> ABRate;
		ABRate *= timeMultiplier;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-DV,AP,AB):" <<endl;
		return false;
	}
	double angle;
	file >> currHeader;
	if(currHeader == "Angle(degrees):"){
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degreesPerHour):" <<endl;
		return false;
	}
	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 2
	GSBp = new RingGrowthFunction(Id, 2, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, CentreX, CentreY, innerR,  outerR, DVRate, APRate,  ABRate, angle);
	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readGrowthType3(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 3, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	int gridX, gridY;
	double*** GrowthMatrix;
	double**  AngleMatrix;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Filename(full-path):"){
		string filepath;
		file >> filepath;
		cerr<<" filename is: "<<filepath<<endl;
		const char* name_growthRates = filepath.c_str();
		ifstream GrowthRateFile;
		GrowthRateFile.open(name_growthRates, ifstream::in);
		if (!(GrowthRateFile.good() && GrowthRateFile.is_open())){
			cerr<<"could not open growth rate file file: "<<name_growthRates<<endl;
			return false;
		}
		//adding the indice of the growth matrix
		cout<<"reading from growth file"<<endl;
		GrowthRateFile >> gridX;
		GrowthRateFile >> gridY;
		float rate;
        //double timeMultiplier = Sim->dt / 3600.0;
        cout<<"initiating growth and angle matrices"<<endl;
		GrowthMatrix = new double**[(const int) gridX];
		AngleMatrix = new double*[(const int) gridX];
		for (int i=0; i<gridX; ++i){
			GrowthMatrix[i] = new double*[(const int) gridY];
			AngleMatrix[i] = new double[(const int) gridY];
			for (int j=0; j<gridY; ++j){
				GrowthMatrix[i][j] = new double[3];
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] = 0.0;
				}
				AngleMatrix[i][j] = 0.0;
			}
		}
		double timeMultiplier = 1.0 / 3600.0; // converting rate per hour to rate per second
		cout<<"reading growth matrix"<<endl;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				for (int k=0; k<3; ++k){
					//cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					GrowthRateFile >> rate;
					//cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					GrowthMatrix[i][j][k] = rate;
					//cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<endl;
				}
			}
		}
		double angle;
		cout<<"reading angle matrix"<<endl;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				GrowthRateFile >> angle;
				AngleMatrix[i][j] = angle; // angles in degrees!
			}
		}
		GrowthRateFile.close();
		for (int i=0; i<gridX; ++i){
			for (int j=0; j<gridY; ++j){
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] *= timeMultiplier;
				}
			}
		}
		//display:
		cout<<"growth matrix: "<<endl;
		for (int i=0; i<gridX; ++i){
			for (int j=0; j<gridY; ++j){
				for (int k=0; k<3; ++k){
					cout<<GrowthMatrix[i][j][k]<<" ";
				}
				cout<<"	";
			}
			cout<<endl;
		}cout<<endl;
		cout<<"angle matrix: "<<endl;
		for (int i=0; i<gridX; ++i){
			for (int j=0; j<gridY; ++j){
				cout<<AngleMatrix[i][j]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Filename(full-path):" <<endl;
		return false;
	}
	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 3
	GSBp = new GridBasedGrowthFunction(Id, 3, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, gridX, gridY, GrowthMatrix, AngleMatrix);
	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readMeshParameters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "MeshInputMode(int-seeDocumentation):"){
		file >> Sim->MeshType;
	}
	else{
		cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: MeshInputMode(seeDocumentation):" <<endl;
		return false;
	}
	if ( Sim->MeshType == 2 ){
		bool Success  = readMeshType2(file);
		if (!Success){
			return false;
		}
	}
	else if ( Sim->MeshType == 4){
		bool Success  = readMeshType4(file);
		if (!Success){
			return false;
		}
	}
	file >> currHeader;
	if(currHeader == "symmetricInX(bool):"){
		file >> Sim->symmetricX;
	}
	else{
		cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: symmetricInX(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "symmetricInY(bool):"){
		file >> Sim->symmetricY;
	}
	else{
		cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: symmetricInY(bool):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readNodeFixingParameters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "FixWithHighViscosity(bool):"){
		file >>Sim->fixWithViscosity;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: FixWithHighViscosity(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixVis(x,y,z):"){
		file >>Sim->fixingViscosity[0];
		file >>Sim->fixingViscosity[1];
		file >>Sim->fixingViscosity[2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: FixVis(x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicalSurfaceFix(bool-x,y,z):"){
		file >>Sim->ApicalNodeFix[0];
		file >>Sim->ApicalNodeFix[1];
		file >>Sim->ApicalNodeFix[2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: ApicalSurfaceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalSurfaceFix(bool-x,y,z):"){
		file >>Sim->BasalNodeFix[0];
		file >>Sim->BasalNodeFix[1];
		file >>Sim->BasalNodeFix[2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: BasalSurfaceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicalCircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[0][0];
		file >>Sim->CircumferentialNodeFix[0][1];
		file >>Sim->CircumferentialNodeFix[0][2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: ApicalCircumferenceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalCircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[1][0];
		file >>Sim->CircumferentialNodeFix[1][1];
		file >>Sim->CircumferentialNodeFix[1][2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: BasalCircumferenceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerApicalCircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[2][0];
		file >>Sim->CircumferentialNodeFix[2][1];
		file >>Sim->CircumferentialNodeFix[2][2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerApicalCircumferenceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerBasalCircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[3][0];
		file >>Sim->CircumferentialNodeFix[3][1];
		file >>Sim->CircumferentialNodeFix[3][2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerBasalCircumferenceFix(bool-x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "CircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[4][0];
		file >>Sim->CircumferentialNodeFix[4][1];
		file >>Sim->CircumferentialNodeFix[4][2];
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: CircumferenceFix(bool-x,y,z):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readManupulationParamters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "AddCurvature(bool):"){
		file >>Sim->addCurvatureToTissue;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: AddCurvature(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "CurvatureDepthAtCentre(double-microns):"){
		file >>Sim->tissueCurvatureDepth;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: CurvatureDepthAtCentre(double-microns):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "AddSoftPeriphery(bool):"){
		file >>Sim->softPeriphery;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: AddSoftPeriphery(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SoftPeripheryRange(double-microns):"){
		file >>Sim->softDebth;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: SoftPeripheryRange(double-microns):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SoftnessFraction(double-fraction):"){
		file >>Sim->softnessFraction;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: SoftnessFraction(double-fraction):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToApicalSurface(bool):"){
		file >>Sim->softPeripheryBooleans[0];
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToApicalSurface(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalSurface(bool):"){
		file >>Sim->softPeripheryBooleans[1];
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToBasalSurface(bool):" <<endl;
		return false;
	}	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
		file >>Sim->softPeripheryBooleans[2];
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
		file >>Sim->softPeripheryBooleans[3];
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	//Reading Random Force Parameters
	file >> currHeader;
	if(currHeader == "AddRandomForce(bool):"){
		file >>Sim->addingRandomForces;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: AddRandomForce(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "RandomForceMean(double):"){
		file >>Sim->randomForceMean;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: RandomForceMean(double):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "RandomForceVar(double):"){
		file >>Sim->randomForceVar;
	}
	else{
		cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: RandomForceVar(double):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readMeshType4(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "MeshFile(full-path):"){
		file >> Sim->inputMeshFileName;
	}
	else{
		cerr<<"Error in reading mesh path, curr string: "<<currHeader<<", should have been: MeshFile(full-path):" <<endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readMeshType2(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "MeshRow(int):"){
		file >> Sim->Row;
	}
	else{
		cerr<<"Error in reading mesh row number, curr string: "<<currHeader<<", should have been: MeshRow(int):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "MeshColumn(int):"){
		file >> Sim->Column;
	}
	else{
		cerr<<"Error in reading mesh column number, curr string: "<<currHeader<<", should have been: MeshColumn(int):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SideLength:"){
		file >> Sim->SideLength;
	}
	else{
		cerr<<"Error in reading side length, curr string: "<<currHeader<<", should have been: SideLength:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "zHeight:"){
		file >> Sim->zHeight;
	}
	else{
		cerr<<"Error in reading z height, curr string: "<<currHeader<<", should have been: zHeight:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicalCircumferenceFixZ:"){
		file >> Sim->ApicalNodeFix[0];
	}
	else{
		cerr<<"Error in reading nodes to fix, curr string: "<<currHeader<<", should have been: ApicalCircumferenceFixZ:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicalCircumferenceFixXY:"){
		file >> Sim->ApicalNodeFix[1];
	}
	else{
		cerr<<"Error in reading nodes to fix, curr string: "<<currHeader<<", should have been: ApicalCircumferenceFixXY:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalCircumferenceFixZ:"){
		file >> Sim->BasalNodeFix[0];
	}
	else{
		cerr<<"Error in reading nodes to fix, curr string: "<<currHeader<<", should have been: BasalCircumferenceFixZ:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalCircumferenceFixXY:"){
		file >> Sim->BasalNodeFix[1];
	}
	else{
		cerr<<"Error in reading nodes to fix, curr string: "<<currHeader<<", should have been: BasalCircumferenceFixXY:" <<endl;
		return false;
	}
	//checking consistency:
	if (Sim->Column>Sim->Row-2){
		Sim->Column = Sim->Row-2;
		Sim->outputFile<<"Too few rows vs. columns, column count cannot be higher than Row-2"<<endl;
		Sim->outputFile<<"Updated to have a mesh: Row: "<<Sim->Row<<" Column: "<<Sim->Column<<endl;
	}
	float aspectratio = Sim->zHeight/Sim->SideLength;
	if ( aspectratio > 10 || aspectratio < 0.01 ){
		Sim->outputFile<<"Warning: The aspect ratio of the shapes are too high or low (aspectratio (z/side): "<<aspectratio<<endl;
	}
	return true;
}

bool ModelInputObject::readPeripodialMembraneParameters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "AddPeripodialMembrane:"){
		file >>Sim->AddPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: AddPeripodialMembrane:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneThickness(fractionOfTissueHeight):"){
		file >>Sim->PeripodialThicnessScale;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneThickness(fractionOfTissueHeight):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneLateralThickness(fractionOfTissueHeight):"){
		file >>Sim->PeripodialLateralThicnessScale;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneLateralThickness(fractionOfTissueHeight):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LumenHeightScale(fractionOfTissueHeight):"){
		file >>Sim->lumenHeightScale;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LumenHeightScale(fractionOfTissueHeight):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneYoungsModulus:"){
		file >>Sim->PeripodialElasticity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneYoungsModulus:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalViscosity:"){
		file >>Sim->PeripodialApicalVisc;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneApicalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalViscosity:"){
		file >>Sim->PeripodialBasalVisc;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneBasalViscosity:" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readLinkerZoneParameters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "BaseOnPeripodialness(bool):"){
		file >>Sim->BaseLinkerZoneParametersOnPeripodialness;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: BaseOnPeripodialness(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalYoungsModulus:"){
		file >>Sim->LinkerZoneApicalElasticity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneApicalYoungsModulus:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalYoungsModulus:"){
		file >>Sim->LinkerZoneBasalYoungsModulus;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneBasalYoungsModulus:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalViscosity:"){
		file >>Sim->LinkerZoneApicalViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneApicalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalViscosity:"){
		file >>Sim->LinkerZoneBasalVisc;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneBasalViscosity:" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readTimeParameters(ifstream& file){
	string currHeader;

	file >> currHeader;
	if(currHeader == "TimeStep(sec):"){
		file >> Sim->dt;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: TimeStep(sec):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "SimulationLength(sec):"){
		float timeInSec;
		file >> timeInSec;
		Sim->SimLength = timeInSec/Sim->dt;
	}
	else{
		cerr<<"Error in reading simulation length, curr string: "<<currHeader<<" should have been: SimulationLength(sec)::" <<endl;
		return false;
	}

	cout<<"Simulation time step	: "<<Sim->dt<<endl;
	cout<<"Simulation Length	: "<<Sim->SimLength<<endl;
	return true;
}

bool ModelInputObject::readPysicalProperties(ifstream& file){
	cout<<"reading physical parameters"<<endl;
	string currHeader;

	file >> currHeader;
	if(currHeader == "YoungsModulusApical:"){
		file >> Sim->EApical;
	}
	else{
		cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusApical:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusBasal:"){
		file >> Sim->EBasal;
	}
	else{
		cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusBasal:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusMid:"){
		file >> Sim->EMid;
	}
	else{
		cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusMid:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[0];
	}
	else{
		cerr<<"Error in reading Young's modulus noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "PoissonsRatio:"){
		file >> Sim->poisson;
	}
	else{
		cerr<<"Error in reading Poisson's ratio, curr string: "<<currHeader<<" should have been: PoissonsRatio:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[1];
	}
	else{
		cerr<<"Error in reading Poisson's Ratio noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "ApicalViscosity:"){
		file >> Sim->ApicalVisc;
	}
	else{
		cerr<<"Error in reading Apical Viscosity, curr string: "<<currHeader<<" should have been: ApicalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[2];
	}
	else{
		cerr<<"Error in reading Viscosity noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalViscosity:"){
		file >> Sim->BasalVisc;
	}
	else{
		cerr<<"Error in reading Basal Viscosity, curr string: "<<currHeader<<" should have been: BasalViscosity:" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readSaveOptions(ifstream& file){
	cout<<"reading save options"<<endl;
	string currHeader;
	file >> currHeader;
	if(currHeader == "SaveImages(bool):"){
		file >> Sim->saveImages;
	}
	else{
		cerr<<"Error in reading image saving option, current string: "<<currHeader<<" should have been: SaveImages(bool):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "SaveData(bool):"){
		file >> Sim->saveData;
	}
	else{
		cerr<<"Error in reading simulation data saving option, current string: "<<currHeader<<" should have been: SaveData(bool):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "ImageSaveInterval(sec):"){
		float timeInSec;
		file >> timeInSec;
		Sim->imageSaveInterval = timeInSec/Sim->dt;
        //the image save interval cannot be smaller than time step!, if this is the case,
        //dataSaveInterval will be zero! Check and correct if necessary
        if (Sim->imageSaveInterval < 1 ){
          Sim->imageSaveInterval =1;
        }
	}
	else{
		cerr<<"Error in reading image save interval, current string: "<<currHeader<<" should have been: ImageSaveInterval(sec)" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "DataSaveInterval(sec):"){
		float timeInSec;
		file >> timeInSec;
		Sim->dataSaveInterval = timeInSec/Sim->dt;
        //the data save interval cannot be smaller than time step!, if this is the case,
        //dataSaveInterval will be zero! Check and correct if necessary
        if (Sim->dataSaveInterval < 1 ){
          Sim->dataSaveInterval =1;
        }
        cout<<"dataSaveInterval as read from file: "<<Sim->dataSaveInterval<<endl;
	}
	else{
		cerr<<"Error in reading data save interval, current string: "<<currHeader<<" should have been: DataSaveInterval(sec)" <<endl;
		return false;
	}
	//cout<<"Sim->saveImages: "<<Sim->saveImages<<" Sim->saveData: "<<Sim->saveData<<" "<<" datainterval: "<<Sim->dataSaveInterval<<" imageinterval: "<<Sim->dataSaveInterval<<endl;
	return true;
}

bool ModelInputObject::readShapeChangeOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	int n;
	cout<<"entered read shape chenge options, current header: "<<currHeader<<endl;
	if(currHeader == "NumberofShapeChangeFunctions(int):"){
		file >> n;
		Sim->nShapeChangeFunctions = n;
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: NumberofShapeChangeFunctions(int):" <<endl;
		return false;
	}
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		//cout<<"inside the loop, read shape change options, current header: "<<currHeader<<endl;
		if(currHeader == "ShapeChangeFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
			cerr<<"Error in reading shpae change type, curr string: "<<currHeader<<", should have been: ShapeChangeFunctionType(int-seeDocumentation):" <<endl;
			return false;
		}
		if (type == 1){
			bool success = readShapeChangeType1(file);
			if (!success){
				return false;
			}
		}
		else{
			cerr<<"Error in reading shape change type, please enter a valid type: {1}, current type: "<<type<<endl;
			return false;
		}
	}
	return true;
}

bool ModelInputObject::readPlasticDeformationOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsPlasticDeformation(bool):"){
		file >> Sim->thereIsPlasticDeformation;
	}
	else{
		cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: ThereIsPlasticDeformation(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "VolumeConserved(bool):"){
		file >> Sim->volumeConservedInPlasticDeformation;
	}
	else{
		cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: VolumeConserved(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DeformationRate(FractionPerHour):"){
		file >> Sim->plasticDeformationRate;
	}
	else{
		cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: DeformationRate(FractionPerHour):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readMyosinOptions(ifstream& file){
	string currHeader;
	//regardless of how many myosin functions I have, I need the diffusion constant and the force per myosin molecule
	file >> currHeader;
	cout<<"entered read myosin concentration options, current header: "<<currHeader<<endl;
	if(currHeader == "MyosinDiffusionRate(double-1/sec):"){
		file >> Sim->kMyo ;
	}
	else{
		cerr<<"Error in reading myosin diffusion rate, curr string: "<<currHeader<<", should have been: MyosinDiffusionRate(double-1/sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ForcePerMyosinMolecule(double):"){
		file >> Sim->forcePerMyoMolecule ;
	}
	else{
		cerr<<"Error in reading force per myosin molecule, curr string: "<<currHeader<<", should have been: ForcePerMyosinMolecule(double):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ThereIsMyosinFeedback(bool):"){
		file >> Sim->thereIsMyosinFeedback;
	}
	else{
		cerr<<"Error in reading force per myosin molecule, curr string: "<<currHeader<<", should have been: ThereIsMyosinFeedback(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MyosinFeedbackCap:"){
		file >> Sim->MyosinFeedbackCap;
	}
	else{
		cerr<<"Error in reading force per myosin molecule, curr string: "<<currHeader<<", should have been: MyosinFeedbackCap:" <<endl;
		return false;
	}
	//Reading number of myosin funcitons as input
	file >> currHeader;
	int n;
	if(currHeader == "NumberofMyosinFunctions(int):"){
		file >> n;
		Sim->nMyosinFunctions = n;
	}
	else{
		cerr<<"Error in reading myosin concentration options, curr string: "<<currHeader<<", should have been: NumberofMyosinFunctions(int):" <<endl;
		return false;
	}
	//now I can read the myosin parameters
	for (int i = 0; i<n; ++i){
		readGridBasedMyosinFunction(file);
		//cout<<"inside the loop, read shape change options, current header: "<<currHeader<<endl;
	}
	return true;
}


bool ModelInputObject::readGridBasedMyosinFunction(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read myosin, current header: "<<currHeader<<endl;
	float initialtimeInSec;
	int initTime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool isApical = true;
	bool isPolarised = false;
	bool manualStripes = false;
	double stripeSize1,stripeSize2, initialPoint, endPoint, manualcEq, tetha;
	int gridX, gridY;
	double cEq;
	double tet;
	double** cEqMatrix;
	double** angleMatrix;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtimeInSec;
		initTime = initialtimeInSec/Sim->dt;
	}
	else{
		cerr<<"Error in reading myosin stimuli options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "isApical(bool):"){
			file >> isApical;
			cout<<"isApical: "<<isApical<<endl;

	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: isApical(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "isPolarised(bool):"){
			file >> isPolarised;
			cout<<"isPolarised: "<<isPolarised<<endl;
	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: isPolarised(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ManualStripes(bool):"){
		file >> manualStripes;
	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: ManualStripes(bool):" <<endl;
		return false;
	}
	if (manualStripes){
		//stripes are manual, I need to read gap sizes, concentration, and x-y limits.
		file >> currHeader;
		if(currHeader == "StripeSizes(micron-0forNoGap):"){
			file >> stripeSize1;
			file >> stripeSize2;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: StripeSizes(micron-0forNoGap):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "LowerEndPoint(micron):"){
			file >> initialPoint;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: LowerEndPoint(micron):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "UpperEndPoint(micron):"){
			file >> endPoint;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: UpperEndPoint(micron):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "EquilibriumMyosinLevel(double):"){
			file >> manualcEq;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: EquilibriumMyosinLevel(double):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "MyosinOrientation(degrees):"){
			file >> tetha;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: MyosinOrientation(degrees):" <<endl;
			return false;
		}
	}
	else{
		//the stripes are not manual, I need to read a grid
		file >> currHeader;
		if(currHeader == "EquilibriumConcentrationFilename(full-path):"){
			string filepath;
			file >> filepath;
			cerr<<" filename is: "<<filepath<<endl;
			const char* name_cEq = filepath.c_str();
			ifstream cEqFile;
			cEqFile.open(name_cEq, ifstream::in);
			if (!(cEqFile.good() && cEqFile.is_open())){
				cerr<<"could not open equilibrium myosin concentrations file: "<<name_cEq<<endl;
				return false;
			}
			//adding the indice of the growth matrix
			cout<<"reading from growth file"<<endl;
			cEqFile >> gridX;
			cEqFile >> gridY;
			cout<<"constructing equilibrium myosin matrix"<<endl;
			cEqMatrix = new double*[(const int) gridX];
			for (int i=0; i<gridX; ++i){
				cEqMatrix[i] = new double[(const int) gridY];
				for (int j=0; j<gridY; ++j){
					cEqMatrix[i][j] = 0.0;
				}
			}
			cout<<"reading equilibrium myosin matrix"<<endl;
			for (int j=gridY-1; j>-1; --j){
				for (int i=0; i<gridX; ++i){
					//cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					cEqFile >> cEq;
					//cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					cEqMatrix[i][j] = cEq;
					//cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<endl;
				}
			}
			cEqFile.close();
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: EquilibriumConcentrationFilename(full-path):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "OrientationAngleFilename(full-path):"){
			string filepath;
			file >> filepath;
			cerr<<" filename is: "<<filepath<<endl;
			const char* name_angle = filepath.c_str();
			ifstream angleFile;
			angleFile.open(name_angle, ifstream::in);
			if (!(angleFile.good() && angleFile.is_open())){
				cerr<<"could not open orientation angle file: "<<name_angle<<endl;
				return false;
			}
			//adding the indice of the growth matrix
			cout<<"reading from orientation angle file"<<endl;
			angleFile >> gridX;
			angleFile >> gridY;
			cout<<"constructing myosin orientation matrix"<<endl;
			angleMatrix = new double*[(const int) gridX];
			for (int i=0; i<gridX; ++i){
				angleMatrix[i] = new double[(const int) gridY];
				for (int j=0; j<gridY; ++j){
					angleMatrix[i][j] = 0.0;
				}
			}
			cout<<"reading myosin orientation matrix"<<endl;
			for (int j=gridY-1; j>-1; --j){
				for (int i=0; i<gridX; ++i){
					//cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					angleFile >> tet;
					angleMatrix[i][j] =  tet;
				}
			}
			angleFile.close();
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: OrientationAngleFilename(full-path):" <<endl;
			return false;
		}
	}
	MyosinFunction* MFp;
	int Id = Sim->myosinFunctions.size();
	if(manualStripes){
		MFp = new MyosinFunction(Id, isApical, isPolarised, initTime, applyToColumnarLayer, applyToPeripodialMembrane, stripeSize1,stripeSize2, initialPoint, endPoint, manualcEq,tetha);

	}
	else{
		MFp = new MyosinFunction(Id, isApical, isPolarised, initTime, applyToColumnarLayer, applyToPeripodialMembrane, gridX, gridY, cEqMatrix, angleMatrix);

	}
	Sim->myosinFunctions.push_back(MFp);
	return true;
}

bool ModelInputObject::readShapeChangeType1(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read shape chenage type 1, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	float Rate;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
		//Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		cerr<<"Error in reading shape chenage options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> Rate;
		Rate  *= timeMultiplier;
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<endl;
		return false;
	}
	cerr<<"Shape change of type : 1, from time: "<<initialtime<<" to "<<finaltime<<" applicable to (C, P) "<<applyToColumnarLayer <<" "<<applyToPeripodialMembrane<<" Rate: "<<Rate<<endl;
	GrowthFunctionBase* GSBp;
	int Id = Sim->ShapeChangeFunctions.size();
	//type is 1
	GSBp = new UniformShapeChangeFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, 1, Rate);
	Sim->ShapeChangeFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readShapeChangeType2(ifstream& file){
	return true;
}

bool ModelInputObject::readStretcherSetup(ifstream& file){
	cout<<"reading stretcher options"<<endl;
	string currHeader;
	file >> currHeader;
	if(currHeader == "StretcherAttached(bool):"){
		bool stretcherAttached;
		file >> stretcherAttached;
		//cerr<<"stretcherAttached "<<stretcherAttached<<endl;
		Sim->stretcherAttached = stretcherAttached;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: StretcherAttached(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ClampedOnDV(bool):"){
		file >> Sim->DVClamp;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: ClampedOnDV(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InitialTime(sec):"){
		double inittime;
		file >> inittime;
		Sim->StretchInitialStep = inittime/Sim->dt;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		double endtime;
		file >> endtime;
		Sim->StretchEndStep = endtime/Sim->dt;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMin:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMin = ClampPos;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMin:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMax:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMax = ClampPos;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMax:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxStrain:"){
		double MaxStrain;
		file >> MaxStrain;
		Sim->StretchStrain = MaxStrain;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMax:" <<endl;
		return false;
	}
	//cout<<"StretcherAttached "<<Sim->stretcherAttached<<"InitialStep "<<Sim->StretchInitialStep<<" EndStep: "<<Sim->StretchEndStep<<" ClapmPos: "<<Sim->StretchMin<<" "<<Sim->StretchMax<<" strain: "<<Sim->StretchStrain<<endl;
	return true;
}

bool ModelInputObject::readPipetteSetup(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "PipetteAspitarionActive(bool):"){
		bool PipetteSuction;
		file >> PipetteSuction;
		//cerr<<"stretcherAttached "<<stretcherAttached<<endl;
		Sim->PipetteSuction = PipetteSuction;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: PipetteAspitarionActive(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InitialTime(sec):"){
		double inittime;
		file >> inittime;
		Sim->PipetteInitialStep = inittime/Sim->dt;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		double endtime;
		file >> endtime;
		Sim->PipetteEndStep = endtime/Sim->dt;
        if (Sim->PipetteEndStep<1){
            Sim->PipetteEndStep =1;
        }
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicalSuciton(bool-will_set_up_basal_suction_if_false):"){
		bool ApicalSuction;
		file >> ApicalSuction;
		Sim->ApicalSuction = ApicalSuction;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: ApicalSuciton(bool-will_set_up_basal_suction_if_false):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre_Position(x,y,z):"){
		double dummy;
		for (int j=0;j<3;++j){
			file >> dummy;
			Sim->pipetteCentre[j] = dummy;
		}
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Centre Position(x,y,z):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Radius(micron):"){
		double dummy;
		file >> dummy;
		Sim->pipetteRadius = dummy;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Radius(micron):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Effect_Debth(micron):"){
		double dummy;
		file >> dummy;
		Sim->pipetteDepth = dummy;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Effect_Debth(micron):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Suction_Pressure(x,y,z-unit):"){
		double dummy;
		for (int j=0;j<3;++j){
			file >> dummy;
			Sim->SuctionPressure[j] = dummy;
		}
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Suction_Pressure(x,y,z-unit):" <<endl;
		return false;
	}
	return true;
}
