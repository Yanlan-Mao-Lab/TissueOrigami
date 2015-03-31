/*
 * ModelInputObject.cpp
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#include "ModelInputObject.h"
#include "Simulation.h"
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
	bool Success = true;
	ifstream parametersFile;
	//const char* name_parametersFile = parameterFileName.c_str();
	//parametersFile.open(name_parametersFile, ifstream::in);
	parametersFile.open(parameterFileName, ifstream::in);
	Success = checkFileStatus(parametersFile,parameterFileName);
	if (!Success){
		return Success;
	}
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
			if(currParameterHeader == "InputMeshParameters:"){
				Success  = readMeshParameters(parametersFile);
			}
			else if(currParameterHeader == "PeripodiumParameters:"){
				Success  = readPeripodiumParameters(parametersFile);
			}
			else if(currParameterHeader == "TimeParameters:"){
				Success  = readTimeParameters(parametersFile);
			}
			else if(currParameterHeader == "PysicalProperties:"){
				Success  = readPysicalProperties(parametersFile);
			}
			else if (currParameterHeader == "SaveOptions:"){
				Success  = readSaveOptions(parametersFile);
			}
			else if (currParameterHeader == "GrowthOptions:"){
				Success  = readGrowthOptions(parametersFile);
			}
			else if(currParameterHeader == "ShapeChangeOptions:"){
				Success  = readShapeChangeOptions(parametersFile);
			}
			else if(currParameterHeader == "Stretcher:"){
				Success  = readStretcherSetup(parametersFile);
			}
			else {
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
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		cout<<"inside the loop, read growth options, current header: "<<currHeader<<endl;
		if(currHeader == "GrowthFunctionType(int-seeDocumentation):"){
			file >> type;
			Sim->GrowthFunctionTypes.push_back(type);
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
	return true;
}

bool ModelInputObject::readGrowthType1(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 1, current header: "<<currHeader<<endl;
	if(currHeader == "InitialTime(sec):"){
		float initialtime;
		file >> initialtime;
		Sim->GrowthParameters.push_back(initialtime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		float finaltime;
		file >> finaltime;
		Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		float maxRate;
		//reading the rate in x,y,z dimensions from input file
		for (int i=0;i<3; ++i){
			file >> maxRate;
			//converting the rate input as fraction per hour, into fraction per time step:
			maxRate /= (60*60/Sim->dt);
			cout<<"maxRate: "<<maxRate<<endl;
			Sim->GrowthParameters.push_back(maxRate);
		}
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readGrowthType2(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 2, current header: "<<currHeader<<endl;
	if(currHeader == "InitialTime(sec):"){
		float initialtime;
		file >> initialtime;
		Sim->GrowthParameters.push_back(initialtime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		float finaltime;
		file >> finaltime;
		Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre:"){
		float x,y;
		file >> x;
		file >> y;
		Sim->GrowthParameters.push_back(x);
		Sim->GrowthParameters.push_back(y);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Centre:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InnerRadius:"){
		float innerR;
		file >> innerR;
		Sim->GrowthParameters.push_back(innerR);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "OuterRadius:"){
		float outerR;
		file >> outerR;
		Sim->GrowthParameters.push_back(outerR);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		float maxRate;
		//reading the rate in x,y,z dimensions from input file
		for (int i=0;i<3; ++i){
			file >> maxRate;
			//converting the rate input as fraction per hour, into fraction per time step:
			maxRate /= (60*60/Sim->dt);
			cout<<"maxRate: "<<maxRate<<endl;
			Sim->GrowthParameters.push_back(maxRate);
		}
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-DV,AP,AB):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readGrowthType3(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 3, current header: "<<currHeader<<endl;
	if(currHeader == "InitialTime(sec):"){
		float initialtime;
		file >> initialtime;
		Sim->GrowthParameters.push_back(initialtime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		float finaltime;
		file >> finaltime;
		Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
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
		Sim->GrowthParameters.push_back(Sim->GrowthMatrices.size());
		cout<<"reading from growth file"<<endl;
		int gridX, gridY;
		GrowthRateFile >> gridX;
		Sim->GrowthParameters.push_back(gridX);
		GrowthRateFile >> gridY;
		Sim->GrowthParameters.push_back(gridY);
		float rate;
		double timeMultiplier = Sim->dt / 3600.0;
		cout<<"constructing growth matrix"<<endl;
		double*** GrowthMatrix;
		GrowthMatrix = new double**[(const int) gridX];
		for (int i=0; i<gridX; ++i){
			GrowthMatrix[i] = new double*[(const int) gridY];
			for (int j=0; j<gridY; ++j){
				GrowthMatrix[i][j] = new double[3];
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] = 0.0;
				}
			}
		}
		cout<<"reading growth matrix"<<endl;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				for (int k=0; k<3; ++k){
					cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					GrowthRateFile >> rate;
					cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					GrowthMatrix[i][j][k] = rate;
					cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<endl;
				}
			}
		}
		GrowthRateFile.close();
		cout<<"growth matrix: "<<endl;
		for (int i=0; i<gridX; ++i){
			for (int j=0; j<gridY; ++j){
				for (int k=0; k<3; ++k){
					cout<<GrowthMatrix[i][j][k]<<" ";
					GrowthMatrix[i][j][k] *= timeMultiplier;
				}
				cout<<"	";
			}
			cout<<endl;
		}cout<<endl;
		/*GrowthMatrix[0][0][0] = 1.00;
		GrowthMatrix[0][0][1] = 0.50;
		GrowthMatrix[0][0][2] = 0.00;
		GrowthMatrix[0][1][0] = 1.00;
		GrowthMatrix[0][1][1] = 0.50;
		GrowthMatrix[0][1][2] = 0.00;
		GrowthMatrix[0][2][0] = 1.00;
		GrowthMatrix[0][2][1] = 0.50;
		GrowthMatrix[0][2][2] = 0.00;

		GrowthMatrix[1][0][0] = 1.00;
		GrowthMatrix[1][0][1] = 0.75;
		GrowthMatrix[1][0][2] = 0.00;
		GrowthMatrix[1][1][0] = 1.00;
		GrowthMatrix[1][1][1] = 0.625;
		GrowthMatrix[1][1][2] = 0.00;
		GrowthMatrix[1][2][0] = 1.00;
		GrowthMatrix[1][2][1] = 0.50;
		GrowthMatrix[1][2][2] = 0.00;

		GrowthMatrix[2][0][0] = 1.00;
		GrowthMatrix[2][0][1] = 1.00;
		GrowthMatrix[2][0][2] = 0.00;
		GrowthMatrix[2][1][0] = 1.00;
		GrowthMatrix[2][1][1] = 0.75;
		GrowthMatrix[2][1][2] = 0.00;
		GrowthMatrix[2][2][0] = 1.00;
		GrowthMatrix[2][2][1] = 0.50;
		GrowthMatrix[2][2][2] = 0.00;
		for (int i=0; i<3; ++i){
			for (int j=0; j<3; ++j){
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] *= timeMultiplier;
				}
			}
		}*/
		Sim->GrowthMatrices.push_back(GrowthMatrix);

	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Filename(full-path):" <<endl;
		return false;
	}
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
	return true;
}



bool ModelInputObject::readMeshType4(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "MeshFile(full-path):"){
		file >> Sim->inputMeshFileName;
	}
	else{
		cerr<<"Error in reading mesh row number, curr string: "<<currHeader<<", should have been: MeshFile(full-path):" <<endl;
		return false;
	}
	//checking consistency:
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

bool ModelInputObject::readPeripodiumParameters(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "AddPeripodium:"){
		file >>Sim->AddPeripodium;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: AddPeripodium:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodiumType:"){
		file >>Sim->PeripodiumType;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodiumYoungsModulus:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodiumYoungsModulus:"){
		file >>Sim->PeripodiumElasticity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodiumYoungsModulus:" <<endl;
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
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: NumberofShapeChangeFunctions(int):" <<endl;
		return false;
	}
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		cout<<"inside the loop, read shape change options, current header: "<<currHeader<<endl;
		if(currHeader == "ShapeChangeFunctionType(int-seeDocumentation):"){
			file >> type;
			Sim->ShapeChangeFunctionTypes.push_back(type);
		}
		else{
			cerr<<"Error in reading growth type, curr string: "<<currHeader<<", should have been: ShapeChangeFunctionType(int-seeDocumentation):" <<endl;
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

bool ModelInputObject::readShapeChangeType1(ifstream& file){
	string currHeader;
	file >> currHeader;
	cout<<"entered read shape change type 1, current header: "<<currHeader<<endl;
	if(currHeader == "InitialTime(sec):"){
		float initialtime;
		file >> initialtime;
		Sim->ShapeChangeParameters.push_back(initialtime);
	}
	else{
		cerr<<"Error in reading shape chenage options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		float finaltime;
		file >> finaltime;
		Sim->ShapeChangeParameters.push_back(finaltime);
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre:"){
		float x,y;
		file >> x;
		file >> y;
		Sim->ShapeChangeParameters.push_back(x);
		Sim->ShapeChangeParameters.push_back(y);
	}
	else{
		cerr<<"Error in reading shape chenge options, curr string: "<<currHeader<<", should have been: Centre:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InnerRadius:"){
		float innerR;
		file >> innerR;
		Sim->ShapeChangeParameters.push_back(innerR);
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "OuterRadius:"){
		float outerR;
		file >> outerR;
		Sim->ShapeChangeParameters.push_back(outerR);
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "AxisOfEffect(DV,AP,AB,DVAP,DVAB,APAB):"){
		string axisOfEffect;
		file >>  axisOfEffect;
		int axis;
		if ( axisOfEffect == "DV"){
			axis = 0;
		}
		else if ( axisOfEffect == "AP"){
			axis = 1;
		}
		else if ( axisOfEffect == "AB"){
			axis = 2;
		}
		else if ( axisOfEffect == "DVAP"){
			axis = 3;
		}
		else if ( axisOfEffect == "DVAB"){
			axis = 4;
		}
		else if ( axisOfEffect == "APAB"){
			axis = 5;
		}
		else {
			cerr<<"Error in reading shape change axis, the options are {DV,AP,AB,DVAP,DVAB,APAB}, current input: "<<axisOfEffect<<endl;
			return false;
		}
		Sim->ShapeChangeParameters.push_back(axis);
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: AxisOfEffect(x,y,z,xy,xz,yz):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour):"){
		float maxRate;
		//reading the rate in x,y,z dimensions from input file
		file >> maxRate;
		//converting the rate input as fraction per hour, into fraction per time step:
		maxRate /= (60*60/Sim->dt);
		cout<<"maxRate: "<<maxRate<<endl;
		Sim->ShapeChangeParameters.push_back(maxRate);
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<endl;
		return false;
	}
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
