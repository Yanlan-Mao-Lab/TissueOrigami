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
};

ModelInputObject::~ModelInputObject(){
};

bool ModelInputObject::readParameters(){
	bool Success = true;
	ifstream parametersFile;
	const char* name_parametersFile = parameterFileName.c_str();
	parametersFile.open(name_parametersFile, ifstream::in);
	Success = checkFileStatus(parametersFile,parameterFileName);
	if (!Success){
		return Success;
	}
	while (!parametersFile.eof()){
		string currline;
		getline(parametersFile,currline);
		if (!parametersFile.eof()){
			if (currline.empty()){
				continue;
			}
			istringstream currSStrem(currline);
			string currParameterHeader;
			currSStrem >> currParameterHeader;

			if(currParameterHeader == "TimeParameters:"){
				Success  = readTimeParameters(parametersFile);
			}
			else if(currParameterHeader == "PysicalProperties:"){
				Success  = readPysicalProperties(parametersFile);
			}
			else if (currParameterHeader == "SaveOptions:"){
				Success  = readSaveOptions(parametersFile);
			}
			else {
				cerr<<"Unidentified parameter input line: "<<endl;
				cerr<<"		"<<currParameterHeader<<endl;
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
	return true;
}

bool ModelInputObject::readSaveOptions(ifstream& file){
	cout<<"reading save options"<<endl;
	return true;
}
