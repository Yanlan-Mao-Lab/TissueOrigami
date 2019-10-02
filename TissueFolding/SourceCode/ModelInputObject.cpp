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
			else if(currParameterHeader == "NodeBindingOptions:"){
				/**
				 * Inputs relating to binding the displacement the nodes of the tissue through the private function ModelInputObject#readNodeBindingParameters
				 */
				Success  = readNodeBindingParameters(parametersFile);
			}
			else if(currParameterHeader == "ExternalViscositySetup:"){
				/**
				 * Inputs relating to the external viscosity felt by the tissue through private function ModelInputObject#readExternalViscosityParameters
				 */
				Success  = readExternalViscosityParameters(parametersFile);
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
			else if(currParameterHeader == "ECM_Perturbation:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readECMPerturbation(parametersFile);
			}
			else if(currParameterHeader == "Stiffness_Perturbation:"){
				/**
				 * Setting perturbations to tissue stiffness, current setup includes changing Young's modulus for apical, basal or whole tissue at a given time point.
				 */
				Success  = readStiffnessPerturbation(parametersFile);
			}
			else if(currParameterHeader == "ApicoBasalVolumeRedistributionOptions:"){
				/**
				 * Setting volume redistribution bewtten apical and basal sides of the tissue, cells redistributing their volume more apically or basally.
				 */
				Success  = readApicoBasalVolumeRedistribution(parametersFile);
			}
			else if(currParameterHeader == "Marker_Ellipses:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readMarkerEllipseBandOptions(parametersFile);
			}
			else if(currParameterHeader == "Cell_Migration:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readCellMigrationOptions(parametersFile);
			}
			else if(currParameterHeader == "ExplicitECMOptions:"){
				/**
				 * Setting perturbations to ECM, current setup includes softening of apical or basal ECM at a given time point.
				 */
				Success  = readExplicitECMOptions(parametersFile);
			}
			else if(currParameterHeader == "AdhesionOptions:"){
				/**
				 * Setting the flag to activate adhesion in the model
				 */
				Success  = readAdhesionOptions(parametersFile);
			}
			else if(currParameterHeader == "NodeCollapseOptions:"){
				/**
				 * Setting the flag to activate node collapsing for elements dangerously close to flipping
				 */
				Success  = readNodeCollapseOptions(parametersFile);
			}
			else if(currParameterHeader == "ExplicitActinOptions:"){
				/**
				 * Setting explicit actin options, The actin layer will not grow in z.
				 */
				Success  = readExplicitActinOptions(parametersFile);
			}
			else if(currParameterHeader == "ColumnViseVolumeConservationOptions:"){
				/**
				 * Setting volume redistributions options, The volume will be redistributed within a column of elements, starting from the apical surface.
				 */
				Success  =readColumnViseVolumeConservationOptions(parametersFile);
			}
			else if(currParameterHeader == "MutationOptions:"){
				/**
				 * Setting mutation options.
				 */
				Success  = readMutationOptions(parametersFile);
			}
			else if(currParameterHeader == "zShellOptions:"){
				/**
				 * Setting z enclosement options.
				 */
				Success  = readEnclosementOptions(parametersFile);
			}
			else if (currParameterHeader == "ArtificialRelaxationOptions:"){
				Success  = readartificialRelaxationOptions(parametersFile);;
			}
			else if(currParameterHeader == "LumenOptions:"){
				/**
				 * Inputs relating to lumen options.
				 */
				Success  = readLumenOptions(parametersFile);
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
	/**
	 * This function will check the status of the input file that was provided. It will
	 * return true if the file is opened successfully, and false otherwise.
	 *
	 */
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
	/**
	 * This function will read the growth options. The generic growth related parameters will be read here, followed by
	 * growth type specific parameters read in specialised functions.
	 * The first parameter to read will be the number of growth options.
	 */
	string currHeader;
	file >> currHeader;
	int n;
	if(currHeader == "NumberofGrowthFunctions(int):"){
		file >> n;
		Sim->nGrowthFunctions = n;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: NumberofGrowthFunctions(int):" <<endl;
		return false;
	}
	/**
	 * The second parameter will identify if the relative positions of elements that identify their growth inputs are pinned to the initial
	 * mesh or is updated at each time step through the simulation.
	 * There will be an option to update the pinned values next.
	 */
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
	/**
	 * Then the following set of parameters defines at which time points the pinning is updated, and the timing of those updates, in seconds.
	 */
	file >> currHeader;
	if(currHeader == "PinningUpdateTimes(number-times(sec)):"){
		file >> Sim->nGrowthPinning;
		if ( Sim->nGrowthPinning > 0){
			Sim->growthPinUpdateTime = new int[(const int) Sim->nGrowthPinning];
			Sim->growthPinUpdateBools = new bool[(const int) Sim->nGrowthPinning];

			for (int j=0; j<Sim->nGrowthPinning; ++j){
				file >> Sim->growthPinUpdateTime[j];
				Sim->growthPinUpdateBools[j] = false;
			}
		}
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: PinningUpdateTimes(number-times(sec)):" <<endl;
		return false;
	}
	/**
	 * The final generic growth parameter will define the interpolation type of growth value read from the grid, given the relative position
	 * of the mesh. Option 0 is a step function, that applies the exact growth rate that is given for the grid point the element finds itself in.
	 * Option 1 carries out a linear interpolation between the nearest 4 corners.
	 *
	 */
	file >> currHeader;
	if (currHeader == "GridGrowthsInterpolationType(0=step,1=linear):"){
		file >> Sim->gridGrowthsInterpolationType;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: GridGrowthsInterpolationType(0=step,1=linear):" <<endl;
		return false;
	}
	/**
	 * Next, the individual specialised growth options are read. The first parameter will read the growth type, and a related specialised function will read
	 * the rest. Growth option types:
	 * 	1: Uniform
	 * 	2: Ring
	 * 	3: From input file grid.
	 */
	for (int i = 0; i<n; ++i){
		file >> currHeader;
		int type;
		if(currHeader == "GrowthFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
			cerr<<"Error in reading growth type, curr string: "<<currHeader<<", should have been: GrowthFunctionType(int-seeDocumentation):" <<endl;
			return false;
		}
		bool Success = true;
		if (type == 1){
			/***
			 * Reading uniform growth options
			 */
			Success = readGrowthType1(file);
		}
		else if (type == 2){
			/***
			 * Reading ring type growth options
			 */
			Success = readGrowthType2(file);
		}
		else if (type == 3){
			/***
			 * Reading grid type growth options
			 */
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
	std::cout<<"Finalised reading growth options"<<endl;
	return true;
}

bool ModelInputObject::readGrowthType1(ifstream& file){
	/**
	 * This funtion reads in uniform growth options.
	 * The format with placeholder values:
	 *   InitialTime(sec): 100
	 *   FinalTime(sec): 200
	 *   ApplyToColumnarLayer(bool): 1
	 *   ApplyToPeripodialMembrane(bool): 0
	 *   ApplyToBasalECM(bool): 0
	 *   ApplyToLateralECM(bool): 1
	 *   MaxValue(fractionPerHour-DV,AP,AB): 0.002 0.001 0.000
	 *   Angle(degrees): 45.0
	 */
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 1, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToBasalECM = false;
	bool applyToLateralECM = false;
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
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		/**
		 * The three values required for the MaxValue represent growth in dorsal-ventral, anterior-posterior and apical-basal axes, respectively.
		 * The growth rate is input per hour but converted to per second for ease in calculation, as the time step of simulation is recorded
		 * in seconds.
		 */
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
		/**
		 * The growth orientation angle is input in degrees, but converted to radians for further calculations.
		 */
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degrees):" <<endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 1
	GSBp = new UniformGrowthFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, DVRate, APRate,  ABRate, angle);
	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readGrowthType2(ifstream& file){
	/**
	 * This function reads in ring shaped growth options. This implements a growth gradient.
	 * Here, the growth is max at the inner radius of the ring
	 * and decays to zero at the outer radius.
	 * The format with placeholder values:
	 *   InitialTime(sec): 100
	 *   FinalTime(sec): 200
	 *   ApplyToColumnarLayer(bool): 1
	 *   ApplyToPeripodialMembrane(bool): 0
	 *   ApplyToBasalECM(bool): 0
	 *   ApplyToLateralECM(bool): 1
	 *   Centre: 5 3
	 *   InnerRadius: 0
	 *   OuterRadius: 10
	 *   MaxValue(fractionPerHour-DV,AP,AB): 0.002 0.001 0.000
	 *   Angle(degrees): 45.0
	 */
	string currHeader;
	file >> currHeader;
	cout<<"entered read growth type 2, current header: "<<currHeader<<endl;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToBasalECM = false;
	bool applyToLateralECM = false;
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
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Centre:"){
		/**
		 * The centre of the ring is given in plane of the tissue with x 7 y coordinates in absolute microns (not relative units).
		 */
		file >> CentreX;
		file >> CentreY;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Centre:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InnerRadius:"){
		/**
		 * The inner radius gives the inner radius of the growth ring, given in absolute micrometers, not relative coordinates.
		 * A value of zero would implement circular growth with max value on the origin.
		 * Note that the inner radius will be the point of maximum growth rate value.
		 */
		file >> innerR;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "OuterRadius:"){
		/**
		 * The outer radius gives the inner radius of the growth ring, , given in absolute micrometers, not relative coordinates.
		 * This is the radius where the implemented growth rate will decay to zero. The gradient is linear.
		 */
		file >> outerR;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour-DV,AP,AB):"){
		/**
		 * The three values required for the MaxValue represent growth in dorsal-ventral, anterior-posterior and apical-basal axes, respectively.
		 * The growth rate is input per hour but converted to per second for ease in calculation, as the time step of simulation is recorded
		 * in seconds.
		 */
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
		/**
		 * The growth orientation angle is input in degrees, but converted to radians for further calculations.
		 */
		file >> angle;
		angle *= M_PI/180.0; 	 // converting to radians
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degrees):" <<endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 2
	GSBp = new RingGrowthFunction(Id, 2, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, CentreX, CentreY, innerR,  outerR, DVRate, APRate,  ABRate, angle);
	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readGrowthType3(ifstream& file){
	/**
	 * This function reads in the grid based growth options.
	 * The format with placeholder values:
	 *   InitialTime(sec): 100
	 *   FinalTime(sec): 200
	 *   ApplyToColumnarLayer(bool): 1
	 *   ApplyToPeripodialMembrane(bool): 0
	 *   ApplyToBasalECM(bool): 0
	 *   ApplyToLateralECM(bool): 1
	 *   Filename(full-path): /path/toWhere/GrowhtIs/myGrowthFile
	 *   zRange: 0 0.5
	 */
	string currHeader;
	file >> currHeader;
	float initialtime;
	float finaltime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool applyToLateralECM = false;
	bool applyToBasalECM = false;
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
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Filename(full-path):"){
		/**
		 * The file path should be a readable path from where the executable is lounced. It can handle relative paths, but full paths are encouraged.
		 */
		string filepath;
		file >> filepath;
		cerr<<" filename is: "<<filepath<<endl;
		const char* name_growthRates = filepath.c_str();
		ifstream GrowthRateFile;
		GrowthRateFile.open(name_growthRates, ifstream::in);
		if (!(GrowthRateFile.good() && GrowthRateFile.is_open())){
			/**
			 * Once the file is attempted to be open, it is checked to make sure it is actually open.
			 * Most common case of a file error is that the path is wrong.
			 * If you checked the file name, check it 5 more times before assuming the error is something else.
			 * No, seriously, it is the file path that is wrong, you have a typo.
			 */
			cerr<<"could not open growth rate file file: "<<name_growthRates<<endl;
			return false;
		}
		//adding the indice of the growth matrix
		/**
		 * The growth grid will contain the size of the grid matrix as the first two
		 * parameters, first these are read.
		 */
		std::cout<<"reading from growth file"<<endl;
		GrowthRateFile >> gridX;
		GrowthRateFile >> gridY;
		float rate;
        cout<<"initiating growth and angle matrices"<<endl;
        /**
         * Then the stacked arrays are generated for the growth rate and the growth orientation matrices, if I have time, I will convert these to vectors, which are
         * much more cleaner then new double arrays. Growth matrix has three layers, x & y coordinates on the bounding box and the 3D growth.
         * Orientation angles have 2 laters, x & y directions only, followed by a single angle value.
         * All initiated as zeros.
         */
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
		/**
		 * Then the growth rates and angles are read in.
		 */
		cout<<"reading growth matrix"<<endl;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				for (int k=0; k<3; ++k){
					//std::cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					GrowthRateFile >> rate;
					//std::cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					GrowthMatrix[i][j][k] = rate;
					//std::cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<endl;
				}
			}
		}
		double angle;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				GrowthRateFile >> angle;
				AngleMatrix[i][j] = angle; // angles in degrees!
			}
		}
		/**
		 * Now the reading from the file is over, I will close the file. All necessary content is copied over to the memory.
		 * I will convert growth rates to be per second for consistency with time step.
		 */
		GrowthRateFile.close();
		for (int i=0; i<gridX; ++i){
			for (int j=0; j<gridY; ++j){
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] *= timeMultiplier;
				}
			}
		}
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Filename(full-path):" <<endl;
		return false;
	}
	/**
	 * The normalised z range of a grid based growth function defined to which height of the tissue the growth function should be applied in.
	 * Note that this is on top of the booleans specifying which components of the tissue the growth function will be applied to.
	 * The z starts from the basal surface of columnar (z=0 is basal columnar surface) and z=1.0 is either the apical surface of the columnar layer
	 * or the basal surface of the peripodial layer (if there is one).
	 */
	double zMin, zMax;
	file >> currHeader;
	if(currHeader == "zRange:"){
		file >> zMin;
		file >> zMax;
	}
	else{
		cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: zRange:" <<endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 3
	GSBp = new GridBasedGrowthFunction(Id, 3, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, gridX, gridY, GrowthMatrix, AngleMatrix);
	GSBp->zMin = zMin;
	GSBp->zMax = zMax;

	Sim->GrowthFunctions.push_back(GSBp);
	return true;
}

bool ModelInputObject::readMeshParameters(ifstream& file){
	/**
	 * This function reads in the mesh parameters, first by obtaining the type of the mesh, and then
	 * calling the specialised mesh reading functions.
	 */
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
		/**
		 * If the mesh type is 2, it will be writing a simple diamond shaped mesh. This was used only in test scenarios.
		 */
		bool Success  = readMeshType2(file);
		if (!Success){
			return false;
		}
	}
	else if ( Sim->MeshType == 4){
		/**
		 * If the mesh type is 4, it will read in the mesh from an input file. This is the main path of
		 * doing things. The specialised function for reading the input mesh file is ModelInputObject#readMeshType4
		 */
		bool Success  = readMeshType4(file);
		if (!Success){
			return false;
		}
	}
	file >> currHeader;
	/**
	 * Once the mesh input is specified, the function reads in the symmetricity options Simulation#symmetricX and Simulation#symmetricY.
	 * If the mesh is symmetric in a specified axis, the code will identify the nodes facing that surface, and fix their position on that given axis.
	 * The selection is relatively crude, and depends on absolute positions wrt to the origin.
	 */
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

bool ModelInputObject::readExternalViscosityParameters(ifstream& file){
	/**
	 * This function reads in the external viscosity parameters, and the form is:
	 * ExternalViscositySetup:
  	 *   ExtendToWholeTissue: 0
  	 *   DiscProperApicalExternalViscosity: 16000.0
  	 *   DiscProperBasalExternalViscosity: 10.0
  	 *   PeripodialMembraneApicalExternalViscosity: 0.0
  	 *   PeripodialMembraneBasalExternalViscosity: 0.0
  	 *   LinkerZoneApicalExternalViscosity: 0.0
  	 *   LinkerZoneBasalExternalViscosity: 0.0
	 *
	 * The identifier line was read to reach this function. The parameter ExtendToWholeTissue will set Simulation#extendExternalViscosityToInnerTissue.
	 * The external viscosity is applied wrt to the movement velocity of the node and its exposed surface.
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ExtendToWholeTissue:"){
		file >>Sim->extendExternalViscosityToInnerTissue;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperApicalExternalViscosity:" <<endl;
		return false;
	}
	/**
	 * The next set of parameters will define the disc proper and peripodial external viscosities. The
	 * linker zone is a specific tissue type that is used when the peripodial and columnar layers are
	 * linked with an explicit layer of mesh elements.
	 */
	file >> currHeader;
	if(currHeader == "DiscProperApicalExternalViscosity:"){
		file >>Sim->externalViscosityDPApical;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperApicalExternalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DiscProperBasalExternalViscosity:"){
		file >>Sim->externalViscosityDPBasal;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperBasalExternalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalExternalViscosity:"){
		file >>Sim->externalViscosityPMApical;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: PeripodialMembraneApicalExternalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalExternalViscosity:"){
		file >>Sim->externalViscosityPMBasal;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: PeripodialMembraneBasalExternalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalExternalViscosity:"){
		file >>Sim->externalViscosityLZApical;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerZoneApicalExternalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalExternalViscosity:"){
		file >>Sim->externalViscosityLZBasal;
	}
	else{
		cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerZoneBasalExternalViscosity:" <<endl;
		return false;
	}
	return true;
}

void ModelInputObject::printErrorMessage(string currentInput, string sourceFuction, string expectedInput){
	cerr<<"Error in reading "<<sourceFuction<<" current input: "<<currentInput<<", should have been: "<<expectedInput<<endl;
}

bool ModelInputObject::readNodeBindingParameters(ifstream& file){
	/**
	 * This function gives the option to bind the degrees of freedom of nodes.
	 * The sample would be:
	 *
	 * NodeBindingOptions:
	 *  bindCircumferenceXYToBasal(bool): 1
  	 *  bindEllipses(bool,nFunctions): 1 2
	 *   bindEllipseBases(bool-x,y,z,nEllipses,ellipseIds): 1 1 0 2 0 1
	 *   bindEllipseBases(bool-x,y,z,nEllipses,ellipseIds): 0 0 1 3 0 1 2
	 * The first boolean binds the circumference of the tissue such that the x&y coordinates of
	 * each column of nodes are bound together, and the circumference cannot bend.
	 * The second boolean gives the option to bind the degrees of freedoms of specified ellipse band
	 * Ids. For ellipse bands, the binding is only at the basal surface.
	 * To help visualise the working of the code, this was initially coded to mimic the attachemt of
	 * the tisue basal surface to the underlying trachea.
	 * The ellipse bands can be explicitly specified in the input, or selected to be emergent with
	 * fold initiation, through the function ModelInputObject#readMarkerEllipseBandOptions.
	 *
	 * Following options provide the selected number of ellipse binding options,
	 * specified to be 2 in the above example (bindEllipses(bool,nFunctions): 1 [2]).
	 * Each function specified the axes to be bound, number of ellipse band ids that are to be bound, and the ids themselves.
	 * such that the second function specifies only the binding of the z degrees of freedom for 3 specified bands: 0 1 and 2.
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "bindCircumferenceXYToBasal(bool):"){
		file >>Sim->thereIsCircumferenceXYBinding;
	}
	else{
		printErrorMessage(currHeader,"NodeBinding options","bindCircumferenceXYToBasal(bool):");
		return false;
	}
	bool ellipseBinding =false;
	int nEllipseFunctions = 0;
	file >> currHeader;
	if(currHeader == "bindEllipses(bool,nFunctions):"){
		file >> ellipseBinding;
		file >> nEllipseFunctions;
	}
	else{
		printErrorMessage(currHeader,"NodeBinding options","bindEllipses(bool,nFunctions):");
		return false;
	}
	for (int i=0; i<nEllipseFunctions; ++i){
		file >> currHeader;
		if(currHeader == "bindEllipseBases(bool-x,y,z,nEllipses,ellipseIds):"){
			vector<bool> boundAxes;
			for (int j=0; j<3; ++j){
				bool axisBound = false;
				file >>axisBound;
				boundAxes.push_back(axisBound);
			}
			Sim->ellipseBasesAreBoundOnAxis.push_back(boundAxes);
			int nEllipses = 0;
			file >>nEllipses;
			vector<int> ellipseIds;
			for (int j=0; j<nEllipses; ++j){
				int ellipseId = -1;
				file >>ellipseId;
				ellipseIds.push_back(ellipseId);
			}
			Sim->ellipseIdsForBaseAxisBinding.push_back(ellipseIds);
		}
		else{
			printErrorMessage(currHeader,"NodeBinding options","bindEllipseBases(bool-x,y,z,nEllipses,ellipseIds):");
			return false;
		}
	}
	return true;
}


bool ModelInputObject::readNodeFixingParameters(ifstream& file){
	/**
	 * This function sets fixing of the nodes. There are two options to limit the
	 * degrees of freedom of specified nodes. A high external viscosity can be applied for
	 * fixing the node, the value specified for x, y & z axes in "FixingViscosity(x,y,z):".
	 * For each of the specified node groups, if the Fix(...)ExtVisc(bool) term is true, the
	 * given external viscosity will be applied. If the boolean is false, the node will be rigidly fixed in position.
	 * Each specified node group is given the option to fix each x,y & z direction independently.
	 * Apical surface will encompass all apical surface including the peripodial if available.
	 * Same with Basal surface. The circumference refers to all the nodes at the boundary of the tissue, excluding
	 * those at symmetricity boundaries (in effect, these are not circumferential).
	 *
	 * Yhe NotumFix option fixes the basal surface only, and the extent of notum is specified in the xFracMin,xFracMax
	 * options, fraction covering the normalised tissue length, notup tip being 0 and pouch tip being 1, at
	 * the DV midline (y-symmetricity border).
	 *
	 * The example would have the syntax:
	 *
	 * NodeFixingOptions:
	 * FixingViscosity(x,y,z): 0   0  32000
  	 * ApicSurfaceFix(bool-x,y,z):   0 0 0   FixApicalExtVisc(bool): 0
  	 * BasalSurfaceFix(bool-x,y,z):  0 0 0   FixBasalExtVisc(bool):  0
  	 * CircumferenceFix(bool-x,y,z): 0 0 0   FixCircWithExtVisc(bool): 0
  	 * ApicCircumFix(bool-x,y,z):    0 0 0   FixApicCircWithExtVisc(bool):  0
  	 * BasalCircumFix(bool-x,y,z):   0 0 0   FixBasalCircWithExtVisc(bool): 0
  	 * LinkerApicCircumFix(bool-x,y,z):  0 0 0  FixLinkerApicCircWithExtVisc(bool):  0
  	 * LinkerBasalCircumFix(bool-x,y,z): 0 0 0  FixLinkerBasalCircWithExtVisc(bool): 0
  	 * NotumFix(bool-x,y,z,double-xFracMin,xFracMax): 0 0 0 -0.1 0.5  FixNotumExtVisc(bool): 0
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "FixingViscosity(x,y,z):"){
		file >>Sim->fixingExternalViscosity[0];
		file >>Sim->fixingExternalViscosity[1];
		file >>Sim->fixingExternalViscosity[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixingViscosity(x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicSurfaceFix(bool-x,y,z):"){
		file >>Sim->ApicalNodeFix[0];
		file >>Sim->ApicalNodeFix[1];
		file >>Sim->ApicalNodeFix[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","ApicSurfaceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixApicalExtVisc(bool):"){
		file >>Sim->ApicalNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixApicalExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalSurfaceFix(bool-x,y,z):"){
		file >>Sim->BasalNodeFix[0];
		file >>Sim->BasalNodeFix[1];
		file >>Sim->BasalNodeFix[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","BasalSurfaceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixBasalExtVisc(bool):"){
		file >>Sim->BasalNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixBasalExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "CircumferenceFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[4][0];
		file >>Sim->CircumferentialNodeFix[4][1];
		file >>Sim->CircumferentialNodeFix[4][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","CircumferenceFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[4];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApicCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[0][0];
		file >>Sim->CircumferentialNodeFix[0][1];
		file >>Sim->CircumferentialNodeFix[0][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","ApicCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixApicCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[0];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixApicCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[1][0];
		file >>Sim->CircumferentialNodeFix[1][1];
		file >>Sim->CircumferentialNodeFix[1][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","BasalCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixBasalCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[1];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixBasalCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerApicCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[2][0];
		file >>Sim->CircumferentialNodeFix[2][1];
		file >>Sim->CircumferentialNodeFix[2][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","LinkerApicCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixLinkerApicCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixLinkerApicCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerBasalCircumFix(bool-x,y,z):"){
		file >>Sim->CircumferentialNodeFix[3][0];
		file >>Sim->CircumferentialNodeFix[3][1];
		file >>Sim->CircumferentialNodeFix[3][2];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","LinkerBasalCircumFix(bool-x,y,z):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixLinkerBasalCircWithExtVisc(bool):"){
		file >>Sim->CircumferentialNodeFixWithHighExternalViscosity[3];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixLinkerBasalCircWithExtVisc(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "NotumFix(bool-x,y,z,double-xFracMin,xFracMax):"){
		file >>Sim->NotumNodeFix[0];
		file >>Sim->NotumNodeFix[1];
		file >>Sim->NotumNodeFix[2];
		file >>Sim->notumFixingRange[0];
		file >>Sim->notumFixingRange[1];
	}
	else{
		printErrorMessage(currHeader,"Fixing options","NotumFix(bool-x,y,z,double-xFracMin,xFracMax):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FixNotumExtVisc(bool):"){
		file >>Sim->NotumNodeFixWithExternalViscosity;
	}
	else{
		printErrorMessage(currHeader,"Fixing options","FixNotumExtVisc(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readManupulationParamters(ifstream& file){
	/**
	 * This option group allows for an ensemble of additional
	 * manipulations on the simulation options, from stiffness perturbations to random forces.
	 * I find it highly unlikely you will use these, probably will remove them soon.
	 *
	 * Manipulations:
	 * AddCurvature(bool): 1
	 * CurvatureDepthAtCentre(double-microns): 2.0
	 * AddSoftPeriphery(bool): 0
	 * SoftPeripheryRange(double-microns): 30.0
	 * SoftnessFraction(double-fraction): 0.1
	 * ApplyToApicalSurface(bool): 1
	 * ApplyToBasalSurface(bool): 0
	 * ApplyToColumnarLayer(bool): 1
	 * ApplyToPeripodialMembrane(bool): 1
	 * AddRandomForce(bool): 0
	 * RandomForceMean(double): 0.0
	 * RandomForceVar(double): 1E-5
	 *
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "AddCurvature(bool):"){
		file >>Sim->addCurvatureToTissue;
	}
	else{
		printErrorMessage(currHeader,"Manipulation options","AddCurvature(bool):");
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
		file >>Sim->softDepth;
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
	/**
	 * The specific function reading the input mesh file name to Simulation#inputMeshFileName
	 * so that the mesh can be read in.
	 *
	 * MeshFile(full-path): <path>
	 */
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
	/**
	 * The options for letting the model write a simple single layered diamond shaped mesh are read in
	 * this function. The user identifies the rows and columns of the deisred mesh, side length of each equilateral triangle
	 * forming the prisms in micrometers and the z height of the prism in micrometers.
	 *
	 * MeshRow(int): 6
	 * MeshColumn(int): 4
	 * SideLength: 1
	 * zHeight: 2
	 */
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
	/**
	 * This section gives the used the option to add a peripodial membrane to the tissue.
	 * The thickness will define the z height of the peripodial layer, and the lateral thickness will
	 * define the thickness of the laterall element layer to be used in peripodial membrane definition.
	 * The lumen height scale will give the empty space height at the centre of the tissue. All length scales
	 * are fractions with respect to the columnar tissue height.
	 *
	 * PeripodialMembraneParameters:
  	 *   AddPeripodialMembrane: 0
  	 *   PeripodialMembraneThickness(fractionOfTissueHeight): 0.45
  	 *   PeripodialMembraneLateralThickness(fractionOfTissueHeight): 0.5
  	 *   LumenHeightScale(fractionOfTissueHeight): 0.25
  	 *   PeripodialMembraneYoungsModulus: 1000.0
  	 *   PeripodialMembraneApicalViscosity: 0.0
  	 *   PeripodialMembraneBasalViscosity: 0.0
	 */
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
		file >>Sim->peripodialApicalViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneApicalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalViscosity:"){
		file >>Sim->peripodialBasalViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneBasalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneMidlineViscosity:"){
		file >>Sim->peripodialMidlineViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneMidlineViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "AdherePeripodialToColumnarInZ(bool):"){
		file >>Sim->adherePeripodialToColumnar;
	}
	else{
		printErrorMessage(currHeader,"Peripodial membrane parameters","AdherePeripodialToColumnarInZ(bool):");
		return false;
	}

	return true;
}

bool ModelInputObject::readLinkerZoneParameters(ifstream& file){
	/**
	 * Most likely will be deleted.
	 */
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
		file >>Sim->linkerZoneApicalViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneApicalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalViscosity:"){
		file >>Sim->linkerZoneBasalViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneBasalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneMidlineViscosity:"){
		file >>Sim->linkerZoneMidlineViscosity;
	}
	else{
		cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneMidlineViscosity:" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readTimeParameters(ifstream& file){
	/**
	 * The time scale parameters of the simulation, giving the
	 * time step and total simulation lenght in seconds.
	 *
	 *  TimeParameters:
  	 *   TimeStep(sec): 120
  	 *   SimulationLength(sec):  147600
	 *
	 */
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
		file >> Sim->SimLength;
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
	/**
	 * The physical properties of the columnar layer are set in this input parameter section.
	 * The Young modulus can be specified independently for apical and basal layers as well as the
	 * remaining tissue in between. There is an option to add random noise to tissue stiffness.
	 * The poisson ratio is uniform for the entire tissue, except for the ECM elements, where it is set to be
	 * 0 (ECM is modelled as compressible - with multiple new safety checks implemented since the addition of this rule, it may be worth to try if it can be relaxed).
	 * The internal viscosity can be set for apical and basal layers independently, the mid line will interpolate.
	 *
	 * PysicalProperties:
	 *   YoungsModulusApical: 1000.0   YoungsModulusBasal: 1000.0      YoungsModulusMid: 1000.0        Noise(%-int): 0
	 *   PoissonsRatio: 0.3            Noise(%-int): 0
	 *   ApicalViscosity: 0.0  Noise(%-int): 0
	 *   BasalViscosity: 0.0
	 */
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
		file >> Sim->discProperApicalViscosity;
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
		file >> Sim->discProperBasalViscosity;
	}
	else{
		cerr<<"Error in reading Basal Viscosity, curr string: "<<currHeader<<" should have been: BasalViscosity:" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MidLineViscosity:"){
		file >> Sim->discProperMidlineViscosity;
	}
	else{
		cerr<<"Error in reading Midline Viscosity, curr string: "<<currHeader<<" should have been: MidLineViscosity:" <<endl;
		return false;
	}

	return true;
}

bool ModelInputObject::readSaveOptions(ifstream& file){
	/**
	 * The save options are specified below. The first two booleans declare the choice to save the
	 * images and data. the second two options identify the frequency (in sec) of saving independently for both.
	 * Since there is a display only version of the model, in most large scale simulations, the options would
	 * be set as below, with the data saved during simulaton and images not saved. For the iterface-free version of
	 * the code, the option to save images will have no functionality.
	 *
	 * SaveOptions:
  	 *   SaveImages(bool): 0
  	 *   SaveData(bool):   1
  	 *   ImageSaveInterval(sec): 1
  	 *   DataSaveInterval(sec):  1800
	 *
	 */
	//cout<<"reading save options"<<endl;
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
	/**
	 * Will change soon.
	 */
	string currHeader;
	file >> currHeader;
	int n;
	//cout<<"entered read shape change options, current header: "<<currHeader<<endl;
	if(currHeader == "NumberofShapeChangeFunctions(int):"){
		file >> n;
		Sim->nShapeChangeFunctions = n;
	}
	else{
		cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: NumberofShapeChangeFunctions(int):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ShapeChangeStartsBelowECMLevel(fraction):"){
			file >> Sim->shapeChangeECMLimit;
		}
		else{
			printErrorMessage(currHeader,"read shape change options","ShapeChangeStartsBelowECMLevel(fraction)");
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
		else if (type ==2){
			bool success = readShapeChangeType2(file);
			if (!success){
				return false;
			}
		}
		else{
			cerr<<"Error in reading shape change type, please enter a valid type: {1},{2} current type: "<<type<<endl;
			return false;
		}
	}
	return true;
}

bool ModelInputObject::readPlasticDeformationOptions(ifstream& file){
	/**
	 *
	 * The plastic defprmation - also termed remodelling - allows for the relaxation of
	 * element elastic deformation to permenant plastic deformation, with the selected
	 * rate. There is a boolean to allow for volume conserved or non-conserved plastic deformation.
	 *
	 * ECM elements always have remodelling, and the remodelling rate is specified with the half-life.
	 * These are set under ModelInputObject#readExplicitECMOption
	 *
	 * PlasticDeformationOptions:
	 *   ThereIsPlasticDeformation(bool): 0
	 *   VolumeConserved(bool): 1
	 *   DeformationRate(FractionPerHour): 0.1
	 *
	 */
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
		if(currHeader == "ApplyToColumnarLayer(bool):"){
				file >> Sim->plasticDeformationAppliedToColumnar;
		}
		else{
			cerr<<"Error in reading  plastic deformation options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToPeripodialMembrane(bool):"){
				file >> Sim->plasticDeformationAppliedToPeripodial;;
		}
		else{
			cerr<<"Error in reading  plastic deformation options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<endl;
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
	if(currHeader == "DeformationHalfLife(hour):"){
		file >> Sim->plasticDeformationHalfLife;
		Sim->plasticDeformationHalfLife *= 3600; //converting to seconds.
	}
	else{
		cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: DeformationHalfLife(hour):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "zDeformationLimits(lowerFraction-upperFraction):"){
		file >> Sim->zRemodellingLowerThreshold;
		file >> Sim->zRemodellingUpperThreshold;
	}
	else{
		cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: zDeformationLimits(lowerFraction-upperFraction):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readMyosinOptions(ifstream& file){
	/**
	 * This is not working properly so should not be included.
	 * We have oscillatoy behaviour.
	 */
	string currHeader;
	//regardless of how many myosin functions I have, I need the diffusion constant and the force per myosin molecule
	file >> currHeader;
	//cout<<"entered read myosin concentration options, current header: "<<currHeader<<endl;
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
	/**
	 * Myosin functionality is not working properly so should not be included.
	 * We have oscillatoy behaviour.
	 */
	string currHeader;
	file >> currHeader;
	//cout<<"entered read myosin, current header: "<<currHeader<<endl;
	float initialtimeInSec;
	int initTime;
	bool applyToColumnarLayer = false;
	bool applyToPeripodialMembrane = false;
	bool isApical = true;
	bool isPolarised = false;
	bool isLateral = false;
	bool manualStripes = false;
	bool useEllipses = false;
	double stripeSize1,stripeSize2, initialPoint, endPoint, manualcEq, tetha;
	double ellipseLateralcEq, ellipsecEq;
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

	}
	else{
		cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: isApical(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "isPolarised(bool):"){
			file >> isPolarised;
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
		file >> currHeader;
		if(currHeader == "UseEllipses(bool):"){
			file >> useEllipses;
		}
		else{
			cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: UseEllipses(bool):" <<endl;
			return false;
		}
		if (useEllipses){
			file >> currHeader;
			if(currHeader == "isLateral(bool):"){
				file >> isLateral;
			}
			else{
				cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: isLateral(bool):" <<endl;
				return false;
			}
			file >> currHeader;
			if(currHeader == "EquilibriumMyosinLevel(double):"){
				file >> ellipsecEq;
			}
			else{
				cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: EquilibriumMyosinLevel(double):" <<endl;
				return false;
			}
			file >> currHeader;
			if(currHeader == "myosinAppliedToEllipses(number,[ellipseId][ellipseId]):"){
				file >> Sim->numberOfMyosinAppliedEllipseBands;
				double ellipseBandId;
				for (int aa=0; aa<Sim->numberOfMyosinAppliedEllipseBands; ++aa){
					file >>ellipseBandId;
					Sim->myosinEllipseBandIds.push_back(ellipseBandId);
					if (ellipseBandId >= 100){
						Sim->thereIsEmergentEllipseMarking = true;
					}
				}
			}
			else{
				cerr<<"Error in reading  myosin stimuli options, curr string: "<<currHeader<<", should have been: myosinAppliedToEllipses(number,[ellipseId][ellipseId]):" <<endl;
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
				cEqFile >> gridX;
				cEqFile >> gridY;
				cEqMatrix = new double*[(const int) gridX];
				for (int i=0; i<gridX; ++i){
					cEqMatrix[i] = new double[(const int) gridY];
					for (int j=0; j<gridY; ++j){
						cEqMatrix[i][j] = 0.0;
					}
				}
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
				angleFile >> gridX;
				angleFile >> gridY;
				angleMatrix = new double*[(const int) gridX];
				for (int i=0; i<gridX; ++i){
					angleMatrix[i] = new double[(const int) gridY];
					for (int j=0; j<gridY; ++j){
						angleMatrix[i][j] = 0.0;
					}
				}
				for (int j=gridY-1; j>-1; --j){
					for (int i=0; i<gridX; ++i){
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
	}
	MyosinFunction* MFp;
	int Id = Sim->myosinFunctions.size();
	if(manualStripes){
		MFp = new MyosinFunction(Id, isApical, isPolarised, initTime, applyToColumnarLayer, applyToPeripodialMembrane, stripeSize1,stripeSize2, initialPoint, endPoint, manualcEq,tetha);

	}
	else if (useEllipses){
		MFp = new MyosinFunction(Id, isApical, isLateral, initTime, applyToColumnarLayer, applyToPeripodialMembrane, ellipseLateralcEq, ellipsecEq);
	}
	else{
		MFp = new MyosinFunction(Id, isApical, isPolarised, initTime, applyToColumnarLayer, applyToPeripodialMembrane, gridX, gridY, cEqMatrix, angleMatrix);

	}
	Sim->myosinFunctions.push_back(MFp);
	return true;
}

bool ModelInputObject::readShapeChangeType1(ifstream& file){
	/**
	 * This will change.
	 */
	string currHeader;
	file >> currHeader;
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
	GSBp = new UniformShapeChangeFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, false /*applyToBasalECM*/, false /*applyToLateralECM*/, 1, Rate);
	Sim->ShapeChangeFunctions.push_back(GSBp);
	return true;
}


bool ModelInputObject::readShapeChangeType2(ifstream& file){
	/**
	 * This will change.
	 */
	string currHeader;
	file >> currHeader;
	float initialtime;
	float finaltime;
	bool applyToBasalECM{false};
	bool applyToLateralECM{false};
	bool applyTissueApical{false};
	bool applyTissueBasal{false};
	bool applyTissueMidline{false};
	bool conserveVolume{false};
	vector <int> markerEllipses;
	double ShapeChangeFractionPerHr;
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","InitialTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
		//Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","FinalTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueApical(bool):"){
			file >> applyTissueApical;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueApical(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueBasal(bool):"){
			file >> applyTissueBasal;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueBasal(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyTissueMidline(bool):"){
			file >> applyTissueMidline;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyTissueMidline(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyToBasalECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ApplyToLateralECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ShapeChangeAppliedToEllipses(number,[ellipseId][ellipseId]):"){
		int n;
		file >> n;
		for (int i=0; i<n; ++i){
			int ellipseId;
			file >> ellipseId;
			markerEllipses.push_back(ellipseId);
		}
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ShapeChangeAppliedToEllipses(number,[ellipseId][ellipseId]):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "xyShapeChange(fractionPerHour):"){
		file >> ShapeChangeFractionPerHr;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","xyShapeChange(fractionPerHour)");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ConserveVolume(bool):"){
		file >> conserveVolume;
	}
	else{
		printErrorMessage(currHeader,"shape change type 2","ConserveVolume(bool)");
		return false;
	}
	GrowthFunctionBase* GSBp;
	int Id = Sim->ShapeChangeFunctions.size();
	//type is 1
	GSBp = new 	markerEllipseBasedShapeChangeFunction(Id, 2, initialtime, finaltime, applyTissueApical, applyTissueBasal, applyTissueMidline, applyToBasalECM, applyToLateralECM, 2, ShapeChangeFractionPerHr, markerEllipses, conserveVolume);
	Sim->ShapeChangeFunctions.push_back(GSBp);
	return true;
}


bool ModelInputObject::readStretcherSetup(ifstream& file){
	/**
	 * The options to attach a stretcher to the tissue are included here.
	 *
	 * The first boolean states if the stretcher will be attached or not, and will set
	 * Simulation#stretcherAttached boolean. The second boolean states the direction of
	 * clamping and sets Simulation#DVClamp. If the value is 1, the tissue is stretched
	 * from the dorsal and ventral tips. If false, the stretch is in the anterior-posterior
	 * direction. The stretching active during the time frame set by InitialTime(sec) and
	 * FinalTime(sec), setting the parameters Simulation#StretcInitialTime and
	 * Simulation#StretchEndTime respectively. The maximum amount of strain in set by  MaxStrain
	 * setting the parameter Simulation#StretchStrain. The maximum is reached linearly during
	 * the stretching period. The scale of clamping is given by the DVClampMin and DVClampMax
	 * parameters (setting Simulation#StretchMin and Simulation#StretchMax, the value being the relative size in the bounding box. In case of AP clamping,
	 * the same scale applies normalised to the AP bounding box length.
	 *
	 * Stretcher:
	 *   StretcherAttached(bool): 1
	 *   ClampedOnDV(bool): 0
	 *   InitialTime(sec): 100
	 *   FinalTime(sec): 200
	 *   DVClampMin: 0.3
	 *   DVClampMax: 0.7
	 *   MaxStrain: 1.5
	 */
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
		Sim->StretchInitialTime = inittime;
	}
	else{
		cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: InitialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		double endtime;
		file >> endtime;
		Sim->StretchEndTime = endtime;
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


bool ModelInputObject::readMarkerEllipseBandOptions(ifstream& file){
	/**
	 * There is a way to draw patterns on the mesh surface, using ellipse bands.
	 * You first set ellipse band number, then give the centres of the ellipses in relative x,
	 * then the inner radia, and outer radia. Positive radia will draw the bands towards the left hand side
	 * of the tissue, whereas negative radia values will draw towards the right-hand-side. The radia are all in relative tissue
	 * length and width, and can be larger than 1 to achieve flat bands.

	 * The parameter set that would perfectly cover the low growth zones of the initial
	 * growth rate are as follows:
	 * Marker_Ellipses:
	 *   numberOfMarkerEllipses(int): 2
	 *   MarkerEllipseXCenters(fractionOfTissueSize): 0.85 0.75
	 *   MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High): 0.22  0.25 0.2  0.25
	 *   MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High): 1.2   1.3   1.2  1.3
    */
	string currHeader;
	file >> currHeader;
	if(currHeader == "numberOfMarkerEllipses(int):"){
		file >> Sim->nMarkerEllipseRanges;
	}
	else{
		cerr<<"Error in reading marker ellipses, curr string: "<<currHeader<<", should have been: numberOfMarkerEllipses(int):" <<endl;
		return false;
	}
	
	file >> currHeader;
	if(currHeader == "MarkerEllipseXCenters(fractionOfTissueSize):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandXCentres.push_back(dummy);
		}
	}
	else{
		cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseXCenters(fractionOfTissueSize):"<<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandR1Ranges.push_back(dummy);
			file >> dummy;	
			Sim->markerEllipseBandR1Ranges.push_back(dummy);
		}
	}
	else{
		cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High):"){
		double dummy;
		for (int i=0; i<Sim->nMarkerEllipseRanges; ++i){
			file >> dummy;	
			Sim->markerEllipseBandR2Ranges.push_back(dummy);
			file >> dummy;	
			Sim->markerEllipseBandR2Ranges.push_back(dummy);
		}
	}
	else{
		cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readApicoBasalVolumeRedistribution(ifstream& file){
	/**
	 * will change
	 *
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "NumberOfVolumeRedistributionFunctions(int):"){
		file >>Sim->nApikobasalVolumeRedistributionFunctions;
	}
	else{
		printErrorMessage(currHeader,"Apicobasal Volume Redistribution","NumberOfVolumeRedistributionFunctions(int):");
		return false;
	}
	for (int i=0; i<Sim->nApikobasalVolumeRedistributionFunctions; ++i){
		file >> currHeader;
		if(currHeader == "timeOfVolumeRedistribution(hr):"){
			double timeInHr;
			file >> timeInHr;
			Sim->apikobasalVolumeRedistributionBeginTimeInSec.push_back(timeInHr*3600);
			file >> timeInHr;
			Sim->apikobasalVolumeRedistributionEndTimeInSec.push_back(timeInHr*3600);
		}
		else{
			printErrorMessage(currHeader,"Apicobasal Volume Redistribution","timeOfVolumeRedistribution(hr):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "shrinksApicalSide(bool):"){
			bool shrinkApical;
			file >> shrinkApical;
			Sim->apikobasalVolumeRedistributionFunctionShrinksApical.push_back(shrinkApical);
		}
		else{
			printErrorMessage(currHeader,"Apicobasal Volume Redistribution","shrinksApicalSide(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "redistributionFractionOver24Hours(double,0-1):"){
			double fraction;
			file >> fraction;
			Sim->apikobasalVolumeRedistributionScales.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"Apicobasal Volume Redistribution","redistributionFractionOver24Hours(double,0-1):");
			return false;
		}

		file >> currHeader;
		if(currHeader == "volumeRedistributionAppliedToEllipses(number,[ellipseId][ellipseId]):"){
			int numberOfVolumeRedistributionAppliesEllipseBands;
			file >>numberOfVolumeRedistributionAppliesEllipseBands;
			Sim->apikobasalVolumeRedistributionFunctionEllipseNumbers.push_back(numberOfVolumeRedistributionAppliesEllipseBands);
			double ellipseBandId;
			Sim->apikobasalVolumeRedistributionFunctionEllipseBandIds.push_back(vector<int>(0));
			for (int aa=0; aa<Sim->apikobasalVolumeRedistributionFunctionEllipseNumbers[i]; ++aa){
				file >>ellipseBandId;
				Sim->apikobasalVolumeRedistributionFunctionEllipseBandIds[i].push_back(ellipseBandId);
				if (ellipseBandId >= 100){
					Sim->thereIsEmergentEllipseMarking = true;
				}
			}
		}
		else{
			printErrorMessage(currHeader,"Apicobasal Volume Redistribution","volumeRedistributionAppliedToEllipses(number,[ellipseId][ellipseId]):");
			return false;
		}
	}
	return true;
}

bool ModelInputObject::readStiffnessPerturbation(ifstream& file){
	/**
	 * Will change
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsStiffnessPerturbation(bool):"){
		file >> Sim->ThereIsStiffnessPerturbation;
	}
	else{
		cerr<<"Error in reading stiffness perturbations, curr string: "<<currHeader<<", should have been: ThereIsStiffnessPerturbation(bool):" <<endl;
		return false;
	}
	int nStiffnessFunctions;
	file >> currHeader;
	if(currHeader == "NumberOfStiffnessPerturbations(int):"){
		file >>nStiffnessFunctions;
	}
	else{
		printErrorMessage(currHeader,"stiffness perturbations","NumberOfStiffnessPerturbations(int):");
		return false;
	}

	for (int i=0; i<nStiffnessFunctions; ++i){
		file >> currHeader;
		if(currHeader == "ApplyToApically(bool):"){
			bool ThereIsApicalStiffnessPerturbation;
			file >>ThereIsApicalStiffnessPerturbation;
			Sim->ThereIsApicalStiffnessPerturbation.push_back(ThereIsApicalStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyToApically(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyBasally(bool):"){
			bool ThereIsBasalStiffnessPerturbation;
			file >>ThereIsBasalStiffnessPerturbation;
			Sim->ThereIsBasalStiffnessPerturbation.push_back(ThereIsBasalStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyBasally(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToWholeTissue(bool):"){
			bool ThereIsWholeTissueStiffnessPerturbation;
			file >>ThereIsWholeTissueStiffnessPerturbation;
			Sim->ThereIsWholeTissueStiffnessPerturbation.push_back(ThereIsWholeTissueStiffnessPerturbation);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","ApplyToWholeTissue(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "Basolateral(bool):"){
			bool basolateral;
			file >>basolateral;
			Sim->ThereIsBasolateralStiffnessPerturbation.push_back(basolateral);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","BasolateralWithApicalRelaxation");
			return false;
		}
		file >> currHeader;
		if(currHeader == "BasolateralWithApicalRelaxation(bool):"){
			bool basolateral;
			file >>basolateral;
			Sim->ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation.push_back(basolateral);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","BasolateralWithApicalRelaxation");
			return false;
		}
		file >> currHeader;
		if(currHeader == "timeOfStiffeningPerturbation(hr):"){
			double timeInHr;
			file >> timeInHr;
			Sim->stiffnessPerturbationBeginTimeInSec.push_back(timeInHr*3600);
			file >> timeInHr;
			Sim->stiffnessPerturbationEndTimeInSec.push_back(timeInHr*3600);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","timeOfStiffeningPerturbation(hr):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessPerturbationAppliedToEllipses(number,[ellipseId][ellipseId]):"){
			int numberOfStiffnessPerturbationAppliesEllipseBands;
			file >>numberOfStiffnessPerturbationAppliesEllipseBands;
			Sim->numberOfStiffnessPerturbationAppliesEllipseBands.push_back(numberOfStiffnessPerturbationAppliesEllipseBands);
			double ellipseBandId;
			Sim->stiffnessPerturbationEllipseBandIds.push_back(vector<int>(0));
			for (int aa=0; aa<Sim->numberOfStiffnessPerturbationAppliesEllipseBands[i]; ++aa){
				file >>ellipseBandId;
				Sim->stiffnessPerturbationEllipseBandIds[i].push_back(ellipseBandId);
				if (ellipseBandId >= 100){
					Sim->thereIsEmergentEllipseMarking = true;
				}
			}
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","stiffnessPerturbationAppliedToEllipses(number,[ellipseId][ellipseId]):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangedToFractionOfOriginal(double):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->stiffnessChangedToFractionOfOriginal.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"stiffness perturbations","stiffnessChangedToFractionOfOriginal(double):");
			return false;
		}
		Sim->startedStiffnessPerturbation.push_back(false);
	}
	return true;
}

bool ModelInputObject::readECMPerturbation(ifstream& file){
	/**
	 * Will change
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsECMStiffnessChange(bool):"){
		file >> Sim->thereIsECMChange;
	}
	else{
		cerr<<"Error in reading ECM perturbations, curr string: "<<currHeader<<", should have been: ThereIsECMSoftening(bool):" <<endl;
		return false;
	}

	int nECMFunctions;
	file >> currHeader;
	if(currHeader == "NumberOfECMPerturbations(int):"){
		file >>nECMFunctions;
	}
	else{
		printErrorMessage(currHeader,"ECM perturbations","NumberOfECMPerturbations(int):");
		return false;
	}
	for (int i=0; i<nECMFunctions; ++i){
		Sim->changedECM.push_back(false);
		file >> currHeader;
		if(currHeader == "ApplyToApicalECM(bool):"){
			bool changeApicalECM;
			file >> changeApicalECM;
			Sim->changeApicalECM.push_back(changeApicalECM);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ApplyToApicalECM(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToBasalECM(bool):"){
			bool changeBasalECM;
			file >> changeBasalECM;
			Sim->changeBasalECM.push_back(changeBasalECM);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ApplyToBasalECM(bool):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "AppliedElementsAreEmergent(bool):"){
			bool emergentApplication;
			file >> emergentApplication;
			Sim->ECMChangeTypeIsEmergent.push_back(emergentApplication);
			if (emergentApplication){
				Sim->thereIsEmergentEllipseMarking = true;
			}
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","AppliedElementsAreEmergent(bool):");
			return false;
		}

		file >> currHeader;
		if(currHeader == "timeOfStiffnessChange(hr):"){
			double timeInHr;
			file >> timeInHr;
			Sim->ECMChangeBeginTimeInSec.push_back(timeInHr*3600);
			file >> timeInHr;
			Sim->ECMChangeEndTimeInSec.push_back(timeInHr*3600);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","timeOfStiffnessChange(hr):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangeAppliedToEllipses(number,[ellipseId][ellipseId]):"){
			int numberOfECMChangeEllipseBands;
			file >> numberOfECMChangeEllipseBands;
			Sim->numberOfECMChangeEllipseBands.push_back(numberOfECMChangeEllipseBands);
			int ellipseBandId;
			Sim->ECMChangeEllipseBandIds.push_back(vector<int>(0));
			for (int aa=0; aa<Sim->numberOfECMChangeEllipseBands[i]; ++aa){
				file >>ellipseBandId;
				Sim->ECMChangeEllipseBandIds[i].push_back(ellipseBandId);
				if (ellipseBandId >= 100){
					Sim->thereIsEmergentEllipseMarking = true;
				}
			}
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","stiffnessChangeAppliedToEllipses(number,[ellipseId][ellipseId]):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "stiffnessChangeFraction(double(0-1.0)):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMStiffnessChangeFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","stiffnessChangeFraction(double(0-1.0)):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ECMRenewalHalfLifeTargetFraction(double(0-1.0)):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMRenewalHalfLifeTargetFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ECMRenewalHalfLifeTargetFraction(double(0-1.0)):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "ECMViscosityChangeFraction(double):"){
			double fraction;
			file >> fraction;
			if (fraction <=0.0) {
				fraction = 0.0001;
			}
			Sim->ECMViscosityChangeFraction.push_back(fraction);
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","ECMViscosityChangeFraction(double):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changeNotumECM(time,fraction):"){
			file >> Sim->notumECMChangeInitTime;
			Sim->notumECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->notumECMChangeEndTime;
			Sim->notumECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->notumECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changeHingeECM(time,fraction):"){
			file >> Sim->hingeECMChangeInitTime;
			Sim->hingeECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->hingeECMChangeEndTime;
			Sim->hingeECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->hingeECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
		file >> currHeader;
		if(currHeader == "changePouchECM(time,fraction):"){
			file >> Sim->pouchECMChangeInitTime;
			Sim->pouchECMChangeInitTime*=3600; //convet to hrs
			file >> Sim->pouchECMChangeEndTime;
			Sim->pouchECMChangeEndTime*=3600; //convet to hrs
			file >> Sim->pouchECMChangeFraction;
		}
		else{
			printErrorMessage(currHeader,"ECM perturbations","changeNotumECM(time,fraction):");
			return false;
		}
	}
	return true;
}

bool ModelInputObject::readCellMigrationOptions(ifstream& file){
	/**
	 * Will delete
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsCellMigration(bool):"){
		file >> Sim->thereIsCellMigration;
	}
	else{
		cerr<<"Error in reading cell migration options: "<<currHeader<<", should have been: ThereIsCellMigration(bool):" <<endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readExplicitActinOptions(ifstream& file){
	/**
	 * This setting will declare the apical surface of the columnar layer to be
	 * the actin layer. The decleration does not change the physical properties of the actin layer.
	 * It changes how the layer is treated in growth functions, the actin later
	 *
	 * ThereIsExplicitActin(bool): 1
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsExplicitActin(bool):"){
		file >> Sim->thereIsExplicitActin;
	}
	else{
		cerr<<"Error in reading explicit actin options: "<<currHeader<<", should have been: ThereIsExplicitActin(bool):" <<endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readColumnViseVolumeConservationOptions(ifstream& file){
	/**
	 * The boolean sets the parameter Simulation#conservingColumnVolumes
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsColumnViseVolumeConservation(bool):"){
		file >> Sim->conservingColumnVolumes;
	}
	else{
		cerr<<"Error in reading column-vise volume conservation options: "<<currHeader<<", should have been: ThereIsColumnViseVolumeConservation(bool):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readLumenOptions(ifstream& file){
	/**
	 * The first boolean states the
	 *
	 * LumenOptions:
  	 *   thereIsLumen(bool): 1
  	 *   LumenBulkModulus(Pa): 32000
  	 *   LumenGrowthRate(foldPer24hr): 2.0
	 *
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "thereIsLumen(bool):"){
		file >> Sim->thereIsExplicitLumen;
	}
	else{
		printErrorMessage(currHeader,"Lumen options","thereIsLumen(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LumenBulkModulus(Pa):"){
		file >> Sim->lumenBulkModulus;
	}
	else{
		printErrorMessage(currHeader,"Lumen options","LumenBulkModulus(Pa)");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LumenGrowthRate(foldPer24hr):"){
		file >> Sim->lumenGrowthFold;
	}
	else{
		printErrorMessage(currHeader,"Lumen options","LumenGrowthRate(foldPer24hr):");
		return false;
	}
	return true;
}

bool ModelInputObject::readartificialRelaxationOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsArtificaialRelaxation(bool):"){
		file >> Sim->thereIsArtificaialRelaxation;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","ThereIsArtificaialRelaxation(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ArtificialRelaxationTime(sec):"){
		file >> Sim->artificialRelaxationTime;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","ArtificialRelaxationTime(sec):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "relaxECM(bool):"){
		file >> Sim->relaxECMInArtificialRelaxation;
	}
	else{
		printErrorMessage(currHeader,"artificial relaxation options","relaxECM(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readEnclosementOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "thereIsEnclosementOfTheTissue(bool):"){
		file >> Sim->encloseTissueBetweenSurfaces;
	}
	else{
		cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: thereIsEnclosementOfTheTissue(bool):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "initialLimits(lowerBound,upperBound):"){
		file >> Sim->initialZEnclosementBoundaries[0];
		file >> Sim->initialZEnclosementBoundaries[1];
	}
	else{
		cerr<<"Error in reading enclosementoptions: "<<currHeader<<", should have been: initialLimits(lowerBound,upperBound):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalLimits(lowerBound,upperBound):"){
		file >> Sim->finalZEnclosementBoundaries[0];
		file >> Sim->finalZEnclosementBoundaries[1];
	}
	else{
		cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: finalLimits(lowerBound,upperBound):" <<endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "initialTime(sec):"){
		file >> Sim->initialTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
		cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: initialTime(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalTime(sec):"){
		file >> Sim->finalTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
		cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: finalTime(sec):" <<endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readMutationOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "numberOfClones(int):"){
		file >> Sim->numberOfClones;
	}
	else{
		printErrorMessage(currHeader,"Mutation options","numberOfClones(int):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):"){
		double relativeX, relativeY, micronRadius, growthRateORFold,useAbsoluteGrowthRate;
		for (int i=0; i<Sim->numberOfClones; ++i){
			file >> relativeX;
			file >> relativeY;
			file >> micronRadius;
			file >> useAbsoluteGrowthRate;
			file >> growthRateORFold;
			Sim->cloneInformationX.push_back(relativeX);
			Sim->cloneInformationY.push_back(relativeY);
			Sim->cloneInformationR.push_back(micronRadius);
			Sim->cloneInformationUsingAbsolueGrowth.push_back(useAbsoluteGrowthRate);
			Sim->cloneInformationGrowth.push_back(growthRateORFold);
		}
	}
	else{
		printErrorMessage(currHeader,"Mutation options","cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):");
		return false;
	}
	return true;
}

bool ModelInputObject::readAdhesionOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsAdhesion(bool):"){
		file >> Sim->thereIsAdhesion;
	}
	else{
		printErrorMessage(currHeader,"adhesion options","ThereIsAdhesion(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "CollapseNodesOnAdhesion(bool):"){
		file >> Sim->collapseNodesOnAdhesion;
	}
	else{
		printErrorMessage(currHeader,"adhesion options","CollapseNodesOnAdhesion(bool):");
		return false;
	}
	return true;
}

bool ModelInputObject::readNodeCollapseOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsNodeCollapse(bool):"){
		file >> Sim->thereNodeCollapsing;
	}
	else{
		printErrorMessage(currHeader,"node collapse options","ThereIsNodeCollapse(bool):");
		return false;
	}
	return true;
}


bool ModelInputObject::readExplicitECMOptions(ifstream& file){
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsExplicitECM(bool):"){
		file >> Sim->thereIsExplicitECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ThereIsExplicitECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "AddLateralECM(bool):"){
		file >> Sim->addLateralECMManually;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","AddLateralECM(bool):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "LateralECMThickness(microns):"){
		file >> Sim->lateralECMThickness;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","LateralECMThickness(micron):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMRemodellingHalfLife(hour):"){
		file >> Sim->ECMRenawalHalfLife;
		Sim->ECMRenawalHalfLife *= 3600; //converting to seconds.
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMRemodellingHalfLife(hour):");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMColumnarYoungsModulus:"){
		file >> Sim->EColumnarECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMColumnarYoungsModulus:");
		return false;
	}
	file >> currHeader;
	if(currHeader == "ECMPeripodialYoungsModulus:"){
		file >> Sim->EPeripodialECM;
	}
	else{
		printErrorMessage(currHeader,"explicit ECM options","ECMPeripodialYoungsModulus:");
		return false;
	}
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
	if(currHeader == "NumberOfPressureStages(int):"){
		file >> Sim->nPipetteSuctionSteps;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  NumberOfPressureStages(int):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InitiationTimes(sec):"){
		for (int i=0;i<Sim->nPipetteSuctionSteps;++i){
			double pressureInitiationTime;
			file >> pressureInitiationTime;
			Sim->pipetteSuctionTimes.push_back(pressureInitiationTime);
		}
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  InitiationTimes(sec):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pressures(Pa):"){
		for (int i=0;i<Sim->nPipetteSuctionSteps;++i){
			double pipetePressure;
			file >> pipetePressure;
			Sim->pipetteSuctionPressures.push_back(pipetePressure);
		}
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  Pressures(Pa):" <<endl;
		return false;
	}
	if (Sim->nPipetteSuctionSteps>0){
		Sim->PipetteInitialStep = Sim->pipetteSuctionTimes[0]/Sim->dt;
	}
	file >> currHeader;
	if(currHeader == "ApicalSuction(bool-will_set_up_basal_suction_if_false):"){
		file >> Sim->ApicalSuction;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: ApicalSuction(bool-will_set_up_basal_suction_if_false):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "TissueStuck(bool-will_fix_the_opposite_surface_in_z):"){
		file >> Sim->TissueStuckOnGlassDuringPipetteAspiration;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: TissueStuck(bool-will_fix_the_opposite_surface_in_z):" <<endl;
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
	if(currHeader == "Pipette_InnerRadius(micron):"){
		file >> Sim->pipetteInnerRadius;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_InnerRadius(micron):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_OuterRadius(micron):"){
		double pippetOuterRad;
		file >> pippetOuterRad;
		Sim->pipetteThickness = pippetOuterRad - Sim->pipetteInnerRadius;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_OuterRadius(micron):" <<endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Effect_Depth(micron):"){
		double dummy;
		file >> dummy;
		Sim->pipetteDepth = dummy;
	}
	else{
		cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Effect_Depth(micron):" <<endl;
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
