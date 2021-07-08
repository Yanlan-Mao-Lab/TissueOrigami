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
#include <memory>

using namespace std;

ModelInputObject::ModelInputObject(){
	meshFileName = "NA";
};

ModelInputObject::~ModelInputObject(){
	//std::cout<<"called the destructor for ModelInputObject class"<<std::endl;
	//std::cout<<"finalised the destructor for ModelInputObject class"<<std::endl;
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
            else if(currParameterHeader == "YoungsModulusTimeseriesGrids:"){
                /**
                 * Timeseries of the Young's modulus of the tissue  through the private function ModelInputObject#readYoungsModulusTimeseriesGrids
                 */
                //std::cout<<"reading YoungsModulusTimeseriesGrids:"<<std::endl;
                Success  = readYoungsModulusTimeseriesGrids(parametersFile);
             }
             else if(currParameterHeader == "TypeOfCoordinateSystem:"){
                /**
                 * Inputs defining the type of coordinate system
                 */
                Success  = readTypeOfCoordinateSystem(parametersFile);
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
				std::cerr<<"Unidentified parameter input line: "<<std::endl;
				std::cerr<<"		"<<currParameterHeader<<std::endl;
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
		std::cout<<"Cannot open parameter input file, "<<fileName<<std::endl;
		return false;
	}
	if (!file.good()) {
		std::cout<<"File does not exist, "<<fileName<<std::endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readGrowthOptions(ifstream& file){
	/**
	 * This function will read the growth options. The generic growth related parameters will be read here, followed by
	 * growth type specific parameters read in specialised functions.
	 * The first parameter to read will be the number of growth options.
     *
     *
     *GrowthOptions:
     *   NumberofGrowthFunctions(int): 5
     *   GridGrowthsPinnedOnInitialMesh(bool): 1
     *   PinningUpdateTimes(number-times(sec)):  2 57600 115200
     *   GridGrowthsInterpolationType(0=step,1=linear): 1
     *
     *   then continue with the each growth option inputs
	 */
	string currHeader;
	file >> currHeader;
	int n;
	if(currHeader == "NumberofGrowthFunctions(int):"){
		file >> n;
		Sim->nGrowthFunctions = n;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: NumberofGrowthFunctions(int):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: GridGrowthsPinnedOnInitialMesh(bool):" <<std::endl;
		return false;
	}
	/**
	 * Then the following set of parameters defines at which time points the pinning is updated, and the timing of those updates, in seconds.
	 */
	file >> currHeader;
	if(currHeader == "PinningUpdateTimes(number-times(sec)):"){
		file >> Sim->nGrowthPinning;
		if ( Sim->nGrowthPinning > 0){
			for (int j=0; j<Sim->nGrowthPinning; ++j){
				int currentUpdateTime;
				file >> currentUpdateTime;
				Sim->growthPinUpdateTime.push_back(currentUpdateTime);
				Sim->growthPinUpdateBools.push_back(false);
			}
		}
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: PinningUpdateTimes(number-times(sec)):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: GridGrowthsInterpolationType(0=step,1=linear):" <<std::endl;
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
			std::cerr<<"Error in reading growth type, curr string: "<<currHeader<<", should have been: GrowthFunctionType(int-seeDocumentation):" <<std::endl;
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
			std::cerr<<"Error in reading growth type, please enter a valid type: {1, 2, 3}, current type: "<<type<<std::endl;
			return false;
		}
		if (!Success){
			return false;
		}
	}
	std::cout<<"Finalised reading growth options"<<std::endl;
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
	std::cout<<"entered read growth type 1, current header: "<<currHeader<<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degrees):" <<std::endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	//GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 1
	//GSBp = new UniformGrowthFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, DVRate, APRate,  ABRate, angle);
	//Sim->GrowthFunctions.push_back(GSBp);
    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<UniformGrowthFunction>(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, DVRate, APRate,  ABRate, angle);
    Sim->GrowthFunctions.push_back(std::move(GSBp));
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
	std::cout<<"entered read growth type 2, current header: "<<currHeader<<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Centre:" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Radius:" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-DV,AP,AB):" <<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Angle(degrees):" <<std::endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	//GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 2
	//GSBp = new RingGrowthFunction(Id, 2, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, CentreX, CentreY, innerR,  outerR, DVRate, APRate,  ABRate, angle);
	//Sim->GrowthFunctions.push_back(GSBp);

    std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<RingGrowthFunction>(Id, 2, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, CentreX, CentreY, innerR,  outerR, DVRate, APRate,  ABRate, angle);
    Sim->GrowthFunctions.push_back(std::move(GSBp));
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
	std::vector<std::vector<std::array<double,3>>> GrowthMatrix; //[grid_i][grid_j][x,y,z]
	std::vector<std::vector<double>>  AngleMatrix; //[grid_i][grid_j][tetha]
	if(currHeader == "InitialTime(sec):"){
		file >> initialtime;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalECM(bool):"){
			file >> applyToBasalECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToBasalECM(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToLateralECM(bool):"){
			file >> applyToLateralECM;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: applyToLateralECM(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Filename(full-path):"){
		/**
		 * The file path should be a readable path from where the executable is lounced. It can handle relative paths, but full paths are encouraged.
		 */
		string filepath;
		file >> filepath;
		std::cerr<<" filename is: "<<filepath<<std::endl;
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
			std::cerr<<"could not open growth rate file file: "<<name_growthRates<<std::endl;
			return false;
		}
		//adding the indice of the growth matrix
		/**
		 * The growth grid will contain the size of the grid matrix as the first two
		 * parameters, first these are read.
		 */
		std::cout<<"reading from growth file"<<std::endl;
		GrowthRateFile >> gridX;
		GrowthRateFile >> gridY;
		float rate;
        std::cout<<"initiating growth and angle matrices"<<std::endl;
        /**
         * Then the stacked arrays are generated for the growth rate and the growth orientation matrices, if I have time, I will convert these to vectors, which are
         * much more cleaner then new double arrays. Growth matrix has three layers, x & y coordinates on the bounding box and the 3D growth.
         * Orientation angles have 2 laters, x & y directions only, followed by a single angle value.
         * All initiated as zeros.
         */
        for (int i= 0; i<gridX; ++i){
            std::vector<std::array<double,3>> tmpGridMatrixY(gridY,std::array<double,3>{0.0});
            GrowthMatrix.push_back(tmpGridMatrixY);
            std::vector<double> tmpShearY(gridY,0.0);
            AngleMatrix.push_back(tmpShearY);
		}
		double timeMultiplier = 1.0 / 3600.0; // converting rate per hour to rate per second
		/**
		 * Then the growth rates and angles are read in.
		 */
		std::cout<<"reading growth matrix"<<std::endl;
		for (int j=gridY-1; j>-1; --j){
			for (int i=0; i<gridX; ++i){
				for (int k=0; k<3; ++k){
					//std::cout<<"i :"<<i<<" j: "<<j<<" k: "<<k<<" ";
					GrowthRateFile >> rate;
					//std::cout<<"rate: "<<rate<<" ";
					//GrowthMatrix[i][j][k] = rate*timeMultiplier;
					GrowthMatrix[i][j][k] = rate;
					//std::cout<<"matrix value: "<<GrowthMatrix[i][j][k]<<std::endl;
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: Filename(full-path):" <<std::endl;
		return false;
	}
	/**
	 * The normalised z range of a grid based growth function defined to which height of the tissue the growth function should be applied in.
	 * Note that this is on top of the booleans specifying which components of the tissue the growth function will be applied to.
         * In the code, the z starts from the basal surface of columnar. However, the zRange that is reported here is the relative z and is starts from top, so z=1.0 for basal columnar surface and z=0.0 is either the apical surface of the columnar layer
         * or the basal surface of the peripodial layer (if there is one).
	 */
	double zMin, zMax;
	file >> currHeader;
	if(currHeader == "zRange:"){
		file >> zMin;
		file >> zMax;
	}
	else{
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: zRange:" <<std::endl;
		return false;
	}
	/**
	 * Once all input parameters are gathered, a new growth function object is initiated with the constructor of UniformGrowthFunction# class.
	 * The pointer to the generated object is recorded in the vector of Growth function objects that the simulation object keeps, Simulation#GrowthFunctions.
	 */
	//GrowthFunctionBase* GSBp;
	int Id = Sim->GrowthFunctions.size();
	//type is 3
	//GSBp = new GridBasedGrowthFunction(Id, 3, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, gridX, gridY, GrowthMatrix, AngleMatrix);
	//GSBp->zMin = zMin;
	//GSBp->zMax = zMax;

	//Sim->GrowthFunctions.push_back(GSBp);
	std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<GridBasedGrowthFunction>(Id, 3, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, applyToBasalECM, applyToLateralECM, gridX, gridY, GrowthMatrix, AngleMatrix);
	GSBp->zMin = zMin;
	GSBp->zMax = zMax;
	Sim->GrowthFunctions.push_back(std::move(GSBp));
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
		std::cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: MeshInputMode(seeDocumentation):" <<std::endl;
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
		std::cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: symmetricInX(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "symmetricInY(bool):"){
		file >> Sim->symmetricY;
	}
	else{
		std::cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: symmetricInY(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "symmetricInZ(bool):"){
		file >> Sim->symmetricZ;
	}
	else{
		std::cerr<<"Error in reading mesh type, curr string: "<<currHeader<<", should have been: symmetricInZ(bool):" <<std::endl;
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
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperApicalExternalViscosity:" <<std::endl;
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
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperApicalExternalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DiscProperBasalExternalViscosity:"){
		file >>Sim->externalViscosityDPBasal;
	}
	else{
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: DiscProperBasalExternalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalExternalViscosity:"){
		file >>Sim->externalViscosityPMApical;
	}
	else{
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: PeripodialMembraneApicalExternalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalExternalViscosity:"){
		file >>Sim->externalViscosityPMBasal;
	}
	else{
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: PeripodialMembraneBasalExternalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalExternalViscosity:"){
		file >>Sim->externalViscosityLZApical;
	}
	else{
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerZoneApicalExternalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalExternalViscosity:"){
		file >>Sim->externalViscosityLZBasal;
	}
	else{
		std::cerr<<"Error in reading Fixing option, curr string: "<<currHeader<<", should have been: LinkerZoneBasalExternalViscosity:" <<std::endl;
		return false;
	}
	return true;
}

void ModelInputObject::printErrorMessage(string currentInput, string sourceFuction, string expectedInput){
	std::cerr<<"Error in reading "<<sourceFuction<<" current input: "<<currentInput<<", should have been: "<<expectedInput<<std::endl;
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
     *   FixingViscosity(x,y,z): 0   0  32000
     *   ApicSurfaceFix(bool-x,y,z):   0 0 0   FixApicalExtVisc(bool): 0
     *   BasalSurfaceFix(bool-x,y,z):  0 0 0   FixBasalExtVisc(bool):  0
     *   CircumferenceFix(bool-x,y,z): 0 0 0   FixCircWithExtVisc(bool): 0
     *   ApicCircumFix(bool-x,y,z):    0 0 0   FixApicCircWithExtVisc(bool):  0
     *   BasalCircumFix(bool-x,y,z):   0 0 0   FixBasalCircWithExtVisc(bool): 0
     *   LinkerApicCircumFix(bool-x,y,z):  0 0 0  FixLinkerApicCircWithExtVisc(bool):  0
     *   LinkerBasalCircumFix(bool-x,y,z): 0 0 0  FixLinkerBasalCircWithExtVisc(bool): 0
     *   NotumFix(bool-x,y,z,double-xFracMin,xFracMax): 0 0 0 -0.1 0.5  FixNotumExtVisc(bool): 0
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
     * manipulations on the simulation options, from stiffness perturbations to random forces. If you are running a simulation with peripodial, you need to copy this into the model input file.
	 * I find it highly unlikely you will use these, probably will remove them soon.
	 *
	 * Manipulations:
     *   AddCurvature(bool): 1
     *   CurvatureDepthAtCentre(double-microns): 2.0
     *   AddSoftPeriphery(bool): 0
     *   SoftPeripheryRange(double-microns): 30.0
     *   SoftnessFraction(double-fraction): 0.1
     *   ApplyToApicalSurface(bool): 1
     *   ApplyToBasalSurface(bool): 0
     *   ApplyToColumnarLayer(bool): 1
     *   ApplyToPeripodialMembrane(bool): 1
     *   AddRandomForce(bool): 0
     *   RandomForceMean(double): 0.0
     *   RandomForceVar(double): 1E-5
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
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: CurvatureDepthAtCentre(double-microns):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "AddSoftPeriphery(bool):"){
		file >>Sim->softPeriphery;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: AddSoftPeriphery(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SoftPeripheryRange(double-microns):"){
		file >>Sim->softDepth;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: SoftPeripheryRange(double-microns):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SoftnessFraction(double-fraction):"){
		file >>Sim->softnessFraction;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: SoftnessFraction(double-fraction):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToApicalSurface(bool):"){
		file >>Sim->softPeripheryBooleans[0];
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToApicalSurface(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToBasalSurface(bool):"){
		file >>Sim->softPeripheryBooleans[1];
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToBasalSurface(bool):" <<std::endl;
		return false;
	}	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
		file >>Sim->softPeripheryBooleans[2];
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
		return false;
	}	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
		file >>Sim->softPeripheryBooleans[3];
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
		return false;
	}
	//Reading Random Force Parameters
	file >> currHeader;
	if(currHeader == "AddRandomForce(bool):"){
		file >>Sim->addingRandomForces;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: AddRandomForce(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "RandomForceMean(double):"){
		file >>Sim->randomForceMean;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: RandomForceMean(double):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "RandomForceVar(double):"){
		file >>Sim->randomForceVar;
	}
	else{
		std::cerr<<"Error in reading manipulations options, curr string: "<<currHeader<<", should have been: RandomForceVar(double):" <<std::endl;
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
		std::cerr<<"Error in reading mesh path, curr string: "<<currHeader<<", should have been: MeshFile(full-path):" <<std::endl;
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
		std::cerr<<"Error in reading mesh row number, curr string: "<<currHeader<<", should have been: MeshRow(int):" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "MeshColumn(int):"){
		file >> Sim->Column;
	}
	else{
		std::cerr<<"Error in reading mesh column number, curr string: "<<currHeader<<", should have been: MeshColumn(int):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "SideLength:"){
		file >> Sim->SideLength;
	}
	else{
		std::cerr<<"Error in reading side length, curr string: "<<currHeader<<", should have been: SideLength:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "zHeight:"){
		file >> Sim->zHeight;
	}
	else{
		std::cerr<<"Error in reading z height, curr string: "<<currHeader<<", should have been: zHeight:" <<std::endl;
		return false;
	}
	//checking consistency:
	if (Sim->Column>Sim->Row-2){
		Sim->Column = Sim->Row-2;
		Sim->outputFile<<"Too few rows vs. columns, column count cannot be higher than Row-2"<<std::endl;
		Sim->outputFile<<"Updated to have a mesh: Row: "<<Sim->Row<<" Column: "<<Sim->Column<<std::endl;
	}
	float aspectratio = Sim->zHeight/Sim->SideLength;
	if ( aspectratio > 10 || aspectratio < 0.01 ){
		Sim->outputFile<<"Warning: The aspect ratio of the shapes are too high or low (aspectratio (z/side): "<<aspectratio<<std::endl;
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
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: AddPeripodialMembrane:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneThickness(fractionOfTissueHeight):"){
		file >>Sim->PeripodialThicnessScale;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneThickness(fractionOfTissueHeight):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneLateralThickness(fractionOfTissueHeight):"){
		file >>Sim->PeripodialLateralThicnessScale;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneLateralThickness(fractionOfTissueHeight):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LumenHeightScale(fractionOfTissueHeight):"){
		file >>Sim->lumenHeightScale;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LumenHeightScale(fractionOfTissueHeight):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneYoungsModulus:"){
		file >>Sim->PeripodialElasticity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneYoungsModulus:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneApicalViscosity:"){
		file >>Sim->peripodialApicalViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneApicalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneBasalViscosity:"){
		file >>Sim->peripodialBasalViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneBasalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "PeripodialMembraneMidlineViscosity:"){
		file >>Sim->peripodialMidlineViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: PeripodialMembraneMidlineViscosity:" <<std::endl;
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
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: BaseOnPeripodialness(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalYoungsModulus:"){
		file >>Sim->LinkerZoneApicalElasticity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneApicalYoungsModulus:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalYoungsModulus:"){
		file >>Sim->LinkerZoneBasalYoungsModulus;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneBasalYoungsModulus:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneApicalViscosity:"){
		file >>Sim->linkerZoneApicalViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneApicalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneBasalViscosity:"){
		file >>Sim->linkerZoneBasalViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneBasalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "LinkerZoneMidlineViscosity:"){
		file >>Sim->linkerZoneMidlineViscosity;
	}
	else{
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: LinkerZoneMidlineViscosity:" <<std::endl;
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
		std::cerr<<"Error in reading time step, curr string: "<<currHeader<<" should have been: TimeStep(sec):" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "SimulationLength(sec):"){
		file >> Sim->SimLength;
	}
	else{
		std::cerr<<"Error in reading simulation length, curr string: "<<currHeader<<" should have been: SimulationLength(sec)::" <<std::endl;
		return false;
	}

	std::cout<<"Simulation time step	: "<<Sim->dt<<std::endl;
	std::cout<<"Simulation Length	: "<<Sim->SimLength<<std::endl;
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
		std::cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusApical:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusBasal:"){
		file >> Sim->EBasal;
	}
	else{
		std::cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusBasal:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "YoungsModulusMid:"){
		file >> Sim->EMid;
	}
	else{
		std::cerr<<"Error in reading Young's modulus, curr string: "<<currHeader<<" should have been: YoungsModulusMid:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[0];
	}
	else{
		std::cerr<<"Error in reading Young's modulus noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "PoissonsRatio:"){
		file >> Sim->poisson;
	}
	else{
		std::cerr<<"Error in reading Poisson's ratio, curr string: "<<currHeader<<" should have been: PoissonsRatio:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[1];
	}
	else{
		std::cerr<<"Error in reading Poisson's Ratio noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "ApicalViscosity:"){
		file >> Sim->discProperApicalViscosity;
	}
	else{
		std::cerr<<"Error in reading Apical Viscosity, curr string: "<<currHeader<<" should have been: ApicalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Noise(%-int):"){
		file >> Sim->noiseOnPysProp[2];
	}
	else{
		std::cerr<<"Error in reading Viscosity noise, curr string: "<<currHeader<<" should have been: Noise(%-int)" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "BasalViscosity:"){
		file >> Sim->discProperBasalViscosity;
	}
	else{
		std::cerr<<"Error in reading Basal Viscosity, curr string: "<<currHeader<<" should have been: BasalViscosity:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MidLineViscosity:"){
		file >> Sim->discProperMidlineViscosity;
	}
	else{
		std::cerr<<"Error in reading Midline Viscosity, curr string: "<<currHeader<<" should have been: MidLineViscosity:" <<std::endl;
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
	//std::cout<<"reading save options"<<std::endl;
	string currHeader;
	file >> currHeader;
	if(currHeader == "SaveImages(bool):"){
		file >> Sim->saveImages;
	}
	else{
		std::cerr<<"Error in reading image saving option, current string: "<<currHeader<<" should have been: SaveImages(bool):" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "SaveData(bool):"){
		file >> Sim->saveData;
	}
	else{
		std::cerr<<"Error in reading simulation data saving option, current string: "<<currHeader<<" should have been: SaveData(bool):" <<std::endl;
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
		std::cerr<<"Error in reading image save interval, current string: "<<currHeader<<" should have been: ImageSaveInterval(sec)" <<std::endl;
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
        std::cout<<"dataSaveInterval as read from file: "<<Sim->dataSaveInterval<<std::endl;
	}
	else{
		std::cerr<<"Error in reading data save interval, current string: "<<currHeader<<" should have been: DataSaveInterval(sec)" <<std::endl;
		return false;
	}
	//std::cout<<"Sim->saveImages: "<<Sim->saveImages<<" Sim->saveData: "<<Sim->saveData<<" "<<" datainterval: "<<Sim->dataSaveInterval<<" imageinterval: "<<Sim->dataSaveInterval<<std::endl;
	return true;
}

bool ModelInputObject::readYoungsModulusTimeseriesGrids(ifstream& file){

    size_t noTimeseriesInputs;
    string currHeader;
    file >> currHeader;
    if(currHeader == "noTimeseriesInputs(size_t):"){
        file >> noTimeseriesInputs;
    }
    else{
        cerr<<"Error in reading number of timeseries inputs: "<<currHeader<<" should have been: noTimeseriesInputs(size_t):" <<endl;
        return false;
    }
    for (size_t i=0;i<noTimeseriesInputs;++i){
        size_t TimePointNoInput;
        bool YoungsModulusMultiplierChangeRateIsLinearInput;
        std::array<bool,7> WhereToApplyInput;
        std::vector<double> WhenToApplyInput;
        std::vector<std::string> nameListInput;

        file >> currHeader;
        if(currHeader == "TimePointNoInput(size_t):"){
            file >> TimePointNoInput;
        }
        else{
            cerr<<"Error in reading number of timepoints in timeseries: "<<currHeader<<" should have been: TimePointNoInput(size_t):" <<endl;
            return false;
        }

        file >> currHeader;
        if(currHeader == "YoungsModulusMultiplierChangeRateIsLinearInput(bool):"){
            file >> YoungsModulusMultiplierChangeRateIsLinearInput;
        }
        else{
            printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","YoungsModulusMultiplierChangeRateIsLinearInput(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToColumnarLayer(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[0]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToColumnarLayer(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToPeripodialMembrane(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[1]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToPeripodialMembrane(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToApicalLayer(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[2]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToApicalLayer(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToMidLayer(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[3]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToMidLayer(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToBasalLayer(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[4]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToBasalLayer(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToPeripodialECM(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[5]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToPeripodialECM(bool):");
            return false;
        }

        file >> currHeader;
        if(currHeader == "applyToColumnarECM(bool):"){
            double currApplyToLayer;
            file >> currApplyToLayer;
            WhereToApplyInput[6]=currApplyToLayer;
        }
        else{
             printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","applyToColumnarECM(bool):");
            return false;
        }

        for (size_t j=0;j<TimePointNoInput;++j){
            file >> currHeader;
            if(currHeader == "Filename(full-path):"){
                std::string currnameListInput;
                file >> currnameListInput;
                nameListInput.push_back(currnameListInput);
            }
            else{
                printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","Filename(full-path):");
                return false;
            }

            file >> currHeader;
            if(currHeader == "WhenToApplyInput(sec):"){
                double currWhenToApplyInput;
                file >> currWhenToApplyInput;
                WhenToApplyInput.push_back(currWhenToApplyInput);
            }
            else{
                 printErrorMessage(currHeader,"readYoungsModulusTimeseriesGrids","WhenToApplyInput(sec):");
                return false;
            }
        }
       std::unique_ptr<YoungsModulusModifier> ym01 = std::make_unique<YoungsModulusModifier>(TimePointNoInput,nameListInput,WhenToApplyInput,WhereToApplyInput,YoungsModulusMultiplierChangeRateIsLinearInput);
       Sim->AllYoungsModulusModifiers.push_back(std::move(ym01));
    }
    return true;
}

bool ModelInputObject::readTypeOfCoordinateSystem(ifstream& file){
    /**
      * Sample:
      * TypeOfCoordinateSystem:
      *     UseXYCoordinatesforRelativePositions(bool): 1
      *     UsePolarCoordinatesforRelativePositions(bool): 0
    **/

    string currHeader;
    file >> currHeader;
    if(currHeader == "UseXYCoordinatesforRelativePositions(bool):"){
        file >> Sim->UseXYCoordinatesforRelativePositions;
    }
    else{
        printErrorMessage(currHeader,"defineTypeOfCoordinateSystem","UseXYCoordinatesforRelativePositions(bool):");
        return false;
    }

    file >> currHeader;
    if(currHeader == "UsePolarCoordinatesforRelativePositions(bool):"){
        file >> Sim->UsePolarCoordinatesforRelativePositions;
    }
    else{
        printErrorMessage(currHeader,"defineTypeOfCoordinateSystem","UsePolarCoordinatesforRelativePositions(bool):");
        return false;
    }

    file >> currHeader;
    if(currHeader == "UseCylindricalCoordintesforRelativePositions(bool):"){
        file >> Sim->UseCylindricalCoordintesforRelativePositions;
    }
    else{
        printErrorMessage(currHeader,"defineTypeOfCoordinateSystem","UseCylindricalCoordintesforRelativePositions(bool):");
        return false;
    }
    return true;
}

bool ModelInputObject::readShapeChangeOptions(ifstream& file){
	/**
     * Thiss function will change the shapes of the elements within the givven time range.
     * For applying to  elements of ellipse amrker 102 (the emergent folds), you should make the initial time 0
     * or any specific time if you want to restrict how early the shape change starts (sometimes folds can initiatw as early as 14-15 hrs (62 -63 hr AEL)
     *The shape chenage does not stop. It will continue till the end point in time.
     *If you are using emergent ellipses (100 or 102) this will create a considerable difference between elemetns that became part of a fol at different
     *times.
     *
     *How can you fix this? Give a new attribute initialTimeShapeChangeStarted to each element.
     * Check the time to stop as the difference between initial and final times, against the new
     * attribute: IF I said apply from time =6hr to 8 hrs, that is 2 hours of shape change.
     * If the shape change of a particular element started at 7 hrs, it should keep on for 2 hours,
     * until 9hrs.
     * YOU HAVE TO SAVE THAT ATTRIBUTE, as it will cause problems while conitnuing from saves!
     *
     * ShapeChangeOptions:
     *   NumberofShapeChangeFunctions(int): 1
     *   ShapeChangeStartsBelowECMLevel(fraction): 1.0

     *   ShapeChangeFunctionType(int-seeDocumentation): 2
     *   InitialTime(sec):  60000
     *   FinalTime(sec):  4200
     *   ApplyTissueApical(bool): 0
     *   ApplyTissueBasal(bool): 0
     *   ApplyTissueMidline(bool): 0
     *   ApplyToBasalECM(bool): 0
     *   ApplyToLateralECM(bool): 0
     *   ShapeChangeAppliedToEllipses(number,[ellipseId][ellipseId]): 1 102
     *   xyShapeChange(fractionPerHour): 0.00699
     *   ConserveVolume(bool): 1
     *
	 */
	string currHeader;
	file >> currHeader;
	int n;
	//std::cout<<"entered read shape change options, current header: "<<currHeader<<std::endl;
	if(currHeader == "NumberofShapeChangeFunctions(int):"){
		file >> n;
		Sim->nShapeChangeFunctions = n;
	}
	else{
		std::cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: NumberofShapeChangeFunctions(int):" <<std::endl;
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
		//std::cout<<"inside the loop, read shape change options, current header: "<<currHeader<<std::endl;
		if(currHeader == "ShapeChangeFunctionType(int-seeDocumentation):"){
			file >> type;
		}
		else{
			std::cerr<<"Error in reading shpae change type, curr string: "<<currHeader<<", should have been: ShapeChangeFunctionType(int-seeDocumentation):" <<std::endl;
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
			std::cerr<<"Error in reading shape change type, please enter a valid type: {1},{2} current type: "<<type<<std::endl;
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
		std::cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: ThereIsPlasticDeformation(bool):" <<std::endl;
		return false;
	}


	file >> currHeader;
		if(currHeader == "ApplyToColumnarLayer(bool):"){
				file >> Sim->plasticDeformationAppliedToColumnar;
		}
		else{
			std::cerr<<"Error in reading  plastic deformation options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
			return false;
		}
		file >> currHeader;
		if(currHeader == "ApplyToPeripodialMembrane(bool):"){
				file >> Sim->plasticDeformationAppliedToPeripodial;;
		}
		else{
			std::cerr<<"Error in reading  plastic deformation options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
			return false;
		}

	file >> currHeader;
	if(currHeader == "VolumeConserved(bool):"){
		file >> Sim->volumeConservedInPlasticDeformation;
	}
	else{
		std::cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: VolumeConserved(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DeformationHalfLife(hour):"){
		file >> Sim->plasticDeformationHalfLife;
		Sim->plasticDeformationHalfLife *= 3600; //converting to seconds.
	}
	else{
		std::cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: DeformationHalfLife(hour):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "zDeformationLimits(lowerFraction-upperFraction):"){
		file >> Sim->zRemodellingLowerThreshold;
		file >> Sim->zRemodellingUpperThreshold;
	}
	else{
		std::cerr<<"Error in reading plastic deformation options, curr string: "<<currHeader<<", should have been: zDeformationLimits(lowerFraction-upperFraction):" <<std::endl;
		return false;
	}
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
		std::cerr<<"Error in reading growth options, curr string: "<<currHeader<<", should have been: InitialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		file >> finaltime;
		//Sim->GrowthParameters.push_back(finaltime);
	}
	else{
		std::cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: FinalTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToColumnarLayer(bool):"){
			file >> applyToColumnarLayer;
	}
	else{
		std::cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: ApplyToColumnarLayer(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ApplyToPeripodialMembrane(bool):"){
			file >> applyToPeripodialMembrane;
	}
	else{
		std::cerr<<"Error in reading shape chenage options, curr string: "<<currHeader<<", should have been: ApplyToPeripodialMembrane(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxValue(fractionPerHour):"){
		double timeMultiplier = 1.0 / 3600.0;
		file >> Rate;
		Rate  *= timeMultiplier;
	}
	else{
		std::cerr<<"Error in reading shape change options, curr string: "<<currHeader<<", should have been: MaxValue(fractionPerHour-xyz):" <<std::endl;
		return false;
	}
	std::cerr<<"Shape change of type : 1, from time: "<<initialtime<<" to "<<finaltime<<" applicable to (C, P) "<<applyToColumnarLayer <<" "<<applyToPeripodialMembrane<<" Rate: "<<Rate<<std::endl;
	//GrowthFunctionBase* GSBp;
	//int Id = Sim->ShapeChangeFunctions.size();
	//type is 1
	//GSBp = new UniformShapeChangeFunction(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, false /*applyToBasalECM*/, false /*applyToLateralECM*/, 1, Rate);
	//Sim->ShapeChangeFunctions.push_back(GSBp);

    int Id = Sim->ShapeChangeFunctions.size();
    //type is 1
    unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<UniformShapeChangeFunction>(Id, 1, initialtime, finaltime, applyToColumnarLayer, applyToPeripodialMembrane, false /*applyToBasalECM*/, false /*applyToLateralECM*/, 1, Rate);
    Sim->ShapeChangeFunctions.push_back(std::move(GSBp));
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
	//GrowthFunctionBase* GSBp;
	//int Id = Sim->ShapeChangeFunctions.size();
	//type is 1
	//GSBp = new 	markerEllipseBasedShapeChangeFunction(Id, 2, initialtime, finaltime, applyTissueApical, applyTissueBasal, applyTissueMidline, applyToBasalECM, applyToLateralECM, 2, ShapeChangeFractionPerHr, markerEllipses, conserveVolume);
	//Sim->ShapeChangeFunctions.push_back(GSBp);

	int Id = Sim->ShapeChangeFunctions.size();
	//type is 2
	std::unique_ptr<GrowthFunctionBase> GSBp = std::make_unique<markerEllipseBasedShapeChangeFunction>(Id, 2, initialtime, finaltime, applyTissueApical, applyTissueBasal, applyTissueMidline, applyToBasalECM, applyToLateralECM, 2, ShapeChangeFractionPerHr, markerEllipses, conserveVolume);
	Sim->ShapeChangeFunctions.push_back(std::move(GSBp));

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
     * I would advise running the stretcher simulations with very small time steps, for short periods, with no external viscosity.
     * Then you can take the end point of the simulation, and continue with an actual larger time step simulation to observe the
     * dynamics of the ststem (such as ECM) post stretch. I would be cautious whie interpretting the actual dynamics ofthe stretching event,
     * and think carefully aboutwhat external and internal viscosities mean if I am to do so.
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
		//std::cerr<<"stretcherAttached "<<stretcherAttached<<std::endl;
		Sim->stretcherAttached = stretcherAttached;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: StretcherAttached(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "ClampedOnDV(bool):"){
		file >> Sim->DVClamp;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: ClampedOnDV(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "InitialTime(sec):"){
		double inittime;
		file >> inittime;
		Sim->StretchInitialTime = inittime;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: InitialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "FinalTime(sec):"){
		double endtime;
		file >> endtime;
		Sim->StretchEndTime = endtime;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: FinalTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMin:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMin = ClampPos;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMin:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "DVClampMax:"){
		double ClampPos;
		file >> ClampPos;
		Sim->StretchMax = ClampPos;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMax:" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "MaxStrain:"){
		double MaxStrain;
		file >> MaxStrain;
		Sim->StretchStrain = MaxStrain;
	}
	else{
		std::cerr<<"Error in reading stretcher setup: "<<currHeader<<", should have been: DVClampMax:" <<std::endl;
		return false;
	}
	//std::cout<<"StretcherAttached "<<Sim->stretcherAttached<<"InitialStep "<<Sim->StretchInitialStep<<" EndStep: "<<Sim->StretchEndStep<<" ClapmPos: "<<Sim->StretchMin<<" "<<Sim->StretchMax<<" strain: "<<Sim->StretchStrain<<std::endl;
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
		std::cerr<<"Error in reading marker ellipses, curr string: "<<currHeader<<", should have been: numberOfMarkerEllipses(int):" <<std::endl;
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
		std::cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseXCenters(fractionOfTissueSize):"<<std::endl;
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
		std::cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseBandR1Ranges(fractionOfTissueSize-x1Low-x1High):" <<std::endl;
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
		std::cerr<<"Error in reading  marker ellipses, curr string: "<<currHeader<<", should have been: MarkerEllipseBandR2Ranges(fractionOfTissueSize-y2Low-y2High):" <<std::endl;
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
     * We did not check if this works with the time series approach.
     *
     *
     * Stiffness_Perturbation:
     *   ThereIsStiffnessPerturbation(bool): 1
     *   NumberOfStiffnessPerturbations(int): 1

     *   ApplyToApically(bool): 1
     *   ApplyBasally(bool): 0
     *   ApplyToWholeTissue(bool): 0
     *   Basolateral(bool): 0
     *   BasolateralWithApicalRelaxation(bool): 0
     *   timeOfStiffeningPerturbation(hr): 30 32
     *   stiffnessPerturbationAppliedToEllipses(number,[ellipseId][ellipseId]): 1 0
     *   stiffnessChangedToFractionOfOriginal(double):  0.75
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsStiffnessPerturbation(bool):"){
		file >> Sim->ThereIsStiffnessPerturbation;
	}
	else{
		std::cerr<<"Error in reading stiffness perturbations, curr string: "<<currHeader<<", should have been: ThereIsStiffnessPerturbation(bool):" <<std::endl;
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
     * ECM change is applied in different ways.
     * First, it is explicitly applied to specified regions (Notum,hinge, pouch).
     * Then on regions specified with ellipses. This are usually the folds.
     * Non-emergent changes are applied between sepcieid timepoint (e.g. timeOfStiffnessChange).
     * For emergent changes, the code will check when the emergent change is applicable and the change will be applied any time after that.
     * The change will continue to be applied until the ChangeFraction is reached. Then it will stop.
     * Emergent ellipse id assignment is done through Simulation#checkForEmergentEllipseFormation.
     * Ellipse ids: 100 is apical node collapse. 101 is basal node collapse and 102 is when the ECM change has reached a threshold.
     * Emergent ECM change is not applied to ellipse id 101.
     * For ellipse id 100, emergent ECM change is applied. When a certain threshold is reach, ellipse id will change to 102 and other changes (e.g. cell shortening) start to emerge as well.
     *
     *ECM_Perturbation:
     *   ThereIsECMStiffnessChange(bool): 0
     *   NumberOfECMPerturbations(int): 1

     *   ApplyToApicalECM(bool): 0
     *    ApplyToBasalECM(bool): 0
     *    AppliedElementsAreEmergent(bool): 1
     *    timeOfStiffnessChange(hr): 20 26
     *    stiffnessChangeAppliedToEllipses(number,[ellipseId][ellipseId]): 1 100
     *    stiffnessChangeFraction(double(0-1.0)):  0.5
     *    ECMRenewalHalfLifeTargetFraction(double(0-1.0)): 1.0
     *    ECMViscosityChangeFraction(double): 0.5
     *    changeNotumECM(time,fraction): 1000 1000 1.0
     *    changeHingeECM(time,fraction): 1000 1000 1.0
     *    changePouchECM(time,fraction): 1000 1000 1.0
     *
	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsECMStiffnessChange(bool):"){
		file >> Sim->thereIsECMChange;
	}
	else{
		std::cerr<<"Error in reading ECM perturbations, curr string: "<<currHeader<<", should have been: ThereIsECMSoftening(bool):" <<std::endl;
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
			for (size_t aa=0; aa<Sim->numberOfECMChangeEllipseBands[i]; ++aa){
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
		std::cerr<<"Error in reading explicit actin options: "<<currHeader<<", should have been: ThereIsExplicitActin(bool):" <<std::endl;
		return false;
	}
	return true;
}


bool ModelInputObject::readColumnViseVolumeConservationOptions(ifstream& file){
	/**
     * The boolean sets the parameter Simulation#conservingColumnVolumes.
     * The functionality is working, but not thoroughly tested.
     *
     *ColumnViseVolumeConservationOptions:
     *  ThereIsColumnViseVolumeConservation(bool): 1

	 */
	string currHeader;
	file >> currHeader;
	if(currHeader == "ThereIsColumnViseVolumeConservation(bool):"){
		file >> Sim->conservingColumnVolumes;
	}
	else{
		std::cerr<<"Error in reading column-vise volume conservation options: "<<currHeader<<", should have been: ThereIsColumnViseVolumeConservation(bool):" <<std::endl;
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
    /**
     * ArtificialRelaxationOptions:
     *   ThereIsArtificaialRelaxation(bool): 1
     *   ArtificialRelaxationTime(sec): 64800
     *   relaxECM(bool): 0
     *
     */
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
    /**
     * zShellOptions:
     *   thereIsEnclosementOfTheTissue(bool): 1
     *   initialLimits(lowerBound,upperBound): -12.50  25.00
     *   finalLimits(lowerBound,upperBound):   -12.50  25.00
     *   initialTime(sec):     0
     *   finalTime(sec):       1800
     */
	string currHeader;
	file >> currHeader;
	if(currHeader == "thereIsEnclosementOfTheTissue(bool):"){
		file >> Sim->encloseTissueBetweenSurfaces;
	}
	else{
		std::cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: thereIsEnclosementOfTheTissue(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "initialLimits(lowerBound,upperBound):"){
		file >> Sim->initialZEnclosementBoundaries[0];
		file >> Sim->initialZEnclosementBoundaries[1];
	}
	else{
		std::cerr<<"Error in reading enclosementoptions: "<<currHeader<<", should have been: initialLimits(lowerBound,upperBound):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalLimits(lowerBound,upperBound):"){
		file >> Sim->finalZEnclosementBoundaries[0];
		file >> Sim->finalZEnclosementBoundaries[1];
	}
	else{
		std::cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: finalLimits(lowerBound,upperBound):" <<std::endl;
		return false;
	}

	file >> currHeader;
	if(currHeader == "initialTime(sec):"){
		file >> Sim->initialTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
		std::cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: initialTime(sec):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "finalTime(sec):"){
		file >> Sim->finalTimeToEncloseTissueBetweenSurfacesSec;
	}
	else{
		std::cerr<<"Error in reading enclosement options: "<<currHeader<<", should have been: finalTime(sec):" <<std::endl;
		return false;
	}
	return true;
}

bool ModelInputObject::readMutationOptions(ifstream& file){
    /**
     * MutationOptions:
     *   numberOfClones(int): 10
     *   cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):
     *  0.8    0.725 3   1 0.1324
     *  0.775  0.75  3   1 0.1324
     *  0.75   0.75  3   1 0.1324
     *  0.7    0.75  3   1 0.1324
     *  0.65   0.75  3   1 0.1324
     *  0.8    0.725  3  1 0.1324
     *  0.775  0.7   3  1 0.1324
     *  0.75   0.7   3  1 0.1324
     *  0.7    0.7   3  1 0.1324
     *  0.65   0.7   3  1 0.1324
     */
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
		for (size_t i=0; i<Sim->numberOfClones; ++i){
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
    /**
     * I believe this should be kept as is in almost all simulations that will be run. I cannot
     * think of a scenario that it should be off.
     *
     * AdhesionOptions:
     *   ThereIsAdhesion(bool): 1
     *   CollapseNodesOnAdhesion(bool): 1
     */
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
    /**
     * This will collapse nodes that are too close such that they can cause element flips,
     * or when nodes are adhered. This should be active unless strictly needed to be off for
     * some reason.
     *
     * NodeCollapseOptions:
     *   ThereIsNodeCollapse(bool): 1
     *
     */
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
    /**
     * DO not add lateral ECM to spherical setups. Be wary when running x & z symmetrixc setups.
     * They have not been thoroughly tested.
     *
     * ExplicitECMOptions:
     *   ThereIsExplicitECM(bool): 1
     *   AddLateralECM(bool): 0
     *   LateralECMThickness(microns): 0.2
     *   ECMRemodellingHalfLife(hour): 8.0
     *   ECMColumnarYoungsModulus:  1600
     *   ECMPeripodialYoungsModulus:  100
     */
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
		//std::cerr<<"stretcherAttached "<<stretcherAttached<<std::endl;
		Sim->PipetteSuction = PipetteSuction;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: PipetteAspitarionActive(bool):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "NumberOfPressureStages(int):"){
		file >> Sim->nPipetteSuctionSteps;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  NumberOfPressureStages(int):" <<std::endl;
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
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  InitiationTimes(sec):" <<std::endl;
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
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been:  Pressures(Pa):" <<std::endl;
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
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: ApicalSuction(bool-will_set_up_basal_suction_if_false):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "TissueStuck(bool-will_fix_the_opposite_surface_in_z):"){
		file >> Sim->TissueStuckOnGlassDuringPipetteAspiration;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: TissueStuck(bool-will_fix_the_opposite_surface_in_z):" <<std::endl;
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
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Centre Position(x,y,z):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_InnerRadius(micron):"){
		file >> Sim->pipetteInnerRadius;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_InnerRadius(micron):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_OuterRadius(micron):"){
		double pippetOuterRad;
		file >> pippetOuterRad;
		Sim->pipetteThickness = pippetOuterRad - Sim->pipetteInnerRadius;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_OuterRadius(micron):" <<std::endl;
		return false;
	}
	file >> currHeader;
	if(currHeader == "Pipette_Effect_Depth(micron):"){
		double dummy;
		file >> dummy;
		Sim->pipetteDepth = dummy;
	}
	else{
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Effect_Depth(micron):" <<std::endl;
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
		std::cerr<<"Error in reading pipette aspiration setup: "<<currHeader<<", should have been: Pipette_Suction_Pressure(x,y,z-unit):" <<std::endl;
		return false;
	}
	return true;
}
