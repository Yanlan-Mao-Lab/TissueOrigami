/*
 * ModelInputObject.h
 *
 *  Created on: 5 Aug 2014
 *      Author: melda
 */

#ifndef MODELINPUTOBJECT_H_
#define MODELINPUTOBJECT_H_

#include <stdio.h>
#include <iostream>
#include <memory>
//#include "Simulation.h"

class Simulation;
using namespace std;

/*! ModelInputObject class */
class ModelInputObject {
private:
	void printErrorMessage(string currentInput, string sourceFuction, string expectedInput); //< This function writes the error message
	bool checkFileStatus(ifstream &file, string fileName); 	///< This function checks the health of given file.
	bool readPysicalProperties(ifstream &file);				///< This function reads the physical parameters of the tissue from file.
	bool readSaveOptions(ifstream &file);					///< This function reads the save options of the simulation from file.
	bool readMeshParameters(ifstream& file);				///< This function reads the mesh structure related parameters of the tissue from file.
	bool readPeripodialMembraneParameters(ifstream& file);	///< This function reads the peripodial membrane related parameters of the tissue from file.
	bool readLinkerZoneParameters(ifstream& file);			///< This function reads the linker zone related parameters of the tissue from file.
	bool readExternalViscosityParameters(ifstream& file);	///< This function reads the external viscosity setup for the whole tissue.
	bool readNodeFixingParameters(ifstream& file);			///< This function reads the inputs relating to fixing the nodes of the tissue from file.
	bool readNodeBindingParameters(ifstream& file);			///< This function reads the inputs relating to binding the nodes to each other for boundary conditions
	bool readManupulationParamters(ifstream& file);			///< This function reads the inputs relating to manipulations made on the tissue after mesh read.
	bool readTimeParameters(ifstream &file);				///< This function reads the time related parameters of the simulation from file.
	bool readMeshType2(ifstream& file);						///< This function reads the mesh structure details for a mesh input as columns and rows of prisms.
	bool readMeshType4(ifstream& file);						///< This function reads the mesh structure details for a mesh given as a pre-assembled mesh in a separate input file.
	bool readGrowthOptions(ifstream& file);					///< This function reads the growth functions and related parameters of the simulation from file.
	bool readGrowthType1(ifstream& file);					///< This function reads the uniform growth parameters from file (UniformGrowthFunction).
	bool readGrowthType2(ifstream& file);					///< This function reads the ring growth parameters from file (RingGrowthFunction).
	bool readGrowthType3(ifstream& file);					///< This function reads the grid based growth parameters from file (GridBasedGrowthFunction). It will utilise a separate input file storing the growth rates.
	bool readShapeChangeOptions(ifstream& file);			///< This function reads the active shape change functions and related parameters of the simulation from file.
	bool readMyosinOptions(ifstream& file);					///< This function reads the active equilibrium myosin concentration stimultion functions and related parameters of the simulation from file.
	bool readPlasticDeformationOptions(ifstream& file);		///< This function reads the parametrs for plastic deformation, as a response to strains and stresses in the tissue.
	bool readGridBasedMyosinFunction(ifstream& file);		///< This function reads the myosin parameters from file (GridBased). It will utilise a separate input file storing the equilibrium myosin levels and orientations.
	bool readShapeChangeType1(ifstream& file);				///< This function reads the shape change  parameters from file (UniformShapeChange).
	bool readShapeChangeType2(ifstream& file);				///< This function reads the shape change  parameters from file (markerEllipseBased).
	bool readMarkerEllipseBandOptions(ifstream& file);		///< This function reads the ellipse band definition parameters
	bool readShapeChangeType3(ifstream& file);				///< This function reads the shape change  parameters from file (GridBasedShapeChange).
	bool readStretcherSetup(ifstream& file);				///< This function reads the stretcher experimental setup parameters of the simulation from file.
	bool readPipetteSetup(ifstream& file);					///< This function reads the pipette aspiration experimental setup parameters of the simulation from file.
	bool readECMPerturbation(ifstream& file);				///< This function reads the ECM perturbation options parameters
	bool readStiffnessPerturbation(ifstream& file);			///< This function reads the stiffness perturbation options parameters
	bool readApicoBasalVolumeRedistribution(ifstream& file);///< This funciton reads the options for apical-basal volume redistribution
	bool readExplicitECMOptions(ifstream& file);			///< This function reads parameters relating to the definitin of an explicit ECM layer.
	bool readAdhesionOptions(ifstream& file);				///< This function reads the adhesion options
	bool readNodeCollapseOptions(ifstream& file);			///< This function reads the options for collapsing the nodes
	bool readExplicitActinOptions(ifstream& file);			///< This function reads the options for explicit apical actin definition
	bool readColumnViseVolumeConservationOptions(ifstream& file); ///< This function reads the options for column-wise volume conservations (volume per columnar cell rather than element)
	bool readMutationOptions(ifstream& file);				///< This funciton reads the mutation options
	bool readartificialRelaxationOptions(ifstream& file);	///< This function reads the artificial stress relaxation options
	bool readLumenOptions(ifstream& file);					///< This function reads the lumen options

	bool readEnclosementOptions(ifstream& file);

public:
	Simulation* Sim;				///< The pointer to the simulation object, for which the parameters are being read from the modelinput file.
	const char*  parameterFileName;	///< The name (including path) of the file containing the model input parameters.
	string meshFileName;			///< The name (including path) of the file containing the input mesh file for tissue geometry.
	ModelInputObject();				///< The constructor of the ModelInputObject.
	~ModelInputObject();

	bool readParameters();			///< This is the main funciton reading the parameters from file
};

#endif /* MODELINPUTOBJECT_H_ */
