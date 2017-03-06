#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
/**
 *  The node class
 *  */
class Node{
private:

public:
	Node(int id, int dim, double* pos, int tissuePos, int tissueType);	//<Constructer of the node
	~Node();
	bool 			*FixedPos;				///< The boolean array stating if the node's position is fixed in any direction, format: [x y z], true = fixed
	int 			Id;						///< The unique identification number of the node
	int 			nDim;					///< The number of dimensions of the node, (2 or 3)
	double 			*Position;				///< The pointer to the position array of the node. The array itself is declared within the constructor, depending on nDim
	double 			*previousStepPosition;	///< The pointer to the previous step position array of the node. The array itself is declared within the constructor, depending on nDim.
	double 			*RKPosition;			///< The pointer to the position array for position during a Runge-Kutta step array of the node. The array itself is declared within the constructor, depending on nDim
	//double 			**Velocity;				///< The pointer(**) to the velocities of the node for each Runge-Kutta step. The final calculated velocity is stored in Velocity[0]
	bool			externalViscositySetInFixing[3];	///< The boolean array stating if the external viscosity of any axis has been set in node fixing options. The node fixing is carried out before the physical parameter settings in most cases. The boolean check is carried out not to overwrite the existing set viscosity in normal viscosity assignment.
	double 			externalViscosity[3];				///< External viscosity of the node, defined by its placement within the tissue. This can be defined as an external adhesion, ECM remodelling, or any other form of viscosity.
	double 			ECMViscosityReductionPerHour[3];
	double			baseExternalViscosity[3];			///< External viscosity of the node, before any remodelling;
	double			displacement;						///< the displacement of the node from previous time step;
	int 			tissuePlacement;		///< The tissue placement is 0 for basal nodes, 1 for apical nodes, and 2 for middle range
	int 			tissueType;		 		///< The tissue type is 0 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
	bool 			atCircumference;		///< Boolean defining if the node is at the circumference of the columnar layer of the tissue.
	double 			mass;					///< The mass of the node, calculated via the elements that use the node as a vertex
	//double 			surface;				///< The surface of the node, calculated via apical or basal elements. Lateral surfaces are not included
	double			viscositySurface;		///< The surface of the node, calculated for application of external viscosity via surface. It is positive for a surface that is to feel external viscosity, zero othervise.
	double  		zProjectedArea;         ///< The surface of the node, as projected in Z, calculated from apical or pasal surfces of elements, lateral surfaces are not included.
   	vector <int> 	immediateNeigs;				///< The list of Id's for immediate neighbours of the node, i.e. the nodes that are shared by the elements that utilise the owner of this node.
	vector <int> 	connectedElementIds;		///< The list of Id's for elements that are utilising this node.
	vector <double>	connectedElementWeights;	///< The list of weights (normalised mass) for elements that are utilising this node, order is linked to Node#connectedElementIds.
	int 			symmetricEquivalentId;		///< The id of the node that this node will move symmetrically in y, if the tissue is symmetric and only half is simulated.
	bool			hasLateralElementOwner;		///< The boolean stating if any lateral element uses this node
	bool			atSymmetricityBorder;		///< The boolean stating if the node is at the border of symmetricity
	bool		insideEllipseBand;
	int		coveringEllipseBandId;	
	//bool 			allOwnersAblated;			///< Boolean stating if the node is ablated. The node is ablated if all the elements making use of the node are ablated.
	//void setExternalViscosity(double ApicalVisc,double BasalVisc, double PeripodialApicalVisc, double PeripodialBasalVisc);///< The function to set the viscosity of the node.
	void setExternalViscosity(double ApicalVisc,double BasalVisc, bool extendExternalViscosityToInnerTissue);///< The function to set the viscosity of the node.
	bool checkIfNeighbour(int IdToCheck); 				///< The function to check if the node with input Id (IdToCheck) is an immediate neighbour of the owner node
	bool checkIfNodeHasPacking();						///< The function to check if the node is eligible for packing.
	void getCurrentPosition(double* pos);				///< return the current position of the node
	void displayConnectedElementIds();					///< This function will print out a list of connected element Id's
	void displayConnectedElementWeights();				///< This function will print out the weights of the connected elements, in the order of  Id s given in connectedElementIds
	void addToImmediateNeigs(int newNodeId);
	void addToConnectedElements(int newElementId, double volumePerNode);
	void removeFromConnectedElements(int newElementId, double volumePerNode);
	double 	getDisplacement();
	void 	updatePreviousPosition();
	void updateECMVisocityWithDeformationRate(double ECMChangeFraction, double averageDisplacement);
	int getId();
};
#endif
