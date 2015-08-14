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
	double 			*RKPosition;			///< The pointer to the position array for position during a Runge-Kutta step array of the node. The array itself is declared within the constructor, depending on nDim
	double 			**Velocity;				///< The pointer(**) to the velocities of the node for each Runge-Kutta step. The final calculated velocity is stored in Velocity[0]
	double 			Viscosity;				///< Viscosity of the node, defined by its placement within the tissue
	int 			tissuePlacement;		///< The tissue placement is 0 for basal nodes, 1 for apical nodes, and 2 for middle range
	int 			tissueType;		 		///< The tissue type is 0 for columnar layer and 1 for peripodial membrane
	bool 			atCircumference;		///< Boolean defining if the node is at the circumference of the columnar layer of the tissue.
	double 			mass;					///< The mass of the node, calculated via the elements that use the node as a vertex
	double 			surface;				///< The surface of the node, calculated via apical or basal elements. Lateral surfaces are not included
    double  		zProjectedArea;         ///< The surface of the node, as projected in Z, calculated from apical or pasal surfces of elements, lateral surfaces are not included.
    bool 			atPeripodialCircumference;						///< Boolean defining if the node is at the circumference of the peripodial membrane of the tissue.
	vector <int> 	AssociatedNodesDueToPeripodialMembrane;			///< The list of columnar Node Id's that are associated with this node of the peripodial membrane circumference. The forces from the peripodial membrane are distributed to the columnar layer, hence anchoring the peripodial membrane to the columnar layer at the circumference, via this list.
	vector <double>	AssociatedNodeWeightsDueToPeripodialMembrane;	///< The list of columnar Node weights (normalised masses) that are associated with this node of the peripodial membrane circumference. The forces from the peripodial membrane are distributed to the columnar layer,hence anchoring the peripodial membrane to the columnar layer at the circumference, as weighted by this list, order is linked to Node#AssociatedNodesDueToPeripodialMembrane.
	int 			LinkedPeripodialNodeId;							///< The peripodial node Id a columnar layer node is linked to. The peripodial node of the Id equal to this variable will distribute its forces on its associated nodes, one of which is the owner of this parameter.
	vector <int> 	immediateNeigs;				///< The list of Id's for immediate neighbours of the node, i.e. the nodes that are shared by the elements that utilise the owner of this parameter.
	vector <int> 	connectedElementIds;		///< The list of Id's for elements that are utilising this node.
	vector <double>	connectedElementWeights;	///< The list of weights (normalised mass) for elements that are utilising this node, order is linked to Node#connectedElementIds.

	void setViscosity(double ApicalVisc,double BasalVisc, double PeripodialViscosity); ///< The function to set the viscosity of the node.
	bool checkIfNeighbour(int IdToCheck); 				///< The function to check if the node with input Id (IdToCheck) is an immediate neighbour of the owner node
	bool checkIfNodeHasPacking();						///< The function to check if the node is eligible for packing.
	void getCurrentRKPosition(int RKId, double* pos);	///< return the current position of the node
	void displayConnectedElementIds();					///< This function will print out a list of connected element Id's
	void displayConnectedElementWeights();				///< This function will print out the weights of the connected elements, in the order of  Id s given in connectedElementIds
};
#endif
