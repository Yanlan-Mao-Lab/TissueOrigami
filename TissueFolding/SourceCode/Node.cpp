#include "Node.h"
#include <math.h>
using namespace std;

Node::Node(int id, int dim, double* pos, int tissuePos, int tissueType){
	Id = id;
	nDim = dim;
	const int n = nDim;
	Position = new double[n];
	RKPosition = new double[n];
	initialPosition = new double[n];
	for (int i=0; i<n; ++i){
		Position[i] = pos[i];
		RKPosition[i] = pos[i];
		initialPosition[i] = pos[i];
	}
	FixedPos = new bool[3];
	for (int i=0; i<3; ++i){
		FixedPos[i] = false;
		externalViscositySetInFixing[i] = false;
		externalViscosity[i] = 0.0;
		ECMViscosityChangePerHour[i] = 0;
	}
	tissuePlacement = tissuePos;
	this->tissueType = tissueType;
	atCircumference = false;
	mass = 0.0;
	//surface = 0.0;
	viscositySurface=0.0;
    zProjectedArea = 0.0;
    symmetricEquivalentId = -1000;
    hasLateralElementOwner = false;
    atSymmetricityBorder = false;
    insideEllipseBand = false;
    coveringEllipseBandId = -1;
    allOwnersECMMimicing = false;
   // allOwnersAblated = false;
}

Node::~Node(){
	delete[] Position;
	delete[] RKPosition;
	delete[] initialPosition;
	delete[] FixedPos;
}

void Node::setExternalViscosity(double ApicalVisc,double BasalVisc, bool extendExternalViscosityToInnerTissue){
	/**
	 *  This node will take in the apical and basal external viscosities of the tissue as inputs, respectively.
	 *  The external viscosity of the node will be assigned via its Node#tissuePlacement and Node#tissueType. On the columnar layer, nodes that are in the mid-zone of the tissue (neither on the
	 *  apical nor on the basal surface, will take the average of the two values.
	 */
	if (tissuePlacement ==0){
		for (int i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = BasalVisc;
			}
		}
	}
	else if (tissuePlacement ==  1){
		for (int i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = ApicalVisc;
			}
		}
	}
	else if (atCircumference){
		//circumferential nodes are equal to the minimum of apical and basal values
		double minV = ApicalVisc;
		if (BasalVisc < ApicalVisc){minV = BasalVisc;};
		for (int i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = minV;
			}
		}
	}
	else if (tissuePlacement == 2 || tissuePlacement == 3){
		if (extendExternalViscosityToInnerTissue){
			//middle or lateral node are equal to the minimum of apical and basal values
			double minV = ApicalVisc;
			if (BasalVisc < ApicalVisc){minV = BasalVisc;};
			for (int i=0; i<3; ++i){
				if (!externalViscositySetInFixing[i]){
					externalViscosity[i] = minV; //(apicalV + basalV) /2.0;
				}
			}
		}
	}
}

bool Node::checkIfNeighbour(int IdToCheck){
	/**
	 *  The function will return true if the node with the unique Node#Id equal to "IdToCheck" is an immediate neighbour of the current node.
	 *  The search will be done through the list Node#immediateNeigs
	 *
	 */
	vector<int>::iterator itInt;
	for(itInt = immediateNeigs.begin(); itInt < immediateNeigs.end(); ++itInt){
		if ((*itInt) == IdToCheck){
			return true;
		}
	}
	return false;
}

bool Node::checkIfNodeHasPacking(){
	/**
	 *  The function will return true if the node is eligible for packing calculation. This packing will ensure volume exclusion, and is calculated via
	 *  function Simulation#calculatePacking. It is not necessary to calculate packing under the following conditions:
	 *  1) The node is at the middle of the columnar layer, the packing should have stopped any other node/element coming close enough to this node, as
	 *  they would need to penetrate through the apical or basal surface of the tissue to reach this node.
	 */
	if (mass == 0){ //IF the node does not have any mass, then means it is ablated, and it should not pack
		return false;
	}
	if (hasLateralElementOwner){ //if the node is owned by any lateral element connecitng peripodial to columnar layers, then it is not affected by packing
		return false;
	}
	if (tissuePlacement == 0 || tissuePlacement == 1){	//Node is apical or basal)
		return true;
	}
	return false;
}

void Node::getCurrentPosition(double* pos){
	/**
	 *  The function will return the current position of the owner node.
	 */
	pos[0] = Position[0];
	pos[1] = Position[1];
	pos[2] = Position[2];
}

void Node::displayConnectedElementIds(){
	/**
	 *  The function will display on screen the list of unique Node#Id s for the elements utilising this node.
	 *  These elements are listed in Node#connectedElementIds.
	 */
	int n = connectedElementIds.size();
	cout<<"	Connected Element Ids: ";
	for (int i=0; i<n ; ++i){
		cout<<connectedElementIds[i]<<"	";
	}
	cout<<endl;
}

void Node::displayConnectedElementWeights(){
	/**
	 *  The function will display on screen output the list of weights (normalised masses) of the connected elements.
	 *  These weights are stored in Node#connectedElementWeights, and the order is linked to the list: Node#connectedElementIds.
	 */
	int n = connectedElementWeights.size();
	cout<<"	Connected Element weights: ";
	for (int i=0; i<n ; ++i){
		cout<<connectedElementWeights[i]<<"	";
	}
	cout<<endl;
}

void  Node::addToImmediateNeigs(int newNodeId){
	/**
	 *  The function will add the input node id to the vector of immediate neighbours of the current node (Node#immediateNeigs).
	 */
	immediateNeigs.push_back(newNodeId);
}

void Node::addToConnectedElements(int newElementId, double volumePerNode){
	/**
	 *  This function adds the input newElementId (int) to the list of elements connected by this node, updating the mass,
	 *  and weights of mass per connected element in the process.
	 *
	 *  First the mass before addition of the new element is recorded. Then the Node#mass
	 *  is updated.
	 */
	double oldMass = mass;
	mass += volumePerNode;
	/**
	 * Each of the already recorded weights (in Node#connectedElementWeights) of the connected elements (in Node#connectedElementIds)
	 * will be updated with the scale newMass / oldMass.
	 */
	double scaler = mass/oldMass;
	int n = connectedElementIds.size();
	for (int j=0;j<n;++j){
		connectedElementWeights[j] /= scaler;
	}
	/**
	 * Then the new element and its corresponding id will be added to the lists of the node, Node#connectedElementIds, and Node#connectedElementWeights, respectively.
	 */
	connectedElementIds.push_back(newElementId);
	connectedElementWeights.push_back(volumePerNode/mass);
}

void Node::removeFromConnectedElements(int ElementId, double volumePerNode){
	/**
	 *  This function removes the input newElementId (int) from the list of elements connected by this node,
	 *  updating the mass, and weights of mass per connected element in the process.
	 *
	 *  First the mass before addition of the new element is recorded. Then the Node#mass
	 *  is updated.
	 */
	double oldMass = mass;
	mass -= volumePerNode;
	double scaler = mass/oldMass;
	/**
	 *  Each of the already recorded weights (in Node#connectedElementWeights) of the connected elements (in Node#connectedElementIds)
	 *  will be updated with the scale newMass / oldMass. The index of the element to be deleted on vector Node#connectedElementIds
	 *  is obtained in the process.
	 */
	int n = connectedElementIds.size();
	int indextToBeDeleted = 0;
	for (int j=0;j<n;++j){
		connectedElementWeights[j] /= scaler;
		if (connectedElementIds[j] == ElementId){
			indextToBeDeleted = j;
		}
	}
	/**
	 *  Then the element with the obtained index is removed from the lists of connected
	 *  element ids and weights. For efficiency, the element to be removed is swapped
	 *  with the last element of each vector, and the vector popped back.
	 */
	vector<int>::iterator itElementId = connectedElementIds.begin();
	itElementId += indextToBeDeleted;
	if (itElementId != connectedElementIds.end()) {
	  using std::swap;
	  // swap the one to be removed with the last element
	  // and remove the item at the end of the container
	  // to prevent moving all items after '5' by one
	  swap(*itElementId, connectedElementIds.back());
	  connectedElementIds.pop_back();
	}
	vector<double>::iterator itElementWeight = connectedElementWeights.begin();
	itElementWeight +=indextToBeDeleted;
	if (itElementWeight != connectedElementWeights.end()) {
	  using std::swap;
	  // swap the one to be removed with the last element
	  // and remove the item at the end of the container
	  // to prevent moving all items after '5' by one
	  swap(*itElementWeight, connectedElementWeights.back());
	  connectedElementWeights.pop_back();
	}

	n = connectedElementIds.size();
}

int  Node::getId(){
	/**
	 * The function returns the Node#Id of the node.
	 */
	return Id;
}
