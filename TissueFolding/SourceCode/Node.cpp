#include "Node.h"
#include <math.h>
using namespace std;

Node::Node(int id, int dim, double* pos, int tissuePos, int tissueType){
	Id = id;
	nDim = dim;
	const int n = nDim;
	Position = new double[n];
	RKPosition = new double[n];
	previousStepPosition = new double[n];
	for (int i=0; i<n; ++i){
		Position[i] = pos[i];
		RKPosition[i] = pos[i];
		previousStepPosition[i] = pos[i];
	}
	FixedPos = new bool[3];
	for (int i=0; i<3; ++i){
		FixedPos[i] = false;
		externalViscositySetInFixing[i] = false;
		externalViscosity[i] = 0.0;
		ECMViscosityReductionPerHour[i] = 0;
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
   // allOwnersAblated = false;
}

Node::~Node(){
	delete[] Position;
	delete[] RKPosition;
	delete[] previousStepPosition;
	delete[] FixedPos;
}

double Node::getDisplacement(){
	displacement = 0;
	for (int i=0; i<nDim; ++i){
		double d = previousStepPosition[i] - Position[i];
		displacement += d*d;
	}
	displacement = pow(displacement,0.5);
	return displacement;
}

void Node::updatePreviousPosition(){
	for (int i=0; i<nDim; ++i){
		previousStepPosition[i] = Position[i];
	}
}

void Node::updateECMVisocityWithDeformationRate(double ECMChangeFraction, double averageDisplacement){
	double multiplier = (1.0+ECMChangeFraction);
	if (displacement < averageDisplacement){
		multiplier = 1/multiplier;
	}
	externalViscosity[0] *= multiplier;
	externalViscosity[1] *= multiplier;
	externalViscosity[2] *= multiplier;
	for (int i=0; i<nDim; ++i){
		externalViscosity[i] *= multiplier;
		if (externalViscosity[i] < baseExternalViscosity[i]){
			externalViscosity[i] = baseExternalViscosity[i];
		}
	}
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
				baseExternalViscosity[i] = externalViscosity[i];
			}
		}
	}
	else if (tissuePlacement ==  1){
		for (int i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = ApicalVisc;
				baseExternalViscosity[i] = externalViscosity[i];
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
				baseExternalViscosity[i] = externalViscosity[i];
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
					baseExternalViscosity[i] = externalViscosity[i];
				}
			}
		}
	}
}

//void Node::setExternalViscosity(double ApicalVisc,double BasalVisc, double PeripodialApicalVisc, double PeripodialBasalVisc){
	/**
	 *  This node will take in the apical columnar layer, basal columnar layer, and peripodial membrane viscosities of the tissue as inputs, respectively.
	 *  The external viscosity of the node will be assigned via its Node#tissuePlacement and Node#tissueType. On the columnar layer, nodes that are in the mid-zone of the tissue (neither on the
	 *  apical nor on the basal surface, will take the average of the two values.
	 */
	/*double apicalV, basalV;
	if (tissueType == 0){		//Columnar layer  node
		apicalV = ApicalVisc;
		basalV  = BasalVisc;
	}
	else if (tissueType == 1) {	//Peripodial Membrane node
		apicalV = PeripodialApicalVisc;
		basalV  = PeripodialBasalVisc;
	}
	else { //linker zone
		apicalV = 0.5*(ApicalVisc + PeripodialApicalVisc);
		basalV  = 0.5*(BasalVisc + PeripodialBasalVisc);
	}
	if (tissuePlacement ==0){
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = basalV;}
		if (!externalViscositySetInFixing[1]){externalViscosity[1] = basalV;}
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = basalV;}
	}
	else if (tissuePlacement ==1){
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = apicalV;}
		if (!externalViscositySetInFixing[1]){externalViscosity[1] = apicalV;}
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = apicalV;}
	}
	else if (tissuePlacement == 2 || tissuePlacement == 3){
		//middle or lateral node are equal to the minimum of apical and basal values
		double minV = apicalV;
		if (basalV < apicalV){minV = basalV;};
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = minV;}//(apicalV + basalV) /2.0;
		if (!externalViscositySetInFixing[2]){externalViscosity[1] = minV;}//(apicalV + basalV) /2.0;
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = minV;}//(apicalV + basalV) /2.0;
	}
}*/

bool Node::checkIfNeighbour(int IdToCheck){
	/**
	 *  The function will return true if the node with the unique Node#Id equal to "IdToCheck" is an immediate neighbour of the current node.
	 *  The search will be done through the list Node#immediateNeigs
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
	 *
	 */
	if (mass == 0){ //IF the node does not have any mass, ot means it is ablated, and it should not pack
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
	 *  The function will display on screen the list of weights (normalised masses) of the connected elements.
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
	immediateNeigs.push_back(newNodeId);
}

void Node::addToConnectedElements(int newElementId, double volumePerNode){
	//
	double oldMass = mass;
	mass += volumePerNode;
	double scaler = mass/oldMass;
	int n = connectedElementIds.size();
	for (int j=0;j<n;++j){
		connectedElementWeights[j] /= scaler;
	}
	connectedElementIds.push_back(newElementId);
	connectedElementWeights.push_back(volumePerNode/mass);
}

void Node::removeFromConnectedElements(int ElementId, double volumePerNode){
	double oldMass = mass;
	mass -= volumePerNode;
	double scaler = mass/oldMass;
	int n = connectedElementIds.size();
	int indextToBeDeleted = 0;
	//cout<<" connected element ids: ";
	for (int j=0;j<n;++j){
		connectedElementWeights[j] /= scaler;
		if (connectedElementIds[j] == ElementId){
			indextToBeDeleted = j;
		}
		//cout<<" "<<connectedElementIds[j]<<" ";
	}
	//cout<<endl;
	//cout<<"Id to be deleted: "<<ElementId<<" index to be deleted: "<<indextToBeDeleted<<endl;
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
