#include "Node.h"

using namespace std;

Node::Node(int id, int dim, double* pos, int tissuePos, int tissueType){
	Id = id;
	nDim = dim;
	const int n = nDim;
	//4 RK steps
	Velocity = new double*[4];
	for (int i=0; i<4; ++i){
		//n dimensions
		Velocity[i] = new double[n];
		for (int j=0; j<n; ++j){
			Velocity[i][j] = 0.0;
		}
	}
	Position = new double[n];
	RKPosition = new double[n];
	for (int i=0; i<n; ++i){
		Position[i] = pos[i];
		RKPosition[i] = pos[i];
	}
	FixedPos = new bool[3];
	for (int i=0; i<3; ++i){
		FixedPos[i] = false;
		externalViscositySetInFixing[i] = false;
	}
	externalViscosity[0] = -10.0;
	externalViscosity[1] = -10.0;
	externalViscosity[2] = -10.0;
	tissuePlacement = tissuePos;
	this->tissueType = tissueType;
	atCircumference = false;
	mass = 0.0;
	surface = 0.0;
    zProjectedArea = 0.0;
    symmetricEquivalentId = -1000;
    hasLateralElementOwner = false;
    atSymmetricityBorder = false;
   // allOwnersAblated = false;
}

Node::~Node(){
	//cout<<"called the destructor for node class"<<endl;
	delete[] Position;
	delete[] RKPosition;
    /*for (int i=0; i<4; ++i){
        for (int j=0; j<3; ++j){
            cout<<"Velocity["<<i<<"]["<<j<<"] of node "<<Id<<": "<<Velocity[i][j]<<endl;
        }
    }*/
	for (int i=0; i<4; ++i){
		delete[] Velocity[i];
	}    
	delete[] Velocity;
	delete[] FixedPos;
	//cout<<"finalised the destructor for node class"<<endl;
}

void Node::setExternalViscosity(double ApicalVisc,double BasalVisc){
	/**
	 *  This node will take in the apical and basal external viscosities of the tissue as inputs, respectively.
	 *  The external viscosity of the node will be assigned via its Node#tissuePlacement and Node#tissueType. On the columnar layer, nodes that are in the mid-zone of the tissue (neither on the
	 *  apical nor on the basal surface, will take the average of the two values.
	 */
	if (tissuePlacement ==0){
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = BasalVisc;}
		if (!externalViscositySetInFixing[1]){externalViscosity[1] = BasalVisc;}
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = BasalVisc;}
	}
	else if (tissuePlacement ==1){
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = ApicalVisc;}
		if (!externalViscositySetInFixing[1]){externalViscosity[1] = ApicalVisc;}
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = ApicalVisc;}
	}
	else if (tissuePlacement == 2 || tissuePlacement == 3){
		//middle or lateral node are equal to the minimum of apical and basal values
		double minV = ApicalVisc;
		if (BasalVisc < ApicalVisc){minV = BasalVisc;};
		if (!externalViscositySetInFixing[0]){externalViscosity[0] = minV;}//(apicalV + basalV) /2.0;
		if (!externalViscositySetInFixing[2]){externalViscosity[1] = minV;}//(apicalV + basalV) /2.0;
		if (!externalViscositySetInFixing[2]){externalViscosity[2] = minV;}//(apicalV + basalV) /2.0;
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
