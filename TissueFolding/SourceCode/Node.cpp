#include "Node.h"
#include <math.h>
#include <algorithm>

using namespace std;

Node::Node(int id, int dim, std::array<double,3> &pos, int tissuePos, int tissueType){
	Id = id;
	nDim = dim;
	FixedPos = {false,false,false};
	externalViscositySetInFixing = {false,false,false};
	externalViscosity = {0.0, 0.0, 0.0};
	externalViscosity = {0.0, 0.0, 0.0};
	initialExternalViscosity = {0.0, 0.0, 0.0};
	ECMViscosityChangePerHour = {0.0, 0.0, 0.0};
	Position = {0.0, 0.0, 0.0};
	RKPosition = {0.0, 0.0, 0.0};
	initialPosition = {0.0, 0.0, 0.0};

	size_t n = nDim;
	for (size_t i=0; i<n; ++i){
		Position[i] = pos[i];
		RKPosition[i] = pos[i];
		initialPosition[i] = pos[i];
	}
	tissuePlacement = tissuePos;
	this->tissueType = tissueType;
	atCircumference = false;
	mass = 0.0;
	displacement = 0;
	viscositySurface=0.0;
    zProjectedArea = 0.0;
    symmetricEquivalentId = -1000;
    hasLateralElementOwner = false;
    atSymmetricityBorder = false;
    insideEllipseBand = false;
    coveringEllipseBandId = -1;
    allOwnersECMMimicing = false;
    slaveTo[0] = -1;
    slaveTo[1] = -1;
    slaveTo[2] = -1;
    isMaster[0] = false;
    isMaster[1] = false;
    isMaster[2] = false;
    adheredTo = -1;
    maximumExternalViscosity[0] = 100000000000;
    maximumExternalViscosity[1] = 100000000000;
    maximumExternalViscosity[2] = 100000000000;
    minimumExternalViscosity[0] = 0;
    minimumExternalViscosity[1] = 0;
    minimumExternalViscosity[2] = 0;
    attachedToPeripodial = false;
    onFoldInitiation = false;
    checkOwnersforEllipseAsignment = false;
    positionUpdateOngoing = false;
    positionUpdateCounter = 0;
    facingLumen = false;
}

Node::~Node(){
}

void Node::setExternalViscosity(double ApicalVisc,double BasalVisc, bool extendExternalViscosityToInnerTissue){
	/**
	 *  This function will take in the apical and basal external viscosities of the tissue as inputs, respectively.
	 *  The external viscosity of the node will be assigned via its Node#tissuePlacement and Node#tissueType. On the columnar layer, nodes that are in the mid-zone of the tissue (neither on the
	 *  apical nor on the basal surface, will take the average of the two values.
	 *
	 */
	if (tissuePlacement ==0){
		for (size_t i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = BasalVisc;
			}
		}
	}
	else if (tissuePlacement ==  1){
		for (size_t i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = ApicalVisc;
			}
		}
	}
	else if (atCircumference){
		//circumferential nodes are equal to the minimum of apical and basal values
		double minV = ApicalVisc;
		if (BasalVisc < ApicalVisc){minV = BasalVisc;};
		for (size_t i=0; i<3; ++i){
			if (!externalViscositySetInFixing[i]){
				externalViscosity[i] = minV;
			}
		}
	}
	else if (tissuePlacement == 2 || tissuePlacement == 3 || tissuePlacement == 4){
		if (extendExternalViscosityToInnerTissue){
			//middle or lateral node are equal to the minimum of apical and basal values
			double minV = ApicalVisc;
			if (BasalVisc < ApicalVisc){minV = BasalVisc;};
			for (size_t i=0; i<3; ++i){
				if (!externalViscositySetInFixing[i]){
					externalViscosity[i] = minV; //(apicalV + basalV) /2.0;
				}
			}
		}
	}
	//externalViscosity[1] *= 0.5;
	for (int i=0; i<3; ++i){
		initialExternalViscosity[i] = externalViscosity[i];
	}

}

bool Node::checkIfNeighbour(int IdToCheck){
	/**
	 *  The function will return true if the node with the unique Node#Id equal to "IdToCheck" is an immediate neighbour of the current node.
	 *  The search will be done through the list Node#immediateNeigs
	 *
	 */
	//vector<int>::iterator itInt;
	for(auto immediateNeigId : immediateNeigs){
		if (immediateNeigId == IdToCheck){
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

std::array<double,3> Node::getCurrentPosition(){
	/**
	 *  The function will return the current position of the owner node.
	 *
	 */
	std::array<double,3> pos = Position;
	return pos;
}

void Node::displayPosition(){
	cout<<"Node "<<Id<<" position: "<<Position[0]<<" "<<Position[1]<<" "<<Position[2]<<endl;
}

void Node::displayConnectedElementIds(){
	/**
	 *  The function will display on screen the list of unique Node#Id s for the elements utilising this node.
	 *  These elements are listed in Node#connectedElementIds.
	 *
	 */
	cout<<"	Connected Element Ids: ";
	for (auto connectedId: connectedElementIds){
		cout<<connectedId<<"	";
	}
	cout<<endl;
}

void Node::displayConnectedElementWeights(){
	/**
	 *  The function will display on screen output the list of weights (normalised masses) of the connected elements.
	 *  These weights are stored in Node#connectedElementWeights, and the order is linked to the list: Node#connectedElementIds.
	 *
	 */
	cout<<"	Connected Element weights: ";
	for (auto connectedWeight : connectedElementWeights){
		cout<<connectedWeight<<"	";
	}
	cout<<endl;
}

void  Node::addToImmediateNeigs(int newNodeId){
	/**
	 *  The function will add the input node id to the vector of immediate neighbours of the current node (Node#immediateNeigs).
	 *
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
	 *
	 */
	double scaler = mass/oldMass;
	for (auto connectedWeight : connectedElementWeights){
		connectedWeight /= scaler;
	}
	/**
	 * Then the new element and its corresponding id will be added to the lists of the node, Node#connectedElementIds, and Node#connectedElementWeights, respectively.
	 *
	 */
	connectedElementIds.push_back(newElementId);
	connectedElementWeights.push_back(volumePerNode/mass);
}

bool Node::isMyNeig(int nodeId){
	/**
	 *  This function will take in a node id, and check if the node is an immediate neighbour of itself.
	 *  The check is done through its Node#immediateNeigs list. The function will return true if the
	 *  input id is a neighbour, and false otherwise.
	 *
	 */
	if (find(immediateNeigs.begin(), immediateNeigs.end(),nodeId)!=immediateNeigs.end()){
		return true;
	}
	return false;
}

bool Node::isNeigWithMyCollapsedNodes(int NodeId, const std::vector <std::unique_ptr<Node>>& Nodes){
	/**
	 *  This function will take in a node id, and the address of the list of node pointers in the simulation.
	 *  It will check if the input node id is an immediate neighbour of any of the nodes it has been collapsed with.
	 *  The check is done through  Node#immediateNeigs lists of all the nodes on its Node#collapsedWith list.
	 *  The function will return true if the input id is a neighbour of any node this node has collapsed with, and false otherwise.
	 *
	 */
	vector<int>::iterator itCollapsedWithNodeIndex;
	for (itCollapsedWithNodeIndex = collapsedWith.begin();itCollapsedWithNodeIndex<collapsedWith.end();++itCollapsedWithNodeIndex){
		if (find(Nodes[(*itCollapsedWithNodeIndex)]->immediateNeigs.begin(), Nodes[(*itCollapsedWithNodeIndex)]->immediateNeigs.end(),NodeId)!=Nodes[(*itCollapsedWithNodeIndex)]->immediateNeigs.end()){
			//the collapsed node of master is neig to slave
			return true;
		}
	}
	return false;
}

//version without element checking, and adhesion
void Node::collapseOnNode(const std::vector <std::unique_ptr<Node>>& Nodes, int masterNodeId){
	/**
	 *  This function will collapse this Node with the input node id (second parameter).
	 *  Collapse will bind the two nodes' all degrees of freedoms to each other, and bring them to the same position,
	 *  mid point of their positions at the time of the function call.
	 *  The function will take in the address of the list of node pointers in the simulation, and a node id.
	 *
	 */
	bool debugDisplay = true;
	//add the nodes master have collapsed with on me:
	vector <int> newCollapseList;
	newCollapseList.push_back(Id);
	newCollapseList.push_back(masterNodeId);
	//take the number of collaspsed nodes before merge
	int nThis = collapsedWith.size();
	int nMaster = Nodes[masterNodeId]->collapsedWith.size();
	for (auto thisNodeAlreadycollapsedWith : collapsedWith){
		newCollapseList.push_back(thisNodeAlreadycollapsedWith);
	}
	for (auto masterNodeAlreadycollapsedWith : Nodes[masterNodeId]->collapsedWith){
		newCollapseList.push_back(masterNodeAlreadycollapsedWith);
	}
	//remove duplicates:
	sort( newCollapseList.begin(), newCollapseList.end() );
	newCollapseList.erase( unique( newCollapseList.begin(), newCollapseList.end() ), newCollapseList.end() );
	size_t nNew = newCollapseList.size();

	//calculate the average weighted pos:
	if(nMaster == 0){nMaster = 1;};
	if(nThis == 0){nThis = 1;};
	double avrPos[3] = {0,0,0};
	avrPos[0] = (nThis*Position[0] + nMaster*Nodes[masterNodeId]->Position[0])/(nThis+nMaster);
	avrPos[1] = (nThis*Position[1] + nMaster*Nodes[masterNodeId]->Position[1])/(nThis+nMaster);
	avrPos[2] = (nThis*Position[2] + nMaster*Nodes[masterNodeId]->Position[2])/(nThis+nMaster);
	//If there is a fixed position in any of the nodes on the list, update to that pos:
	bool fix[3] = {false, false, false};
	for (size_t i = 0; i<nNew; ++i){
		for (size_t j=0; j<3; ++j){
			if(Nodes[newCollapseList[i]]->FixedPos[j]){
				avrPos[j] = Nodes[newCollapseList[i]]->Position[j];
				fix[j] = true;
			}
		}
	}

	//update adhesion:
	int baseNode = 0;
	int adhesionPoint = -1;
	for (size_t i = 0; i<nNew; ++i){
		if(Nodes[newCollapseList[i]]->adheredTo > -1){
			adhesionPoint = Nodes[newCollapseList[i]]->adheredTo;
			baseNode = newCollapseList[i];
			break;
		}
	}
	if (adhesionPoint > -1){
		//at least on of the nodes is adheres. make all of them adhered:
		for (size_t i = 0; i<nNew; ++i){
			if(newCollapseList[i] == adhesionPoint){
				//the adhered node is bound to be on this list, I keep that one adhered to the base
				Nodes[newCollapseList[i]]->adheredTo = baseNode;
			}
			else{
				//all other nodes on the list are adhered to the selected adheredNode
				Nodes[newCollapseList[i]]->adheredTo = adhesionPoint;
			}
		}
	}

	if (debugDisplay){
		cout<<"in single element collapse, masterId: "<<masterNodeId<<" slaveId: "<<Id<<endl;
		cout<<"  nMaster "<<nMaster<<" nThis "<<nThis<<endl;
		cout<<"  avrPos:     "<<avrPos[0]<<" "<<avrPos[1]<<" "<<avrPos[2]<<endl;
		cout<<"  Position:   "<<Position[0]<<" "<<Position[1]<<" "<<Position[2]<<endl;
		cout<<"  Master pos: "<<Nodes[masterNodeId]->Position[0]<<" "<<Nodes[masterNodeId]->Position[1]<<" "<<Nodes[masterNodeId]->Position[2]<<endl;
		cout<<"  collapsed with list before update: ";
		for(size_t i=0; i<collapsedWith.size(); ++i){
			cout<<collapsedWith[i]<<" ";
		}
		cout<<endl;
		cout<<"  masters collapsed list before update: ";
		for(size_t i=0; i<Nodes[masterNodeId]->collapsedWith.size(); ++i){
			cout<<Nodes[masterNodeId]->collapsedWith[i]<<" ";
		}
		cout<<endl;
	}
	//update collapsed list of all the members to be the same, update positions:
	nNew = newCollapseList.size();
	for (size_t i = 0; i<nNew; ++i){
		Nodes[newCollapseList[i]]->collapsedWith = newCollapseList;
		for (int j=0; j<3; ++j){
			if(!Nodes[newCollapseList[i]]->FixedPos[j]){
				Nodes[newCollapseList[i]]->Position[j] = avrPos[j];
				Nodes[newCollapseList[i]]->FixedPos[j] = fix[j];
			}
		}
	}
	if (debugDisplay){
		cout<<"  collapsed with list after update: ";
		for(size_t i=0; i<collapsedWith.size(); ++i){
			cout<<collapsedWith[i]<<" ";
		}
		cout<<endl;
	}
}

void Node::clearDuplicatesFromCollapseList(){
	//remove duplicates:
	sort( collapsedWith.begin(), collapsedWith.end() );
	collapsedWith.erase( unique( collapsedWith.begin(), collapsedWith.end() ), collapsedWith.end() );
}

void Node::getNewCollapseListAndAveragePos(vector<int> &newCollapseList, double* avrPos, bool* fix, const std::vector <std::unique_ptr<Node>>& Nodes, int masterNodeId){
	for (size_t j=0; j<3; ++j){
		avrPos[j] = 0;
		fix[j] = false;
	}

	newCollapseList.push_back(Id);
	newCollapseList.push_back(masterNodeId);
	size_t nThis = collapsedWith.size();
	for (auto collapsedNodeId : collapsedWith){
		newCollapseList.push_back(collapsedNodeId);
	}
	size_t nMaster = Nodes[masterNodeId]->collapsedWith.size();
	for (auto collapseNodeIdOnMAster : Nodes[masterNodeId]->collapsedWith){
		newCollapseList.push_back(collapseNodeIdOnMAster);
	}
	//remove duplicates:
	sort( newCollapseList.begin(), newCollapseList.end() );
	newCollapseList.erase( unique( newCollapseList.begin(), newCollapseList.end() ), newCollapseList.end() );
	size_t nNew = newCollapseList.size();

	//calculate the average weighted pos:
	if(nMaster == 0){nMaster = 1;};
	if(nThis == 0){nThis = 1;};
	avrPos[0] = (nThis*Position[0] + nMaster*Nodes[masterNodeId]->Position[0])/(nThis+nMaster);
	avrPos[1] = (nThis*Position[1] + nMaster*Nodes[masterNodeId]->Position[1])/(nThis+nMaster);
	avrPos[2] = (nThis*Position[2] + nMaster*Nodes[masterNodeId]->Position[2])/(nThis+nMaster);
	//If there is a fixed position in any of the nodes on the list, update to that pos:
	for (size_t i = 0; i<nNew; ++i){
		for (size_t j=0; j<3; ++j){
			if(Nodes[newCollapseList[i]]->FixedPos[j]){
				avrPos[j] = Nodes[newCollapseList[i]]->Position[j];
				fix[j] = true;
			}
		}
	}
}


//version that will check element collapse:
void Node::collapseOnNodeInStages( vector<int> &newCollapseList, double* avrPos, bool* fix, const std::vector <std::unique_ptr<Node>>& Nodes){
	//update adhesion:
	size_t nNew = newCollapseList.size();
	int baseNode = 0;
	int adhesionPoint = -1;
	for (size_t i = 0; i<nNew; ++i){
		if(Nodes[newCollapseList[i]]->adheredTo > -1){
			adhesionPoint = Nodes[newCollapseList[i]]->adheredTo;
			baseNode = newCollapseList[i];
			break;
		}
	}
	if (adhesionPoint > -1){
		//at least on of the nodes is adhered. make all of them adhered:
		for (size_t i = 0; i<nNew; ++i){
			if(newCollapseList[i] == adhesionPoint){
				//the adhered node is bound to be on this list, I keep that one adhered to the base
				Nodes[newCollapseList[i]]->adheredTo = baseNode;
			}
			else{
				//all other nodes on the list are adhered to the selected adheredNode
				Nodes[newCollapseList[i]]->adheredTo = adhesionPoint;
			}
		}
	}
	//update collapsed list of all the members to be the same, update positions:
	nNew = newCollapseList.size();
	for (size_t i = 0; i<nNew; ++i){
		Nodes[newCollapseList[i]]->collapsedWith = newCollapseList;
		for (size_t j=0; j<3; ++j){
			if(!Nodes[newCollapseList[i]]->FixedPos[j]){
				Nodes[newCollapseList[i]]->positionUpdateOngoing = true;
			}
		}
		Nodes[newCollapseList[i]]->updatePositionTowardsPoint(avrPos,fix);
	}
}

void Node::collapseOnNode(vector<int> &newCollapseList, double* avrPos, bool* fix, const std::vector <std::unique_ptr<Node>>& Nodes, int masterNodeId){
	bool debugDisplay = true;
	//add the nodes master have collapsed with on me:
	//update adhesion:
	size_t nNew = newCollapseList.size();
	int baseNode = 0;
	int adhesionPoint = -1;
	for (size_t i = 0; i<nNew; ++i){
		if(Nodes[newCollapseList[i]]->adheredTo > -1){
			adhesionPoint = Nodes[newCollapseList[i]]->adheredTo;
			baseNode = newCollapseList[i];
			break;
		}
	}
	if (adhesionPoint > -1){
		//at least on of the nodes is adhered. make all of them adhered:
		for (size_t i = 0; i<nNew; ++i){
			if(newCollapseList[i] == adhesionPoint){
				//the adhered node is bound to be on this list, I keep that one adhered to the base
				Nodes[newCollapseList[i]]->adheredTo = baseNode;
			}
			else{
				//all other nodes on the list are adhered to the selected adheredNode
				Nodes[newCollapseList[i]]->adheredTo = adhesionPoint;
			}
		}
	}

	if (debugDisplay){
		cout<<"masterId: "<<masterNodeId<<" slaveId: "<<Id<<endl;
		cout<<"  avrPos:     "<<avrPos[0]<<" "<<avrPos[1]<<" "<<avrPos[2]<<endl;
		cout<<"  Position:   "<<Position[0]<<" "<<Position[1]<<" "<<Position[2]<<endl;
		cout<<"  Master pos: "<<Nodes[masterNodeId]->Position[0]<<" "<<Nodes[masterNodeId]->Position[1]<<" "<<Nodes[masterNodeId]->Position[2]<<endl;
		cout<<"  collapsed with list before update: ";
		for(size_t i=0; i<collapsedWith.size(); ++i){
			cout<<collapsedWith[i]<<" ";
		}
		cout<<endl;
		cout<<"  masters collapsed list before update: ";
		for(size_t i=0; i<Nodes[masterNodeId]->collapsedWith.size(); ++i){
			cout<<Nodes[masterNodeId]->collapsedWith[i]<<" ";
		}
		cout<<endl;
	}
	//update collapsed list of all the members to be the same, update positions:
	nNew = newCollapseList.size();
	for (size_t i = 0; i<nNew; ++i){
		Nodes[newCollapseList[i]]->collapsedWith = newCollapseList;
		for (size_t j=0; j<3; ++j){
			if(!Nodes[newCollapseList[i]]->FixedPos[j]){
				Nodes[newCollapseList[i]]->Position[j] = avrPos[j];
				Nodes[newCollapseList[i]]->FixedPos[j] = fix[j];
			}
		}
	}
	if (debugDisplay){
		cout<<"  collapsed with list after update: ";
		cout<<"  Position:   "<<Position[0]<<" "<<Position[1]<<" "<<Position[2]<<endl;
		cout<<"  Master pos: "<<Nodes[masterNodeId]->Position[0]<<" "<<Nodes[masterNodeId]->Position[1]<<" "<<Nodes[masterNodeId]->Position[2]<<endl;
		for(size_t i=0; i<collapsedWith.size(); ++i){
			cout<<collapsedWith[i]<<" ";
		}
		cout<<endl;
	}
}

void Node::updatePositionTowardsPoint(double* avrPos,bool* fix){
	positionUpdateCounter ++;
	for (size_t j=0; j<3; ++j){
		if(!FixedPos[j]){
			double d =  avrPos[j]- Position[j];
			if (positionUpdateCounter>2){
				Position[j]  = avrPos[j];
			}
			else{
				//collapse at three time steps
				d /= 2.0;
				Position[j] += d;
			}
			FixedPos[j] = fix[j];
		}
	}
	if (positionUpdateCounter>2){
		positionUpdateOngoing = false;
		positionUpdateCounter = 0;
	}
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
	size_t n = connectedElementIds.size();
	size_t indextToBeDeleted = 0;
	for (size_t j=0;j<n;++j){
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
	std::vector<int>::iterator itElementId = connectedElementIds.begin();
	itElementId += indextToBeDeleted;
	if (itElementId != connectedElementIds.end()) {
	  using std::swap;
	  // swap the one to be removed with the last element
	  // and remove the item at the end of the container
	  // to prevent moving all items after '5' by one
	  swap(*itElementId, connectedElementIds.back());
	  connectedElementIds.pop_back();
	}
	std::vector<double>::iterator itElementWeight = connectedElementWeights.begin();
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

bool Node::isECMChangeAppliedToNode(bool changeApicalECM, bool changeBasalECM, vector<int> &ECMChangeEllipseBandIds, size_t numberOfECMChangeEllipseBands){
	if (allOwnersECMMimicing){
		if  (    (changeApicalECM && tissuePlacement == 1 )
			  || (changeBasalECM  && tissuePlacement == 0 )
					  //|| (changeStiffnessBasalECM  && tissuePlacement == 2 )
			){
			if(insideEllipseBand){
				for (size_t ECMReductionRangeCounter = 0; ECMReductionRangeCounter<numberOfECMChangeEllipseBands; ++ECMReductionRangeCounter){
					if (coveringEllipseBandId == ECMChangeEllipseBandIds[ECMReductionRangeCounter]){
						return true;
					}
				}
			}

		}
	}
	return false;
}

bool Node::haveCommonOwner(Node* nodeSlave){
	vector<int>::iterator itOwnerId;
	itOwnerId = find_first_of ( connectedElementIds.begin(),  connectedElementIds.end(), nodeSlave->connectedElementIds.begin(), nodeSlave->connectedElementIds.end());
	if (itOwnerId!=connectedElementIds.end()){
		return true;
	}
	return false;
}

int Node::getCommonOwnerId(Node* nodeSlave){
	vector<int>::iterator itOwnerId;
	itOwnerId = find_first_of ( connectedElementIds.begin(),  connectedElementIds.end(), nodeSlave->connectedElementIds.begin(), nodeSlave->connectedElementIds.end());
	if (itOwnerId!=connectedElementIds.end()){
		return (*itOwnerId);
	}
	return -1;
}

