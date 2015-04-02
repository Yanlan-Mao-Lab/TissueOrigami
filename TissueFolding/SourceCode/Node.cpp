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
	}
	Viscosity = -10.0;
	tissuePlacement = tissuePos;
	this->tissueType = tissueType;
	atCircumference = false;
	atPeripodiumCircumference = false;
	mass = 0.0;
	surface = 0.0;
	LinkedPeripodiumNodeId = -1;
}

Node::~Node(){
	//cout<<"called the destructor for node class"<<endl;
	delete[] Position;
	delete[] RKPosition;
	for (int i=0; i<4; ++i){
		delete[] Velocity[i];
	}
	delete[] Velocity;
	delete[] FixedPos;
	//cout<<"finalised the destructor for node class"<<endl;
}


void Node::setViscosity(double ApicalVisc,double BasalVisc){
	if (tissuePlacement ==0){
		Viscosity = BasalVisc;
	}
	else if (tissuePlacement ==1){
		Viscosity = ApicalVisc;
	}
	else if (tissuePlacement == 2 || tissuePlacement == 3){
		//middle or lateral node
		Viscosity = (ApicalVisc + BasalVisc) /2.0;
	}
}

bool Node::checkIfNeighbour(int IdToCheck){
	vector<int>::iterator itInt;
	for(itInt = immediateNeigs.begin(); itInt < immediateNeigs.end(); ++itInt){
		if ((*itInt) == IdToCheck){
			return true;
		}
	}
	return false;
}

bool Node::checkIfNodeHasPacking(){
	if (atPeripodiumCircumference){
		return false;
	}
	if (tissuePlacement == 2){	//Node is midline node (neither apical nor basal)
		return false;
	}
	return true;
}

void Node::getCurrentRKPosition(int RKId, double* pos){
	if (RKId == 0){
		pos[0] = Position[0];
		pos[1] = Position[1];
		pos[2] = Position[2];
	}
	else{
		pos[0]= RKPosition[0];
		pos[1]= RKPosition[1];
		pos[2]= RKPosition[2];
	}
}

void Node::displayConnectedElementIds(){
	int n = connectedElementIds.size();
	cout<<"	Connected Element Ids: ";
	for (int i=0; i<n ; ++i){
		cout<<connectedElementIds[i]<<"	";
	}
	cout<<endl;
}

void Node::displayConnectedElementWeights(){
	int n = connectedElementWeights.size();
	cout<<"	Connected Element weights: ";
	for (int i=0; i<n ; ++i){
		cout<<connectedElementWeights[i]<<"	";
	}
	cout<<endl;
}
