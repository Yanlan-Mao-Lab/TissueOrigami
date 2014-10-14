#include "Node.h"

using namespace std;

Node::Node(int id, int dim, double* pos){
	Id = id;
	nDim = dim;
	const int n = nDim;
	Position = new double[n];
	Velocity = new double[n];
	FixedPos = new bool[3];
	for (int i=0; i<n; ++i){
		Position[i] = pos[i];
		Velocity[i] = 0.0;
		FixedPos[i] = false;
	}
	Viscosity = -10.0;
}

Node::~Node(){
	delete[] Position;
	delete[] Velocity;
	delete[] FixedPos;
}
