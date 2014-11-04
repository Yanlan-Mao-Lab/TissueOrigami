#include "Node.h"

using namespace std;

Node::Node(int id, int dim, double* pos){
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
		RKPosition[i] = 0.0;
	}
	FixedPos = new bool[3];
	for (int i=0; i<3; ++i){
		FixedPos[i] = false;
	}
	Viscosity = -10.0;
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
