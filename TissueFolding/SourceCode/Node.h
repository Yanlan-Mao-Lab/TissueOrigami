#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

class Node{
private:

public:
	Node(int id, int dim, double* pos, int tissuePos, int tissueType);
	~Node();
	bool *FixedPos;
	int Id;
	int nDim;
	double *Position;
	double *RKPosition;
	double **Velocity;
	double Viscosity;
	int tissuePlacement;	//1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
	int tissueType;		 	//0 -> columnar layer, 1->peripodium
	bool atCircumference;
	double mass;
	bool atPeripodiumCircumference;
	vector <int> AssociatedNodesDueToPeripodium;
	vector <double> AssociatedNodeWeightsDueToPeripodium;
	vector <int> immediateNeigs;
	vector <int> connectedElementIds;
	vector <double> connectedElementWeights;

	void setViscosity(double ApicalVisc,double BasalVisc);
	bool checkIfNeighbour(int IdToCheck);
	void displayConnectedElementIds();
	void displayConnectedElementWeights();
};
#endif
