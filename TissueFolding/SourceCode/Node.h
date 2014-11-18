#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>

class Node{
private:

public:
	Node(int id, int dim, double* pos, int tissuePos);
	~Node();
	bool *FixedPos;
	int Id;
	int nDim;
	double *Position;
	double *RKPosition;
	double **Velocity;
	double Viscosity;
	int tissuePlacement; //1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
	double mass;
	void setViscosity(double ApicalVisc,double BasalVisc);
};
#endif
