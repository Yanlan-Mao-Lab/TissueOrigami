#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>

class Node{
private:

public:
	Node(int id, int dim, double* pos);
	~Node();
	bool *FixedPos;
	int Id;
	int nDim;
	double *Position;
	double *RKPosition;
	double **Velocity;
	double Viscosity;

};
#endif
