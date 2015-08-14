#ifndef ReferenceShapeBase_H
#define ReferenceShapeBase_H

#include <iostream>
#include <stdio.h>
#include <string>


using namespace std;

class ReferenceShapeBase{
protected:
	int ShapeType;
	int Id;
	int nNodes;


	void setShapeType(string TypeName);
	void setNodeNumber();
public:
	double** 	Positions;
	double		Volume;
	double 		BasalArea;
	double 		height; //slab height for 2D elements, value is -100 for 3D elements
	ReferenceShapeBase(string SyapeType);
	~ReferenceShapeBase();
	int getShapeType();
};

#endif
