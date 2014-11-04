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
	double* 	CurrentNormal;
	double		Volume;

	ReferenceShapeBase(string SyapeType);
	~ReferenceShapeBase();
	int getShapeType();
};

#endif
