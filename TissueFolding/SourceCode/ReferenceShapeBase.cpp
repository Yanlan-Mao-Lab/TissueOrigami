#include "ReferenceShapeBase.h"

using namespace std;

ReferenceShapeBase::ReferenceShapeBase(string ShapeType){
	height = -100;
	setShapeType(ShapeType);
	setNodeNumber();
}

ReferenceShapeBase::~ReferenceShapeBase(){
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
	}
	delete[] Positions;
	delete[] CurrentNormal;
}

void ReferenceShapeBase::setShapeType(string TypeName){
	if (TypeName == "Prism"){
		this->ShapeType = -1;
	}
	else if (TypeName == "PrismLateral"){
		this->ShapeType = -2;
	}
	else if (TypeName == "Tetrahedron"){
		this->ShapeType = -3;
	}
	else if (TypeName == "Triangle"){
		this->ShapeType = -4;
	}
	else{
		this->ShapeType= 100;
	};
}

void ReferenceShapeBase::setNodeNumber(){
	if (ShapeType == -1){
			nNodes = 6;
	}
	if (ShapeType == -2){
			nNodes = 6;
	}
	if (ShapeType == -3){
			nNodes = 4;
	}
	if (ShapeType == -4){
		nNodes = 3;
	}
}
