#include "ReferenceShapeBase.h"

using namespace std;

ReferenceShapeBase::ReferenceShapeBase(string ShapeType){
	setShapeType(ShapeType);
	setNodeNumber();
}

void ReferenceShapeBase::setShapeType(string TypeName){
	if (TypeName == "Prism"){
		this->ShapeType = -1;
	}
	else{
		this->ShapeType= 100;
	};
}

void ReferenceShapeBase::setNodeNumber(){
	if (ShapeType == -1){
			nNodes = 6;
	}
}
