#ifndef GrowthFunctionBase_H
#define GrowthFunctionBase_H

#include <stdio.h>
#include <iostream>
//#include <vector>
using namespace std;

class GrowthFunctionBase{
private:

public:
	GrowthFunctionBase(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane){
		this->Id = id;
		this->Type = type;
		this->initTime = initTime;
		this->endTime = endTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
		//the implementation of tilted elements is not correct as yet. Do not toggle this on as of now
		correctForTiltedElements = false;
	}
	~GrowthFunctionBase(){};

	int Type;
	int Id;
	float initTime;
	float endTime;
	bool applyToColumnarLayer;
	bool applyToPeripodialMembrane;
	bool correctForTiltedElements;


	void 	ParentErrorMessage(string functionName){
		cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	}
	double 	ParentErrorMessage(string functionName, double returnValue){
		cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
		return returnValue;
	}
	int 	ParentErrorMessage(string functionName, int returnValue){
		cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
		return returnValue;
	}
	virtual void		writeSummary(ofstream &saveFileSimulationSummary, double dt){ParentErrorMessage("writeSummary");};
	virtual void 		getCentre(float &centreX, float &centreY){ParentErrorMessage("getCentre");};
	virtual float 		getInnerRadius(){return ParentErrorMessage("getInnerRadius",0.0);};
	virtual float 		getOuterRadius(){return ParentErrorMessage("getOuterRadius",0.0);};
	virtual void 		getGrowthRate(double *maxValue){ParentErrorMessage("getGrowthRate");};
	virtual int			getGridX(){return ParentErrorMessage("getGridX",0);};
	virtual int			getGridY(){return ParentErrorMessage("getGridY",0);};
	virtual double*** 	getGrowthMatrix(){ParentErrorMessage("getGrowthMatrix");double*** a;return a;}
	virtual	double*** 	getShearValuesGrowthMatrix(){ParentErrorMessage("getGrowthMatrix");double*** a;return a;}
	virtual	double 		getGrowthMatrixElement(int i, int j, int k){return ParentErrorMessage("getGrowthMatrixElement",0.0);};
	virtual	double 		getShearValuesGrowthMatrixElement(int i, int j, int k){return ParentErrorMessage("getShearValuesGrowthMatrixElement",0.0);};
	virtual void		setGrowtRate(double ex, double ey, double ez){ParentErrorMessage("setGrowtRate");};
	virtual void		setShearValuesGrowthRate(double exy, double exz, double eyz){ParentErrorMessage("setShearValuesGrowthRate");};
	virtual void		setGrowthMatrixElement(double ex, double ey, double ez, int i, int j){ParentErrorMessage("setGrowtRate");};
	virtual void		setShearValuesGrowthMatrixElement(double exy, double exz, double eyz, int i, int j){ParentErrorMessage("setShearValuesGrowthRate");};

};
#endif

