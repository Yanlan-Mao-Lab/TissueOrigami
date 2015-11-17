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
		/**
		 *  integer id will set GrowthFunctionBase#Id. \n
		 *  integer type will set GrowthFunctionBase#Type. \n
		 *  floats initTime and endTime will set GrowthFunctionBase#initTime and GrowthFunctionBase#endTime respectively. \n
		 *	booleans applyToColumnarLayer and applyToPeripodialMembrane will set GrowthFunctionBase#applyToColumnarLayer and GrowthFunctionBase#applyToPeripodialMembrane, respectively. \n
		 *   */
		this->Id = id;
		this->Type = type;
		this->initTime = initTime;
		this->endTime = endTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
	} ///< The constructor of GrowthFunctionBase. Different growth functions will be derived from this class
	~GrowthFunctionBase(){};

	int Type; 						///< The type of the growth function, 1: uniform growth, 2: Ring shaped growth, 3: Grid based growth, where growth rates read from a separate input file
	int Id;							///< The unique identification number of the growth function
	float initTime;					///< The initiation time of the growth, in seconds
	float endTime;					///< The end time of the growth, in seconds.
	bool applyToColumnarLayer;		///< Boolean stating if the growth should be applied to columnar layer
	bool applyToPeripodialMembrane; ///< Boolean stating if the growth should be applied to peripodial membrane


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
	virtual gsl_matrix* getShearAngleRotationMatrix(){ParentErrorMessage("getShearAngleRotationMatrix");}; // this is used by unoiform growth
	virtual double 		getShearAngle(){ParentErrorMessage("getShearAngle");return 0.0;};
	virtual int			getGridX(){return ParentErrorMessage("getGridX",0);};
	virtual int			getGridY(){return ParentErrorMessage("getGridY",0);};
	virtual double*** 	getGrowthMatrix(){ParentErrorMessage("getGrowthMatrix");double*** a;return a;}
	virtual double** 	getXyShearAngleMatrix(){ParentErrorMessage("getXyShearMatrix");double** a;return a;}
	virtual	double 		getGrowthMatrixElement(int i, int j, int k){return ParentErrorMessage("getGrowthMatrixElement",0.0);};
	virtual	double 		getXyShearAngleMatrixElement(int i, int j){return ParentErrorMessage("getXyShearhMatrixElement",0.0);};
	virtual bool 		isAspectRatioOverOne(int i, int j){return ParentErrorMessage("isAspectRatioOverOne",0);};
	virtual gsl_matrix* getXyShearRotationsMatrixElement(int i, int j){ParentErrorMessage("getShearAngleRotationMatrixElement");}; //this is used by grid based growth

	virtual void		setGrowtRate(double ex, double ey, double ez){ParentErrorMessage("setGrowtRate");};
	virtual void		setGrowthMatrixElement(double ex, double ey, double ez, int i, int j){ParentErrorMessage("setGrowtRate");};

};
#endif

