/*
 * GrowthFunctionTypes.h
 *
 *  Created on: 24 Apr 2015
 *      Author: melda
 */

#ifndef GROWTHFUNCTIONTYPES_H
#define GROWTHFUNCTIONTYPES_H

#include "GrowthFunctionBase.h"
using namespace std;

class UniformGrowthFunction : public GrowthFunctionBase{
private:

public:
	UniformGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double DVGrowth, double APGrowth, double ABGrowth) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane ){
		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
		//this will be necessary for personalised growth rate storing in tilted elements
		ShearValuesGrowthRate[0] = 0.0;
		ShearValuesGrowthRate[1] = 0.0;
		ShearValuesGrowthRate[2] = 0.0;
		correctForTiltedElements = false;
	}
	~UniformGrowthFunction(){};
	double GrowthRate[3]; //The uniform growth rate throughout the tissue, in [ DV axis (x), AP axis (y), and AB axis (z)]
	double ShearValuesGrowthRate[3];

	void getGrowthRate(double* maxValues){
		maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
		maxValues[3] = ShearValuesGrowthRate[0];
		maxValues[4] = ShearValuesGrowthRate[1];
		maxValues[5] = ShearValuesGrowthRate[2];
	}
	void setGrowtRate(double ex, double ey, double ez){
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	};
	void setShearValuesGrowthRate(double exy, double exz, double eyz){
		ShearValuesGrowthRate[0] = exy;
		ShearValuesGrowthRate[1] = exz;
		ShearValuesGrowthRate[2] = eyz;
	};
	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		saveFileSimulationSummary<<"Growth Type:  Uniform (1)"<<endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
		saveFileSimulationSummary<<endl;
	}
};


class RingGrowthFunction : public GrowthFunctionBase{
private:

public:
	RingGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double Cx, double Cy, double innerR, double outerR, double DVGrowth, double APGrowth, double ABGrowth) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane){
		this->innerRadius = innerR;
		this->outerRadius = outerR;
		this->centre[0] = Cx;
		this->centre[1] = Cy;
		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
		//this will be necessary for personalised growth rate storing in tilted elements
		ShearValuesGrowthRate[0] = 0.0;
		ShearValuesGrowthRate[1] = 0.0;
		ShearValuesGrowthRate[2] = 0.0;
	}
	~RingGrowthFunction(){};
	double centre[2];
	double innerRadius;
	double outerRadius;
	double GrowthRate[3]; //The uniform growth rate throughout the tissue, in [ DV axis (x), AP axis (y), and AB axis (z)]
	double ShearValuesGrowthRate[3];
	void 	getCentre(float &centreX, float &centreY){
		centreX = centre[0];
		centreY = centre[1];
	}
	float 	getInnerRadius(){
		return innerRadius;
	}
	float 	getOuterRadius(){
		return outerRadius;
	}
	void getGrowthRate(double* maxValues){
		maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
		maxValues[3] = ShearValuesGrowthRate[0];
		maxValues[4] = ShearValuesGrowthRate[1];
		maxValues[5] = ShearValuesGrowthRate[2];
	}
	void setGrowtRate(double ex, double ey, double ez){
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	};
	void setShearValuesGrowthRate(double exy, double exz, double eyz){
		ShearValuesGrowthRate[0] = exy;
		ShearValuesGrowthRate[1] = exz;
		ShearValuesGrowthRate[2] = eyz;
	};
	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		saveFileSimulationSummary<<"Growth Type:  Ring (2)"<<endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	Centre(micron): ";
		saveFileSimulationSummary<<centre[0];
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<centre[1];
		saveFileSimulationSummary<<"	Inner radius(micron): ";
		saveFileSimulationSummary<<innerRadius;
		saveFileSimulationSummary<<"	Outer radius(micron): ";
		saveFileSimulationSummary<<outerRadius;
		saveFileSimulationSummary<<"	GrowthRate(fraction/hr): ";
		saveFileSimulationSummary<<GrowthRate[0]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[1]/dt*3600.0;
		saveFileSimulationSummary<<"  ";
		saveFileSimulationSummary<<GrowthRate[2]/dt*3600.0;
		saveFileSimulationSummary<<endl;
	}
};


class GridBasedGrowthFunction : public GrowthFunctionBase{
private:

public:
	GridBasedGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, int nX, int nY, double*** GrowthMat) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane){
		this ->nGridX = nX;
		this ->nGridY = nY;
		GrowthMatrix = new double**[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			GrowthMatrix[i] = new double*[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				GrowthMatrix[i][j] = new double[3];
				for (int k=0; k<3; ++k){
					GrowthMatrix[i][j][k] = GrowthMat[i][j][k];
				}
			}
		}
		//this will be necessary for personalised growth rate storing in tilted elements
		ShearValuesGrowthMatrix = new double**[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			ShearValuesGrowthMatrix[i] = new double*[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				ShearValuesGrowthMatrix[i][j] = new double[3];
				for (int k=0; k<3; ++k){
					ShearValuesGrowthMatrix[i][j][k] = 0.0;
				}
			}
		}
	}
	~GridBasedGrowthFunction(){
		for (int i=0; i<nGridX; ++i){
			for (int j=0; j<nGridY; ++j){
				delete[] GrowthMatrix[i][j];
				delete[] ShearValuesGrowthMatrix[i][j];
			}
		}
		for (int i=0; i<nGridX; ++i){
			delete[] GrowthMatrix[i];
			delete[] ShearValuesGrowthMatrix[i];
		}
		delete[] GrowthMatrix;
		delete[] ShearValuesGrowthMatrix;
	}
	int nGridX;
	int nGridY;
	double ***GrowthMatrix;
	double ***ShearValuesGrowthMatrix;

	int getGridX(){
		return nGridX;
	}
	int getGridY(){
		return nGridY;
	}
	double*** getGrowthMatrix(){
		return GrowthMatrix;
	}
	double*** getShearValuesGrowthMatrix(){
		return ShearValuesGrowthMatrix;
	}
	double getGrowthMatrixElement(int i, int j, int k){
		return GrowthMatrix[i][j][k];
	}
	double getShearValuesGrowthMatrixElement(int i, int j, int k){
		return ShearValuesGrowthMatrix[i][j][k];
	}
	void setGrowthMatrixElement(double ex, double ey, double ez, int i, int j){
		GrowthMatrix[i][j][0] = ex;
		GrowthMatrix[i][j][1] = ey;
		GrowthMatrix[i][j][2] = ez;
	};
	void setShearValuesGrowthMatrixElement(double exy, double exz, double eyz, int i, int j){
		ShearValuesGrowthMatrix[i][j][0] = exy;
		ShearValuesGrowthMatrix[i][j][1] = exz;
		ShearValuesGrowthMatrix[i][j][2] = eyz;
	};

	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		saveFileSimulationSummary<<"Growth Type:  growth From File (3)"<<endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	Growth matrix mesh size: ";
		saveFileSimulationSummary<<nGridX <<" "<<nGridY<<endl;
	}
};

#endif /* GROWTHFUNCTIONTYPES_H */
