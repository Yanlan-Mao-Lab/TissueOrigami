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
	double GrowthRate[3]; ///< The double array stating the uniform growth rate throughout the tissue. Growth rate is in 1/sec, format: [ DV axis (x), AP axis (y), and AB axis (z)]

	UniformGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double DVGrowth, double APGrowth, double ABGrowth) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane ){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  doubles DVGrowth, APGrowth and ABGrowth will set the  UniformGrowthFunction#GrowthRate, in the given order.
		 */
		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
	}///< The constructor of UniformGrowthFunction

	~UniformGrowthFunction(){};

	void getGrowthRate(double* maxValues){
		/**
		 *  This function will write the UniformGrowthFunction#GrowthRate of the current growth function to the input double array
		 *  pointer. The double array pointer should be set to point at a double array of size 3 (or higher) before calling the fucntion.
		 */
		maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
	} ///< The function is to get the 3D growth rate of the current growth function.
	void setGrowtRate(double ex, double ey, double ez){
		/**
		 *  This function will set the UniformGrowthFunction#GrowthRate of the current growth function to the input values
		 *  The parameters are in the order [ DV axis (x), AP axis (y), and AB axis (z)].
		 */
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	}///< The function is to set the 3D growth rate of the current growth function.
	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the UniformGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Uniform (1)
		 *			Initial time(sec): UniformGrowthFunction#initTime	FinalTime time(sec): UniformGrowthFunction#endTime	GrowthRate(fraction/hr): DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
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
	}///< The function is to write the growth function summary to simulation summary file
};


class RingGrowthFunction : public GrowthFunctionBase{
private:

public:
	double centre[2];		///< The double array of 2, giving the centre of the ring in micro-meters, format [x, y].
	double innerRadius;		///< The inner radius of the ring, inner boundary of the growth region in micro-meters. The growth will be zero at and inside the inner radius. This value can be set to zero to have circular growth.
	double outerRadius; 	///< The outer radius of the ring, outer boundary of the growth region in micro-meters. The growth will be at the maximum value set by RingGrowthFunction#GrowthRate at the outer radius.
	double GrowthRate[3]; 	///< The maximum growth rate at the RingGrowthFunction#outerRadius, in (1/sec), format: [ DV axis (x), AP axis (y), and AB axis (z)]

	RingGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double Cx, double Cy, double innerR, double outerR, double DVGrowth, double APGrowth, double ABGrowth) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  doubles Cx and Cy will set RingGrowthFunction#centre[0] and RingGrowthFunction#centre[1], respectively. \n
		 *  doubles innerR and outerR will set RingGrowthFunction#innerRadius and RingGrowthFunction#outerRadius, respectively. \n
		 *  doubles DVGrowth, APGrowth and ABGrowth will set the  RingGrowthFunction#GrowthRate, in the given order.
		 */
		this->innerRadius = innerR;
		this->outerRadius = outerR;
		this->centre[0] = Cx;
		this->centre[1] = Cy;
		GrowthRate[0] = DVGrowth;	//x
		GrowthRate[1] = APGrowth;	//y
		GrowthRate[2] = ABGrowth;	//z
	}///< The constructor of RingGrowthFunction
	~RingGrowthFunction(){};
	void 	getCentre(float &centreX, float &centreY){
		centreX = centre[0];
		centreY = centre[1];
	} ///< This function writes RingGrowthFunction#centre into the input float addresses centreX and centreY, in this order.
	float 	getInnerRadius(){
		return innerRadius;
	}///< This function returns RingGrowthFunction#innerRadius.
	float 	getOuterRadius(){
		return outerRadius;
	}///< This function returns RingGrowthFunction#outerRadius.
	void getGrowthRate(double* maxValues){
		/**
		 *  This function will write the RingGrowthFunction#GrowthRate of the current growth function to the input double array
		 *  pointer. The double array pointer should be set to point at a double array of size 3 (or higher) before calling the function.
		 */
		maxValues[0] = GrowthRate[0];
		maxValues[1] = GrowthRate[1];
		maxValues[2] = GrowthRate[2];
	}///< The function is to set the 3D maximum growth rate of the current ring growth function.
	void setGrowtRate(double ex, double ey, double ez){
		/**
		 *  This function will set the RingGrowthFunction#GrowthRate of the current growth function to the input values
		 *  The parameters are in the order [ DV axis (x), AP axis (y), and AB axis (z)].
		 */
		GrowthRate[0] = ex;
		GrowthRate[1] = ey;
		GrowthRate[2] = ez;
	}///< The function is to set the 3D maximum growth rate of the current ring growth function.
	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the RingGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  Ring (2)
		 *			Initial time(sec): RingGrowthFunction#initTime	FinalTime time(sec): RingGrowthFunction#endTime	Centre(micron): RingGrowthFunction#centre[0]  RingGrowthFunction#centre[1]	Inner radius(micron): RingGrowthFunction#innerRadius	Outer radius(micron): RingGrowthFunction#outerRadius	GrowthRate(fraction/hr):DVGrowth(in 1/hr)  APGrowth(in 1/hr) ABGrowth(in 1/hr)
		 */
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
	}///< The function is to write the growth function summary to simulation summary file
};


class GridBasedGrowthFunction : public GrowthFunctionBase{
private:

public:
	int nGridX;	///< The number of grid points that discretise the tissue in x
	int nGridY;	///< The number of grid points that discretise the tissue in y
	double ***GrowthMatrix;	///<The matrix of growth rates in (1/sec). It is a matrix of double triplets for growth rate at each grid point. The dimensions of the matrix are equal to (GridBasedGrowthFunction::nGridX, GridBasedGrowthFunction::nGridY), and set in constructor of the GridBasedGrowthFunction. The triplets store the growth rate in [ DV axis (x), AP axis (y), and AB axis (z)].
	GridBasedGrowthFunction(int id, int type, float initTime, float endTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, int nX, int nY, double*** GrowthMat) : GrowthFunctionBase(id, type, initTime, endTime, applyToColumnarLayer, applyToPeripodialMembrane){
		/**
		 *  The first six parameters will be directed to the parent constructor, GrowthFunctionBase#GrowthFunctionBase. \n
		 *  integers nX and nY will set GridBasedGrowthFunction#nGridX and GridBasedGrowthFunction#nGridY, respectively.  GridBasedGrowthFunction::GrowthMatrix will be initiated to point at a 2 dimensional matrix of double triplets the size(nX, nY). \n
		 *  double*** GrowthMat is the pointer to the 2-dimensional matrix of double triplets, holding the 3D growth rates at each grid point. Values stored in GrowthMat will set the values in GridBasedGrowthFunction::GrowthMatrix.
		 *  The matrix storing the growth rates have been read from an input file through ModelInputObject#readGrowthOptions and related functions therein \n
		 */
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
	} ///< The constructor of GridBasedGrowthFunction
	~GridBasedGrowthFunction(){
		for (int i=0; i<nGridX; ++i){
			for (int j=0; j<nGridY; ++j){
				delete[] GrowthMatrix[i][j];
			}
		}
		for (int i=0; i<nGridX; ++i){
			delete[] GrowthMatrix[i];
		}
		delete[] GrowthMatrix;
	}

	int getGridX(){
		return nGridX;
	}///< This function returns GridBasedGrowthFunction#nGridX.
	int getGridY(){
		return nGridY;
	}///< This function returns GridBasedGrowthFunction#nGridY.
	double*** getGrowthMatrix(){
		return GrowthMatrix;
	}///< This function returns GridBasedGrowthFunction#GrowthMatrix.
	double getGrowthMatrixElement(int i, int j, int k){
		return GrowthMatrix[i][j][k];
	}///< This function returns a the growth rate at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY), for the growth dimension [k] (as in  [ DV axis (x), AP axis (y), and AB axis (z)] ).

	void setGrowthMatrixElement(double ex, double ey, double ez, int i, int j){
		GrowthMatrix[i][j][0] = ex;
		GrowthMatrix[i][j][1] = ey;
		GrowthMatrix[i][j][2] = ez;
	}///< This function sets a the growth rate at grid point [i]\[j\] (in dimensions GridBasedGrowthFunction#nGridX, GridBasedGrowthFunction#nGridY), to the growth rate [ex, ey, ez] in the format [ DV axis (x), AP axis (y), and AB axis (z)].


	void writeSummary(ofstream &saveFileSimulationSummary,double dt){
		/**
		 *  This function will write the GridBasedGrowthFunction details into the simulation summary file, provided as the first input.
		 *  Time step (dt) of the simulation is provided as second input, to report the growth rates per hour.
		 *  The output should look like: \n
		 *			Growth Type:  growth From File (3)
		 *			Initial time(sec): GridBasedGrowthFunction#initTime	FinalTime time(sec): GridBasedGrowthFunction#endTime	Growth matrix mesh size: GridBasedGrowthFunction#nGridX GridBasedGrowthFunction#nGridY
		 */
		saveFileSimulationSummary<<"Growth Type:  growth From File (3)"<<endl;
		saveFileSimulationSummary<<"	Initial time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	FinalTime time(sec): ";
		saveFileSimulationSummary<<endTime;
		saveFileSimulationSummary<<"	Growth matrix mesh size: ";
		saveFileSimulationSummary<<nGridX <<" "<<nGridY<<endl;
	}///< The function is to write the growth function summary to simulation summary file
};

#endif /* GROWTHFUNCTIONTYPES_H */
