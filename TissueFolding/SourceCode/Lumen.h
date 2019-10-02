/*
 * Lumen.h
 *
 *  Created on: 22 Jun 2018
 *      Author: melda
 */

#ifndef LUMEN_H_
#define LUMEN_H_

#include <vector>
//#include <math.h>
#include "ShapeBase.h"

class Lumen{

protected:
   //double getApicalSideLengthAverage();
	int nTriangleSize;	//the number of triangles forming the surface of lumen

public:
	Lumen(vector<ShapeBase*>& Elements,vector<Node*>& Nodes, double lumenBulkModulus, double lumenGrowthFold);
	~Lumen();

	void updateMatrices(vector<Node*>& Nodes);
	void calculateCurrentVolume();
	void calculateResiduals(vector<Node*>& Nodes);
	void calculateLumengFromElementalResiduals(gsl_matrix* g);
	void calculateJacobian();
	void writeLumenJacobianToSystemJacobian(gsl_matrix* K,vector<Node*>& Nodes);
	void growLumen(double currentTimeInSec);


	int Dim ;		//dimensions of the system, should be 3D.
	double rV;	//normalised deviation from ideal volume (V - V0 )/V0
	vector <ShapeBase*> encapsulatingElements;
	vector <int> nodeIdsList;
	double growthRate;	//growth rate per sec
	double initialIdealVolume;
	double currentIdealVolume;
	double bulkModulus;
	double currentVolume;
	gsl_matrix* Kv;
	gsl_matrix* KvNumerical;
	vector <gsl_matrix*> xCap1;
	vector <gsl_matrix*> xCap2;
	vector <gsl_matrix*> xCap3;
	vector <gsl_matrix*> x1;
	vector <gsl_matrix*> x2;
	vector <gsl_matrix*> x3;
	vector <gsl_matrix*> g1;
	vector <gsl_matrix*> g2;
	vector <gsl_matrix*> g3;

	void 	displayMatrix(gsl_matrix* mat, string matname){
	    int m = mat->size1;
	    int n = mat->size2;
	    cout<<matname<<": "<<endl;

	    for (int i =0; i<m; i++){
	        for (int j =0; j<n; j++){
	            cout.precision(4);
	            cout.width(6);
	            cout<<gsl_matrix_get(mat,i,j)<<" ";
	        }
	        cout<<endl;
	    }
	    cout<<endl;
	}

};



#endif /* LUMEN_H_ */
