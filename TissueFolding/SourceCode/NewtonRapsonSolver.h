#ifndef NewtonRapsonSolver_H
#define NewtonRapsonSolver_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "Node.h"
#include "ShapeBase.h"

using namespace std;
/**
 *  The node class
 *  */
class NewtonRapsonSolver{
private:
	int nDim;
	int nNodes;
	double threshold;

public:
	NewtonRapsonSolver(int nDim, int nNodes);
	~NewtonRapsonSolver();

	gsl_matrix* un;
	gsl_matrix* mvisc;
	gsl_matrix* mviscPerDt;
    gsl_matrix* ge;
	gsl_matrix* gvInternal;
	gsl_matrix* gvExternal;
	gsl_matrix* gExt;
	gsl_vector* gSum;
	gsl_matrix* uk;
	gsl_matrix* displacementPerDt;
	gsl_vector* deltaU;
	gsl_matrix* K;
	gsl_matrix* Knumerical;
	bool numericalParametersSet;

	void setMatricesToZeroAtTheBeginningOfIteration(bool thereIsNumericalCalculation);
	void setMatricesToZeroInsideIteration();

	void constructUnMatrix(vector <Node*>& Nodes);
	void initialteUkMatrix();
	void constructLumpedMassExternalViscosityDtMatrix(vector <Node*>& Nodes, double dt);
	void calculateDisplacementMatrix(double dt);
	void calcutateFixedK(vector <Node*>& Nodes);
	void calculateForcesAndJacobianMatrixNR(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double dt, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile);
	void writeForcesTogeAndgvInternal(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double** SystemForces);
	void writeImplicitElementalKToJacobian(vector <ShapeBase*>& Elements);
	void calculateExternalViscousForcesForNR();
	void addImplicitKViscousExternalToJacobian();
	void checkJacobianForAblatedNodes(vector <int> & AblatedNodes);
	void calculateSumOfInternalForces();
	void addExernalForces();

	void solveForDeltaU();
	int  solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables);
	void constructiaForPardiso(int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec);
	void writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a);
	void writeginPardisoFormat(double* b, const int n);
	bool checkConvergenceViaDeltaU();
	bool checkConvergenceViaForce();
	void updateUkInIteration();
	void useNumericalJacobianInIteration();
	void calculateDifferenceBetweenNumericalAndAnalyticalJacobian(vector <Node*>& Nodes, bool displayMAtricesDuringNumericalCalculation);

	void displayMatrix(gsl_matrix* mat, string matname);

};

#endif
