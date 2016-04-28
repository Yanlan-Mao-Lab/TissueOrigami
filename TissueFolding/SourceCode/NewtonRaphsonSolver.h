#ifndef NewtonRaphsonSolver_H
#define NewtonRaphsonSolver_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_linalg.h>
#include "Node.h"
#include "ShapeBase.h"

using namespace std;
/**
 *  The Newton-Raphson solver class
 *  */
class NewtonRaphsonSolver{
private:
	int nDim;			//< Dimension of the space (3D)
	int nNodes;			//< Number of nodes of the system
	double threshold;	//< Convergence threshold for iterations

public:
	NewtonRaphsonSolver(int nDim, int nNodes); 	//<Constructer of the N-R solver
	~NewtonRaphsonSolver();						//<Desturctor of the N-R solver

	gsl_matrix* un;								//< The initial positions of the nodes, as calculated at the end of previous step "n"
	gsl_matrix* mvisc;							//< The matrix containing ( mass * external viscosity ) term for each node, of size (nDim x nNodes , nDim*nNodes). This is a diagonal matrix.
	gsl_matrix* mviscPerDt;						//< The matrix containing ( mass * external viscosity / time step ) term for each node, of size (nDim x nNodes , nDim*nNodes). This is a diagonal matrix.
    gsl_matrix* ge;								//< The matrix containing elastic forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
	gsl_matrix* gvInternal;						//< The matrix containing internal viscous forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
	gsl_matrix* gvExternal;						//< The matrix containing external viscous forces on each node, size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
	gsl_matrix* gExt;							//< The matrix containing external forces on each node (currently includes packing and myosin forces), size (nDim*nNodes,1). Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
	gsl_vector* gSum;							//< The matrix containing sum of NewtonRaphsonSolver#ge, NewtonRaphsonSolver#gvInternal, NewtonRaphsonSolver#gvExternal, NewtonRaphsonSolver#gExt. Organisation is [Node0,x ; Node0,y ; Node0,z; ... ; Noden,x ; Noden,y ; Noden,z]
	gsl_matrix* uk;								//< The matrix storing the position of each node at iteration "k". Initiated in function NewtonRaphsonSolver#initialteUkMatrix at the beginning of each step, and updated by function NewtonRaphsonSolver#updateUkInIteration during the iteartions.
	gsl_matrix* displacementPerDt;				//< The displacement per time step of each node in current iteration "k", from its position at the end of the last time step "n"
	gsl_vector* deltaU;							//< The incremental change in positions as calculated in current iteration, resulting from the imbalance of elastic, viscous and any other external forces acting on each nodes. The solver minimises this value, convergence occurs when all incremental movements for all nodes sufficiently close to zero.
	gsl_matrix* K;								//< The Jacobian matrix, derivative of sum of forces acting on each node with respect to displacements.
	gsl_matrix* Knumerical;						//< The Jacobian calculated by numerical methods
	bool numericalParametersSet;				//< The boolean stating if the numerical Jacobian matrix has been initiated or not. For memory management purposes, this large matrix is not initiated unless necessary.

	void setMatricesToZeroAtTheBeginningOfIteration(bool thereIsNumericalCalculation); 		//< The function setting the calculation matrices to zero at the beginning of each iteration.
	void setMatricesToZeroInsideIteration();												//< The function setting the relevant matrices to zero at each iteration.

	void constructUnMatrix(vector <Node*>& Nodes);											//< This function constructs NewtonRaphsonSolver#un matrix at the beginning of the iterations.
	void initialteUkMatrix();																//< This function initiates NewtonRaphsonSolver#uk matrix at the beginning of the iterations, it is initiated to be equal to NewtonRaphsonSolver#un.
	void constructLumpedMassExternalViscosityDtMatrix(vector <Node*>& Nodes, double dt);	//< This function constructs NewtonRaphsonSolver#mvisc and  NewtonRaphsonSolver#mviscPerDt for external viscosity related calculations.
	void calculateDisplacementMatrix(double dt);											//< This function calculates the displacement of each node in current iteration "k", from their positions at the end of the previous step "n" (NewtonRaphsonSolver#uk - NewtonRaphsonSolver#un)
	void calcutateFixedK(vector <Node*>& Nodes);											//< This function updates the Jacobian to account for nodes  that are fixed in certain dimensions in space, as part of boundary conditions.
	void calculateForcesAndJacobianMatrixNR(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double dt, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile);	//< This function calculates elemental forces and Jacobians, later to be combined in NewtonRaphsonSolver#K and NewtonRaphsonSolver#gSum
	void writeForcesTogeAndgvInternal(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double** SystemForces);	//< This function writes the values of elemental elastic (ShapeBase#ge) and internal viscous forces (ShapeBase#gvInternal) into the system elastic and internal viscous forces, NewtonRaphsonSolver#ge, and NewtonRaphsonSolver#gvInternal, respectively.
	void writeImplicitElementalKToJacobian(vector <ShapeBase*>& Elements);	//< This function writes the elemental values for elastic part of the Jacobian - stiffness matrix - (ShapeBase#TriPointKe) and for viscous part of Jacobian (ShapeBase#TriPointKv) into the system Jacobian NewtonRaphsonSolver#K.
	void calculateExternalViscousForcesForNR();		//< This function calculates the external viscous forces acting on each node, the values are sotred in NewtonRaphsonSolver#gvExternal
	void addImplicitKViscousExternalToJacobian();	//< This function adds the external related terms of the Jacobian to the system Jacobian NewtonRaphsonSolver#K.
	void checkJacobianForAblatedNodes(vector <int> & AblatedNodes);	//< This functions checks the Jacobian to ensure the diagonal terms are non-zero for ablated nodes.
	void calculateSumOfInternalForces();			//< This function adds the ealsticity and viscosity related forces (NewtonRaphsonSolver#ge, NewtonRaphsonSolver#gvInternal, NewtonRaphsonSolver#gvExternal) to sum of forces, NewtonRaphsonSolver#gSum.
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
