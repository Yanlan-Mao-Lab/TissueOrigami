/*
 * NewtonRaphsonSolver.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: melda
 */

#include "NewtonRaphsonSolver.h"
#include <math.h>

//#include "Node.h"
//#include <gsl/gsl_linalg.h>

using namespace std;

NewtonRaphsonSolver::NewtonRaphsonSolver(int dim, int n){
	/**
	*  The function initiates (memory allocates) the necessary matrices for
	*  the numerical solver.
	*
	*/
	threshold = 1E-8;
	nDim = dim;
	nNodes = n;
	externalViscosityVolumeBased = false; //using external surface
	numericalParametersSet = false;
	thereIsLumen = false;

	un = gsl_matrix_calloc(nDim*nNodes,1);
	ge = gsl_matrix_calloc(nDim*nNodes,1);
	gvInternal = gsl_matrix_calloc(nDim*nNodes,1);
	gvExternal = gsl_matrix_calloc(nDim*nNodes,1);
	gExt = gsl_matrix_calloc(nDim*nNodes,1);
	gSum = gsl_vector_calloc(nDim*nNodes);
	uk = gsl_matrix_calloc(nDim*nNodes,1);
	displacementPerDt = gsl_matrix_calloc(nDim*nNodes,1);
	deltaU = gsl_vector_calloc(nDim*nNodes);
	K = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	Knumerical = 0;
	boundNodesWithSlaveMasterDefinition = false;
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(){
	gsl_matrix_free(uk);
	std::cout<<" deleted uk"<<std::endl;
   	gsl_matrix_free(un);
	std::cout<<" deleted un"<<std::endl;
	gsl_matrix_free(ge);
	std::cout<<" deleted ge"<<std::endl;
	gsl_matrix_free(gvExternal);
	std::cout<<" deleted gvExternal"<<std::endl;
	gsl_matrix_free(gvInternal);
	std::cout<<" deleted gvInternal"<<std::endl;
	gsl_matrix_free(gExt);
	std::cout<<" deleted gExt"<<std::endl;
	gsl_vector_free(gSum);
	std::cout<<" deleted gSum"<<std::endl;
	gsl_matrix_free(displacementPerDt);
	std::cout<<" deleted displacementPerDt"<<std::endl;
	gsl_vector_free(deltaU);
	std::cout<<" deleted deltaU"<<std::endl;
	gsl_matrix_free(K);
	std::cout<<" deleted K"<<std::endl;
	if (Knumerical != 0){
		gsl_matrix_free(Knumerical);
	}
	else{
		free(Knumerical);
	}
	std::cout<<" deleted Knumerical"<<std::endl;
}

void NewtonRaphsonSolver::setMatricesToZeroAtTheBeginningOfIteration(bool thereIsNumericalCalculation){
	gsl_matrix_set_zero(un);
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_matrix_set_zero(gExt);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(uk);
	gsl_matrix_set_zero(displacementPerDt);
	gsl_vector_set_zero(deltaU);
	gsl_matrix_set_zero(K);
	if (thereIsNumericalCalculation){
		if (!numericalParametersSet){
			Knumerical = gsl_matrix_calloc(K->size1,K->size2);
			numericalParametersSet = true;
		}
		else{
			gsl_matrix_set_zero(Knumerical);
		}
	}
}

void NewtonRaphsonSolver::setMatricesToZeroInsideIteration(){
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(gExt);
	gsl_matrix_set_zero(K);


}

void NewtonRaphsonSolver::constructUnMatrix(const std::vector<std::unique_ptr<Node>>& Nodes){
    for (size_t i = 0; i<nNodes; ++i ){
        for (size_t j=0; j<nDim; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void NewtonRaphsonSolver::initialteUkMatrix(){
    gsl_matrix_memcpy(uk,un);
}


void NewtonRaphsonSolver::calculateDisplacementMatrix(double dt){
	/**
	* The displacement of nodes per Simulation#dt, for each iteration defined as \f$ \frac{\boldsymbol{u_k} - \boldsymbol{u_n}}{dt} \f$
	*/
	gsl_matrix_memcpy(displacementPerDt,uk);
	gsl_matrix_sub(displacementPerDt,un);
	gsl_matrix_scale(displacementPerDt,1.0/dt);
}

void NewtonRaphsonSolver::calculateBoundKWithSlavesMasterDoF(){
	if (boundNodesWithSlaveMasterDefinition){
        /** The matrix calculation for node binding is in the form: \n
         * \f{eqnarray*}
          \boldsymbol{K}_{bound} & = & \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} + \boldsymbol{\bar{I}}\\
          \boldsymbol{g}_{bound} & = & \boldsymbol{N^{T}}  \boldsymbol{g} \\
          \boldsymbol{N_{ij}} &= & \left\{ \begin{matrix}
                                        1 & if		& i = j \text{ and } i \text{ is not a slave}\\
                                        1 & if 	& i  \text{ is slave to } j \\
                                        0 &		& \text{elsewhere}
                                    \end{matrix} \right. \\
          \boldsymbol{\bar{I}_{ij}} & = &  \left\{ \begin{matrix}
                                             1 & if		& i = j \text{ and } i \text{ is a slave}\\
                                             0 &		& \text{elsewhere}
                                           \end{matrix} \right. \\
         * \f} \n
         * As this actual matrix calculation requires a significant memory allocation, these operations are
         *carried out on a row/column basis, rather than using full matrices. \n
         */
		int totalNumberOfDoF = K->size1; //nDim * nNodes
		for (vector< vector<int> >::iterator itSlaveMasterCouple=slaveMasterList.begin(); itSlaveMasterCouple<slaveMasterList.end(); ++itSlaveMasterCouple){
			/** First all forces on slave to master, set slave forces to zero on NewtonRaphsonSolver#gSum,
			*equivalent of \f$ \boldsymbol{g}_{bound} = \boldsymbol{N^{T}}  \boldsymbol{g} \f$.
			*/
			int slaveIndex = (*itSlaveMasterCouple)[0];
			int masterIndex = (*itSlaveMasterCouple)[1];
			double slaveForce = gsl_vector_get(gSum,slaveIndex);
			double masterForce = gsl_vector_get(gSum,masterIndex);
			masterForce +=  slaveForce;
			gsl_vector_set(gSum,slaveIndex,0);
			gsl_vector_set(gSum,masterIndex,masterForce);
			/** Then start the manipulation of the Jacobian, by adding all elements of the slave degrees of freedom
			* row to master degrees of freedom row on NewtonRaphsonSolver#K, equivalent of operation \f$ \boldsymbol{N^{T}} \boldsymbol{K} \f$.
			*/
             		// format of view: [ matrix*, origin i, origin j, rows, columns ]
			gsl_matrix_view KSlaveRow = gsl_matrix_submatrix (K, slaveIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_view KMasterRow = gsl_matrix_submatrix (K, masterIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_add(&(KMasterRow.matrix),&(KSlaveRow.matrix));
			gsl_matrix_set_zero(&(KSlaveRow.matrix));
			//add all elements of the slave column to master column - N^T K N
			/** Followed by adding all elements of the slave degrees of freedom column to master degrees of freedom
			*  columno n NewtonRaphsonSolver#K, equivalent of operation  (cumulatively) \f$ \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} \f$.
			*/
			gsl_matrix_view KSlaveColumn = gsl_matrix_submatrix (K, 0, slaveIndex, totalNumberOfDoF,1);
			gsl_matrix_view KMasterColumn = gsl_matrix_submatrix (K,0, masterIndex, totalNumberOfDoF,1);
			gsl_matrix_add(&(KMasterColumn.matrix),&(KSlaveColumn.matrix));
			gsl_matrix_set_zero(&(KSlaveColumn.matrix));
			/** Finally, make the diagonal element of NewtonRaphsonSolver#K, \f$ \boldsymbol{K}_{DOFslave,DOFslave} \f$ to unity,
			* equivalent of operation  (cumulatively)
			* \f$ \boldsymbol{N^{T}} \boldsymbol{K} \boldsymbol{N} + \boldsymbol{\bar{I}}\f$.
			*/
			gsl_matrix_set(K,slaveIndex,slaveIndex,1);
			//views do not need to be freed, they are addresses to original matrices
		}
	}
}

void NewtonRaphsonSolver::equateSlaveDisplacementsToMasters(){
	int n = slaveMasterList.size();
	for( int slaveIterator = 0; slaveIterator <n; ++slaveIterator){
		int slaveIndex = slaveMasterList[slaveIterator][0];
		int masterIndex = slaveMasterList[slaveIterator][1];
		double masterDisplacement = gsl_vector_get(deltaU,masterIndex);
		gsl_vector_set(deltaU,slaveIndex,masterDisplacement);
	}
}

void NewtonRaphsonSolver::calcutateFixedK(const std::vector<std::unique_ptr<Node>>& Nodes){
    /** Some degrees of freedom are fixed for some nodes, as defined by the user input boundary conditions.
     * this will be recorded in the 3 dimensional boolean array of each node Node#FixedPos
     * for x,y and z coordinates. In the node has a fixed degree of freedom, then the  sum of
     * elastic and viscous forces recorded on NewtonRaphsonSolver#gSum is made zero. Then in the Jacobian, the
     * diagonal term for the degree of freedom is set to unity, all remaining terms of the column and row of
     * the degree of freedom is set to zero.
     */
    size_t dim = 3;
    size_t Ksize = K->size1;
    for(size_t i=0; i<nNodes; i++){
        for (size_t j=0; j<dim; ++j){
            if (Nodes[i]->FixedPos[j]){
                size_t index1 = i*dim+j;
                gsl_vector_set(gSum,index1,0.0); // making the forces zero
                for (size_t k =0; k<Ksize; ++k){
                    double value =0.0;
                    if (index1 == k ){value =1.0;}
                    gsl_matrix_set(K, index1, k, value);
                    gsl_matrix_set(K, k, index1, value); //K is symmetric;
                }
            }
        }
    }
}

void NewtonRaphsonSolver::calculateForcesAndJacobianMatrixNR(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements, double dt){
	#ifndef DO_NOT_USE_OMP
	/** If DO_NOT_USE_OMP is not defined,I will be using omp. This
	* is necessary as omp is not set up on mac
	*/
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	#endif
	for(std::vector<std::unique_ptr<ShapeBase>>::const_iterator itElement = Elements.begin(); itElement < Elements.end(); ++itElement){
	   	/** The calculation of foces and their derivatives in each element starts with
		 * calculation of forces via ShapeBase#calculateForces. A series of calculations necessary for
		 * Jacobian calculation are obtained at this stage. Then the Elastic and Viscous parts of the elemetnal Jacobian
		 * are calculates through ShapeBase#calculateImplicitKElastic and ShapeBase#calculateImplicitKViscous ,respectively
		  */
          	if (!(*itElement)->IsAblated){
          		(*itElement)->calculateForces(Nodes, displacementPerDt);
		}
          	(*itElement)->calculateImplicitKElastic(); //This is the stiffness matrix, elastic part of the jacobian matrix
          	(*itElement)->calculateImplicitKViscous(displacementPerDt, dt); //This is the viscous part of jacobian matrix
	}
	if (thereIsLumen){
		tissueLumen->updateMatrices(Nodes, Elements);
		tissueLumen->calculateCurrentVolume();
		tissueLumen->calculateResiduals(Nodes, Elements);
		//calculateLumenNumericalJacobian(tissueLumen,Nodes,Elements);
		tissueLumen->calculateJacobian(Elements);
		tissueLumen->writeLumenJacobianToSystemJacobian(K,Nodes);
	}
}


void NewtonRaphsonSolver::updateElementPositions(const std::vector <std::unique_ptr<Node>>&  Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	for(const auto& itElement : Elements){
		itElement->updatePositions(Nodes);
	}
}

void 	NewtonRaphsonSolver::calculateLumenNumericalJacobian(Lumen* tissueLumen, const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	gsl_matrix_set_zero(tissueLumen->Kv);
	double perturbation = 1E-6;
	const int nEle = tissueLumen->encapsulatingElementIds.size();
	//take backup of elemental elastic system forces and reset to zero
	double backupElasticForces [nEle][6][3];
	for (int i=0; i<nEle; ++i){
		for (int j=0; j<6; j++){
			for (int k =0; k<3;++k){
				int eleId = tissueLumen->encapsulatingElementIds[i];
				backupElasticForces[i][j][k] = Elements[eleId]->getElementalElasticForce(j,k);
				Elements[eleId]->setElementalElasticForce(j,k,0.0);
			}
		}
	}
	int  nNode = tissueLumen->nodeIdsList.size();
	gsl_matrix* gOriginal = gsl_matrix_calloc(tissueLumen->Dim*nNode,1);
	double rVover6V0 =  tissueLumen->rV /6.0 / tissueLumen->currentIdealVolume;
	tissueLumen->calculateLumengFromElementalResiduals(gOriginal, Elements);
	gsl_matrix_scale(gOriginal,-1.0*tissueLumen->bulkModulus*rVover6V0);

	gsl_matrix* gPerturbed = gsl_matrix_calloc(tissueLumen->Dim*nNode,1);
	for (int nodeIndex = 0; nodeIndex<nNode; ++nodeIndex){
		int currId = tissueLumen->nodeIdsList[nodeIndex];
		for (int currDim = 0; currDim <3 ; ++currDim){
			if (Nodes[currId]->FixedPos[currDim] ){
				continue;
			}
			//std::cout<<" perturbing node "<<currId<<" at dim: "<<currDim<<std::endl;
			Nodes[currId]->Position[currDim]  += perturbation;
			updateElementPositions(Nodes, Elements);
			tissueLumen->updateMatrices(Nodes, Elements);
			tissueLumen->calculateCurrentVolume();
			double rVover6V0 =  tissueLumen->rV /6.0 / tissueLumen->currentIdealVolume;			
			//reset perturbed g:
			gsl_matrix_set_zero(gPerturbed);
			tissueLumen->calculateResiduals(Nodes, Elements);
			tissueLumen->calculateLumengFromElementalResiduals(gPerturbed, Elements);
			gsl_matrix_scale(gPerturbed,-1.0*tissueLumen->bulkModulus*rVover6V0);
			//writing the perturbation to global jacobian
			gsl_matrix_sub(gPerturbed,gOriginal);
			gsl_matrix_scale(gPerturbed,1.0/perturbation);
			for (int affectedNodeIndex = 0; affectedNodeIndex<nNode; ++affectedNodeIndex){
				for (int affectedDim = 0; affectedDim <3 ; affectedDim++){
					double numericalValueAffectedDim = -1*gsl_matrix_get(gPerturbed,affectedNodeIndex*3+affectedDim,0);
					int affectedNodeId = tissueLumen->nodeIdsList[affectedNodeIndex];
					if (Nodes[affectedNodeId]->FixedPos[affectedDim] ){
						continue;
					}
					int currIndexOnLumenK = nodeIndex*3 + currDim;
					int affectedIndexOnLumenK = affectedNodeIndex*3 + affectedDim;
					double valueOnLumenK = gsl_matrix_get(tissueLumen->Kv,currIndexOnLumenK,affectedIndexOnLumenK);
					double newValue = valueOnLumenK+numericalValueAffectedDim;
					gsl_matrix_set(tissueLumen->Kv,currIndexOnLumenK,affectedIndexOnLumenK,newValue);
					if (currIndexOnLumenK == 0 && affectedIndexOnLumenK == 7){
						std::cout<<"numerical ["<<currIndexOnLumenK<<", "<<affectedIndexOnLumenK<<"]: value to add"<<numericalValueAffectedDim<<" valueOnLumenK "<<valueOnLumenK<<std::endl;
					}
				}
			}
			//end of writing the perturbation
			//revert the perturbation, no need to update positions here as the loop continues
			Nodes[currId]->Position[currDim] -= perturbation;
		}
	}
	//outside the loop, correct positions and all calculated values:
	updateElementPositions(Nodes, Elements);
	tissueLumen->updateMatrices(Nodes, Elements);
	tissueLumen->calculateCurrentVolume();
	//set back the elemental elastic forces:
	for (int i=0; i<nEle; ++i){
		for (int j=0; j<6; j++){
			for (int k =0; k<3;++k){
				int eleId = tissueLumen->encapsulatingElementIds[i];
				Elements[eleId]->setElementalElasticForce(j,k,backupElasticForces[i][j][k]);
			}
		}
	}
	tissueLumen->calculateResiduals(Nodes, Elements);
	//free memory
	gsl_matrix_free(gOriginal);
	gsl_matrix_free(gPerturbed);
}

void NewtonRaphsonSolver::writeForcesTogeAndgvInternal(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements, std::vector<std::array<double,3>>& SystemForces){
    for(const auto&  itElement : Elements){
        if (!itElement->IsAblated){
        	itElement->writeInternalForcesTogeAndgv(ge,gvInternal,SystemForces,Nodes);
        }
    }
}

void NewtonRaphsonSolver::writeImplicitElementalKToJacobian(const std::vector <std::unique_ptr<ShapeBase>>& Elements){
    /** All elemental elastic (ShapeBase#Ke) and internal viscous (ShapeBase#Kv) jacobians are added onto the
     * system Jacobian (NewtonRaphsonSolver#K) in this function.
     */
	for(const auto&  itElement : Elements){
		itElement->writeKviscousToMainKatrix(K);
	}
	for(const auto&  itElement : Elements){
		itElement->writeKelasticToMainKatrix(K);
	}
}

void NewtonRaphsonSolver::calculateExternalViscousForcesForNR(const std::vector<std::unique_ptr<Node>>&  Nodes){
    //the mass is already updated for symmetricity boundary nodes, and the viscous forces will be calculated correctly,
	//as mvisdt is already accounting for the doubling of mass
    //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mvisc, displacementPerDt,0.0, gvExternal);
	for (size_t i = 0; i<nNodes; ++i ){
		for (size_t j=0; j<nDim; ++j){
			/** For all nodes of the system, the external viscous forces will depend on the
			* exposed surface associated with each node (Node#viscositySurface),
			* calculated through ShapeBase::assignViscositySurfaceAreaToNodes. This surface
			* will be muliplied by the local viscosity Node#externalViscosity and the
			* displacement of the node in current iteration k, from its position at the end of
			* previous time step n.
			*/
			if (externalViscosityVolumeBased){
				double massTimesViscosity = Nodes[i]->mass*Nodes[i]->externalViscosity[j];
				double displacementValue = gsl_matrix_get(displacementPerDt,3*i+j,0);
				gsl_matrix_set(gvExternal,3*i+j,0,massTimesViscosity*displacementValue);
			}
			else{
				double surfaceAreaTimesViscosity = Nodes[i]->viscositySurface*Nodes[i]->externalViscosity[j];
				double displacementValue = gsl_matrix_get(displacementPerDt,3*i+j,0);
				gsl_matrix_set(gvExternal,3*i+j,0,surfaceAreaTimesViscosity*displacementValue);
				if (std::isnan(surfaceAreaTimesViscosity)){
					std::cout<<" node: "<<i<<" dimention: "<<j<<" surfaceAreaTimesViscosity is nan: "<<surfaceAreaTimesViscosity<<std::endl;
				}
				if (std::isnan(displacementValue)){
					std::cout<<" node: "<<i<<" dimention: "<<j<<" displacementValue is nan: "<<displacementValue<<std::endl;
				}
			}
		}
	}
	/** Once all the viscous forces are calculated, the force direction will be inverted, as I am interested in the visouc
	* drag applied by the media to the node.
	*/
 	gsl_matrix_scale(gvExternal,-1.0);
}

void NewtonRaphsonSolver::addImplicitKViscousExternalToJacobian(const std::vector<std::unique_ptr<Node>>&  Nodes, double dt){
	/** This function will add the derivatives of external viscous drag forces with respect to
	* nodal displacement onto the system Jacobian.
	*
	*/
	for (size_t i = 0; i<nNodes; ++i ){
		if (externalViscosityVolumeBased){
			double massPerDt = Nodes[i]->mass/dt;
			for (size_t j=0; j<nDim; ++j){
				double curKValue = gsl_matrix_get(K,3*i+j,3*i+j);
				double massTimesViscosityPerDt = massPerDt*Nodes[i]->externalViscosity[j];
				//double matrixValuePerDt = gsl_matrix_get(mvisc,0,3*i+j)/dt;
				gsl_matrix_set(K,3*i+j,3*i+j,massTimesViscosityPerDt+curKValue);
			}
		}
		else{
			double surfaceAreaPerDt = Nodes[i]->viscositySurface/dt;
			for (size_t j=0; j<nDim; ++j){
				double curKValue = gsl_matrix_get(K,3*i+j,3*i+j);
				double surfaceAreaTimesViscosityPerDt = surfaceAreaPerDt*Nodes[i]->externalViscosity[j];
				gsl_matrix_set(K,3*i+j,3*i+j,surfaceAreaTimesViscosityPerDt+curKValue);
			}
		}
	}
}

void NewtonRaphsonSolver::checkJacobianForAblatedNodes(std::vector <int> & AblatedNodes){
	/** If there are ablated nodes, as recorded in the input vector AblatedNodes,
	* these nodes should have identity Jacobians. The Jacobian is cleared accordingly. 
	*/
	int nAblatedNode = AblatedNodes.size();
	for (int a = 0; a<nAblatedNode; ++a){
		int NodeId = AblatedNodes[a]*3;
		for (int aa= 0; aa<3; ++aa){
			double Kdiagonal = gsl_matrix_get(K,NodeId+aa,NodeId+aa);
			if (Kdiagonal == 0){
				gsl_matrix_set(K,NodeId+aa,NodeId+aa,1);
			}
		}
	}
}

void NewtonRaphsonSolver::calculateSumOfInternalForces(){
	for (size_t i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_matrix_get(ge,i,0)+gsl_matrix_get(gvInternal,i,0)+gsl_matrix_get(gvExternal,i,0));
	}
}

void NewtonRaphsonSolver::addExernalForces(){
	/** This function adds the external forces on the system sum of forces, NewtonRaphsonSolver#gSum)
	*/
	for (size_t i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
		double value = gsl_vector_get(gSum,i);
		if (std::isnan(value)){
		      std::cout<<" gSUM is nan at matrix point: "<<i<<std::endl;
		}
	}
}



void NewtonRaphsonSolver::solveForDeltaU(){
	/** This function will solve for new displacements from the system forces and the Jacobian.
	*This requires solving a sparse system of linear equations, and the operation is handled by Pardiso solver.
	*Please refer to the manual of Pardiso to follow the necessary steps in this function/
	*/
	const int nmult  = nDim*nNodes;
	int *ia = new int[nmult+1];
	double *b = new double[nmult];
	vector <int> ja_vec;
	vector <double> a_vec;
	constructiaForPardiso(ia, nmult, ja_vec, a_vec);
	const int nNonzero = ja_vec.size();
	int* ja = new int[nNonzero];
	double* a = new double [nNonzero];
	writeKinPardisoFormat(nNonzero, ja_vec, a_vec, ja, a);
	writeginPardisoFormat(b,nmult);
//	int error = solveWithPardiso(a, b, ia, ja, nmult);
//	if (error != 0){std::cerr<<"Pardiso solver did not return success!!"<<std::endl;}
//	if (boundNodesWithSlaveMasterDefinition){
//		equateSlaveDisplacementsToMasters();
//	}
	delete[] ia;
	delete[] ja;
	delete[] a;
	delete[] b;
}

// PARDISO prototype. //
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, double *, int    *,    int *, int *,   int *, int *,   int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


//int NewtonRaphsonSolver::solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables){

//    // I am copying my libraries to a different location for this to work:
//    // On MAC:
//    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgfortran.3.dylib /usr/local/lib/
//    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgomp.1.dylib /usr/local/lib/
//    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libquadmath.0.dylib /usr/local/lib/
//    // cp libpardiso500-MACOS-X86-64.dylib usr/local/lib
//    //
//    // compilation:
//    // g++ pardiso_sym.cpp -o pardiso_sym  -L./ -L/usr/local/lib -L/usr/lib/  -lpardiso500-MACOS-X86-64 -llapack


//    // On ubuntu,
//    // cp libpardiso500-GNU461-X86-64.so /usr/lib/
//    //
//    // sometimes linux cannot recognise liblapack.so.3gf or liblapack.so.3.0.1 or others like this, are essentially liblapack.so
//    // on ubuntu you can get this solved by installing liblapack-dev:
//    // sudo apt-get install liblapack-dev
//    //
//    // compilation:
//    // gcc test.cpp -o testexe  -L/usr/lib/  -lpardiso500-GNU461-X86-64  -fopenmp  -llapack

//    //
//    // also for each terminal run:
//    // export OMP_NUM_THREADS=1
//    // For mkl this is :
//    // export MKL_PARDISO_OOC_MAX_CORE_SIZE=10000
//    // export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
//    // fo
//    // MSGLVL: the level of verbal output, 0 is no output.

//    int    n = n_variables;
//    int    nnz = ia[n];
//    int    mtype = 11;        /* Real unsymmetric matrix */

//    /* RHS and solution vectors. */
//    int      nrhs = 1;          /* Number of right hand sides. */
//    double   x[n_variables];//, diag[n_variables];
//    /* Internal solver memory pointer pt,                  */
//    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
//    /* or void *pt[64] should be OK on both architectures  */
//    void    *pt[64];

//    /* Pardiso control parameters. */
//    int      iparm[64];
//    double   dparm[64];
//    int      maxfct, mnum, phase, error, msglvl, solver;

//    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

//    /* Number of processors. */
//    int      num_procs;

//    /* Auxiliary variables. */
//    char    *var;
//    int      i;// k;

//    double   ddum;              /* Double dummy */
//    int      idum;              /* Integer dummy. */


///* -------------------------------------------------------------------- */
///* ..  Setup Pardiso control parameters.                                */
///* -------------------------------------------------------------------- */

//    error = 0;
//    solver = 0; /* use sparse direct solver */
//    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

//    if (error != 0)
//    {
//        if (error == -10 )
//           printf("No license file found \n");
//        if (error == -11 )
//           printf("License is expired \n");
//        if (error == -12 )
//           printf("Wrong username or hostname \n");
//         return 1;
//    }
//    else
//        //printf("[PARDISO]: License check was successful ... \n");

//    /* Numbers of processors, value of OMP_NUM_THREADS */
//    var = getenv("OMP_NUM_THREADS");
//    if(var != NULL)
//        sscanf( var, "%d", &num_procs );
//    else {
//        printf("Set environment OMP_NUM_THREADS to 1");
//        exit(1);
//    }
//    iparm[2]  = num_procs;

//    maxfct = 1;		    /* Maximum number of numerical factorizations.  */
//    mnum   = 1;         /* Which factorization to use. */

//    iparm[10] = 0; /* no scaling  */
//    iparm[12] = 0; /* no matching */

//    msglvl = 0;         /* Print statistical information  */
//    error  = 0;         /* Initialize error flag */

///* -------------------------------------------------------------------- */
///* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
///*     notation.                                                        */
///* -------------------------------------------------------------------- */
//    for (i = 0; i < n+1; i++) {
//        ia[i] += 1;
//    }
//    for (i = 0; i < nnz; i++) {
//        ja[i] += 1;
//    }

///* -------------------------------------------------------------------- */
///*  .. pardiso_chk_matrix(...)                                          */
///*     Checks the consistency of the given matrix.                      */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//    bool carryOutDebuggingChecks = false;
//    if (carryOutDebuggingChecks){
//        pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//        if (error != 0) {
//            printf("\nERROR in consistency of matrix: %d", error);
//            exit(1);
//        }
//    }
///* -------------------------------------------------------------------- */
///* ..  pardiso_chkvec(...)                                              */
///*     Checks the given vectors for infinite and NaN values             */
///*     Input parameters (see PARDISO user manual for a description):    */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */

//    if (carryOutDebuggingChecks){
//        pardiso_chkvec (&n, &nrhs, b, &error);
//        if (error != 0) {
//            printf("\nERROR  in right hand side: %d", error);
//            exit(1);
//        }
//    }
///* -------------------------------------------------------------------- */
///* .. pardiso_printstats(...)                                           */
///*    prints information on the matrix to STDOUT.                       */
///*    Use this functionality only for debugging purposes                */
///* -------------------------------------------------------------------- */
//    if (carryOutDebuggingChecks){
//        pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
//        if (error != 0) {
//            printf("\nERROR right hand side: %d", error);
//            exit(1);
//        }
//    }
///* -------------------------------------------------------------------- */
///* ..  Reordering and Symbolic Factorization.  This step also allocates */
///*     all memory that is necessary for the factorization.              */
///* -------------------------------------------------------------------- */
//    phase = 11;
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error, dparm);
////std::cout<<"symbolic factorisation"<<std::endl;
//    if (error != 0) {
//        printf("\nERROR during symbolic factorization: %d", error);
//        exit(1);
//    }
//    //printf("\nReordering completed ... ");
//    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

///* -------------------------------------------------------------------- */
///* ..  Numerical factorization.                                         */
///* -------------------------------------------------------------------- */
//    phase = 22;
////    iparm[32] = 1; /* compute determinant */

//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
////std::cout<<"numerical factorisation"<<std::endl;
//    if (error != 0) {
//        printf("\nERROR during numerical factorization: %d", error);
//        exit(2);
//    }
//    //printf("\nFactorization completed ...\n ");

///* -------------------------------------------------------------------- */
///* ..  Back substitution and iterative refinement.                      */
///* -------------------------------------------------------------------- */
//   /* phase = 33;

//    iparm[7] = 1;       // Max numbers of iterative refinement steps.

//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, b, x, &error,  dparm);

//    if (error != 0) {
//        printf("\nERROR during solution: %d", error);
//        exit(3);
//    }
//    bool displayResult = false;
//    if (displayResult){
//        printf("\nSolve completed ... ");
//        printf("\nThe solution of the system is: ");
//        for (i = 0; i < n; i++) {
//            printf("\n x [%d] = % f", i, x[i] );
//        }
//        printf ("\n");
//    }
//    //Write x into deltaU:
//    for (int i=0; i<n_variables; ++i){
//        gsl_vector_set(deltaU,i,x[i]);
//    }
//    */
///* -------------------------------------------------------------------- */
///* ..  Back substitution with tranposed matrix A^t x=b                  */
///* -------------------------------------------------------------------- */

//	phase = 33;
//	//iparm[4]  = 61;	 /*changing the precision of convergence with pre-conditioning, not sure what it does, I added as trial, but did not change anything */
//	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
//	iparm[11] = 1;       /* Solving with transpose matrix. */

//	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//			 &n, a, ia, ja, &idum, &nrhs,
//			 iparm, &msglvl, b, x, &error,  dparm);

//	if (error != 0) {
//		printf("\nERROR during solution: %d", error);
//		exit(3);
//	}

//	bool displayResult = false;
//	if (displayResult){
//		printf("\nSolve completed ... ");
//		printf("\nThe solution of the system is: ");
//		for (i = 0; i < n; i++) {
//			printf("\n x [%d] = % f", i, x[i] );
//		}
//		printf ("\n");
//	}
//    //Write x into deltaU:
//    for (int i=0; i<n_variables; ++i){
//        gsl_vector_set(deltaU,i,x[i]);
//    }

///* -------------------------------------------------------------------- */
///* ..  Convert matrix back to 0-based C-notation.                       */
///* -------------------------------------------------------------------- */
//    for (i = 0; i < n+1; i++) {
//        ia[i] -= 1;
//    }
//    for (i = 0; i < nnz; i++) {
//        ja[i] -= 1;
//    }

///* -------------------------------------------------------------------- */
///* ..  Termination and release of memory.                               */
///* -------------------------------------------------------------------- */
//    phase = -1;                 /* Release internal memory. */

//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, &ddum, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//    return 0;
//}


void NewtonRaphsonSolver::constructiaForPardiso(int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
    double negThreshold = -1E-13, posThreshold = 1E-13;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=0; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //std::cout<<"writing is for row "<<i<<" column is: "<<j<<std::endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}




void NewtonRaphsonSolver::writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}

void NewtonRaphsonSolver::writeginPardisoFormat(double* b, const int n){
    for (int i=0; i<n; ++i){
        b[i] = gsl_vector_get(gSum,i);
    }
}


bool NewtonRaphsonSolver::checkConvergenceViaDeltaU(){
	/**
	* This function checks the norm of NewtonRaphsonSolver#deltaU, the change in
	*nodal positions from previous iteration to this one. If the norm is below the threshold of convergence
	*NewtonRaphsonSolver#threshold, then the solution has been achieved for the current time step.
	*/
	bool converged = true;

	double d = gsl_blas_dnrm2 (deltaU);

	if (d>threshold){
		converged = false;
		std::cout<<" not  yet converged via du: norm "<<d<<std::endl;
	}
	else{
		std::cout<<"converged with displacement: norm"<<d<<std::endl;
	}
	return converged;
}

bool NewtonRaphsonSolver::checkConvergenceViaForce(){
	/** The system can converge with zero forces as well. This function
	* will check the norm of the sum of all nodal forces against the threshold. This check
	* is not currently used.
	*/
	bool converged = true;
	double d = gsl_blas_dnrm2 (gSum);
	if (d>threshold){
		converged = false;
		std::cout<<" not  yet converged via forces: norm "<<d<<std::endl;
	}
	else{
		std::cout<<"converged with forces: norm"<<d<<std::endl;
	}
	return converged;
}

void NewtonRaphsonSolver::updateUkInIteration(){
	/**
	* This function updates the positions of nodes for next iteration $ \f \boldsymbol{u}_{k+1} $ \f
	*from the positions of current iteration the $ \f \boldsymbol{u}_{k} $ \f and $ \f \boldsymbol{\delta u}_{k} $ \f.
	*/
	int n = uk->size1;
	for (int i=0; i<n;++i){
		double newValue = gsl_matrix_get(uk,i,0)+gsl_vector_get(deltaU,i);
		gsl_matrix_set(uk,i,0,newValue);
	}
}

void NewtonRaphsonSolver::calculateDifferenceBetweenNumericalAndAnalyticalJacobian(const std::vector<std::unique_ptr<Node>>&  Nodes, bool displayMatricesDuringNumericalCalculation){
	/**
	 * This function calculates the difference between the numerical and analytical Jacobians.
	 * It is here for debugging purposes only, it should be used with caution as it displays the difference matrix.
	 */
	//The normal K still includes the values for fixed nodes. Should correct the fixed nodes,
	//then calculate the difference
	calcutateFixedK(Nodes);
	//calculate the difference:
	gsl_matrix* Kdiff = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	double d = 0;
	for (size_t i=0; i<nDim*nNodes; ++i){
		for (size_t j=0; j<nDim*nNodes; ++j){
			double value = gsl_matrix_get(Knumerical,i,j) - gsl_matrix_get(K,i,j);
			gsl_matrix_set(Kdiff,i,j,value);
			d += value*value;
		}
	}
	d = pow(d,0.5);
	if(displayMatricesDuringNumericalCalculation){
		displayMatrix(K,"normalK");
		displayMatrix(Kdiff,"differenceKMatrix");
	}
	std::cout<<"norm of difference between numerical and analytical K: "<<d<<std::endl;
	gsl_matrix_free(Kdiff);
}

void NewtonRaphsonSolver::useNumericalJacobianInIteration(){
	gsl_matrix_memcpy(K,Knumerical);
}

void NewtonRaphsonSolver::displayMatrix(gsl_matrix* mat, string matname){
    int m = mat->size1;
    int n = mat->size2;
    std::cout<<matname<<": "<<std::endl;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            cout.precision(4);
            cout.width(6);
            std::cout<<gsl_matrix_get(mat,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void NewtonRaphsonSolver::displayMatrix(gsl_vector* mat, string matname){
    int m = mat->size;
    std::cout<<matname<<": "<<std::endl;
    for (int i =0; i<m; i++){
		std::cout.precision(4);
		std::cout.width(6);
		std::cout<<gsl_vector_get(mat,i)<<std::endl;
    }
}

bool NewtonRaphsonSolver::checkIfCombinationExists(int dofSlave, int dofMaster){
	/**
	 * The input gives a potential new slave degrees of freedom, and its potential new master.
	 * This function checks if the pair is already recorded.
	 */
	int n= slaveMasterList.size();
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][0] == dofSlave && slaveMasterList[i][1] == dofMaster){
			return false; //continue addition? false
		}
		if(slaveMasterList[i][1] == dofSlave && slaveMasterList[i][0] == dofMaster){
			return false; //continue addition? false
		}
	}
	return true;
}

void NewtonRaphsonSolver::checkMasterUpdate(int& dofMaster, int& masterId){
	/**
	 * The input gives a potential new master degrees of freedom, and its node Id.
	 * This function will check if the master is already slave to other degrees of freedom.
	 * If the master is indeed a slave (recorded in the NewtonRaphsonSolver#slaveMasterList
	 * 2D array's first dimension), then the master degree of freedom, and the master nod eid are
	 * updated to be the original potential master's existing master.
	 */
	//order in NewtonRaphsonSolver#slaveMasterList [slave][master]
	int n= slaveMasterList.size();
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][0] == dofMaster){
			dofMaster = slaveMasterList[i][1];
			int dim = dofMaster % 3;
			masterId = (dofMaster-dim)/3;
			break;
		}
	}
}

void NewtonRaphsonSolver::cleanPeripodialBindingFromMaster(int masterDoF, const std::vector<std::unique_ptr<Node>>&  Nodes){
	//order is [slave][master]
	//std::cout<<" searching for master: "<<masterDoF<<std::endl;
	for (vector< vector<int> >::iterator iter = slaveMasterList.begin(); iter != slaveMasterList.end(); ) {
	    //std::cout<<"checking slaveMasterList of size "<<slaveMasterList.size()<<" slave/master "<< (*iter)[0]<<" / "<<(*iter)[1]<<std::endl;
		bool deletedItem = false;
		if ((*iter)[1] == masterDoF){
	    	//I have found a couple where this is a master, is the slave a peripodiaal node?
	    	int dofSlave = (*iter)[0];
	    	int dim = dofSlave % 3;
	    	int slaveId = (dofSlave-dim)/3;
	    	if (Nodes[slaveId]->tissueType == 1){
	    		//clearing the peripodial slave
	    		//std::cout<<"deleting binding from peripodial node "<<Nodes[slaveId]<<" dim: "<<dim<<std::endl;
	    		Nodes[slaveId]->slaveTo[dim] = -1;
	    		iter = slaveMasterList.erase(iter);
	    		deletedItem = true;
	    	}
	    }
	    if(!deletedItem){
	        ++iter;
	    }
	}

}

bool NewtonRaphsonSolver::checkIfSlaveIsAlreadyMasterOfOthers(int dofSlave, int dofMaster){
	/**
	 * The input gives a potential new slave degrees of freedom, and its potential new master.
	 * This function will check if the slave is already master of other degrees of freedom.
	 * If the slave is indeed a master (recorded in the NewtonRaphsonSolver#slaveMasterList
	 * 2D array's second dimension), then the master of those nodes are assigned to be the
	 * new master. The function informs the caller that it has made changes.
	 */
	int n= slaveMasterList.size();
	bool madeChange = false;
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][1] == dofSlave){
			std::cout<<"making change, slave was master of "<<slaveMasterList[i][0]<<std::endl;
			//proposed slave is already a master, update the master to the proposed master
			slaveMasterList[i][1] = dofMaster;
			madeChange = true;
		}
	}
	return madeChange;
}
