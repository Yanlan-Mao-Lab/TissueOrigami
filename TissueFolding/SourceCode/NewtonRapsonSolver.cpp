/*
 * NewtonRapsonSolver.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: melda
 */

#include "NewtonRapsonSolver.h"
//#include "Node.h"
//#include <gsl/gsl_linalg.h>

using namespace std;

NewtonRapsonSolver::NewtonRapsonSolver(int dim, int n){
    threshold = 1E-10;

	nDim = dim;
	nNodes = n;

	numericalParametersSet = false;
	un = gsl_matrix_calloc(dim*nNodes,1);
	mvisc = gsl_matrix_calloc(dim*nNodes,dim*nNodes);
	mviscPerDt = gsl_matrix_calloc(dim*nNodes,dim*nNodes);
	ge = gsl_matrix_calloc(dim*nNodes,1);
	gvInternal = gsl_matrix_calloc(dim*nNodes,1);
	gvExternal = gsl_matrix_calloc(dim*nNodes,1);
	gExt = gsl_matrix_calloc(dim*nNodes,1);
	gSum = gsl_vector_calloc(dim*nNodes);
	uk = gsl_matrix_calloc(dim*nNodes,1);
	displacementPerDt = gsl_matrix_calloc(dim*nNodes,1);
	deltaU = gsl_vector_calloc(dim*nNodes);
	K = gsl_matrix_calloc(dim*nNodes,dim*nNodes);

}

NewtonRapsonSolver::~NewtonRapsonSolver(){
    gsl_matrix_free(un);
    gsl_matrix_free(mvisc);
    gsl_matrix_free(mviscPerDt);
    gsl_matrix_free(ge);
    gsl_matrix_free(gvInternal);
    gsl_matrix_free(gvExternal);
    gsl_matrix_free(gExt);
    gsl_vector_free(gSum);
    gsl_matrix_free(uk);
    gsl_matrix_free(displacementPerDt);
    gsl_vector_free(deltaU);
    gsl_matrix_free(K);
    gsl_matrix_free(Knumerical);
}

void NewtonRapsonSolver::setMatricesToZeroAtTheBeginningOfIteration(bool thereIsNumericalCalculation){
	gsl_matrix_set_zero(un);
	gsl_matrix_set_zero(mvisc);
	gsl_matrix_set_zero(mviscPerDt);
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

void NewtonRapsonSolver::setMatricesToZeroInsideIteration(){
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(K);
}


void NewtonRapsonSolver::constructUnMatrix(vector <Node*>& Nodes){
    for (int i = 0; i<nNodes; ++i ){
        for (int j=0; j<nDim; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void NewtonRapsonSolver::initialteUkMatrix(){
    gsl_matrix_memcpy(uk,un);
}

void NewtonRapsonSolver::constructLumpedMassExternalViscosityDtMatrix(vector <Node*>& Nodes, double dt){
    for (int i = 0; i<nNodes; ++i ){
        //double matrixValue = Nodes[i]->mass*Nodes[i]->Viscosity / dt;
        for (int j=0; j<nDim; ++j){
        	double matrixValue = Nodes[i]->mass*Nodes[i]->externalViscosity[j];
            gsl_matrix_set(mvisc,3*i+j,3*i+j,matrixValue);
        }
        //cout<<" Node "<<i<<" - mass: "<<Nodes[i]->mass<<" visc: "<<Nodes[i]->Viscosity<<" matrixvalues: "<<gsl_matrix_get(mviscdt,3*i,3*i)<<" "<<gsl_matrix_get(mviscdt,3*i+1,3*i+1)<<" "<<gsl_matrix_get(mviscdt,3*i+2,3*i+2)<<endl;
    }
    gsl_matrix_memcpy(mviscPerDt,mvisc);
    gsl_matrix_scale(mviscPerDt,1/dt);
}

void NewtonRapsonSolver::calculateDisplacementMatrix(double dt){
	gsl_matrix_memcpy(displacementPerDt,uk);
	gsl_matrix_sub(displacementPerDt,un);
	gsl_matrix_scale(displacementPerDt,1.0/dt);

}

void NewtonRapsonSolver::calcutateFixedK(vector <Node*>& Nodes){
    int dim = 3;
    int Ksize = K->size1;
    for(int i=0; i<nNodes; i++){
        for (int j=0; j<dim; ++j){
            //if ( (Xfixed && j == 0) || (Yfixed && j == 1) || (Zfixed && j == 2) ){
            if (Nodes[i]->FixedPos[j]){
                int index1 = i*dim+j;
                gsl_vector_set(gSum,index1,0.0); // making the forces zero
                for (int k =0; k<Ksize; ++k){
                    double value =0.0;
                    if (index1 == k ){value =1.0;}
                    gsl_matrix_set(K, index1, k, value);
                    gsl_matrix_set(K, k, index1, value); //K is symmetric;
                }
            }
        }
    }
}

void NewtonRapsonSolver::calculateForcesAndJacobianMatrixNR(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double dt, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile){
//    omp_set_num_threads(4);
//#pragma omp parallel for
//    for (int a =0; a<2; ++a){
//    	int initialpoint = 0;
//    	int breakpoint =  nElements/2;
//    	if (a == 1){
//    		initialpoint = breakpoint;
//    		breakpoint = nElements;
//    	}
		for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	//for(vector<ShapeBase*>::iterator  itElement=Elements.begin()+ initialpoint; itElement<Elements.begin()+breakpoint; ++itElement){
			//int nThreads = omp_get_thread_num();
			//cout<<" omp threads: "<<nThreads<<endl;
			if (!(*itElement)->IsAblated){
				(*itElement)->calculateForces(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile);
			}
			(*itElement)->calculateImplicitKElastic(); //This is the stiffness matrix, elastic part of the jacobian matrix
			(*itElement)->calculateImplicitKViscous(displacementPerDt, dt); //This is the viscous part of jacobian matrix
		}
//	}
}

void NewtonRapsonSolver::writeForcesTogeAndgvInternal(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double** SystemForces){
    for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->writeInternalForcesTogeAndgv(ge,gvInternal,SystemForces,Nodes);
        }
    }
}

void NewtonRapsonSolver::writeImplicitElementalKToJacobian(vector <ShapeBase*>& Elements){
    //writing all elements K values into big K matrix:
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	//if element is ablated, current elemental K matrix will be identity
    	//(*itElement)->writeKelasticToMainKatrix(K);
    	//activate this for internal viscosity
    	(*itElement)->writeKviscousToMainKatrix(K);
    }
    //Elements[0]->displayMatrix(K,"normalKonlyViscous");
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	//if element is ablated, current elemental K matrix will be identity
    	(*itElement)->writeKelasticToMainKatrix(K);
    	//activate this for internal viscosity
    }
}

void NewtonRapsonSolver::calculateExternalViscousForcesForNR(){
    //the mass is already updated for symmetricity boundary nodes, and the viscous forces will be calculated correctly,
	//as mvisdt is already accounting for the doubling of mass
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mvisc, displacementPerDt,0.0, gvExternal);
    //added this line with sign change
    gsl_matrix_scale(gvExternal,-1.0);
}

void NewtonRapsonSolver::addImplicitKViscousExternalToJacobian(){
    gsl_matrix_add(K,mviscPerDt);
}

void NewtonRapsonSolver::checkJacobianForAblatedNodes(vector <int> & AblatedNodes){
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
void NewtonRapsonSolver::calculateSumOfInternalForces(){
	for (int i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_matrix_get(ge,i,0)+gsl_matrix_get(gvInternal,i,0)+gsl_matrix_get(gvExternal,i,0));
	}
}

void NewtonRapsonSolver::addExernalForces(){
	//cout<<"after random forces addition"<<endl;
	for (int i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
	}
}


void NewtonRapsonSolver::solveForDeltaU(){
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
    int error = solveWithPardiso(a, b, ia, ja, nmult);
    if (error != 0){cerr<<"Pardiso solver did not return success!!"<<endl;}
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
}

#include <math.h>
/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, double *, int    *,    int *, int *,   int *, int *,   int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


int NewtonRapsonSolver::solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables){

    // I am copying my libraries to a different location for this to work:
    // On MAC:
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgfortran.3.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgomp.1.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libquadmath.0.dylib /usr/local/lib/
    // cp libpardiso500-MACOS-X86-64.dylib usr/local/lib
    //
    // compilation:
    // g++ pardiso_sym.cpp -o pardiso_sym  -L./ -L/usr/local/lib -L/usr/lib/  -lpardiso500-MACOS-X86-64 -llapack


    // On ubuntu,
    // cp libpardiso500-GNU461-X86-64.so /usr/lib/
    //
    // sometimes linux cannot recognise liblapack.so.3gf or liblapack.so.3.0.1 or others like this, are essentially liblapack.so
    // on ubuntu you can get this solved by installing liblapack-dev:
    // sudo apt-get install liblapack-dev
    //
    // compilation:
    // gcc test.cpp -o testexe  -L/usr/lib/  -lpardiso500-GNU461-X86-64  -fopenmp  -llapack

    //
    // also for each terminal run:
    // export OMP_NUM_THREADS=1
    // For mkl this is :
    // export MKL_PARDISO_OOC_MAX_CORE_SIZE=10000
    // export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
    //
    // MSGLVL: the level of verbal output, 0 is no output.

    int    n = n_variables;
    int    nnz = ia[n];
    int    mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    double   x[n_variables], diag[n_variables];
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i, k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */

    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        //printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		    /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    iparm[10] = 0; /* no scaling  */
    iparm[12] = 0; /* no matching */

    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    bool carryOutDebuggingChecks = false;
    if (carryOutDebuggingChecks){
        pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
        if (error != 0) {
            printf("\nERROR in consistency of matrix: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    if (carryOutDebuggingChecks){
        pardiso_chkvec (&n, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR  in right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */
    if (carryOutDebuggingChecks){
        pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11;
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//cout<<"symbolic factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */
    phase = 22;
//    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//cout<<"numerical factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    //printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */
   /* phase = 33;

    iparm[7] = 1;       // Max numbers of iterative refinement steps.

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    bool displayResult = false;
    if (displayResult){
        printf("\nSolve completed ... ");
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n");
    }
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }
    */
/* -------------------------------------------------------------------- */
/* ..  Back substitution with tranposed matrix A^t x=b                  */
/* -------------------------------------------------------------------- */

	phase = 33;
	//iparm[4]  = 61;	 /*changing the precision of convergence with pre-conditioning, not sure what it does, I added as trial, but did not change anything */
	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
	iparm[11] = 1;       /* Solving with transpose matrix. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			 &n, a, ia, ja, &idum, &nrhs,
			 iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	bool displayResult = false;
	if (displayResult){
		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		printf ("\n");
	}
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }

/* -------------------------------------------------------------------- */
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    return 0;
}

void NewtonRapsonSolver::constructiaForPardiso(int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
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
                    //cout<<"writing is for row "<<i<<" column is: "<<j<<endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}
void NewtonRapsonSolver::writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}

void NewtonRapsonSolver::writeginPardisoFormat(double* b, const int n){
    for (int i=0; i<n; ++i){
        b[i] = gsl_vector_get(gSum,i);
    }
}

bool NewtonRapsonSolver::checkConvergenceViaDeltaU(){
    bool converged = true;

    double d = gsl_blas_dnrm2 (deltaU);

    if (d>threshold){
        converged = false;
        cout<<" not  yet converged via du: norm "<<d<<endl;
    }
    else{
        cout<<"converged with displacement: norm"<<d<<endl;
    }
    return converged;
}

bool NewtonRapsonSolver::checkConvergenceViaForce(){
    bool converged = true;
    double d = gsl_blas_dnrm2 (gSum);
    if (d>threshold){
        converged = false;
        cout<<" not  yet converged via forces: norm "<<d<<endl;
    }
    else{
        cout<<"converged with forces: norm"<<d<<endl;
    }
    return converged;
}

void NewtonRapsonSolver::updateUkInIteration(){
    int n = uk->size1;
    for (int i=0; i<n;++i){
        double newValue = gsl_matrix_get(uk,i,0)+gsl_vector_get(deltaU,i);
        gsl_matrix_set(uk,i,0,newValue);
    }
}
