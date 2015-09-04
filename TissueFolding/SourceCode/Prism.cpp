#include "Prism.h"
#include "ReferenceShapeBase.h"
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

Prism::Prism(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId){
	//cout<<"constructing prism"<<endl;
	nNodes = 6;
	nDim = 3;
	Id = CurrId;
	ShapeDim = 3;	//3D shape
	NodeIds = new int[6];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
    lambda = E*v /(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);

    D = gsl_matrix_calloc(6,6);
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[6];
	CurrGrowthStrainAddition = new double[6];
	ApicalNormalForPacking =  new double[3];
	BasalNormalForPacking =  new double[3];
	RelativePosInBoundingBox = new double[3];
	for (int i=0; i<3; ++i){
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
		ApicalNormalForPacking[i] = 0;
		BasalNormalForPacking[i] = 0;
		RelativePosInBoundingBox[i] =0;
	}
	columnarGrowthWeight = 1.0;
	peripodialGrowthWeight = 0.0;
	for (int i=0; i<6; ++i){
		CurrGrowthStrainAddition[i] = 0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	CurrGrowthStrainsUpToDate = false;
	IsGrowing = false;
	IsChangingShape = false;
	GrewInThePast = false;
	ChangedShapeInThePast = false;
	ApicalNormalForPackingUpToDate = false;
	BasalNormalForPackingUpToDate = false;
	IsAblated = false;
	IsClippedInDisplay = false;
	capElement = false;
    rotatedGrowth = false;

	setIdentificationColour();
	setShapeType("Prism");
	ReferenceShape = new ReferenceShapeBase("Prism");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);

    ShapeFuncDerivatives = new gsl_matrix*[3];
    ShapeFuncDerStacks = new gsl_matrix*[3];
    InvdXdes = new gsl_matrix*[3];
    detdXdes = new double[3];
    Bmatrices = new gsl_matrix*[3];
    FeMatrices = new gsl_matrix*[3];
    CMatrices = new gsl_matrix*[3];
    detFs = new double[3];
    invJShapeFuncDerStack = new gsl_matrix*[3];
    invJShapeFuncDerStackwithFe  = new gsl_matrix*[3];
    elasticStress = new gsl_matrix*[3];
    for (int i=0; i<3; ++i){
        ShapeFuncDerivatives[i] = gsl_matrix_calloc(nDim, nNodes);
        ShapeFuncDerStacks[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        InvdXdes[i] = gsl_matrix_calloc(nDim, nDim);
        detdXdes[i] = 0.0;
        Bmatrices[i] = gsl_matrix_calloc(nNodes,nDim*nNodes);
        FeMatrices[i] = gsl_matrix_calloc(3,3);
        CMatrices[i] = gsl_matrix_calloc(6,6);
        invJShapeFuncDerStack[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        invJShapeFuncDerStackwithFe[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        elasticStress[i] = gsl_matrix_calloc(3,3);
    }
    Strain = gsl_matrix_calloc(6,1);
    RK1Strain = gsl_matrix_calloc(6,1);
    GrowthStrainsRotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(GrowthStrainsRotMat);
    Fg = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fg);
    InvFg = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFg);
    Fsc = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fsc);
    InvFsc = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFsc);
    TriPointF = gsl_matrix_calloc(3,3);
    TriPointKe = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	RotatedElement = false;    

	CurrShapeChangeToAdd[0] = 0;
	CurrShapeChangeToAdd[1] = 0;
	CurrShapeChangeToAdd[2] = 0;

	VolumePerNode = 0;
    ZProjectedBasalArea=0.0;
    ZProjectedApicalArea=0.0;
}

Prism::~Prism(){
	//cout<<"called the destructor for prism class"<<endl;
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
	}
    delete[] Positions;
	delete[] RelativePosInBoundingBox;
    delete[] NodeIds;
	delete[] IdentifierColour;
	delete[] GrowthRate;
    delete[] ShapeChangeRate;
	delete	 ReferenceShape;

    //freeing matrices allocated
    gsl_matrix_free(D);
    gsl_matrix_free(CoeffMat);
    gsl_matrix_free(Fg);
    gsl_matrix_free(InvFg);
    gsl_matrix_free(Fsc);
    gsl_matrix_free(InvFsc);
    gsl_matrix_free(TriPointF);
    gsl_matrix_free(Strain);
    gsl_matrix_free(RK1Strain);
    gsl_matrix_free(TriPointKe);
    gsl_matrix_free(GrowthStrainsRotMat);
    for (int i=0; i<3; ++i){
        gsl_matrix_free (ShapeFuncDerivatives[i]);
        gsl_matrix_free (ShapeFuncDerStacks[i]);
        gsl_matrix_free (InvdXdes[i]);
        gsl_matrix_free (Bmatrices[i]);
        gsl_matrix_free (FeMatrices[i]);
        gsl_matrix_free (CMatrices[i]);
        gsl_matrix_free (invJShapeFuncDerStack[i]);
        gsl_matrix_free (invJShapeFuncDerStackwithFe[i]);
        gsl_matrix_free (elasticStress[i]);
    }
    delete[] ShapeFuncDerivatives;
    delete[] ShapeFuncDerStacks;
    delete[] InvdXdes;
    delete[] detdXdes;
    delete[] Bmatrices;
    delete[] FeMatrices;
    delete[] CMatrices;
    delete[] detFs;
    delete[] invJShapeFuncDerStack;
    delete[] invJShapeFuncDerStackwithFe;
    delete[] elasticStress;
}

void Prism::setCoeffMat(){
    CoeffMat = gsl_matrix_calloc(6, nDim*nDim);
    gsl_matrix_set(CoeffMat,0,0,1);
    gsl_matrix_set(CoeffMat,1,4,1);
    gsl_matrix_set(CoeffMat,2,8,1);
    gsl_matrix_set(CoeffMat,3,1,1);
    gsl_matrix_set(CoeffMat,3,3,1);
    gsl_matrix_set(CoeffMat,4,5,1);
    gsl_matrix_set(CoeffMat,4,7,1);
    gsl_matrix_set(CoeffMat,5,2,1);
    gsl_matrix_set(CoeffMat,5,6,1);
}

void Prism::checkRotationConsistency3D(){
	//The nodes should be ordered in counter-clock-wise order
	//The view-vector is from apical towards basal
	//(If the system is not rotated, and is in standard coordinate system,
	//I am looking from top, view is (-)ve z;
	//If they are not, correct the order here

	double *vec1   = new double[3];
	double *vec2   = new double[3];
	double *view   = new double[3];
	double *normal = new double[3];
	for (int i= 0; i<nDim; ++i){
		vec1[i] = Positions[1][i] - Positions[0][i];
		vec2[i] = Positions[2][i] - Positions[0][i];
		view[i] = Positions[0][i] - Positions[3][i];
	}
	crossProduct3D(vec1,vec2,normal);
	double  dot = dotProduct3D(view,normal);
	if (dot > 0) {
		cerr<<"prism: "<<Id<<" nodes are ordered clockwise, correcting"<<endl;
		delete[] vec1;
		delete[] vec2;
		delete[] view;
		delete[] normal;
		cout<<"Positions before swap, element: "<<Id<<endl;
		displayPositions();
		//swapping node ids:
		int ids[2] = { NodeIds[1], NodeIds[4]};
		NodeIds[1] = NodeIds[2];
		NodeIds[4] = NodeIds[5];
		NodeIds[2] = ids[0];
		NodeIds[5] = ids[1];
		//swapping positions:
		for (int i = 0; i<nDim; ++i){
			double pos[2] = {Positions[1][i],Positions[4][i]};
			double refpos[2] = {ReferenceShape->Positions[1][i],ReferenceShape->Positions[4][i]};
			Positions[1][i] = Positions[2][i];
			Positions[4][i] = Positions[5][i];
			Positions[2][i] = pos[0];
			Positions[5][i] = pos[1];
			ReferenceShape->Positions[1][i] = ReferenceShape->Positions[2][i];
			ReferenceShape->Positions[4][i] = ReferenceShape->Positions[5][i];
			ReferenceShape->Positions[2][i] = refpos[0];
			ReferenceShape->Positions[5][i] = refpos[1];
		}
		cout<<"Positions after swap, element: "<<Id<<endl;
		displayPositions();
	}
}

void  Prism::calculateBasalNormal(double * normal){
	double * u;
	u = new double[3];
	double * v;
	v = new double[3];
	for (int i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[1][i] - ReferenceShape->Positions[0][i];
		v[i] = ReferenceShape->Positions[2][i] - ReferenceShape->Positions[0][i];
		normal[i] = 0.0;
	}
	cerr<<"		u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<" v: "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
	crossProduct3D(u,v,normal);
	cerr<<"		normal before normalisation: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	normaliseVector3D(normal);
	cerr<<"		normal after normalisation: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	for (int i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[3][i] - ReferenceShape->Positions[0][i];
	}
	cerr<<"		vector to apical: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
	double  dot = dotProduct3D(u,normal);
	cerr<<"		dot product: "<<dot<<endl;
	if (dot<0){
		for (int i=0; i<nDim; ++i){
			normal[i] *=(-1.0);
		}
	}
	cerr<<"		normal after direction correction: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	delete[] v;
	delete[] u;
}

void  Prism::AlignReferenceBaseNormalToZ(){
	//getting the normal of the reference element basal surface
	double * normal;
	normal = new double[3];
	cerr<<"Element: "<<Id<<endl;
	calculateBasalNormal(normal);
	//Now I have the normal pointing towards the apical surface, I need to align it with (+ve)z vector
	double* z = new double[3];
	z[0] = 0;
	z[1] = 0;
	z[2] = 1;
	double c, s;
	calculateRotationAngleSinCos(normal,z,c,s); //align normal onto z
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(normal,z,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
		//calculateRotationAxis(normal,z,rotAx);	//calculating the rotation axis that is perpendicular to both u and v
		constructRotationMatrix(c,s,rotAx,rotMat);
		rotateReferenceElementByRotationMatrix(rotMat);
		delete[] rotAx;
		delete[] rotMat;
	}
	delete[] normal;
	delete[] z;
}

void  Prism::setElasticProperties(double EApical, double EBasal, double EMid, double v){
	this -> E = EMid;
	if (tissuePlacement == 0 ){
		this -> E = EBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> E = EApical;
	}
	this -> v = v; //poisson ratio
	if (v>0.5){v = 0.5;}
	else if (v<0.0){v = 0.0;}

    lambda = E*v/(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);

	double multiplier = E/((1+v)*(1-2*v));
    gsl_matrix_set(D,0,0,  multiplier*(1-v));
    gsl_matrix_set(D,0,1,  multiplier*v);
    gsl_matrix_set(D,0,2,  multiplier*v);
    gsl_matrix_set(D,1,0,  multiplier*v);
    gsl_matrix_set(D,1,1,  multiplier*(1-v));
    gsl_matrix_set(D,1,2,  multiplier*v);
    gsl_matrix_set(D,2,0,  multiplier*v);
    gsl_matrix_set(D,2,1,  multiplier*v);
    gsl_matrix_set(D,2,2,  multiplier*(1-v));
    gsl_matrix_set(D,3,3,  multiplier*(1-2*v)/2);
    gsl_matrix_set(D,4,4,  multiplier*(1-2*v)/2);
    gsl_matrix_set(D,5,5,  multiplier*(1-2*v)/2);
	
    //calculating 4th order tensor D
    double Idouble[3][3] = {{1.0,0.0,0.0} , {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    for (int I = 0; I<nDim; ++I){
        for (int J = 0; J<nDim; ++J){
            for (int K = 0; K<nDim; ++K){
                for (int L = 0; L<nDim; ++L){
                    D81[I][J][K][L] = lambda*Idouble[K][L]*Idouble[I][J] + mu * ( Idouble[I][K]*Idouble[J][L] + Idouble[I][L]*Idouble[J][K] );
                }
            }
        }
    }
    //cout<<" Element: "<<Id<<" E : "<<E<<" v: "<<v<<" lambda: "<<lambda<< " mu: "<<mu<<endl;
}


void Prism::getCurrRelaxedShape(gsl_matrix* CurrRelaxedShape){
	for (int i =0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
            gsl_matrix_set(CurrRelaxedShape,i,j,ReferenceShape->Positions[i][j]);
		}
	}
}

void Prism::setShapeFunctionDerivatives(gsl_matrix* ShapeFuncDer, double eta, double zeta, double nu){
	double alpha  = (1 - zeta)/2;
	double beta = (1 + zeta)/2;
	double lambda = 1-eta-nu;

    gsl_matrix_set(ShapeFuncDer,0,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,0,1,  alpha);
    gsl_matrix_set(ShapeFuncDer,0,2,  0);
    gsl_matrix_set(ShapeFuncDer,0,3, -beta);
    gsl_matrix_set(ShapeFuncDer,0,4,  beta);
    gsl_matrix_set(ShapeFuncDer,0,5, 0);

    gsl_matrix_set(ShapeFuncDer,1,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,1,1,  0);
    gsl_matrix_set(ShapeFuncDer,1,2,  alpha);
    gsl_matrix_set(ShapeFuncDer,1,3, -beta);
    gsl_matrix_set(ShapeFuncDer,1,4,  0);
    gsl_matrix_set(ShapeFuncDer,1,5,  beta);

    gsl_matrix_set(ShapeFuncDer,2,0, -lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,1, -eta/2);
    gsl_matrix_set(ShapeFuncDer,2,2, -nu/2);
    gsl_matrix_set(ShapeFuncDer,2,3,  lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,4,  eta/2);
    gsl_matrix_set(ShapeFuncDer,2,5,  nu/2);
}

void Prism::setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack){
	int n = nNodes;
	int dim = nDim;
	for (int i=0; i<n;++i){
        for (int k=0; k<dim; ++k){
            gsl_matrix_set(ShapeFuncDerStack,k,i*dim, gsl_matrix_get(ShapeFuncDer,k,i));
        }
        //subrange(ShapeFuncDerStack, 0,dim,i*dim,i*dim+1) = subrange(ShapeFuncDer,0,dim,i,i+1);
	}
    for (int k =dim; k<2*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
    for (int k =2*dim; k<3*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
    //subrange(ShapeFuncDerStack, dim,2*dim,1,dim*n) = subrange(ShapeFuncDerStack, 0,dim,0,dim*n-1);
    //subrange(ShapeFuncDerStack, 2*dim,3*dim,1,dim*n) = subrange(ShapeFuncDerStack, dim,2*dim,0,dim*n-1);
}

void Prism::calculateElementShapeFunctionDerivatives(){
    //Shape Function Derivatives 3-point:
    double points[3][3]={{1.0/6.0,1.0/6.0,0.0},{2.0/3.0,1.0/6.0,0.0},{1.0/6.0,2.0/3.0,0.0}};
    //Setting up the current reference shape position matrix:
    gsl_matrix * CurrRelaxedShape = gsl_matrix_calloc(nNodes, nDim);
    gsl_matrix * dXde  = gsl_matrix_calloc(nDim, nDim);
    getCurrRelaxedShape(CurrRelaxedShape);
    for (int iter =0; iter<3;++iter){
        double eta  = points[iter][0];
        double nu   = points[iter][1];
        double zeta = points[iter][2];
        //calculating shape function derivatives:
        setShapeFunctionDerivatives(ShapeFuncDerivatives[iter],eta,zeta,nu);
        setShapeFunctionDerivativeStack(ShapeFuncDerivatives[iter],ShapeFuncDerStacks[iter]);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDerivatives[iter], CurrRelaxedShape, 0.0, dXde);
        gsl_matrix_transpose(dXde);
        detdXdes[iter] = determinant3by3Matrix(dXde);
        bool inverted = InvertMatrix(dXde, InvdXdes[iter]);
        if (!inverted){
            cerr<<"dXde not inverted at point "<<iter<<"!!"<<endl;
        }
    }
}

void Prism::calculateCurrNodalForces(gsl_matrix *currg, gsl_matrix *currF, int pointNo){
    const int n = nNodes;
    const int dim = nDim;
    gsl_matrix *currFscFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix *currFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix *CurrShape = gsl_matrix_alloc(n,dim);

    //Getting the current shape positions matrix:
    getPos(CurrShape);

    //calculating dx/de (Jacobian) and reading out dX/de, shape function derivaties:
    gsl_matrix * ShapeFuncDer = ShapeFuncDerivatives[pointNo];
    gsl_matrix * ShapeFuncDerStack = ShapeFuncDerStacks[pointNo];
    gsl_matrix * InvdXde = InvdXdes[pointNo];
    gsl_matrix * Jacobian = gsl_matrix_calloc(dim, dim);
    gsl_matrix * B = Bmatrices[pointNo];
    gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix * invJShFuncDerSWithFe =invJShapeFuncDerStackwithFe[pointNo];
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
    gsl_matrix_transpose(Jacobian);

    //calculating F:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);


    //calculating Fe:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currF, InvFg, 0.0, currFscFe);	///< Removing growth
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currFscFe, InvFsc, 0.0, currFe);	///< Removing shape change

    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(Fg, "Fg");
    //displayMatrix(currFe, "Fe");
    //displayMatrix(currF, "F");
    gsl_matrix * currFeT = gsl_matrix_alloc(dim, dim);
    gsl_matrix_transpose_memcpy(currFeT,currFe);
    createMatrixCopy(FeMatrices[pointNo], currFe); // storing Fe for use in implicit elastic K calculation.

    //calculating E (E = 1/2 *(Fe^T*Fe-I):
    gsl_matrix * E = calculateEForNodalForces(currFe,currFeT);

    //calculating S: (S = D:E)
    gsl_matrix* S = calculateSForNodalForces(E);

    //calculating stress (stress = detFe^-1 Fe S Fe^T):
    gsl_matrix_set_zero(elasticStress[pointNo]);
    gsl_matrix* compactStress  = calculateCompactStressForNodalForces(currFe,S,currFeT, elasticStress[pointNo]);
    //displayMatrix(E,"E");
    //cout<<"pointNo: "<<pointNo<<endl;
    //displayMatrix(compactStress,"compactStress");
    //displayMatrix(currFe,"Fe");
    //displayMatrix(S,"S");
    //Now from stress, I will calculate nodal forces via B.
    //Calculating the inverse Jacobian stack matrix:
    gsl_matrix* InvJacobianStack = calculateInverseJacobianStackForNodalForces(Jacobian);

    //Calculating currB^T:
    detFs[pointNo] = determinant3by3Matrix(currF);
    gsl_matrix* currBT = calculateBTforNodalForces(InvJacobianStack,ShapeFuncDerStack, B, invJShFuncDerS);

    //Calculate invJShapeFuncDerStackwithFe for K calculation (using F in inverse jacobian calculation rather than Fe):
    calculateInvJShFuncDerSWithFe(currFe, InvdXde, ShapeFuncDerStack, invJShFuncDerSWithFe);

    //calculating nodal forces as B^T compactStress detF
    gsl_matrix_scale(currBT,detFs[pointNo]);
    gsl_matrix_scale(currBT,detdXdes[pointNo]);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currBT, compactStress,0.0, currg);

/*
    gsl_matrix* C = CMatrices[pointNo];
    calculateCMatrix(pointNo);
    gsl_matrix* tmp1 = gsl_matrix_calloc(3,3);
    gsl_matrix* InvFe = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* tmpFeForInversion = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFeForInversion,currFe);
    bool inverted = InvertMatrix(tmpFeForInversion, InvFe);
    gsl_matrix* InvFeT = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix_transpose_memcpy(InvFeT,InvFe);
    gsl_matrix* b = gsl_matrix_calloc(nDim,nDim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvFeT, E,0.0, tmp1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp1, InvFe,0.0, b);
    gsl_matrix* bVoigt = gsl_matrix_calloc(6,1);
    gsl_matrix_set(bVoigt,0,0, gsl_matrix_get(b,0,0));
    gsl_matrix_set(bVoigt,1,0, gsl_matrix_get(b,1,1));
    gsl_matrix_set(bVoigt,2,0, gsl_matrix_get(b,2,2));
    gsl_matrix_set(bVoigt,3,0, 2.0*gsl_matrix_get(b,0,1));
    gsl_matrix_set(bVoigt,4,0, 2.0*gsl_matrix_get(b,2,1));
    gsl_matrix_set(bVoigt,5,0, 2.0*gsl_matrix_get(b,0,2));
    gsl_matrix* sigma2 = gsl_matrix_calloc(6,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, C, bVoigt,0.0, sigma2);

    displayMatrix(sigma2,"sigma2");
    displayMatrix(compactStress,"compactStress");
*/
    //freeing the matrices allocated in this function
    gsl_matrix_free(currFeT);
    gsl_matrix_free(currFscFe);
    gsl_matrix_free(E);
    gsl_matrix_free(S);
    gsl_matrix_free(compactStress);
    gsl_matrix_free(InvJacobianStack);
    gsl_matrix_free(currBT);

    //cout<<"Finished calculate nodel forces"<<endl;
}


void Prism::calculateReferenceVolume(){
	double height = 0.0;
	for (int i = 0; i<3; ++i){
        double d = ReferenceShape->Positions[0][i] - ReferenceShape->Positions[3][i];
		d *= d;
		height +=d;
	}
	height = pow(height,0.5);

	double basesideVec1[3];
	double basesideVec2[3];
	double baseSide1 = 0.0;
	double baseSide2 = 0.0;
	double costet =0.0;
	for (int i = 0; i<3; ++i){
		basesideVec1[i]= ReferenceShape-> Positions[1][i] - ReferenceShape-> Positions[0][i];
		basesideVec2[i]= ReferenceShape-> Positions[2][i] - ReferenceShape-> Positions[0][i];
		costet += basesideVec1[i]*basesideVec2[i];
		baseSide1 += basesideVec1[i] * basesideVec1[i];
		baseSide2 += basesideVec2[i] * basesideVec2[i];
	}
	baseSide1 = pow(baseSide1,0.5);
	baseSide2 = pow(baseSide2,0.5);
	costet /= (baseSide1*baseSide2);
	double sintet = pow((1-costet*costet),0.5);
	double baseArea = baseSide1* baseSide2 * sintet / 2.0;
	ReferenceShape->BasalArea = baseArea;
	ReferenceShape->Volume = height * baseArea;
    GrownVolume = ReferenceShape->Volume;
    VolumePerNode = GrownVolume/nNodes;
	//cout<<"baseSide1: "<<baseSide1<<" baseSide2: "<<baseSide2<<" costet: "<<costet<<" sintet: "<<sintet<<endl;
	//cout<<"basearea: "<<baseArea<<" heignt: "<<	height<<" Volume: "<<ReferenceShape->Volume<<endl;
}

void Prism::checkHealth(){
	double** normals;
	normals =new double*[8];
	for (int i =0;i<8; ++i){
		normals[i] = new double[3];
		normals[i][0] = 0.0;
		normals[i][1] = 0.0;
		normals[i][2] = 0.0;
	}
	calculatePlaneNormals(normals);
	bool elementsAreHealthy = checkNodePlaneConsistency(normals);
	if (!elementsAreHealthy){
		cerr<<" Element not healthy! : "<<Id<<endl;
	}
	for (int i=0; i<8;++i){
		delete[] normals[i];
	}
	delete[] normals;
}

void Prism::calculatePlaneNormals(double** normals){
	//Calculating plane normals:
	//plane 0 -> normal 0 - > normal for nodes 0  1  2
	//plane 1 -> normal 1 - > normal for nodes 0  1  3
	//plane 2 -> normal 2 - > normal for nodes 0  2  3
	//plane 3 -> normal 3 - > normal for nodes 3  4  5
	//plane 4 -> normal 4 - > normal for nodes 3  1  4
	//plane 5 -> normal 5 - > normal for nodes 3  2  5
	//plane 6 -> normal 6 - > normal for nodes 1  2  4
	//plane 7 -> normal 7 - > normal for nodes 2  4  5
	int List[8][3]={{0,1,2},{0,3,1},{0,2,3},{3,5,4},{3,4,1},{3,2,5},{1,4,2},{2,4,5}};
	double *u,*v;
	u = new double[3];
	v = new double[3];
	for (int i=0; i<8; ++i){
		assignNodalVector(u,List[i][0],List[i][1]);
		assignNodalVector(v,List[i][0],List[i][2]);
		crossProduct3D(u,v,normals[i]);
		//cout<<" Id: "<<Id<<" u "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
		//cout<<" Id: "<<Id<<" v "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
		//cout<<" Id: "<<Id<<" normals["<<i<<"] "<<normals[i][0]<<" "<<normals[i][1]<<" "<<normals[i][2]<<endl;
	}
	delete[] u;
	delete[] v;
}

void Prism::assignNodalVector(double* vec, int id0, int id1){
	//vector from NodeId id0 to NodeId id1
	if (id0 <0 || id1 <0 || id0>=nNodes || id1>= nNodes){
		cerr<<"Error in node input in nodal vector assignment!"<<endl;
	}
	for (int i=0; i<nDim; ++i){
        vec[i] = Positions[id1][i] - Positions[id0][i];
	}
}

bool Prism::checkNodePlaneConsistency(double** normals){
	//cout<<"inside check consistency, Id: "<<Id<<endl;
	//List of constricting planes for each node:
	//format is for each node (0->5): [plane1, plane2, plane3, node id for plane1, node id  for plane2&3]
	//int ListOfBorder[6][3] = {{4,1,5},{3,0,5},{3,0,4},{4,1,5},{3,0,5},{3,0,4}};
	int List[6][5] = {{3,6,7,3,2},{3,2,5,4,2},{3,1,4,5,1},{0,6,7,0,2},{0,2,5,1,2},{0,1,4,2,1}};
	bool elementHealthy = true;
	double *u;
	u = new double[3];
	for (int i =0; i<nNodes; ++i){
		assignNodalVector(u,List[i][3],i);
		double dotp[3];
		dotp[0] = dotProduct3D(u,normals[List[i][0]]);
		if (dotp[0]<0){
			cerr <<"The element is not consistent! - top/bottom plane, Id: "<<Id<<endl;
			cerr<<"i: "<<i<<endl;
			cerr<<"normals["<<List[i][0]<<"]: "<<normals[List[i][0]][0]<<" "<<normals[List[i][0]][1]<<" "<<normals[List[i][0]][2]<<endl;
			cerr<<"u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
			cerr<<"dotp: "<<dotp[0]<<endl;
			elementHealthy =  false;
		}
		assignNodalVector(u,List[i][4],i);
		dotp[1] = dotProduct3D(u,normals[List[i][1]]);
		dotp[2] = dotProduct3D(u,normals[List[i][2]]);
		if(dotp[1]<0 ||  dotp[2]<0){
			cerr <<"The element is not consistent! side planes, Id: "<<Id<<endl;
			cerr<<"dot 1: "<<dotp[1]<<" dot2 :"<<dotp[2]<<endl;
			cerr<<"1 : normals["<<List[i][1]<<"]: "<<normals[List[i][1]][0]<<" "<<normals[List[i][1]][1]<<" "<<normals[List[i][1]][2]<<endl;
			cerr<<"2 : normals["<<List[i][2]<<"]: "<<normals[List[i][2]][0]<<" "<<normals[List[i][2]][1]<<" "<<normals[List[i][2]][2]<<endl;
			for (int i=0; i<nNodes;++i){
				for (int j =0; j<nDim; ++j){
					cout<<Positions[i][j]<<"  ";
				}
				cout<<endl;
			}
			elementHealthy =  false;
		}
	}
	delete[] u;
	return elementHealthy;
}

double Prism::getApicalSideLengthAverage(){
	double dx,dy,dz;
	double dsum =0.0;
	int pairs[3][2] = {{3,4},{3,5},{4,5}};
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		dsum += pow((dx*dx + dy*dy+dz*dz),0.5);
	}
	dsum /= 3.0;
	return dsum;
}

void Prism::getApicalTriangles(vector <int> &ApicalTriangles){
	ApicalTriangles.push_back(NodeIds[3]);
	ApicalTriangles.push_back(NodeIds[4]);
	ApicalTriangles.push_back(NodeIds[5]);
}

int Prism::getCorrecpondingApical(int currNodeId){
	if (NodeIds[0] == currNodeId){
		return NodeIds[3];
	}
	if (NodeIds[1] == currNodeId){
		return NodeIds[4];
	}
	if (NodeIds[2] == currNodeId){
		return NodeIds[5];
	}
	return -100;
}

bool Prism::IsThisNodeMyBasal(int currNodeId){
	if (NodeIds[0] == currNodeId || NodeIds[1] == currNodeId || NodeIds[2] == currNodeId ){
		return true;
	}
	return false;
}

double Prism::getElementHeight(){
	double dx = Positions[0][0] - Positions[3][0];
	double dy = Positions[0][1] - Positions[3][1];
	double dz = Positions[0][2] - Positions[3][2];
	return pow((dx*dx + dy*dy + dz*dz),0.5);
}


void Prism::AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, int RKId,  double ***SystemForces,  double ***PackingForces, vector<Node*> &Nodes){
	int Id0, Id1, Id2;
	if (tissueplacement == 0){//basal surface:
		Id0 = 0;
		Id1 = 1;
		Id2 = 2;
	}
	if (tissueplacement == 1){//apical surface:
		Id0 = 3;
		Id1 = 4;
		Id2 = 5;
	}
	double F[3];
	F[0] = Fx / 3.0;
	F[1] = Fy / 3.0;
	F[2] = Fz / 3.0;
	for(int j=0; j<nDim; ++j){
		if (!Nodes[NodeIds[Id0]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[Id0]][j] -= F[j];
            PackingForces[RKId][NodeIds[Id0]][j] -= F[j];
		}
		if (!Nodes[NodeIds[Id1]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[Id1]][j] -= F[j];
			PackingForces[RKId][NodeIds[Id1]][j] -= F[j];
		}
		if (!Nodes[NodeIds[Id2]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[Id2]][j] -= F[j];
			PackingForces[RKId][NodeIds[Id2]][j] -= F[j];
		}
	}
}



void Prism::calculateNormalForPacking(int tissuePlacement){
	double * u;
	u = new double[3];
	double * v;
	v = new double[3];
	if (tissuePlacement == 1 ){ //calculating apical packing
		for (int i=0; i<nDim; ++i){
			u[i] = Positions[4][i] - Positions[3][i];
			v[i] = Positions[5][i] - Positions[3][i];
			ApicalNormalForPacking[i] = 0.0;
		}
		crossProduct3D(u,v,ApicalNormalForPacking);
		normaliseVector3D(ApicalNormalForPacking);
	}
	else if (tissuePlacement == 0 ){ //calculating basal packing
		for (int i=0; i<nDim; ++i){
			u[i] = Positions[1][i] - Positions[0][i];
			v[i] = Positions[2][i] - Positions[0][i];
			BasalNormalForPacking[i] = 0.0;
		}
		crossProduct3D(u,v,BasalNormalForPacking);
		normaliseVector3D(BasalNormalForPacking);
	}
	//cerr<<"	Element: "<<Id<<"	u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<" v: "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
	//cerr<<"		normal before normalisation: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	//cerr<<"		normal after normalisation: "<<normalForPacking[0]<<" "<<normalForPacking[1]<<" "<<normalForPacking[2]<<endl;
	for (int i=0; i<nDim; ++i){
		u[i] = Positions[0][i] - Positions[3][i];
	}
	//cerr<<"		vector to basal: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
	double  dot;
	if (tissuePlacement == 1 ){ //calculating apical packing
		dot = dotProduct3D(u,ApicalNormalForPacking);
		if (dot>0){
			for (int i=0; i<nDim; ++i){
				ApicalNormalForPacking[i] *=(-1.0);
			}
		}
		ApicalNormalForPackingUpToDate = true;
	}
	else if (tissuePlacement == 0 ){ //calculating basal packing
		dot = dotProduct3D(u,BasalNormalForPacking);
		if (dot>0){
			for (int i=0; i<nDim; ++i){
				BasalNormalForPacking[i] *=(-1.0);
			}
		}
		BasalNormalForPackingUpToDate = true;
	}
	delete[] v;
	delete[] u;
}


bool Prism::IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode){
	float threshold = 1000.0;
	if (tissueType == 1){ //element is on the peripodial membrane, use the distance threshold of the peripodial membrane
		threshold = Peripodialthreshold;
	}
	else{
		threshold = Columnarthreshold;
	}
	int initial =0, final = 6;
	//check against all the nodes if tissue type of the node is peripodial, only check against the necessary side if node is on the columnar layer
	if (TissueTypeOfPackingNode == 0){  //node is on the columnar layer;
		if (TissuePlacementOfPackingNode == 0){
			//tissue placement of the node is basal, should check it against basal nodes only
			final = 3;
		}
		else if (TissuePlacementOfPackingNode == 1){
			//tissue placement of the node is apical, should check it against basal nodes only
			initial = 3;
		}
	}
	float dmin = 2.0*threshold;
	float dminNeg = (-2.0)*threshold;
	for (int i=initial; i<final; ++i){
		float dx =100.0, dy = 100.0, dz = 100.0;
		dx = Pos[0]-Positions[i][0];
		dy = Pos[1]-Positions[i][1];
		dz = Pos[2]-Positions[i][2];
		if ((dx >0 && dx < dmin) || (dx <0 && dx >dminNeg)){
			if ((dy >0 && dy < dmin) || (dy <0 && dy >dminNeg)){
				if ((dz >0 && dz < dmin) || (dz <0 && dz >dminNeg)){
					return true;
				}
			}
		}
	}
	return false;
}

void  Prism::getApicalNodePos(double* posCorner){
	posCorner[0] = Positions[3][0];
	posCorner[1] = Positions[3][1];
	posCorner[2] = Positions[3][2];
}

void  Prism::getBasalNodePos(double* posCorner){
	posCorner[0] = Positions[0][0];
	posCorner[1] = Positions[0][1];
	posCorner[2] = Positions[0][2];
}

bool  Prism::IspointInsideTriangle(int tissueplacement,double x, double y,double z){
	bool isInside = false;
	//Using Barycentric coordinates:
	int  E0Index = -1, E1Index = -1, E2Index = -1;
	if (tissueplacement ==0 ){ //checking basal surface
		E0Index = 0;
		E1Index = 1;
		E2Index = 2;
	}
	else if (tissueplacement == 1 ){ //checking apical surface
		E0Index = 3;
		E1Index = 4;
		E2Index = 5;
	}
	double *E0E1 = new double[3];
	double *E0E2 = new double[3];

	//cout<<" surface positions: "<<endl;
	//cout<<"	"<<Positions[E0Index][0]<<" "<<Positions[E0Index][1]<<" "<<Positions[E0Index][2]<<endl;
	//cout<<"	"<<Positions[E1Index][0]<<" "<<Positions[E1Index][1]<<" "<<Positions[E1Index][2]<<endl;
	//cout<<"	"<<Positions[E2Index][0]<<" "<<Positions[E2Index][1]<<" "<<Positions[E2Index][2]<<endl;
	//cout<<" point: "<<x<<" "<<y<<" "<<z<<endl;

	for (int i=0; i<3; ++i){
		E0E1[i]=Positions[E1Index][i] - Positions[E0Index][i];
		E0E2[i]=Positions[E2Index][i] - Positions[E0Index][i];
	}
	double *CrossPMain = new double [3];
	crossProduct3D(E0E1,E0E2,CrossPMain);
	double DoubleArea = calculateMagnitudeVector3D(CrossPMain);
	double alpha =0.0, beta = 0.0, gamma = 0.0;
	double *PE1 = new double[3];
	PE1[0] = Positions[E1Index][0] - x;
	PE1[1] = Positions[E1Index][1] - y;
	PE1[2] = Positions[E1Index][2] - z;
	double *PE2 = new double[3];
	PE2[0] = Positions[E2Index][0] - x;
	PE2[1] = Positions[E2Index][1] - y;
	PE2[2] = Positions[E2Index][2] - z;
	double *CrossP = new double [3];
	crossProduct3D(PE1,PE2,CrossP);
	//the vectors should look at te same direction:
	double dotp = dotProduct3D(CrossP,CrossPMain);
	//cout<<"dotp for alpha:  "<<dotp<<" ";
	if (dotp>0){
		alpha = calculateMagnitudeVector3D(CrossP);
		alpha /= DoubleArea;
		//cout<<" alpha: "<<alpha<<" ";
		if (alpha >-1E-10 && alpha <= 1.0+1E-10){
			double *PE0 = new double[3];
			PE0[0] = Positions[E0Index][0] - x;
			PE0[1] = Positions[E0Index][1] - y;
			PE0[2] = Positions[E0Index][2] - z;
			crossProduct3D(PE2,PE0,CrossP);
			dotp = dotProduct3D(CrossP,CrossPMain);
			//cout<<"dotp for beta:  "<<dotp<<" ";
			if (dotp>0){
				beta = calculateMagnitudeVector3D(CrossP);
				beta /= DoubleArea;
				//cout<<" beta: "<<beta<<" ";
				if (beta >-1E-10 && beta <= 1.0+1E-10){
					gamma = 1 - alpha - beta;
					//cout<<" gamma: "<<gamma<<" ";
					crossProduct3D(PE1,PE0,CrossP);
					if (gamma >-1E-10 && gamma <1.0+1E-10){
						isInside =  true;
					}
				}
			}
			delete[] PE0;
		}
		//cout<<endl;
	}

	delete[] E0E1;
	delete[] E0E2;
	delete[] CrossPMain;
	delete[] CrossP;
	delete[] PE1;
	delete[] PE2;
	return isInside;
}
