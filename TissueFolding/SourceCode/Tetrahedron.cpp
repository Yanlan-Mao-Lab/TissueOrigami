/*
 * Tetrahedron.cpp
 *
 *  Created on: 19 Dec 2014
 *      Author: melda
 */

#include "Tetrahedron.h"
#include "ReferenceShapeBase.h"

using namespace std;

Tetrahedron::Tetrahedron(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId){
	cout<<"constructing tetrahedron"<<endl;
	nNodes = 4;
	nDim = 3;
	Id = CurrId;
	NodeIds = new int[4];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[3];
	CurrGrowthStrainAddition = new double[3];
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	CurrGrowthStrainsUpToDate = false;
	WorldToTissueRotMatUpToDate= false;
	GrowthStrainsRotMatUpToDate= false;
	IsGrowing = false;
	IsChangingShape = false;
	GrewInThePast = false;
	ChangedShapeInThePast = false;

	setIdentificationColour();
	setShapeType("Tetrahedron");
	ReferenceShape = new ReferenceShapeBase("Tetrahedron");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);

	Strain = boost::numeric::ublas::zero_vector<double>(6);
	RK1Strain = boost::numeric::ublas::zero_vector<double>(6);
	StrainTissueMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	PlasticStrain = boost::numeric::ublas::zero_vector<double>(6);
	CurrPlasticStrainsInTissueCoordsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalGrowthStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);

	WorldToTissueRotMat= boost::numeric::ublas::identity_matrix<double>(3,3);
	GrowthStrainsRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
	RotatedElement = false;
	WorldToReferenceRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
	//setting rotation matrices to identity;
	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			int matrixElement;
			if (i==j){
				matrixElement=1;
			}
			else{
				matrixElement=0;
			}
			WorldToReferenceRotMat(i,j) = matrixElement;
		}
	}

	CurrShapeChangeToAdd[0] = 0;
	CurrShapeChangeToAdd[1] = 0;
	CurrShapeChangeToAdd[2] = 0;
	TissueCoordinateSystem = new double[9];
	TissueCoordinateSystem[0]=1.0;
	TissueCoordinateSystem[1]=0.0;
	TissueCoordinateSystem[2]=0.0;
	TissueCoordinateSystem[3]=0.0;
	TissueCoordinateSystem[4]=1.0;
	TissueCoordinateSystem[5]=0.0;
	TissueCoordinateSystem[6]=0.0;
	TissueCoordinateSystem[7]=0.0;
	TissueCoordinateSystem[8]=1.0;
	//RefToTissueRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	//RefToTissueRotMatT = boost::numeric::ublas::zero_matrix<double>(3,3);
	//TissueToWorldRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	//setTissueCoordsRotationsBuffers();
}

Tetrahedron::~Tetrahedron(){
	//cout<<"called the destructor for prism class"<<endl;
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
	}
	delete[] Positions;
	delete[] PositionsAlignedToReference;
	//delete[] PositionsInTissueCoord;
	delete[] NodeIds;
	delete[] IdentifierColour;
	delete[] GrowthRate;
	delete[] ShapeChangeRate;
	delete ReferenceShape;
    //cout<<"finalised the destructor for tetrahedron class"<<endl;
}

void Tetrahedron::setCoeffMat(){
	using namespace boost::numeric::ublas;
	CoeffMat = zero_matrix<int> (6, nDim*nDim);
	CoeffMat(0,0)=1;
	CoeffMat(1,4)=1;
	CoeffMat(2,8)=1;
	CoeffMat(3,1)=1;CoeffMat(3,3)=1;
	CoeffMat(4,5)=1;CoeffMat(4,7)=1;
	CoeffMat(5,2)=1;CoeffMat(5,6)=1;
}


void  Tetrahedron::setElasticProperties(double EApical, double EBasal, double EMid, double v){
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
	using namespace boost::numeric::ublas;
	D = zero_matrix<double>(6,6);
	double multiplier = E/((1+v)*(1-2*v));
	D(0,0)= multiplier*(1-v);	D(0,1)=	multiplier*v;		D(0,2)=	multiplier*v;
	D(1,0)= multiplier*v;		D(1,1)= multiplier*(1-v);	D(1,2)= multiplier*v;
	D(2,0)= multiplier*v;		D(2,1)= multiplier*v;		D(2,2)= multiplier*(1-v);
	D(3,3)= multiplier*(1-2*v)/2;
	D(4,4)= multiplier*(1-2*v)/2;
	D(5,5)= multiplier*(1-2*v)/2;

	//make it behave like a sheet only: (remove all z-response)
	//D(2,0)= 0.0;
	//D(2,1)= 0.0;
	//D(0,2)= 0.0;
	//D(1,2)= 0.0;
	//making it more resistive in Ez
	//D(2,2) *= 3.0;
}


void Tetrahedron::getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape){
	using namespace boost::numeric::ublas;
	for (int i =0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
			CurrRelaxedShape(i,j) = ReferenceShape->Positions[i][j];
		}
	}
}

void Tetrahedron::setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer, double eta, double zeta, double nu){

	ShapeFuncDer(0,0)= -1;
	ShapeFuncDer(0,1)=  1;
	ShapeFuncDer(0,2)=  0;
	ShapeFuncDer(0,3)=	0;

	ShapeFuncDer(1,0)= -1;
	ShapeFuncDer(1,1)=	0;
	ShapeFuncDer(1,2)=  1;
	ShapeFuncDer(1,3)=  0;

	ShapeFuncDer(2,0)= -1;
	ShapeFuncDer(2,1)=  0;
	ShapeFuncDer(2,2)=  0;
	ShapeFuncDer(2,3)=	1;
}

void Tetrahedron::setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack){
	int n = nNodes;
	int dim = nDim;
	for (int i=0; i<n;++i){
		subrange(ShapeFuncDerStack, 0,dim,i*dim,i*dim+1) = subrange(ShapeFuncDer,0,dim,i,i+1);
	}
	subrange(ShapeFuncDerStack, dim,2*dim,1,dim*n) = subrange(ShapeFuncDerStack, 0,dim,0,dim*n-1);
	subrange(ShapeFuncDerStack, 2*dim,3*dim,1,dim*n) = subrange(ShapeFuncDerStack, dim,2*dim,0,dim*n-1);
}

void Tetrahedron::calculateReferenceStiffnessMatrix(){
	const int n = nNodes;
	const int dim = nDim;

	using namespace boost::numeric::ublas;
	//Setting up the current reference shape position matrix:
	matrix<double> CurrRelaxedShape (n, dim);
	getCurrRelaxedShape(CurrRelaxedShape);

	double GaussPoints[4] = {-0.861136312, -0.339981044, 0.339981044, 0.861136312};
	double GaussCoeff[4] = {0.347854845, 0.652145155,  0.652145155, 0.347854845};

	double eta, zeta, nu;
	k  = zero_matrix<double>(dim*n, dim*n);
	B  = zero_matrix<double>(6, dim*n);
	BE = zero_matrix<double>(dim*n,6);
	Bo = zero_matrix<double>(dim*dim,dim*n);
	for (int etaiter = 0; etaiter<4; ++etaiter){
		float EtaLimits[2] = {0,1.0};
		eta = GaussPoints[etaiter];
		eta = (EtaLimits[1]-EtaLimits[0])/2 * eta + (EtaLimits[1]+EtaLimits[0])/2;
		float etaMultiplier = GaussCoeff[etaiter] * (EtaLimits[1]-EtaLimits[0])/2;
		matrix<double> kSumNu  = zero_matrix<double>(dim*n, dim*n);
		matrix<double> BSumNu  = zero_matrix<double>(6, dim*n);
		matrix<double> BESumNu = zero_matrix<double>(dim*n, 6);
		matrix<double> BoSumNu = zero_matrix<double>(dim*dim, dim*n);
		for(int nuiter = 0; nuiter<4; ++nuiter){
			float NuLimits[2] = {0, 1-eta};
			nu = GaussPoints[nuiter];
			nu = (NuLimits[1]-NuLimits[0])/2*nu + (NuLimits[1]+NuLimits[0])/2;
			float nuMultiplier = GaussCoeff[nuiter] * (NuLimits[1]-NuLimits[0])/2;
			matrix<double> kSumZeta  = zero_matrix<double>(dim*n, dim*n);
			matrix<double> BSumZeta  = zero_matrix<double>(6, dim*n);
			matrix<double> BESumZeta = zero_matrix<double>(dim*n, 6);
			matrix<double> BoSumZeta = zero_matrix<double>(dim*dim, dim*n);
			for (int zetaiter = 0; zetaiter<4; ++zetaiter){
				float ZetaLimits[2] = {0, 1-eta-nu};
				zeta = GaussPoints[zetaiter];
				zeta = (ZetaLimits[1]-ZetaLimits[0])/2*zeta + (ZetaLimits[1]+ZetaLimits[0])/2;
				float zetaMultiplier = GaussCoeff[nuiter] * (ZetaLimits[1]-ZetaLimits[0])/2;
				matrix<double> currk  (dim*n, dim*n);
				matrix<double> currB  (6, dim*n);
				matrix<double> currBE (dim*n, 6);
				matrix<double> currBo (dim*dim, dim*n);
				calculateCurrk(currk, currB, currBE, currBo, eta,zeta,nu);
				kSumZeta  = kSumZeta + zetaMultiplier   * currk;
				BSumZeta  = BSumZeta + zetaMultiplier * currB;
				BESumZeta = BESumZeta + zetaMultiplier  * currBE;
				BoSumZeta = BoSumZeta + zetaMultiplier  * currBo;
			}
			kSumNu  = kSumNu  + nuMultiplier * kSumZeta;
			BSumNu  = BSumNu  + nuMultiplier * BSumZeta;
			BESumNu = BESumNu + nuMultiplier * BESumZeta;
			BoSumNu = BoSumNu + nuMultiplier * BoSumZeta;
		}
		k  = k  + etaMultiplier * kSumNu;
		B  = B  + etaMultiplier * BSumNu;
		BE = BE + etaMultiplier * BESumNu;
		Bo = Bo + etaMultiplier * BoSumNu;
	}
	//displayMatrix(k,"k");
}

void Tetrahedron::calculateCurrk(boost::numeric::ublas::matrix<double>& currk, boost::numeric::ublas::matrix<double>& currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double>& currBo, double eta, double zeta, double nu){
	const int n = nNodes;
	const int dim = nDim;

	currB  = boost::numeric::ublas::zero_matrix<double>(6, dim*n);
	currBE = boost::numeric::ublas::zero_matrix<double>(dim*n, 6);
	currk = boost::numeric::ublas::zero_matrix<double>(dim*n, dim*n);
	currBo  = boost::numeric::ublas::zero_matrix<double>(dim*dim, dim*n);
	using namespace boost::numeric::ublas;
	//Setting up the current reference shape position matrix:
	matrix<double> CurrRelaxedShape (n, dim);
	getCurrRelaxedShape(CurrRelaxedShape);

	matrix<double> ShapeFuncDer (dim, n);
	setShapeFunctionDerivatives(ShapeFuncDer,eta,zeta,nu);

	//Generating the shape function derivatives stack:
	int dim2 = dim*dim;
	matrix<double> ShapeFuncDerStack = zero_matrix<double>(dim2, dim*n);
	setShapeFunctionDerivativeStack(ShapeFuncDer,ShapeFuncDerStack);

	matrix<double> Jacobian (dim, dim);
	Jacobian  = zero_matrix<double> (dim,dim);
	//Jacobian =  ShapeFuncDer*CurrRelaxedShape
	boost::numeric::ublas::axpy_prod(ShapeFuncDer,CurrRelaxedShape,Jacobian);

	//Getting the inverse of Jacobian:
	matrix<double> InvJacobian (dim, dim);
	double detJ = determinant3by3Matrix(Jacobian);
	bool inverted = InvertMatrix(Jacobian, InvJacobian);
	if (!inverted){
		cerr<<"Jacobian not inverted!!"<<endl;
	}

	//Generating the inverse Jacobian stack:
	matrix<double> InvJacobianStack = zero_matrix<double>(dim2,dim2);
	for (int i =0; i<dim; i++){
		subrange(InvJacobianStack, i*dim,(i+1)*dim,i*dim,(i+1)*dim) = InvJacobian;
	}

	//Generating currB:
	matrix<double> tmpMat1(6, dim2);
	tmpMat1 = zero_matrix<double>(6, dim2);
	boost::numeric::ublas::axpy_prod(CoeffMat,InvJacobianStack,tmpMat1);
	boost::numeric::ublas::axpy_prod(tmpMat1,ShapeFuncDerStack,currB);
	currBo = ShapeFuncDerStack;

	//Generating currk:
	matrix<double> currBT = trans(currB);
	boost::numeric::ublas::axpy_prod(currBT,detJ*D,currBE);
	boost::numeric::ublas::axpy_prod(currBE,currB,currk);
	//currB = currB*detJ;
	//cout<<"Id: "<<Id<<" detJ: "<<detJ<<endl;
}

void Tetrahedron::calculateReferenceVolume(){
	double *vec1 = new double[3];
	double *vec2 = new double[3];
	double *vec3 = new double[3];
	for (int i=0; i<nDim; ++i){
		vec1[i] = ReferenceShape->Positions[1][i] - ReferenceShape->Positions[0][i];
		vec2[i] = ReferenceShape->Positions[2][i] - ReferenceShape->Positions[0][i];
		vec3[i] = ReferenceShape->Positions[3][i] - ReferenceShape->Positions[0][i];
	}
	double *cross = new double[3];
	crossProduct3D(vec2,vec3,cross);
	ReferenceShape->Volume = (1.0/6.0)* dotProduct3D(vec1,cross);
	if (ReferenceShape->Volume<0){
		ReferenceShape->Volume *=(-1.0);
	}
}

void Tetrahedron::checkHealth(){
	//element heath check function here
}
