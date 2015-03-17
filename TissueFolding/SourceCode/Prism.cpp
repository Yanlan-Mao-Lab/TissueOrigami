#include "Prism.h"
#include "ReferenceShapeBase.h"


using namespace std;

Prism::Prism(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId){
	//cout<<"constructing prism"<<endl;
	nNodes = 6;
	nDim = 3;
	Id = CurrId;
	NodeIds = new int[6];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[3];
	CurrGrowthStrainAddition = new double[3];
	ApicalNormalForPacking =  new double[3];
	BasalNormalForPacking =  new double[3];
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
		ApicalNormalForPacking[i] = 0;
		BasalNormalForPacking[i] = 0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	CurrGrowthStrainsUpToDate = false;
	WorldToTissueRotMatUpToDate= false;
	GrowthStrainsRotMatUpToDate= false;
	IsGrowing = false;
	IsChangingShape = false;
	GrewInThePast = false;
	ChangedShapeInThePast = false;
	ApicalNormalForPackingUpToDate = false;
	BasalNormalForPackingUpToDate = false;
	IsAblated = false;
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
	VolumePerNode = 0;
	//RefToTissueRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	//RefToTissueRotMatT = boost::numeric::ublas::zero_matrix<double>(3,3);
	//TissueToWorldRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	//setTissueCoordsRotationsBuffers();
}

Prism::~Prism(){
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
    //cout<<"finalised the destructor for prism class"<<endl;
}

void Prism::setCoeffMat(){
	using namespace boost::numeric::ublas;
	CoeffMat = zero_matrix<int> (6, nDim*nDim);
	CoeffMat(0,0)=1;
	CoeffMat(1,4)=1;
	CoeffMat(2,8)=1;
	CoeffMat(3,1)=1;CoeffMat(3,3)=1;
	CoeffMat(4,5)=1;CoeffMat(4,7)=1;
	CoeffMat(5,2)=1;CoeffMat(5,6)=1;
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


void Prism::getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape){
	using namespace boost::numeric::ublas;
	for (int i =0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
			CurrRelaxedShape(i,j) = ReferenceShape->Positions[i][j];
		}
	}
}

void Prism::setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer, double eta, double zeta, double nu){
	double alpha  = (1 - zeta)/2;
	double beta = (1 + zeta)/2;
	double lambda = 1-eta-nu;

	ShapeFuncDer(0,0)= -alpha;
	ShapeFuncDer(0,1)=  alpha;
	ShapeFuncDer(0,2)=  0;
	ShapeFuncDer(0,3)= -beta;
	ShapeFuncDer(0,4)=	beta;
	ShapeFuncDer(0,5)=	0;

	ShapeFuncDer(1,0)= -alpha;
	ShapeFuncDer(1,1)=	0;
	ShapeFuncDer(1,2)=  alpha;
	ShapeFuncDer(1,3)= -beta;
	ShapeFuncDer(1,4)=  0;
	ShapeFuncDer(1,5)=  beta;

	ShapeFuncDer(2,0)= -lambda/2;
	ShapeFuncDer(2,1)= -eta/2;
	ShapeFuncDer(2,2)= -nu/2;
	ShapeFuncDer(2,3)=	lambda/2;
	ShapeFuncDer(2,4)=	eta/2;
	ShapeFuncDer(2,5)=	nu/2;
}

void Prism::setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack){
	int n = nNodes;
	int dim = nDim;
	for (int i=0; i<n;++i){
		subrange(ShapeFuncDerStack, 0,dim,i*dim,i*dim+1) = subrange(ShapeFuncDer,0,dim,i,i+1);
	}
	subrange(ShapeFuncDerStack, dim,2*dim,1,dim*n) = subrange(ShapeFuncDerStack, 0,dim,0,dim*n-1);
	subrange(ShapeFuncDerStack, 2*dim,3*dim,1,dim*n) = subrange(ShapeFuncDerStack, dim,2*dim,0,dim*n-1);
}

void Prism::calculateReferenceStiffnessMatrix(){
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
				zeta = GaussPoints[zetaiter];
				matrix<double> currk  (dim*n, dim*n);
				matrix<double> currB  (6, dim*n);
				matrix<double> currBE (dim*n, 6);
				matrix<double> currBo (dim*dim, dim*n);
				calculateCurrk(currk, currB, currBE, currBo, eta,zeta,nu);
				kSumZeta  = kSumZeta + GaussCoeff[zetaiter]   * currk;
				BSumZeta  = BSumZeta + GaussCoeff[zetaiter] * currB;
				BESumZeta = BESumZeta + GaussCoeff[zetaiter]  * currBE;
				BoSumZeta = BoSumZeta + GaussCoeff[zetaiter]  * currBo;
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

void Prism::calculateCurrk(boost::numeric::ublas::matrix<double>& currk, boost::numeric::ublas::matrix<double>& currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double>& currBo, double eta, double zeta, double nu){
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
	double detJ = determinant3by3Matrix(Jacobian);
	//Getting the inverse of Jacobian:
	matrix<double> InvJacobian (dim, dim);
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

void Prism::calculateReferenceVolume(){
	double height = 0.0;
	for (int i = 0; i<3; ++i){
		double d = ReferenceShape-> Positions[0][i] - ReferenceShape-> Positions[3][i];
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
	ReferenceShape->Volume = height * baseArea;
	VolumePerNode = ReferenceShape->Volume/nNodes;
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
		vec[i] = PositionsAlignedToReference[id1][i] - PositionsAlignedToReference[id0][i];
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
