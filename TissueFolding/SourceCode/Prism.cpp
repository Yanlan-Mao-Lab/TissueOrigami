#include "Prism.h"
#include "ReferenceShapeBase.h"


using namespace std;

Prism::Prism(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId){
	cout<<"constructing prism"<<endl;
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
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
	}
	updatedReference = false;
	TissueCoordinateSystemUpToDate = false;
	CurrShapeChangeStrainsUpToDate = false;
	CurrGrowthStrainsUpToDate = false;
	IsGrowing = false;
	IsChangingShape = false;
	GrewInThePast = false;
	ChangedShapeInThePast = false;

	setIdentificationColour();
	setShapeType("Prism");
	ReferenceShape = new ReferenceShapeBase("Prism");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	alingmentTurn=0;
	setNormals();
	setRefShapePosBuffers();


	Strain = boost::numeric::ublas::zero_vector<double>(6);
	StrainTissueMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	PlasticStrain = boost::numeric::ublas::zero_vector<double>(6);
	LocalGrowthStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalShapeChangeStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalPlasticStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrGrowthStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrShapeChangeStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrPlasticStrainsInTissueCoordsMat= boost::numeric::ublas::zero_matrix<double>(3,3);

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
	RefToTissueRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatT = boost::numeric::ublas::zero_matrix<double>(3,3);
	TissueToWorldRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	setTissueCoordsRotationsBuffers();
}

Prism::~Prism(){
	//cout<<"called the destructor for prism class"<<endl;
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
	}
	delete[] Positions;
	delete[] NodeIds;
	delete[] IdentifierColour;
	delete[] GrowthRate;
	delete[] ShapeChangeRate;
	delete[] RefShapePosBottomAlignedBuffer;
	delete[] RefShapePosTopAlignedBuffer;
	delete[] TissueCoordinateSystemTopAlignedBuffer;
	delete[] TissueCoordinateSystemBottomAlignedBuffer;
	delete ReferenceShape;
    //cout<<"finalised the destructor for prism class"<<endl;
}

void Prism::setNormals(){
	const int Dim = nDim;
	CurrentNormal = new double[Dim];
	ReferenceShape->CurrentNormal = new double[Dim];
	for (int i=0; i<Dim; ++i){
		CurrentNormal[i]=0.0;
		ReferenceShape->CurrentNormal[i]=0.0;
	}
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

void  Prism::setElasticProperties(double E, double v){
	this -> E = E;
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

	//making it more resistive in Ez
	D(2,2) *= 3.0;
}

void Prism::setViscosity(double ApicalVisc,double BasalVisc, vector <Node*>& Nodes){
	for (int i=0; i<3; ++i){
		Nodes[NodeIds[i]]->Viscosity = BasalVisc;
	}
	for (int i=3; i<6; ++i){
		Nodes[NodeIds[i]]->Viscosity = ApicalVisc;
	}
}

void Prism::setRefShapePosBuffers(){
	const int n = nNodes;
	const int dim = nDim;
	RefShapePosBottomAlignedBuffer = new double*[n];
	RefShapePosTopAlignedBuffer = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		RefShapePosBottomAlignedBuffer[i] = new double[dim];
		RefShapePosTopAlignedBuffer[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			RefShapePosBottomAlignedBuffer[i][j] = Positions[i][j];
			RefShapePosTopAlignedBuffer[i][j] = Positions[i][j];
		}
	}

}

void Prism::setStiffnessMatrixBuffers(){
	//cout<<"assigning stiffness matrix buffers"<<endl;
	kTopAlignedBuffer = k;
	kBottomAlignedBuffer = k;
	BTopAlignedBuffer = B;
	BBottomAlignedBuffer = B;
	BETopAlignedBuffer = BE;
	BEBottomAlignedBuffer = BE;
}

void Prism::setTissueCoordsRotationsBuffers(){
	calculateTissueCoordinateSystem();
	TissueCoordinateSystemTopAlignedBuffer = new double[9];
	TissueCoordinateSystemBottomAlignedBuffer = new double[9];
	for (int i=0;i<9;++i){
		TissueCoordinateSystemTopAlignedBuffer[i] = TissueCoordinateSystem[i];
		TissueCoordinateSystemBottomAlignedBuffer[i] = TissueCoordinateSystem[i];
	}
	TissueToWorldRotMatTopAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	TissueToWorldRotMatBottomAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	TissueToWorldRotMatTopAlignedBuffer = TissueToWorldRotMat;
	TissueToWorldRotMatBottomAlignedBuffer = TissueToWorldRotMat;

	calculateRotationMatrixReferenceToTissue();
	RefToTissueRotMatTopAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatBottomAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatTTopAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatTBottomAlignedBuffer = boost::numeric::ublas::zero_matrix<double>(3,3);
	RefToTissueRotMatTopAlignedBuffer = RefToTissueRotMat;
	RefToTissueRotMatBottomAlignedBuffer = RefToTissueRotMat;
	RefToTissueRotMatTTopAlignedBuffer = RefToTissueRotMatT;
	RefToTissueRotMatTBottomAlignedBuffer = RefToTissueRotMatT;
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
	for (int etaiter = 0; etaiter<4; ++etaiter){
		float EtaLimits[2] = {0,1.0};
		eta = GaussPoints[etaiter];
		eta = (EtaLimits[1]-EtaLimits[0])/2 * eta + (EtaLimits[1]+EtaLimits[0])/2;
		float etaMultiplier = GaussCoeff[etaiter] * (EtaLimits[1]-EtaLimits[0])/2;
		matrix<double> kSumNu  = zero_matrix<double>(dim*n, dim*n);
		matrix<double> BSumNu  = zero_matrix<double>(6, dim*n);
		matrix<double> BESumNu = zero_matrix<double>(dim*n, 6);
		for(int nuiter = 0; nuiter<4; ++nuiter){
			float NuLimits[2] = {0, 1-eta};
			nu = GaussPoints[nuiter];
			nu = (NuLimits[1]-NuLimits[0])/2*nu + (NuLimits[1]+NuLimits[0])/2;
			float nuMultiplier = GaussCoeff[nuiter] * (NuLimits[1]-NuLimits[0])/2;
			matrix<double> kSumZeta  = zero_matrix<double>(dim*n, dim*n);
			matrix<double> BSumZeta  = zero_matrix<double>(6, dim*n);
			matrix<double> BESumZeta = zero_matrix<double>(dim*n, 6);
			for (int zetaiter = 0; zetaiter<4; ++zetaiter){
				zeta = GaussPoints[zetaiter];
				matrix<double> currk  (dim*n, dim*n);
				matrix<double> currB  (6, dim*n);
				matrix<double> currBE (dim*n, 6);
				calculateCurrk(currk, currB, currBE, eta,zeta,nu);
				kSumZeta  = kSumZeta + GaussCoeff[zetaiter]   * currk;
				BSumZeta  = BSumZeta + GaussCoeff[zetaiter] * currB;
				BESumZeta = BESumZeta + GaussCoeff[zetaiter]  * currBE;
			}
			kSumNu  = kSumNu  + nuMultiplier * kSumZeta;
			BSumNu  = BSumNu  + nuMultiplier * BSumZeta;
			BESumNu = BESumNu + nuMultiplier * BESumZeta;
		}
		k  = k  + etaMultiplier * kSumNu;
		B  = B  + etaMultiplier * BSumNu;
		BE = BE + etaMultiplier * BESumNu;
	}
	//displayMatrix(k,"k");
}

void Prism::calculateCurrk(boost::numeric::ublas::matrix<double>& currk, boost::numeric::ublas::matrix<double>& currB, boost::numeric::ublas::matrix<double>& currBE, double eta, double zeta, double nu){
	const int n = nNodes;
	const int dim = nDim;

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
	//Jacobian =  ShapeFuncDer*CurrRelaxedShape
	boost::numeric::ublas::axpy_prod(ShapeFuncDer,CurrRelaxedShape,Jacobian);

	//Getting the inverse of Jacobian:
	matrix<double> InvJacobian (dim, dim);
	double detJ;
	bool inverted = InvertMatrix(Jacobian, InvJacobian, detJ);
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

	boost::numeric::ublas::axpy_prod(CoeffMat,InvJacobianStack,tmpMat1);
	boost::numeric::ublas::axpy_prod(tmpMat1,ShapeFuncDerStack,currB);

	//Generating currk:
	matrix<double> currBT = trans(currB);

	boost::numeric::ublas::axpy_prod(currBT,detJ*D,currBE);
	boost::numeric::ublas::axpy_prod(currBE,currB,currk);

	//display the matrices:
	//displayMatrix(ShapeFuncDer,"ShapeFuncDerStack");
	//displayMatrix(ShapeFuncDerStack,"ShapeFuncDerStack");
	//displayMatrix(Jacobian,"Jacobian");
	//displayMatrix(InvJacobian,"InvJacobian");
	//displayMatrix(InvJacobianStack,"InvJacobianStack");
	//displayMatrix(currB,"currB");
	//displayMatrix(currk,"currk");
}

void Prism::calculateNormals(){
	if (alingmentTurn == 0){
		calculateNormalToBottom();
		calculateReferenceNormalToBottom();
	}
	else{
		calculateNormalToTop();
		calculateReferenceNormalToTop();
	}
}

void Prism::updateAlignmentTurn(){
	if (alingmentTurn == 0){
		alingmentTurn = 1;
	}
	else{
		alingmentTurn = 0;
	}
}

void Prism::updateReferenceShapeBaseFromBuffer(){
	//cout<<"updating from buffer"<<endl;
	if (alingmentTurn == 0){
		//the alignment is to the bottom of the prism
		//I will take the reference shape and the k & B from buffer
		//maybe I would not need to do the update:
		for (int i = 0 ; i < nNodes; ++i ){
			for (int j=0; j < nDim; ++j){
				RefShapePosTopAlignedBuffer[i][j] = ReferenceShape->Positions[i][j];
				ReferenceShape->Positions[i][j] = RefShapePosBottomAlignedBuffer[i][j];
			}
		}
		kTopAlignedBuffer = k;
		BTopAlignedBuffer = B;
		BETopAlignedBuffer = BE;
		k = kBottomAlignedBuffer;
		B = BBottomAlignedBuffer;
		BE = BEBottomAlignedBuffer;

		TissueToWorldRotMatTopAlignedBuffer = TissueToWorldRotMat;
		RefToTissueRotMatTopAlignedBuffer = RefToTissueRotMat;
		RefToTissueRotMatTTopAlignedBuffer = RefToTissueRotMatT;
		for (int i=0; i<9;++i){
			TissueCoordinateSystemTopAlignedBuffer[i] = TissueCoordinateSystem[i];
		}
		CurrShapeChangeStrainsUpToDate = false;
		CurrGrowthStrainsUpToDate = false;

		TissueToWorldRotMat = TissueToWorldRotMatBottomAlignedBuffer;
		RefToTissueRotMat = RefToTissueRotMatBottomAlignedBuffer;
		RefToTissueRotMatT = RefToTissueRotMatTBottomAlignedBuffer;
		for (int i=0; i<9;++i){
			TissueCoordinateSystem[i] = TissueCoordinateSystemBottomAlignedBuffer[i];
		}
	}
	else {
		for (int i = 0 ; i < nNodes; ++i ){
			for (int j=0; j < nDim; ++j){
				RefShapePosBottomAlignedBuffer[i][j] = ReferenceShape->Positions[i][j];
				ReferenceShape->Positions[i][j] = RefShapePosTopAlignedBuffer[i][j];
			}
		}
		kBottomAlignedBuffer = k;
		BBottomAlignedBuffer = B;
		BEBottomAlignedBuffer = BE;
		k = kTopAlignedBuffer;
		B = BTopAlignedBuffer;
		BE = BETopAlignedBuffer;

		TissueToWorldRotMatBottomAlignedBuffer = TissueToWorldRotMat;
		RefToTissueRotMatBottomAlignedBuffer = RefToTissueRotMat;
		RefToTissueRotMatTBottomAlignedBuffer = RefToTissueRotMatT;
		for (int i=0; i<9;++i){
			TissueCoordinateSystemBottomAlignedBuffer[i] = TissueCoordinateSystem[i];
		}
		CurrShapeChangeStrainsUpToDate = false;
		CurrGrowthStrainsUpToDate = false;

		TissueToWorldRotMat = TissueToWorldRotMatTopAlignedBuffer;
		RefToTissueRotMat = RefToTissueRotMatTopAlignedBuffer;
		RefToTissueRotMatT = RefToTissueRotMatTTopAlignedBuffer;
		for (int i=0; i<9;++i){
			TissueCoordinateSystem[i] = TissueCoordinateSystemTopAlignedBuffer[i];
		};
	}
	//cout<<"finalised updating from buffer"<<endl;
}

void Prism::resetBuffersAfterGrowth(){
	for (int i = 0 ; i < nNodes; ++i ){
		for (int j=0; j < nDim; ++j){
			RefShapePosBottomAlignedBuffer[i][j] = ReferenceShape->Positions[i][j];
			RefShapePosTopAlignedBuffer[i][j] =ReferenceShape->Positions[i][j];
		}
	}
}

void Prism::calculateNormalToBottom(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = Positions[1][0] - Positions[0][0];
	vec1[1] = Positions[1][1] - Positions[0][1];
	vec1[2] = Positions[1][2] - Positions[0][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = Positions[2][0] - Positions[0][0];
	vec2[1] = Positions[2][1] - Positions[0][1];
	vec2[2] = Positions[2][2] - Positions[0][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	//now I have the normal, but the normal must look into the prism
	//For the bottom case, the vector should have a less than 90 degree angle
	//with the vector towards the top. lets say from node 0 to 3:
	double* vec3;
	vec3 = new double[3];
	vec3[0] = Positions[3][0] - Positions[0][0];
	vec3[1] = Positions[3][1] - Positions[0][1];
	vec3[2] = Positions[3][2] - Positions[0][2];
	double dotp = dotProduct3D(vec3, CurrentNormal);
	if(dotp < 0){
		CurrentNormal[0] *= (-1.0);
		CurrentNormal[1] *= (-1.0);
		CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(CurrentNormal);
	//cout<<"Normal calculated: "<<CurrentNormal[0]<<" "<<CurrentNormal[1]<<" "<<CurrentNormal[2]<<endl;
}

void Prism::calculateReferenceNormalToBottom(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = ReferenceShape->Positions[1][0] - ReferenceShape->Positions[0][0];
	vec1[1] = ReferenceShape->Positions[1][1] - ReferenceShape->Positions[0][1];
	vec1[2] = ReferenceShape->Positions[1][2] - ReferenceShape->Positions[0][2];
	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = ReferenceShape->Positions[2][0] - ReferenceShape->Positions[0][0];
	vec2[1] = ReferenceShape->Positions[2][1] - ReferenceShape->Positions[0][1];
	vec2[2] = ReferenceShape->Positions[2][2] - ReferenceShape->Positions[0][2];

	crossProduct3D(vec1,vec2,ReferenceShape->CurrentNormal);
	//now I have the normal, but the normal must look into the prism
	//For the bottom case, the vector should have a less than 90 degree angle
	//with the vector towards the top. lets say from node 0 to 3:
	double* vec3;
	vec3 = new double[3];
	vec3[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
	vec3[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
	vec3[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	double dotp = dotProduct3D(vec3, ReferenceShape->CurrentNormal);
	if(dotp < 0 ){
		ReferenceShape->CurrentNormal[0] *= (-1.0);
		ReferenceShape->CurrentNormal[1] *= (-1.0);
		ReferenceShape->CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(ReferenceShape->CurrentNormal);
}

void Prism::calculateNormalToTop(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = Positions[4][0] - Positions[3][0];
	vec1[1] = Positions[4][1] - Positions[3][1];
	vec1[2] = Positions[4][2] - Positions[3][2];

	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = Positions[5][0] - Positions[3][0];
	vec2[1] = Positions[5][1] - Positions[3][1];
	vec2[2] = Positions[5][2] - Positions[3][2];

	//normal vector to the triangle
	crossProduct3D(vec1,vec2,CurrentNormal);
	double* vec3;
	vec3 = new double[3];
	vec3[0] = Positions[0][0] - Positions[3][0];
	vec3[1] = Positions[0][1] - Positions[3][1];
	vec3[2] = Positions[0][2] - Positions[3][2];
	double dotp = dotProduct3D(vec3, CurrentNormal);
	if(dotp < 0 ){
		CurrentNormal[0] *= (-1.0);
		CurrentNormal[1] *= (-1.0);
		CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(CurrentNormal);
	//cout<<"Normal calculated: "<<CurrentNormal[0]<<" "<<CurrentNormal[1]<<" "<<CurrentNormal[2]<<endl;
}

void Prism::calculateReferenceNormalToTop(){
	//vector from point 0 to point 1:
	double* vec1;
	vec1 = new double[3];
	vec1[0] = ReferenceShape->Positions[4][0] - ReferenceShape->Positions[3][0];
	vec1[1] = ReferenceShape->Positions[4][1] - ReferenceShape->Positions[3][1];
	vec1[2] = ReferenceShape->Positions[4][2] - ReferenceShape->Positions[3][2];
	//vector from point 0 to point 2:
	double* vec2;
	vec2 = new double[3];
	vec2[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[4][0];
	vec2[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[4][1];
	vec2[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[4][2];

	crossProduct3D(vec1,vec2,ReferenceShape->CurrentNormal);
	double* vec3;
	vec3 = new double[3];
	vec3[0] = ReferenceShape->Positions[0][0] - ReferenceShape->Positions[3][0];
	vec3[1] = ReferenceShape->Positions[0][1] - ReferenceShape->Positions[3][1];
	vec3[2] = ReferenceShape->Positions[0][2] - ReferenceShape->Positions[3][2];
	double dotp = dotProduct3D(vec3, ReferenceShape->CurrentNormal);
	if(dotp < 0 ){
		ReferenceShape->CurrentNormal[0] *= (-1.0);
		ReferenceShape->CurrentNormal[1] *= (-1.0);
		ReferenceShape->CurrentNormal[2] *= (-1.0);
	}
	normaliseVector3D(ReferenceShape->CurrentNormal);
}

void Prism::getCurrentAlignmentSides(double* RefSide, double* ShapeSide){
	//side vector for alignment:
	if (alingmentTurn == 0){
		ShapeSide[0] = Positions[1][0] - Positions[0][0];
		ShapeSide[1] = Positions[1][1] - Positions[0][1];
		ShapeSide[2] = Positions[1][2] - Positions[0][2];
		RefSide[0] = ReferenceShape->Positions[1][0] - ReferenceShape->Positions[0][0];
		RefSide[1] = ReferenceShape->Positions[1][1] - ReferenceShape->Positions[0][1];
		RefSide[2] = ReferenceShape->Positions[1][2] - ReferenceShape->Positions[0][2];
	}
	else{
		ShapeSide[0] = Positions[4][0] - Positions[3][0];
		ShapeSide[1] = Positions[4][1] - Positions[3][1];
		ShapeSide[2] = Positions[4][2] - Positions[3][2];
		RefSide[0] = ReferenceShape->Positions[4][0] - ReferenceShape->Positions[3][0];
		RefSide[1] = ReferenceShape->Positions[4][1] - ReferenceShape->Positions[3][1];
		RefSide[2] = ReferenceShape->Positions[4][2] - ReferenceShape->Positions[3][2];
	}
	normaliseVector3D(RefSide);
	normaliseVector3D(ShapeSide);
}

void Prism::getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace){
	double* Ref2ndSide;
	Ref2ndSide = new double[3];
	double* Shape2ndSide;
	Shape2ndSide = new double[3];
	//side vector for alignment:
	if (alingmentTurn == 0){
		Shape2ndSide[0] = Positions[2][0] - Positions[0][0];
		Shape2ndSide[1] = Positions[2][1] - Positions[0][1];
		Shape2ndSide[2] = Positions[2][2] - Positions[0][2];
		Ref2ndSide[0] = ReferenceShape->Positions[2][0] - ReferenceShape->Positions[0][0];
		Ref2ndSide[1] = ReferenceShape->Positions[2][1] - ReferenceShape->Positions[0][1];
		Ref2ndSide[2] = ReferenceShape->Positions[2][2] - ReferenceShape->Positions[0][2];
	}
	else{
		Shape2ndSide[0] = Positions[5][0] - Positions[3][0];
		Shape2ndSide[1] = Positions[5][1] - Positions[3][1];
		Shape2ndSide[2] = Positions[5][2] - Positions[3][2];
		Ref2ndSide[0] = ReferenceShape->Positions[5][0] - ReferenceShape->Positions[3][0];
		Ref2ndSide[1] = ReferenceShape->Positions[5][1] - ReferenceShape->Positions[3][1];
		Ref2ndSide[2] = ReferenceShape->Positions[5][2] - ReferenceShape->Positions[3][2];
	}
	crossProduct3D(ShapeSide,Shape2ndSide,ShapeFace);
	crossProduct3D(RefSide,Ref2ndSide,RefFace);
	cout<<"Shape2ndSide: "<<Shape2ndSide[0]<<" "<<Shape2ndSide[1]<<" "<<Shape2ndSide[2]<<endl;
	cout<<"Ref2ndSide: "<<Ref2ndSide[0]<<" "<<Ref2ndSide[1]<<" "<<Ref2ndSide[2]<<endl;
	cout<<"ShapeFace: "<<ShapeFace[0]<<" "<<ShapeFace[1]<<" "<<ShapeFace[2]<<endl;
	cout<<"RefFace: "<<RefFace[0]<<" "<<RefFace[1]<<" "<<RefFace[2]<<endl;
}

void Prism::calculateZVecForTissueCoordAlignment(double* u){
	u[0] = ReferenceShape->Positions[3][0] - ReferenceShape->Positions[0][0];
	u[1] = ReferenceShape->Positions[3][1] - ReferenceShape->Positions[0][1];
	u[2] = ReferenceShape->Positions[3][2] - ReferenceShape->Positions[0][2];
	normaliseVector3D(u);
}

void Prism::calculateXVecForTissueCoordAlignment(double* u ){
	u[0] = ReferenceShape->Positions[1][0] - ReferenceShape->Positions[0][0];
	u[1] = ReferenceShape->Positions[1][1] - ReferenceShape->Positions[0][1];
	u[2] = ReferenceShape->Positions[1][2] - ReferenceShape->Positions[0][2];
	normaliseVector3D(u);
}
