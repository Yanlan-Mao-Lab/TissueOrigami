/*
 * Triangle.cpp
 *
 *  Created on: 6 Jan 2015
 *      Author: melda
 */

#include "Triangle.h"
#include "ReferenceShapeBase.h"

using namespace std;

Triangle::Triangle(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId, double h){
	//cout<<"constructing triangle: "<<CurrId<<endl;
	nNodes = 3;
	nDim = 3;	//the triangle has its  nodes in 3D, but the calculations will only use 2 dimensions.
	Id = CurrId;
	ShapeDim = 2;	//2D shape
	NodeIds = new int[3];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
	slabHeight = h;
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[3];
	CurrGrowthStrainAddition = new double[3];
	ApicalNormalForPacking =  new double[3];
	BasalNormalForPacking =  new double[3];
	RelativePosInBoundingBox = new double[3];
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] = 0;
		ApicalNormalForPacking[i] = 0.0;
		BasalNormalForPacking[i] = 0.0;
		RelativePosInBoundingBox[i] = 0.0;
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
	IsClippedInDisplay = false;
	capElement = false;
	tiltedElement = false;
	setIdentificationColour();
	setShapeType("Triangle");
	ReferenceShape = new ReferenceShapeBase("Triangle");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	ReferenceShape->height = slabHeight;
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);

	Strain = boost::numeric::ublas::zero_vector<double>(6);
	RK1Strain = boost::numeric::ublas::zero_vector<double>(6);
	StrainTissueMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	PlasticStrain = boost::numeric::ublas::zero_vector<double>(6);
	CurrPlasticStrainsInTissueCoordsMat = boost::numeric::ublas::zero_matrix<double>(3,3);
	LocalGrowthStrainsMat = boost::numeric::ublas::zero_matrix<double>(3,3);

	xGrowthScaling  = boost::numeric::ublas::zero_matrix<double>(3,3);
	yGrowthScaling  = boost::numeric::ublas::zero_matrix<double>(3,3);
	zGrowthScaling  = boost::numeric::ublas::zero_matrix<double>(3,3);

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

	normalCrossOrder[0] = 1;
	normalCrossOrder[1] = 2;
	VolumePerNode = 0;
	//apicalZDir = +1.0;
	//displayPositions();
	//cout<<"finalised construction"<<endl;
}

Triangle::~Triangle(){
	//cout<<"called the destructor for triangle class"<<endl;
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
	}
	delete[] Positions;
	delete[] PositionsAlignedToReference;
	delete[] RelativePosInBoundingBox;
	//delete[] PositionsInTissueCoord;
	delete[] NodeIds;
	delete[] IdentifierColour;
	delete[] GrowthRate;
	delete[] ShapeChangeRate;
	delete ReferenceShape;
    //cout<<"finalised the destructor for tetrahedron class"<<endl;
}


void Triangle::setCoeffMat(){
	using namespace boost::numeric::ublas;
	CoeffMat = zero_matrix<int> (3, (nDim-1.0)*(nDim-1.0));
	CoeffMat(0,0)=1;
	CoeffMat(1,3)=1;
	CoeffMat(2,1)=1;CoeffMat(2,2)=1;
}

void  Triangle::setElasticProperties(double EApical, double EBasal, double EMid, double v){
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
	D = zero_matrix<double>(3,3);
	double multiplier = E/((1+v)*(1-2*v));
	D(0,0)= multiplier*(1-v);	D(0,1)=	multiplier*v;
	D(1,0)= multiplier*v;		D(1,1)= multiplier*(1-v);
	D(2,2)= multiplier*(1-2*v)/2;
}

void Triangle::getCurrRelaxedShape(boost::numeric::ublas::matrix<double> & CurrRelaxedShape){
	using namespace boost::numeric::ublas;
	for (int i =0; i<nNodes; ++i){
		for (int j=0; j<nDim-1; ++j){ //return the x&y coordinates only, the reference shape is flat on z
			CurrRelaxedShape(i,j) = ReferenceShape->Positions[i][j];
		}
	}
}


void Triangle::setShapeFunctionDerivatives(boost::numeric::ublas::matrix<double> &ShapeFuncDer, double eta, double nu){
	ShapeFuncDer(0,0)= -1.0;
	ShapeFuncDer(0,1)=  1.0;
	ShapeFuncDer(0,2)=  0.0;

	ShapeFuncDer(1,0)= -1.0;
	ShapeFuncDer(1,1)=  0.0;
	ShapeFuncDer(1,2)=  1.0;
}

void Triangle::setShapeFunctionDerivativeStack(boost::numeric::ublas::matrix<double> &ShapeFuncDer,boost::numeric::ublas::matrix<double> &ShapeFuncDerStack){
	int n = nNodes;
	int dim = nDim-1;	//calculation in 2D
	for (int i=0; i<n;++i){
		subrange(ShapeFuncDerStack, 0,dim,i*dim,i*dim+1) = subrange(ShapeFuncDer,0,dim,i,i+1);
	}
	subrange(ShapeFuncDerStack, dim,2*dim,1,dim*n) = subrange(ShapeFuncDerStack, 0,dim,0,dim*n-1);
}

void Triangle::calculateReferenceStiffnessMatrix(){
	using namespace boost::numeric::ublas;
	const int n = nNodes;
	const int dim = nDim-1;	//calculation in 2D
	double Area = ReferenceShape->Volume/ReferenceShape->height;
	double height = ReferenceShape->height;
	matrix<double> ShapeFuncDer (dim, n);
	setShapeFunctionDerivatives(ShapeFuncDer,1.0,1.0);
	//Generating the shape function derivatives stack:
	int dim2 = dim*dim;
	matrix<double> ShapeFuncDerStack = zero_matrix<double>(dim2, dim*n);
	setShapeFunctionDerivativeStack(ShapeFuncDer,ShapeFuncDerStack);
	double y23 = (1.0/2.0/Area) * (ReferenceShape->Positions[1][1]-ReferenceShape->Positions[2][1]);
	double y31 = (1.0/2.0/Area) * (ReferenceShape->Positions[2][1]-ReferenceShape->Positions[0][1]);
	double y12 = (1.0/2.0/Area) * (ReferenceShape->Positions[0][1]-ReferenceShape->Positions[1][1]);
	double x32 = (1.0/2.0/Area) * (ReferenceShape->Positions[2][0]-ReferenceShape->Positions[1][0]);
	double x13 = (1.0/2.0/Area) * (ReferenceShape->Positions[0][0]-ReferenceShape->Positions[2][0]);
	double x21 = (1.0/2.0/Area) * (ReferenceShape->Positions[1][0]-ReferenceShape->Positions[0][0]);
	matrix<double> manualk  = zero_matrix<double>(dim*n, dim*n);
	matrix<double> manualB  = zero_matrix<double>(3, dim*n);
	manualB(0,0)=y23; manualB(0,1)=0.0; manualB(0,2)=y31; manualB(0,3)=0.0; manualB(0,4)=y12;  manualB(0,5)=0.0;
	manualB(1,0)=0.0; manualB(1,1)=x32; manualB(1,2)=0.0; manualB(1,3)=x13; manualB(1,4)=0.0;  manualB(1,5)=x21;
	manualB(2,0)=x32; manualB(2,1)=y23; manualB(2,2)=x13; manualB(2,3)=y31; manualB(2,4)=x21;  manualB(2,5)=y12;
	/*if (this->Id == 1166){
		cout<<"Element : "<<Id<<endl;
		cout<<"Reference Shape Positions: "<<endl;
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				cout<<" 	"<<ReferenceShape->Positions[i][j]<<" ";
			}
			cout<<endl;
		}
		cout<<"Area: "<<Area<<endl;
		cout<<"Volume: "<<ReferenceShape->Volume<<endl;
		cout<<"Height: "<<height<<endl;
		cout<<"x13: "<<x13<<" "<<x13*Area*2.0<<endl;
		cout<<"x21: "<<x21<<" "<<x21*Area*2.0<<endl;
		cout<<"x32: "<<x32<<" "<<x32*Area*2.0<<endl;
		cout<<"y12: "<<y12<<" "<<y12*Area*2.0<<endl;
		cout<<"y23: "<<y23<<" "<<y23*Area*2.0<<endl;
		cout<<"y31: "<<y31<<" "<<y31*Area*2.0<<endl;
		displayMatrix(manualB, "manulaB");
	}*/
	matrix<double> manualBE  = zero_matrix<double>(dim*n, 3);
	matrix<double> manualBT = trans(manualB);
	boost::numeric::ublas::axpy_prod(manualBT,height*Area*D,manualBE);
	boost::numeric::ublas::axpy_prod(manualBE,manualB,manualk);
	k = manualk;
	B = manualB;
	BE = manualBE;
	Bo = ShapeFuncDerStack;
	//cout<<"finished stiffness matrix of triangle"<<endl;
}
/*void Triangle::calculateReferenceStiffnessMatrix(){
	const int n = nNodes;
	const int dim = nDim-1;	//calculation in 2D

	using namespace boost::numeric::ublas;
	//Setting up the current reference shape position matrix:
	matrix<double> CurrRelaxedShape (n, dim);
	getCurrRelaxedShape(CurrRelaxedShape);

	double GaussPoints[4] = {-0.861136312, -0.339981044, 0.339981044, 0.861136312};
	double GaussCoeff[4] = {0.347854845, 0.652145155,  0.652145155, 0.347854845};

	double eta, nu, zeta;
	k  = zero_matrix<double>(dim*n, dim*n);
	B  = zero_matrix<double>(3, dim*n);
	BE = zero_matrix<double>(dim*n,3);
	Bo = zero_matrix<double>(dim*dim,dim*n);

	for (int etaiter = 0; etaiter<4; ++etaiter){
		float EtaLimits[2] = {0,1.0};
		eta = GaussPoints[etaiter];
		eta = (EtaLimits[1]-EtaLimits[0])/2 * eta + (EtaLimits[1]+EtaLimits[0])/2;
		float etaMultiplier = GaussCoeff[etaiter] * (EtaLimits[1]-EtaLimits[0])/2;
		matrix<double> kSumNu  = zero_matrix<double>(dim*n, dim*n);
		matrix<double> BSumNu  = zero_matrix<double>(3, dim*n);
		matrix<double> BESumNu = zero_matrix<double>(dim*n, 3);
		matrix<double> BoSumNu = zero_matrix<double>(dim*dim, dim*n);
		for(int nuiter = 0; nuiter<4; ++nuiter){
			//float NuLimits[2] = {0, 1-eta};
			float NuLimits[2] = {0, 1};

			nu = GaussPoints[nuiter];
			nu = (NuLimits[1]-NuLimits[0])/2*nu + (NuLimits[1]+NuLimits[0])/2;
			float nuMultiplier = GaussCoeff[nuiter] * (NuLimits[1]-NuLimits[0])/2;
			matrix<double> currk  (dim*n, dim*n);
			matrix<double> currB  (3, dim*n);
			matrix<double> currBE (dim*n, 3);
			matrix<double> currBo (dim*dim, dim*n);
			calculateCurrk(currk, currB, currBE, currBo, eta, nu);
			kSumNu  = kSumNu  + nuMultiplier * currk;
			BSumNu  = BSumNu  + nuMultiplier * currB;
			BESumNu = BESumNu + nuMultiplier * currBE;
			BoSumNu = BoSumNu + nuMultiplier * currBo;
			//displayMatrix(currk,"currk");
			//displayMatrix(currB,"currB");
			//displayMatrix(currBE,"currBE");
			//displayMatrix(currBo,"currBo");
		}
		k  = k  + etaMultiplier * kSumNu;
		B  = B  + etaMultiplier * BSumNu;
		BE = BE + etaMultiplier * BESumNu;
		Bo = Bo + etaMultiplier * BoSumNu;

		//displayMatrix(kSumNu,"kSumNu");
		//displayMatrix(BSumNu,"BSumNu");
		//displayMatrix(BESumNu,"BESumNu");
		//displayMatrix(BoSumNu,"BoSumNu");
	}
	double A = ReferenceShape->Volume;
	double h=5.0;
	double y23 = (1.0/2.0/A) * (ReferenceShape->Positions[1][1]-ReferenceShape->Positions[2][1]);
	double y31 = (1.0/2.0/A) * (ReferenceShape->Positions[2][1]-ReferenceShape->Positions[0][1]);
	double y12 = (1.0/2.0/A) * (ReferenceShape->Positions[0][1]-ReferenceShape->Positions[1][1]);
	double x32 = (1.0/2.0/A) * (ReferenceShape->Positions[2][0]-ReferenceShape->Positions[1][0]);
	double x13 = (1.0/2.0/A) * (ReferenceShape->Positions[0][0]-ReferenceShape->Positions[2][0]);
	double x21 = (1.0/2.0/A) * (ReferenceShape->Positions[1][0]-ReferenceShape->Positions[0][0]);
	matrix<double> manualk  = zero_matrix<double>(dim*n, dim*n);
	matrix<double> manualB  = zero_matrix<double>(3, dim*n);
	manualB(0,0)=y23; manualB(0,1)=0.0; manualB(0,2)=y31; manualB(0,3)=0.0; manualB(0,4)=y12;  manualB(0,5)=0.0;
	manualB(1,0)=0.0; manualB(1,1)=x32; manualB(1,2)=0.0; manualB(1,3)=x13; manualB(1,4)=0.0;  manualB(1,5)=x21;
	manualB(2,0)=x32; manualB(2,1)=y23; manualB(2,2)=x13; manualB(2,3)=y31; manualB(2,4)=x21;  manualB(2,5)=y12;

	matrix<double> manualBE  = zero_matrix<double>(dim*n, 3);
	matrix<double> manualBT = trans(manualB);
	boost::numeric::ublas::axpy_prod(manualBT,h/4.0/A*D,manualBE);
	boost::numeric::ublas::axpy_prod(manualBE,manualB,manualk);
	displayMatrix(k,"k");
	displayMatrix(manualk,"manualk");
	displayMatrix(B,"B");
	displayMatrix(manualB,"manualB");
	displayMatrix(BE,"BE");
	displayMatrix(manualBE,"manualBE");
	displayMatrix(Bo,"Bo");


}*/

void Triangle::calculateCurrk(boost::numeric::ublas::matrix<double>& currk, boost::numeric::ublas::matrix<double>& currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double>& currBo, double eta, double nu){
	const int n = nNodes;
	const int dim = nDim - 1; //calculation in 2D
	currB  = boost::numeric::ublas::zero_matrix<double>(3, dim*n);
	currBE = boost::numeric::ublas::zero_matrix<double>(dim*n, 3);
	currk = boost::numeric::ublas::zero_matrix<double>(dim*n, dim*n);
	currBo  = boost::numeric::ublas::zero_matrix<double>(dim*dim, dim*n);
	using namespace boost::numeric::ublas;
	//Setting up the current reference shape position matrix:
	matrix<double> CurrRelaxedShape (n, dim);
	getCurrRelaxedShape(CurrRelaxedShape);

	matrix<double> ShapeFuncDer (dim, n);
	setShapeFunctionDerivatives(ShapeFuncDer,eta,nu);
	//Generating the shape function derivatives stack:
	int dim2 = dim*dim;
	matrix<double> ShapeFuncDerStack = zero_matrix<double>(dim2, dim*n);
	setShapeFunctionDerivativeStack(ShapeFuncDer,ShapeFuncDerStack);
	matrix<double> Jacobian (dim, dim);
	Jacobian  = zero_matrix<double> (dim,dim);
	boost::numeric::ublas::axpy_prod(ShapeFuncDer,CurrRelaxedShape,Jacobian);
	//Getting the inverse of Jacobian:
	matrix<double> InvJacobian (dim, dim);
	double detJ = determinant2by2Matrix(Jacobian);
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
	matrix<double> tmpMat1(3, dim2);
	tmpMat1 = zero_matrix<double>(3, dim2);
	boost::numeric::ublas::axpy_prod(CoeffMat,InvJacobianStack,tmpMat1);
	boost::numeric::ublas::axpy_prod(tmpMat1,ShapeFuncDerStack,currB);
	currBo = ShapeFuncDerStack;
	//Generating currk:
	matrix<double> currBT = trans(currB);
	//This line is changed but not tested!!!!!
	boost::numeric::ublas::axpy_prod(currBT,ReferenceShape->height*2*detJ*D,currBE);
	boost::numeric::ublas::axpy_prod(currBE,currB,currk);
	//displayMatrix(currB,"CurrB");
	//currB = currB*detJ;
}

void Triangle::calculateReferenceVolume(){
	double x1 = ReferenceShape->Positions[0][0];
	double y1 = ReferenceShape->Positions[0][1];
	double z1 = ReferenceShape->Positions[0][2];
	double x2 = ReferenceShape->Positions[1][0];
	double y2 = ReferenceShape->Positions[1][1];
	double z2 = ReferenceShape->Positions[1][2];
	double x3 = ReferenceShape->Positions[2][0];
	double y3 = ReferenceShape->Positions[2][1];
	double z3 = ReferenceShape->Positions[2][2];
	double* v1 = new double[3];
	v1[0] = x2-x1;
	v1[1] = y2-y1;
	v1[2] = z2-z1;
	double* v2 = new double [3];
	v2[0] = x3-x1;
	v2[1] = y3-y1;
	v2[2] = z3-z1;
	double* cross = new double[3];
	crossProduct3D(v1,v2,cross);
	double area = 0.5 * pow(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2],0.5) ;
	ReferenceShape->BasalArea = area;
	ReferenceShape->Volume =  ReferenceShape->height * area ;
	/*if (Id == 1166){
		cout<<"Id "<<Id<<" , area: "<<area<<" height: "<< ReferenceShape->height<<" volume: "<<ReferenceShape->Volume <<endl;
		cout<<"Id "<<Id<<" , v1: "<<v1[0]<<" "<<v1[1]<<" "<<v1[2]<<" v2: "<<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<" cross: "<<cross[0]<<" "<<cross[1]<<" "<<cross[2]<<endl;
		cout<<"pnt1: "<<x1<<" "<<y1<<" pnt2: "<<x2<<" "<<y2<<" pnt3: "<<x3<<" "<<endl;
		displayPositions();
		displayReferencePositions();
	}*/
	//ReferenceShape->Volume = 0.5 * ReferenceShape->height * (x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3);
	if (ReferenceShape->Volume<0){
		ReferenceShape->Volume *=(-1.0);
	}
	if (ReferenceShape->BasalArea<0){
		ReferenceShape->BasalArea *=(-1.0);
	}
	VolumePerNode = ReferenceShape->Volume/nNodes;
	//delete[] v1;
	//delete[] v2;
	//delete[] cross;
}

void Triangle::checkHealth(){
}

void Triangle::AlignReferenceApicalNormalToZ(double* SystemCentre){
	calculateApicalNormalCrossOrder(SystemCentre);
	//having the order of vector calculation, I need to get the normal of the apical side for the reference triangle
	double* vec0;
	double* vec1;
	vec0 = new double[3];	//vector from node 0 to node 1 or 2, the selection is done so that the cross-product points towards apical direction (towards the lumen - between peripodial membrane and columnar layer)
	vec1 = new double[3]; 	//vector from node 0 to (same as above)
	for (int i = 0; i<nDim; ++i){
		vec0[i] = ReferenceShape->Positions[normalCrossOrder[0]][i] - ReferenceShape->Positions[0][i];
		vec1[i] = ReferenceShape->Positions[normalCrossOrder[1]][i] - ReferenceShape->Positions[0][i];
	}
	double* normal;
	normal = new double[3];
	crossProduct3D(vec0,vec1,normal);
	normaliseVector3D(normal);

	//then rotate the reference to have the vector pointing towards the lumen, to align to (+)ve z;
	double* z;
	z = new double[3];
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = +1;

	double c, s;
	calculateRotationAngleSinCos(normal,z,c,s);  //align normal to z
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(normal,z,rotAx,c);	//calculating the rotation axis that is perpendicular to both normal and z
		constructRotationMatrix(c,s,rotAx,rotMat);	//calculating the rotation matrix
		//You need to carry out a rotation:
		rotateReferenceElementByRotationMatrix(rotMat);
		delete[] rotAx;
		delete[] rotMat;
	}
	//then re-calculate the stiffness matrix.
	delete[] vec0;
	delete[] vec1;
	delete[] normal;
	delete[] z;
}

void Triangle::calculateApicalNormalCrossOrder(double* SystemCentre){
	double* vec0;
	double* vec1;
	vec0 = new double[3];	//vector from node 0 to 1
	vec1 = new double[3]; 	//vector from node 0 to 2
	for (int i = 0; i<nDim; ++i){
		vec0[i] = Positions[1][i] - Positions[0][i];
		vec1[i] = Positions[2][i] - Positions[0][i];
	}
	double* normal;
	normal = new double[3];
	crossProduct3D(vec0,vec1,normal);
	normaliseVector3D(normal);
	double* ElementCentre = getCentre();
	double* vecToCentre = new double[3];
	double targetPoint[3];
	if(tissueType == 1) {  // If tissue type is peropodial membrane, then point to tissue centre
		targetPoint[0] = SystemCentre[0];
		targetPoint[1] = SystemCentre[1];
		targetPoint[2] = SystemCentre[2];
	}
	else{ // If these are columnar layer elements that are 2D, then aim for hegither than the system centre, to ensure pointing at apical layer
		targetPoint[0] = SystemCentre[0]+slabHeight/2;
		targetPoint[1] = SystemCentre[1]+slabHeight/2;
		targetPoint[2] = SystemCentre[2]+slabHeight/2;
	}
	for (int i = 0; i<nDim; ++i){
		vecToCentre[i] = targetPoint[i] - ElementCentre[i];
	}
	/*if (vecToCentre[2]>=0){
		apicalZDir = +1.0;
	}
	else{
		apicalZDir = -1.0;
	}*/
	normaliseVector3D(vecToCentre);
	double dotP = dotProduct3D(normal, vecToCentre);
	if (dotP < 0){
		//the normal is pointing towards the "basal side" that is outside the tissue,
		//not into the lumen between columnar and peripodial layers
		//the normal pointing towards the apical surface is constructed with vec0 X vec1, the order should be reversed.
		//The triangle has nodes ordered in clock-wise order!
		//Need to correct - swapping 1 & 2:
		int NodeId = NodeIds[1];
		double Pos[3]= {Positions[1][0], Positions[1][1],Positions[1][2]};
		double RefPos[3]= {ReferenceShape->Positions[1][0], ReferenceShape->Positions[1][1], ReferenceShape->Positions[1][2]};
		NodeIds[1] = NodeIds[2];
		NodeIds[2] = NodeId;
		for (int i=0;i<nDim;i++){
			Positions[1][i] = Positions[2][i];
			Positions[2][i]= Pos[i];
			ReferenceShape->Positions[1][i] = ReferenceShape->Positions[2][i];
			ReferenceShape->Positions[2][i]= RefPos[i];
		}
	}
		/*normalCrossOrder[0] =2;
		normalCrossOrder[1] =1;

	}
	else{
		normalCrossOrder[0] =1;
		normalCrossOrder[1] =2;
	}*/
	delete[] vec0;
	delete[] vec1;
	delete[] normal;
	delete[] vecToCentre;
	delete[] ElementCentre;

}

void 	Triangle::correctFor2DAlignment(){
	double* vec0;
	double* vec1;
	vec0 = new double[3];	//vector from node 0 to node 1 or 2, the selection is done so that the cross-product points towards apical direction
	vec1 = new double[3]; 	//vector from node 0 to (same as above)
	for (int i = 0; i<nDim; ++i){
		vec0[i] = PositionsAlignedToReference[normalCrossOrder[0]][i] - PositionsAlignedToReference[0][i];
		vec1[i] = PositionsAlignedToReference[normalCrossOrder[1]][i] - PositionsAlignedToReference[0][i];
	}
	double* normal;
	normal = new double[3];
	crossProduct3D(vec0,vec1,normal);
	normaliseVector3D(normal);
	/*//This normal should be in aligned with the direction of z towards apical side - towards the lumen between columnar and paripodial membrane (which is the same vector of the reference element as I have already aligned it!
	double* z;
	z = new double[3];
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = apicalZDir;
	*/
	//This normal should be in aligned with the direction of z towards apical side - towards the lumen between columnar and paripodial membrane (which is the same vector of the reference element as I have already aligned it!
	double* z;
	z = new double[3];
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = +1;
	double c, s;
	calculateRotationAngleSinCos(normal,z,c,s);  //align normal to z
	if (c<0.9998){
		//cout<<"Element: "<<Id<<" correction for 2D alignment, c: "<<c<<endl;
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(normal,z,rotAx,c);	//calculating the rotation axis that is perpendicular to both normal and z
		constructRotationMatrix(c,s,rotAx,rotMat);	//calculating the rotation matrix
		//cout<<"Element: "<<Id<<" rotMat param - c: "<<c<<" s: "<<s<<" normal: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" rotAx: "<<rotAx[0]<<" "<<rotAx[1]<<" "<<rotAx[2]<<endl;
		//You need to carry out a rotation:
		double u[3];
		for (int i=0; i<nNodes; ++i){
			u[0] = PositionsAlignedToReference[i][0];
			u[1] = PositionsAlignedToReference[i][1];
			u[2] = PositionsAlignedToReference[i][2];
			rotateVectorByRotationMatrix(u,rotMat);
			PositionsAlignedToReference[i][0] = u[0];
			PositionsAlignedToReference[i][1] = u[1];
			PositionsAlignedToReference[i][2] = u[2];
		}
		int counter = 0;
		boost::numeric::ublas::matrix<double>CurrentRotMat(3,3);
		for (int i=0; i<3; ++i){
			for (int j=0; j<3; ++j){
				CurrentRotMat(i,j) = rotMat[counter];
				counter++;
			}
		}
		boost::numeric::ublas::matrix<double>tmpMat(3,3);
		tmpMat = boost::numeric::ublas::zero_matrix<double>(3,3);
		boost::numeric::ublas::axpy_prod(CurrentRotMat,WorldToReferenceRotMat, tmpMat);
		WorldToReferenceRotMat = tmpMat;
		delete[] rotAx;
		delete[] rotMat;
	}
	//then re-calculate the stiffness matrix.
	delete[] vec0;
	delete[] vec1;
	delete[] normal;
	delete[] z;
}

double Triangle::getApicalSideLengthAverage(){
	double dx,dy,dz;
	double dsum =0.0;
	int pairs[3][2] = {{0,1},{0,2},{1,2}};
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		dsum += pow((dx*dx + dy*dy + dz*dz),0.5);
	}
	dsum /= 3.0;
	return dsum;
}

double Triangle::getElementHeight(){
	return slabHeight;
}

void Triangle::AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, int RKId,  double ***SystemForces, double ***PackingForces, vector<Node*> &Nodes){
	double F[3];
	F[0] = Fx / 3.0;
	F[1] = Fy / 3.0;
	F[2] = Fz / 3.0;
	for(int j=0; j<nDim; ++j){
		if (!Nodes[NodeIds[0]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[0]][j]  -= F[j];
			PackingForces[RKId][NodeIds[0]][j] -= F[j];
		}
		if (!Nodes[NodeIds[1]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[1]][j]  -= F[j];
			PackingForces[RKId][NodeIds[1]][j] -= F[j];
		}
		if (!Nodes[NodeIds[2]]->FixedPos[j]){
			SystemForces[RKId][NodeIds[2]][j]  -= F[j];
			PackingForces[RKId][NodeIds[2]][j] -= F[j];
		}
	}
}

void Triangle::calculateNormalForPacking(int tissuePlacement){
	double * u;
	u = new double[3];
	double * v;
	v = new double[3];
	for (int i=0; i<nDim; ++i){
		u[i] = Positions[1][i] - Positions[0][i];
		v[i] = Positions[2][i] - Positions[0][i];
		ApicalNormalForPacking[i] = 0.0;
	}
	//cerr<<"		u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<" v: "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
	crossProduct3D(u,v,ApicalNormalForPacking);
	//cerr<<"		normal before normalisation: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	normaliseVector3D(ApicalNormalForPacking);
	ApicalNormalForPackingUpToDate = true;
	BasalNormalForPacking[0] = (-1.0)*ApicalNormalForPacking[0];
	BasalNormalForPacking[1] = (-1.0)*ApicalNormalForPacking[1];
	BasalNormalForPacking[2] = (-1.0)*ApicalNormalForPacking[2];
	BasalNormalForPackingUpToDate = true;
	delete[] v;
	delete[] u;

}

bool Triangle::IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode){
	float threshold = 1000.0;
	if (tissueType == 1){ //element is on the peripodial membrane, use the distance threshold of the peripodial membrane
		threshold = Peripodialthreshold;
	}
	else{
		threshold = Columnarthreshold;
	}
	float dmin = 2.0*threshold;
	float dminNeg = (-2.0)*threshold;
	for (int i=0; i<3; ++i){
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

void  Triangle::getApicalNodePos(double* posCorner){
	posCorner[0] = Positions[0][0];
	posCorner[1] = Positions[0][1];
	posCorner[2] = Positions[0][2];
}

void  Triangle::getBasalNodePos(double* posCorner){
	posCorner[0] = Positions[0][0];
	posCorner[1] = Positions[0][1];
	posCorner[2] = Positions[0][2];
}

bool Triangle::IspointInsideTriangle(int tissueplacement,double x, double y,double z){
	bool isInside = false;
	//cout<<"Called inside tri check for triangle"<<endl;
	int  E0Index, E1Index, E2Index;
	E0Index = 0;
	E1Index = 1;
	E2Index = 2;
	double *E0E1 = new double[3];
	double *E0E2 = new double[3];
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
	double dotp = dotProduct3D(CrossP,CrossPMain);
	//cout<<"dotp for alpha:  "<<dotp<<" ";
	if (dotp>0){
		alpha = calculateMagnitudeVector3D(CrossP);
		alpha /= DoubleArea;
		//cout<<" alpha: "<<alpha<<" ";
		if (alpha >0 && alpha <= 1.0+1E-10){
			double *PE0 = new double[3];
			PE0[0] = Positions[E0Index][0] - x;
			PE0[1] = Positions[E0Index][1] - y;
			PE0[2] = Positions[E0Index][2] - z;
			crossProduct3D(PE2,PE0,CrossP);
			dotp = dotProduct3D(CrossP,CrossPMain);
			if (dotp>0){
				beta = calculateMagnitudeVector3D(CrossP);
				beta /= DoubleArea;
				delete[] PE0;
				if (beta >0 && beta <= 1.0+1E-10){
					gamma = 1 - alpha - beta;
					if (gamma >0 && gamma <1.0){
						isInside = true;
					}
				}
			}
		}
	}
	delete[] E0E1;
	delete[] E0E2;
	delete[] CrossPMain;
	delete[] CrossP;
	delete[] PE1;
	delete[] PE2;
	return isInside;
}
