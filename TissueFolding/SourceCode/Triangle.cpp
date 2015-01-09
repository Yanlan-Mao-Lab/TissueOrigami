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
	cout<<"constructing triangle"<<endl;
	nNodes = 3;
	nDim = 3;	//the triangle has its  nodes in 3D, but the calculations will only use 2 dimensions.
	Id = CurrId;
	NodeIds = new int[3];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[3];
	CurrGrowthStrainAddition = new double[3];
	for (int i=0; i<3; ++i){
		CurrGrowthStrainAddition[i] = 0;
		GrowthRate[i] = 0;
		ShapeChangeRate[i] = 0;
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
	setShapeType("Triangle");
	ReferenceShape = new ReferenceShapeBase("Triangle");
	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
	setCoeffMat();
	ReferenceShape->height = 5.0;
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

	normalCrossOrder[0] = 1;
	normalCrossOrder[1] = 2;
	cout<<"finalised construction"<<endl;
}

Triangle::~Triangle(){
	//cout<<"called the destructor for triangle class"<<endl;
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
	matrix<double> manualBE  = zero_matrix<double>(dim*n, 3);
	matrix<double> manualBT = trans(manualB);
	boost::numeric::ublas::axpy_prod(manualBT,height*Area*D,manualBE);
	boost::numeric::ublas::axpy_prod(manualBE,manualB,manualk);
	k = manualk;
	B = manualB;
	BE = manualBE;
	Bo = ShapeFuncDerStack;
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
	double x2 = ReferenceShape->Positions[1][0];
	double y2 = ReferenceShape->Positions[1][1];
	double x3 = ReferenceShape->Positions[2][0];
	double y3 = ReferenceShape->Positions[2][1];
	ReferenceShape->Volume = 0.5 * ReferenceShape->height * (x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3);
	if (ReferenceShape->Volume<0){
		ReferenceShape->Volume *=(-1.0);
	}
}

void Triangle::checkHealth(){
}

void Triangle::AlignReferenceApicalNormalToZ(double* SystemCentre){
	calculateApicalNormalCrossOrder(SystemCentre);
	//having the order of vector calculation, I need to get the normal of the apical side for the reference triangle
	double* vec0;
	double* vec1;
	vec0 = new double[3];	//vector from node 0 to node 1 or 2, the selection is done so that the cross-product points towards apical direction
	vec1 = new double[3]; 	//vector from node 0 to (same as above)
	for (int i = 0; i<nDim; ++i){
		vec0[i] = ReferenceShape->Positions[normalCrossOrder[0]][i] - ReferenceShape->Positions[0][i];
		vec1[i] = ReferenceShape->Positions[normalCrossOrder[1]][i] - ReferenceShape->Positions[0][i];
	}
	double* normal;
	normal = new double[3];
	crossProduct3D(vec0,vec1,normal);
	normaliseVector3D(normal);
	//then rotate the reference to have this vector pointing towards (+ve) z;
	double* z;
	z = new double[3];
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = -1.0;
	double c, s;
	calculateRotationAngleSinCos(normal,z,c,s);  //align normal to z
	cout<<"vec0: "<<vec0[0]<<" "<<vec0[1]<<" "<<vec0[2]<<" vec1: "<<vec1[0]<<" "<<vec1[1]<<" "<<vec1[2]<<endl;
	cout<<"normal: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<" sin: "<<s<<" cos: "<<c<<endl;
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(normal,z,rotAx);	//calculating the rotation axis that is perpendicular to both normal and z
		constructRotationMatrix(c,s,rotAx,rotMat);	//calculating the rotation matrix
		rotateReferenceElementByRotationMatrix(rotMat);
		//You need to carry out a rotation:
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
	for (int i = 0; i<nDim; ++i){
		vecToCentre[i] = SystemCentre[i] - ElementCentre[i];
	}
	normaliseVector3D(vecToCentre);

	double dotP = dotProduct3D(normal, vecToCentre);
	if (dotP < 0){
		//the normal is pointing towards the "basal side" that is outside the tissue,
		//not into the lumen between columnar and peripodial layers
		//the normal pointing towards the apical surface is constructed with vec0 X vec1, the order should be reversed.
		normalCrossOrder[0] =2;
		normalCrossOrder[1] =1;
	}
	else{
		normalCrossOrder[0] =1;
		normalCrossOrder[1] =2;
	}
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
	//This normal should be in aligned with (-)ve z; (which is the same vector of the reference element as I have already aligned it!
	double* z;
	z = new double[3];
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = -1.0;
	double c, s;
	calculateRotationAngleSinCos(normal,z,c,s);  //align normal to z
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(normal,z,rotAx);	//calculating the rotation axis that is perpendicular to both normal and z
		constructRotationMatrix(c,s,rotAx,rotMat);	//calculating the rotation matrix
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
