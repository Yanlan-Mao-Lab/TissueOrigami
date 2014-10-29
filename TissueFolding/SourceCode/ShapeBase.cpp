#include "ShapeBase.h"
#include "Node.h"
#include <sstream>


void 	ShapeBase::ParentErrorMessage(){
	cerr<<"You are calling the function from a parent here, check declaration is via pointers"<<endl;
}

void	ShapeBase::setShapeType(string TypeName){
	if (TypeName == "Prism"){
		this->ShapeType = 1;
	}
	else if (TypeName == "PrismLateral"){
		this->ShapeType = 2;
	}
	else{
		this->ShapeType= -100;
	};
}

void	ShapeBase::setIdentificationColour(){
	IdentifierColour[2] = Id % 255;
	int a = (Id - IdentifierColour[2]) / 255;
	IdentifierColour[1] = ( a ) % 255;
	if (a>255){
		IdentifierColour[0] = (a - IdentifierColour[1]) / 255;
	}
	else{
		IdentifierColour[0] = 0;
	}
	//cout<<"IdentifierColour: "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}

int		ShapeBase::getShapeType(){
	return ShapeType;
}

int 	ShapeBase::getId(){
	return Id;
}

int 	ShapeBase::getNodeNumber(){
	return nNodes;
}

int* 	ShapeBase::getNodeIds(){
	return NodeIds;
}

int 	ShapeBase::getDim(){
	return nDim;
}

string 	ShapeBase::getName(){
	string name;
	if (ShapeType == 1){
		name = "Prism";
	}
	else if (ShapeType == 2){
		name = "PrismLateral";
	}
	else{
		name = "Unknown";
	}
	stringstream inter;
	inter.fill('0');
	inter.width(4);
	inter<<Id;
	name = name + inter.str();
	return name;
}

double** ShapeBase::getReferencePos(){
	return ReferenceShape->Positions;
}

double*	 ShapeBase::getReferenceNormal(){
	return ReferenceShape->CurrentNormal;
}

double 	ShapeBase::getYoungModulus(){
	return E;
}

double 	ShapeBase::getPoissonRatio(){
	return v;
}

double* ShapeBase::getGrowthRate(){
	return GrowthRate;
}

double* ShapeBase::getShapeChangeRate(){
	return ShapeChangeRate;
}

double* ShapeBase::getCentre(){
	double* d = new double[3];
	d[0]= 0.0; d[1]= 0.0; d[2]=0.0;
	for (int i = 0; i<nNodes; ++i ){
		for (int j = 0; j< nDim; ++j){
			d[j] += Positions[i][j];
		}
	}
	d[0] /= nNodes; d[1] /= nNodes; d[2] /= nNodes;
	return d;
}

void 	ShapeBase::readNodeIds(int* inpNodeIds){
	for (int i=0; i<nNodes; ++i){
		this->NodeIds[i] = inpNodeIds[i];
	}
}

void 	ShapeBase::displayName(){
	cout<<"Type: "<<this->ShapeType<<" Id: "<<this->Id<<endl;
}

void 	ShapeBase::setPositionMatrix(vector<Node*>& Nodes){
	const int n = nNodes;
	const int dim = nDim;
	Positions = new double*[n];
	PositionsAlignedToReference = new double*[n];
	PositionsInTissueCoord = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		Positions[i] = new double[dim];
		PositionsAlignedToReference[i] = new double[dim];
		PositionsInTissueCoord[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
			PositionsInTissueCoord[i][j] = 0.0;
		}
	}

}

void 	ShapeBase::setReferencePositionMatrix(){
	const int n = nNodes;
	const int dim = nDim;
	ReferenceShape -> Positions = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		ReferenceShape -> Positions[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			ReferenceShape -> Positions[i][j] = Positions[i][j];
		}
	}
}

void 	ShapeBase::updateShapeFromSave(ifstream& file){
	updateNodeIdsFromSave(file);
	updateReferencePositionMatrixFromSave(file);
}

void 	ShapeBase::updateNodeIdsFromSave(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;
		NodeIds[i] = savedId;
	}
}

void 	ShapeBase::updateReferencePositionMatrixFromSave(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			ReferenceShape -> Positions[i][j] = savedPos;
			cout<<"savedPos: "<<savedPos<<endl;
		}
	}
}

void 	ShapeBase::displayPositions(){
	for (int i=0; i<nNodes;++i){
		for (int j =0; j<nDim; ++j){
			cout<<Positions[i][j]<<"  ";
		}
		cout<<endl;
	}
}

int*	ShapeBase::getIdentifierColour(){
	return IdentifierColour;
}

void 	ShapeBase::getStrain(int type, float &StrainMag){
	updateTissueCoordStrain();
	StrainMag = 0.0;
	if (type == 0){
		//this is the average strain
		for (int i=0; i<3; ++i){
			StrainMag += ( StrainTissueMat(i,i) ) ;
		}
		StrainMag /= 3;
	}
	else if (type == 1){
		StrainMag = ( StrainTissueMat(0,0) );
	}
	else if (type == 2){
		StrainMag = ( StrainTissueMat(1,1) );
	}
	else if (type == 3){
		StrainMag = ( StrainTissueMat(2,2) );
	}
	else{
		return;
	}
}

void 	ShapeBase::getPlasticStrain(int type, float &StrainMag){
	updateTissueCoordPlasticStrain(); //this is not updated within the force calculation, as I only need it for display purposes
	StrainMag = 0.0;
	if (type == 0){
		for (int i=0; i<3; ++i){
			StrainMag += CurrPlasticStrainsInTissueCoordsMat(i,i);
		}
		StrainMag /= 3;
	}
	else if (type == 1){
		StrainMag = CurrPlasticStrainsInTissueCoordsMat(0,0);
	}
	else if (type == 2){
		StrainMag = CurrPlasticStrainsInTissueCoordsMat(1,1);
	}
	else if (type == 3){
		StrainMag = CurrPlasticStrainsInTissueCoordsMat(2,2);
	}
	else{
		return;
	}
}

void 	ShapeBase::getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag){
	PysPropMag = 0.0;
	if (type == 0){
		PysPropMag = Nodes[NodeIds[NodeNo]] -> Viscosity;
	}
}

void 	ShapeBase::getPysProp(int type, float &PysPropMag){
	PysPropMag = 0.0;
	if (type ==1){
		PysPropMag = getYoungModulus();
	}
	else if (type == 2 ){
		PysPropMag = getPoissonRatio();
	}
	else if (type ==3){
		double* growth;
		growth = getGrowthRate();
		for (int i =0 ; i< nDim ; ++i){
			PysPropMag += growth[i];
		}
		PysPropMag /= nDim;
	}
	else if (type ==4){
		double* shapechange;
		shapechange = getShapeChangeRate();
		PysPropMag = shapechange[2];
	}
}

void 	ShapeBase::displayIdentifierColour(){
	cout <<" IdentifierColour:  "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}

void 	ShapeBase::resetCurrStepGrowthData(){
	for (int i=0;i<3;++i){
		CurrGrowthStrainAddition[i]  = 0.0;
	}
	CurrGrowthStrainsUpToDate = false;
	IsGrowing = false;
}

void 	ShapeBase::resetCurrStepShapeChangeData(){
	for (int i=0;i<3;++i){
		CurrShapeChangeToAdd[i] = 0.0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	IsChangingShape = false;
}

void 	ShapeBase::updateGrowthToAdd(double* growthscale){
	IsGrowing = true;
	GrewInThePast = true;
	for (int i=0;i<3;++i){
			CurrGrowthStrainAddition[i]  = CurrGrowthStrainAddition[i] + growthscale[i];
	}
}

void 	ShapeBase::growShape(){
	calculateGrowthInLocalCoordinates(CurrGrowthStrainAddition);
	// The growth scale is the growth of this element should feel,
	// on the coordinates of the tissue, on Dorsal-Ventral, Anterior-Posterior, and Apical-Basal axes.
	// I need to get the orientation of this tissue coordinate system first.
	// This coordinate system, and all associated rotation matrices wiill change only when the
	// reference prism rotates, and it could have been updated in shape change functions up till this point.
	// There fore, I will re-calculate only if necessary
}

void 	ShapeBase::updateTissueCoordStrain(){
	boost::numeric::ublas::matrix<double> tmpMat1(3,3);
	boost::numeric::ublas::matrix<double> StrainMat(3,3);
	StrainMat(0,0) = Strain(0);
	StrainMat(1,1) = Strain(1);
	StrainMat(2,2) = Strain(2);
	StrainMat(1,0) = Strain(3);
	StrainMat(1,2) = Strain(4);
	StrainMat(2,0) = Strain(5);
	StrainMat(0,1) = Strain(3);
	StrainMat(2,1) = Strain(4);
	StrainMat(0,2) = Strain(5);
	axpy_prod(trans(WorldToReferenceRotMat),StrainMat,tmpMat1);
	axpy_prod(tmpMat1,WorldToReferenceRotMat,StrainTissueMat);

}

void 	ShapeBase::updateTissueCoordPlasticStrain(){
	boost::numeric::ublas::matrix<double> tmpMat1(3,3);
	boost::numeric::ublas::matrix<double> StrainMat(3,3);
	StrainMat(0,0) = PlasticStrain(0);
	StrainMat(1,1) = PlasticStrain(1);
	StrainMat(2,2) = PlasticStrain(2);
	StrainMat(1,0) = PlasticStrain(3);
	StrainMat(1,2) = PlasticStrain(4);
	StrainMat(2,0) = PlasticStrain(5);
	StrainMat(0,1) = PlasticStrain(3);
	StrainMat(2,1) = PlasticStrain(4);
	StrainMat(0,2) = PlasticStrain(5);
	axpy_prod(trans(WorldToReferenceRotMat),StrainMat,tmpMat1);
	axpy_prod(tmpMat1,WorldToReferenceRotMat,CurrPlasticStrainsInTissueCoordsMat);
}

void 	ShapeBase::calculateGrowthInLocalCoordinates(double * strainsToAdd){
	//get the (+)ve z of the reference aligned onto current coordinate system (use transpose of rotation matrix you already have for alignment)
	boost::numeric::ublas::matrix<double>  ReferenceCoordSystem;
	boost::numeric::ublas::matrix<double>  tmpMat1;
	ReferenceCoordSystem = boost::numeric::ublas::identity_matrix<double>(3,3);
	tmpMat1 = boost::numeric::ublas::zero_matrix<double>(3,3);
	//tmp1 will be reference coordinate system on current shape
	axpy_prod(trans(WorldToReferenceRotMat),ReferenceCoordSystem,tmpMat1);
	//define (+)z vector in world coordinates.
	double* u = new double[3];
	u[0] = 0;
	u[1] = 0;
	u[2] = 1;
	//get the (+)ve z of the element on current coordinates
	double* v = new double[3];
	v[0]=tmpMat1(0,2);
	v[1]=tmpMat1(1,2);
	v[2]=tmpMat1(2,2);
	//get the rotation matrix to align current normal onto (+) ve z (R world to tissue)
	double c, s;
	calculateRotationAngleSinCos(u,v,c,s); //align u onto v
	boost::numeric::ublas::matrix<double> WorldToTissueRotMat;
	WorldToTissueRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(u,v,rotAx);	//calculating the rotation axis that is perpendicular to both u and v
		constructRotationMatrix(c,s,rotAx,rotMat);
		WorldToTissueRotMat(0,0)=rotMat[0];
		WorldToTissueRotMat(0,1)=rotMat[1];
		WorldToTissueRotMat(0,2)=rotMat[2];
		WorldToTissueRotMat(1,0)=rotMat[3];
		WorldToTissueRotMat(1,1)=rotMat[4];
		WorldToTissueRotMat(1,2)=rotMat[5];
		WorldToTissueRotMat(2,0)=rotMat[6];
		WorldToTissueRotMat(2,1)=rotMat[7];
		WorldToTissueRotMat(2,2)=rotMat[8];
		//rotate the (+)ve x &y of the reference aligned to current, with the rotation matrix R world to tissue
		axpy_prod(trans(WorldToTissueRotMat),tmpMat1,ReferenceCoordSystem);
	}
	else{
		ReferenceCoordSystem = tmpMat1;
	}
	boost::numeric::ublas::matrix<double> TissueCoordSystem;
	TissueCoordSystem =boost::numeric::ublas::identity_matrix<double>(3,3);
	tmpMat1 = boost::numeric::ublas::identity_matrix<double>(3,3);
	axpy_prod(WorldToTissueRotMat,tmpMat1,TissueCoordSystem);
	TissueCoordinateSystem[0] = TissueCoordSystem(0,0);
	TissueCoordinateSystem[1] = TissueCoordSystem(1,0);
	TissueCoordinateSystem[2] = TissueCoordSystem(2,0);
	TissueCoordinateSystem[3] = TissueCoordSystem(0,1);
	TissueCoordinateSystem[4] = TissueCoordSystem(1,1);
	TissueCoordinateSystem[5] = TissueCoordSystem(2,1);
	TissueCoordinateSystem[6] = TissueCoordSystem(0,2);
	TissueCoordinateSystem[7] = TissueCoordSystem(1,2);
	TissueCoordinateSystem[8] = TissueCoordSystem(2,2);
	//Now this rotated vectors x, y, and (+)ve z defines the coordinate system you want the growth strains in. Get the rotation matrix, and calculate strains.
	u[0] = 1;
	u[1] = 0;
	u[2] = 0;
	v[0]=ReferenceCoordSystem(0,0);
	v[1]=ReferenceCoordSystem(1,0);
	v[2]=ReferenceCoordSystem(2,0);
	calculateRotationAngleSinCos(u,v,c,s);  //align u to v: align the growth in tissue to growth in local
	boost::numeric::ublas::matrix<double>  CurrLocalGrowthToAdd;
	CurrLocalGrowthToAdd= boost::numeric::ublas::zero_matrix<double>(3,3);
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(u,v,rotAx);	//calculating the rotation axis that is perpendicular to both u and v
		constructRotationMatrix(c,s,rotAx,rotMat);
		boost::numeric::ublas::matrix<double> GrowthStrainsRotMat;
		GrowthStrainsRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
		GrowthStrainsRotMat(0,0)=rotMat[0];
		GrowthStrainsRotMat(0,1)=rotMat[1];
		GrowthStrainsRotMat(0,2)=rotMat[2];
		GrowthStrainsRotMat(1,0)=rotMat[3];
		GrowthStrainsRotMat(1,1)=rotMat[4];
		GrowthStrainsRotMat(1,2)=rotMat[5];
		GrowthStrainsRotMat(2,0)=rotMat[6];
		GrowthStrainsRotMat(2,1)=rotMat[7];
		GrowthStrainsRotMat(2,2)=rotMat[8];
		boost::numeric::ublas::matrix<double>  CurrGrowthToAddTissue;
		CurrGrowthToAddTissue = boost::numeric::ublas::zero_matrix<double>(3,3);
		tmpMat1 =  boost::numeric::ublas::zero_matrix<double>(3,3);
		//cout<<"constructing CurrGrowthToAddTissue"<<endl;
		CurrGrowthToAddTissue(0,0)= strainsToAdd[0];
		CurrGrowthToAddTissue(1,1)= strainsToAdd[1];
		CurrGrowthToAddTissue(2,2)= strainsToAdd[2];
		boost::numeric::ublas::axpy_prod(GrowthStrainsRotMat,CurrGrowthToAddTissue,tmpMat1);
		boost::numeric::ublas::axpy_prod(tmpMat1,trans(GrowthStrainsRotMat),CurrLocalGrowthToAdd);
	}
	else{
		CurrLocalGrowthToAdd(0,0)= strainsToAdd[0];
		CurrLocalGrowthToAdd(1,1)= strainsToAdd[1];
		CurrLocalGrowthToAdd(2,2)= strainsToAdd[2];
	}
	LocalGrowthStrainsMat(0,0) = ( (1.0 + LocalGrowthStrainsMat(0,0)) * (1.0 + CurrLocalGrowthToAdd(0,0)) ) - 1.0;
	LocalGrowthStrainsMat(1,1) = ( (1.0 + LocalGrowthStrainsMat(1,1)) * (1.0 + CurrLocalGrowthToAdd(1,1)) ) - 1.0;
	LocalGrowthStrainsMat(2,2) = ( (1.0 + LocalGrowthStrainsMat(2,2)) * (1.0 + CurrLocalGrowthToAdd(2,2)) ) - 1.0;
}

void 	ShapeBase::updatePositionsAlignedToReferenceWithBuffers(){
	boost::numeric::ublas::matrix<double> PositionsMat(nDim,nNodes);
	boost::numeric::ublas::matrix<double> tmpMat1(nDim,nNodes);
	boost::numeric::ublas::matrix<double> PositionsOnReferenceMat(nDim,nNodes);
	for (int i=0;i<nNodes; ++i){
		for (int j=0;j<nDim; ++j){
			PositionsMat(j,i)=Positions[i][j];
		}
	}
	axpy_prod(WorldToReferenceRotMat,PositionsMat,PositionsOnReferenceMat);
	for (int i=0;i<nNodes; ++i){
		for (int j=0;j<nDim; ++j){
			PositionsAlignedToReference[i][j] = PositionsOnReferenceMat(j,i);
		}
	}
}

void 	ShapeBase::bringShapePositionsToOrigin(double** RefNormalised, double* refCentre){
	for (int j = 0; j<nDim; ++j){
		refCentre[j]=0.0;
		for (int i = 0; i<nNodes; ++i){
			RefNormalised[i][j]=0.0;
		}
	}
	double centre[3] = {0.0,0.0,0.0};
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			centre[j] +=  PositionsAlignedToReference[i][j];
			refCentre[j] +=  ReferenceShape ->Positions[i][j];
		}
	}
	for (int j = 0; j<nDim; ++j){
		centre[j] /= nNodes;
		refCentre[j] /= nNodes;
	}
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			 PositionsAlignedToReference[i][j] -= centre[j];
			 RefNormalised[i][j] = ReferenceShape ->Positions[i][j]-refCentre[j];
		}
	}
}

bool 	ShapeBase::calculateAlignmentScore(double** RefNormalised){
	bool needAlignment=true;
	double threshold = 1E-10;
	double sum = 0.0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			double d = RefNormalised[i][j] -  PositionsAlignedToReference[i][j];
			d *= d;
			sum += d;
		}
	}
	if (sum < threshold){
		needAlignment = false;
	}
	return needAlignment;
}

bool 	ShapeBase::calculateAlignmentRotationMatrix(double** RefNormalised, double* rotMat){
	//Calculating the covarience matrix:
	boost::numeric::ublas::matrix<double>S(3,3);
	boost::numeric::ublas::matrix<double>X(nDim,nNodes);
	boost::numeric::ublas::matrix<double>YT(nNodes,nDim);
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			X(j,i)= PositionsAlignedToReference[i][j];
			YT(i,j)= RefNormalised[i][j];
		}
	}
	boost::numeric::ublas::axpy_prod(X,YT,S);
	gsl_matrix * Sgsl = gsl_matrix_alloc (3, 3);
	gsl_matrix * V = gsl_matrix_alloc (3, 3);
	gsl_matrix * R = gsl_matrix_alloc (3, 3);
	gsl_vector * Sig = gsl_vector_alloc (3);
	gsl_vector * workspace = gsl_vector_alloc (3);
	for (int i = 0; i < 3; i++){
	    for (int j = 0; j < 3; j++){
	    	gsl_matrix_set (Sgsl, i, j, S(i,j));
	    }
	}
	//Singular Value Decomposition of covariance matrix S
	int a  =  gsl_linalg_SV_decomp (Sgsl, V, Sig, workspace);
	boost::numeric::ublas::matrix<double>Vublas(3,3);
	boost::numeric::ublas::matrix<double>UT(3,3);
	boost::numeric::ublas::matrix<double>CurrentRotMat(3,3);
	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			UT(i,j) =  gsl_matrix_get(Sgsl,j,i);
			Vublas(i,j) =  gsl_matrix_get(V,i,j);
		}
	}
	boost::numeric::ublas::axpy_prod(Vublas,UT, CurrentRotMat);
	int counter =0;
	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			rotMat[counter] = CurrentRotMat(i,j);
			counter++;
		}
	}
	double det = determinant3by3Matrix(rotMat);
	if (det<0){
		cout<<"Error! Flipped element, Id: "<<Id<<endl;
		//writing a matrix that is identity, except for the last diagonal element, which is det
		boost::numeric::ublas::matrix<double> weights(3,3);
		weights = boost::numeric::ublas::identity_matrix<double>(3,3);
		boost::numeric::ublas::matrix<double> tempMat(3,3);
		S(2,2)=det;
		//calculating the new rotation matrix, doing my best to recover flipped elements
		boost::numeric::ublas::axpy_prod(Vublas,S, tempMat);
		boost::numeric::ublas::axpy_prod(tempMat,UT, CurrentRotMat);
		int counter =0;
		for (int i=0; i<3; ++i){
			for (int j=0; j<3; ++j){
				rotMat[counter] = CurrentRotMat(i,j);
				counter++;
			}
		}
	}
	gsl_matrix_free (Sgsl);
	gsl_matrix_free (V);
	gsl_matrix_free (R);
	gsl_vector_free (Sig);
	gsl_vector_free (workspace);
	//Now I need to check if there sis only numerical error accumulationg on rotMat, or there is an actual rotation (above 1 degrees):
	double threshold = 0.002; //this is sine 1 degrees
	for (int i=0;i<3;++i){
		for (int j=0;j<3;++j){
			if(i != j){
				if (CurrentRotMat(i,j)>threshold || CurrentRotMat(i,j)< (-1.0*threshold)) {
					return true;
				}
			}
		}
	}
	return false; //none of the off - diagonal terms of the matrix are above the threshold, the current rotation is only niumerical error.

}

double 	ShapeBase::determinant3by3Matrix(double* rotMat){
	double det =0.0;
	det =  rotMat[0]*(rotMat[4]*rotMat[8]-rotMat[5]*rotMat[7]);
	det -= rotMat[1]*(rotMat[3]*rotMat[8]-rotMat[5]*rotMat[6]);
	det += rotMat[2]*(rotMat[3]*rotMat[7]-rotMat[4]*rotMat[6]);
	return det;
}

void 	ShapeBase::alignElementOnReference(){
	updatePositionsAlignedToReferenceWithBuffers();
	const int n = nNodes;
	const int dim = nDim;
	double* refCentre = new double[dim];
	double** RefNormalised = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		RefNormalised[i] = new double[dim];
	}
	bringShapePositionsToOrigin(RefNormalised,refCentre);
	bool needAlignment = calculateAlignmentScore(RefNormalised);
	double* rotMat;
	rotMat = new double[9];
	if (needAlignment){
		bool calculateRotation = calculateAlignmentRotationMatrix(RefNormalised, rotMat);
		if (calculateRotation){
			double u[3];
			for (int i=0; i<nNodes; ++i){
				u[0] = PositionsAlignedToReference[i][0];
				u[1] = PositionsAlignedToReference[i][1];
				u[2] = PositionsAlignedToReference[i][2];
				rotateVectorByRotationMatrix(u,rotMat);
				PositionsAlignedToReference[i][0] = u[0] + refCentre[0];
				PositionsAlignedToReference[i][1] = u[1] + refCentre[1];
				PositionsAlignedToReference[i][2] = u[2] + refCentre[2];
			}
			int counter = 0;
			boost::numeric::ublas::matrix<double>CurrentRotMat(3,3);
			for (int i=0; i<3; ++i){
				for (int j=0; j<3; ++j){
					//WorldToReferenceNormalRotMat(i,j) =  rotMat[counter];
					CurrentRotMat(i,j) = rotMat[counter];
					counter++;
				}
			}
			boost::numeric::ublas::matrix<double>tmpMat(3,3);
			boost::numeric::ublas::axpy_prod(CurrentRotMat,WorldToReferenceRotMat, tmpMat);
			WorldToReferenceRotMat = tmpMat;
		}
	}
	delete[] RefNormalised;
	delete[] rotMat;
	delete 	 refCentre;
}
/*
bool	ShapeBase::areSidesFacingSameDirection(double* RefSide, double* ShapeSide){
	double* RefFace;
	RefFace = new double[3];
	double* ShapeFace;
	ShapeFace = new double[3];
	getCurrentAlignmentFaces(RefSide, ShapeSide, RefFace, ShapeFace);
	double dotp = dotProduct3D(RefFace,ShapeFace);
	if ( dotp < 0){
		return false;
	}
	else{
		return true;
	}
	delete[] RefFace;
	delete[] ShapeFace;
}
*/
void	ShapeBase::calculateRotationAngleSinCos(double* u, double* v, double& c, double& s){
	//aligning u onto v:
	c = dotProduct3D(u,v);
	if (c > 1.0){
		c = 1.0;
		s = 0.0;

	}
	else if( c<-1.0){
		c = -1.0;
		s = 0.0;
	}
	else{
		double tet = acos(c);
		s = sin(tet);
	}
}

void	ShapeBase::calculateRotationAxis(double* u, double* v,double* rotAx){
	//aligning u onto v:
	crossProduct3D(u,v,rotAx);
	normaliseVector3D(rotAx);
}

void	ShapeBase::constructRotationMatrix(double c, double s, double* rotAx, double* rotMat){
	rotMat[0] = c + rotAx[0]*rotAx[0]*(1 - c);
	rotMat[1] = rotAx[0]*rotAx[1]*(1 - c) - rotAx[2]*s;
	rotMat[2] = rotAx[0]*rotAx[2]*(1 - c) + rotAx[1]*s;

	rotMat[3] = rotAx[1]*rotAx[0]*(1 - c) + rotAx[2]*s;
	rotMat[4] = c + rotAx[1]*rotAx[1]*(1 - c);
	rotMat[5] = rotAx[1]*rotAx[2]*(1 - c) - rotAx[0]*s;

	rotMat[6] = rotAx[2]*rotAx[0]*(1 - c) - rotAx[1]*s;
	rotMat[7] = rotAx[2]*rotAx[1]*(1 - c) + rotAx[0]*s;
	rotMat[8] = c + rotAx[2]*rotAx[2]*(1 - c);
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,double* rotMat){
	double x = rotMat[0]*u[0]+rotMat[1]*u[1]+rotMat[2]*u[2];
	double y = rotMat[3]*u[0]+rotMat[4]*u[1]+rotMat[5]*u[2];
	double z = rotMat[6]*u[0]+rotMat[7]*u[1]+rotMat[8]*u[2];
	u[0] = x;
	u[1] = y;
	u[2] = z;
}

void	ShapeBase::calculateForces(int RKid, double **SystemForces, vector <Node*>& Nodes){
	//cout<<"calculating forces"<<endl;
	const int nMult = nNodes*nDim;

	using namespace boost::numeric::ublas;
	boost::numeric::ublas::vector<double> displacement(nMult);
	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]-ReferenceShape->Positions[i][j];
			counter++;
		}
	}

	Forces = zero_vector<double>(nMult);
	Strain = zero_vector<double>(6);
	boost::numeric::ublas::vector<double> PlasticStrainForces = zero_vector<double>(nMult);
	//displayMatrix(BE,"BE");

	//cout<<"multiplying to get forces"<<endl;
	boost::numeric::ublas::axpy_prod(k,displacement,Forces);
	//cout<<"multiplying to get strain"<<endl;
	boost::numeric::ublas::axpy_prod(B,displacement,Strain);
	PlasticStrain(0)= LocalGrowthStrainsMat(0,0);
	PlasticStrain(1)= LocalGrowthStrainsMat(1,1);
	PlasticStrain(2)= LocalGrowthStrainsMat(2,2);
	PlasticStrain(3)= LocalGrowthStrainsMat(1,0);
	PlasticStrain(4)= LocalGrowthStrainsMat(1,2);
	PlasticStrain(5)= LocalGrowthStrainsMat(2,0);
	//displayMatrix(PlasticStrain,"PlasticStrain");
	boost::numeric::ublas::axpy_prod(BE,PlasticStrain,PlasticStrainForces);
	Forces = Forces - PlasticStrainForces;
	//Now I have the forces in tissue coordinate system, I need the forces in world coordinates:
	boost::numeric::ublas::matrix<double>forcesInReferenceCoordsMat(nDim,nNodes);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			forcesInReferenceCoordsMat(j,i)= Forces(counter);
			counter++;
		}
	}
	boost::numeric::ublas::matrix<double>forcesInWorldT(nDim,nNodes);
	boost::numeric::ublas::axpy_prod(trans(WorldToReferenceRotMat),forcesInReferenceCoordsMat,forcesInWorldT);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			Forces(counter) = forcesInWorldT(j,i);
			counter++;
		}
	}
	//Now put the forces in world coordinates into system forces
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			if (!Nodes[NodeIds[i]]->FixedPos[j]){
				SystemForces[NodeIds[i]][j] = SystemForces[NodeIds[i]][j] - Forces(counter);
			}
			counter++;
		}
	}
}

void	ShapeBase::updatePositions(vector<Node*>& Nodes){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
				Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
		}
	}
}

void 	ShapeBase::setGrowthRate(double x, double y, double z){
	GrowthRate[0] = x;
	GrowthRate[1] = y;
	GrowthRate[2] = z;
}

void 	ShapeBase::setShapeChangeRate(double x, double y, double z){
	ShapeChangeRate[0] = x;
	ShapeChangeRate[1] = y;
	ShapeChangeRate[2] = z;
}

void 	ShapeBase::updateGrowthRate(double scalex, double scaley, double scalez){
	//This value is stored as per hour:
	GrowthRate[0] += scalex;
	GrowthRate[1] += scaley;
	GrowthRate[2] += scalez;
}
/*
void 	ShapeBase::updateShapeChangeRate(double scale, int axis){
	int compansatingaxes[2];
	if (axis == 0){
		compansatingaxes[0] = 1;
		compansatingaxes[1] = 2;
	}
	else if (axis == 1){
		compansatingaxes[0] = 0;
		compansatingaxes[1] = 2;
	}
	else if (axis == 2){
		compansatingaxes[0] = 0;
		compansatingaxes[1] = 1;
	}
	float effect = ( 1.0 + ShapeChangeRate[axis] ) * scale;
	float compansation = -0.5 * effect;
	ShapeChangeRate[axis] += effect;
	ShapeChangeRate[compansatingaxes[0]] += compansation;
	ShapeChangeRate[compansatingaxes[2]] += compansation;
}
*/
bool 	ShapeBase::InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse, double& det){
	//Matrix inversion routine.
	//Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<double> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	det = 1.0;
	for(unsigned int i = 0; i < A.size1(); i++) {
		det *= A(i,i); // multiply by elements on diagonal
	    det = det * determinant_sign( pm );
	}
	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

int 	ShapeBase::determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}

void	ShapeBase::crossProduct3D(double* u, double* v, double* cross){
	cross[0] = u[1]*v[2] - u[2]*v[1];
	cross[1] = u[2]*v[0] - u[0]*v[2];
	cross[2] = u[0]*v[1] - u[1]*v[0];
}

void	ShapeBase::normaliseVector3D(double* v){
	double mag2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
		double mag = pow(mag2,0.5);
		v[0] /= mag;
		v[1] /= mag;
		v[2] /= mag;
	}
}

double 	ShapeBase::dotProduct3D(double* u, double* v){
	double dot = 0;
	dot = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	return dot;
}

void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname){
	int m = mat.size1();
	int n = mat.size2();
	cout<<matname<<": "<<endl;

	for (int i =0; i<m; i++){
		for (int j =0; j<n; j++){
			cout.precision(4);
			cout.width(6);
			cout<<mat(i,j)<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname){
	int m = mat.size1();
	int n = mat.size2();
	cout<<matname<<": "<<endl;

	for (int i =0; i<m; i++){
		for (int j =0; j<n; j++){
			cout.precision(4);
			cout.width(6);
			cout<<mat(i,j)<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

void	ShapeBase::displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname){
	int m = vec.size();
	cout<<matname<<": "<<endl;
	for (int i =0; i<m; i++){
		cout.precision(4);
		cout.width(6);
		cout<<vec(i)<<" ";
	}
	cout<<endl;
}
