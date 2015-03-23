#include "ShapeBase.h"
#include "Node.h"
#include <sstream>


void 	ShapeBase::ParentErrorMessage(string functionName){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
}

bool 	ShapeBase::ParentErrorMessage(string functionName, bool returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

double 	ShapeBase::ParentErrorMessage(string functionName, double returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

int 	ShapeBase::ParentErrorMessage(string functionName, int returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

void	ShapeBase::setShapeType(string TypeName){
	//cout<<"inside set shape type"<<endl;
	if (TypeName == "Prism"){
		this->ShapeType = 1;
	}
	else if (TypeName == "PrismLateral"){
		this->ShapeType = 2;
	}
	else if (TypeName == "Tetrahedron"){
		this->ShapeType = 3;
	}
	else if (TypeName == "Triangle"){
		//cout<<"set shape type to triangle"<<endl;
		this->ShapeType = 4;
	}
	else{
		this->ShapeType= -100;
	};
	//cout<<"finalised set shape type"<<endl;
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
	else if (ShapeType == 3){
		name = "Tetrahedron";
	}
	else if (ShapeType == 4){
		name = "Triangle";
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
			PositionsAlignedToReference[i][j] = 0.0;
		}
	}
}

void 	ShapeBase::setTissuePlacement(vector<Node*>& Nodes){
	bool hasApicalNode = false;
	bool hasBasalNode = false;
	bool hasLateralNode = false;
	for (int i = 0; i<nNodes; ++i){
		if (Nodes[NodeIds[i]]->tissuePlacement == 1){
			hasApicalNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissuePlacement == 0){
			hasBasalNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissuePlacement == 3){
			hasLateralNode = true;
		}
	}
	if (hasLateralNode){
		tissuePlacement = 3;
	}
	else{
		if (hasApicalNode){
			if (hasBasalNode){
				//the element spans through the whole tissue, the mid-line value should be used
				tissuePlacement = 2;
			}
			else{
				//the element has only apical and midline nodes, it is apical
				tissuePlacement = 1;
			}
		}
		else if (hasBasalNode){
			//the element only has basal and mid-line nodes, it is basal
			tissuePlacement = 0;
		}
		else{
			//the element has only mid-line nodes, it is mid-line
			tissuePlacement = 2;
		}
	}
}



void 	ShapeBase::setTissueType(vector<Node*>& Nodes){
	bool hasColumnarNode = false;
	bool hasPeripodiumNode = false;
	for (int i = 0; i<nNodes; ++i){
		if (Nodes[NodeIds[i]]->tissueType == 0){
			hasColumnarNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissueType == 1){
			hasPeripodiumNode = true;
		}
	}
	if (hasPeripodiumNode){
		tissueType = 1;
	}
	else if (hasColumnarNode){
		//ASK PERIPODIUM FIRST, SOME PERIPODIUM ELEMENTS CAN HAVE COLUMNAR NODES, NO COLUMNAR ELEMENT SHOULD HAVE A PERIPODIUM NODE
		tissueType = 0;
	}
	else {
		cerr<<"Element is not placed into tissue correctly, Id: "<<Id<<endl;
	}
	//cout<<"Element : "<<Id<<" hasColumnarNode: "<<hasColumnarNode<<" hasPeripodiumNode "<<hasPeripodiumNode<<" tissueType: "<<tissueType<<endl;
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
	file >> IsAblated;
	updateNodeIdsFromSave(file);
	updateReferencePositionMatrixFromSave(file);
	//displayName();
	//displayPositions();
	//displayReferencePositions();
	//displayNodeIds();
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
			//cout<<"savedPos: "<<savedPos<<endl;
		}
	}
}

void 	ShapeBase::updateReferencePositionMatrixFromMeshInput(ifstream& file){
	updateReferencePositionMatrixFromSave(file);
}

void ShapeBase::updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes){
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);
}

void 	ShapeBase::displayNodeIds(){
	for (int i=0; i<nNodes;++i){
			cout<<NodeIds[i]<<"  ";
		cout<<endl;
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

void 	ShapeBase::displayReferencePositions(){
	for (int i=0; i<nNodes;++i){
		for (int j =0; j<nDim; ++j){
			cout<<ReferenceShape ->Positions[i][j]<<"  ";
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
		for (int i=0; i<6; ++i){
			for (int j=i; j<3; ++j){
				if (i==j || StrainTissueMat(i,j)>0 ){
					StrainMag += ( StrainTissueMat(i,j) ) ;
				}
				else{
					StrainMag -= ( StrainTissueMat(i,j) ) ;
				}
			}
		}
		StrainMag /= 6;
		//if (Id == 300 || Id == 301){
		//	cout<<"Element: "<<Id<<endl;
		//	displayMatrix(Strain,"Strain");
		//	displayMatrix(StrainTissueMat,"StrainTissueMat");
		//}
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
			for (int j=i; j<3; ++j){
				if (i==j || CurrPlasticStrainsInTissueCoordsMat(i,j)>0 ){
					StrainMag += ( CurrPlasticStrainsInTissueCoordsMat(i,j) ) ;
				}
				else{
					StrainMag -= ( CurrPlasticStrainsInTissueCoordsMat(i,j) ) ;
				}
			}
		}
		StrainMag /= 6;
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
			CurrGrowthStrainAddition[i]  += growthscale[i];
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
	tmpMat1 = boost::numeric::ublas::zero_matrix<double> (3,3);
	StrainMat = boost::numeric::ublas::zero_matrix<double> (3,3);
	StrainMat(0,0) = Strain(0);
	StrainMat(1,1) = Strain(1);
	StrainMat(2,2) = Strain(2);
	StrainMat(1,0) = Strain(3);
	StrainMat(1,2) = Strain(4);
	StrainMat(2,0) = Strain(5);
	StrainMat(0,1) = Strain(3);
	StrainMat(2,1) = Strain(4);
	StrainMat(0,2) = Strain(5);
	if(!GrowthStrainsRotMatUpToDate){
		//updating the rotation matrix:
		double* RefCoords;
		RefCoords = new double[9];
		calculateReferenceCoordSysAlignedToTissue(RefCoords);
		//Now this rotated vectors x, y, and (+)ve z defines the coordinate system you want the growth strains in. Get the rotation matrix, and calculate strains.
		//get (+)x vector in reference coordinates aligned to tissue in z.
		double* v = new double[3];
		v[0]=RefCoords[0];
		v[1]=RefCoords[1];
		v[2]=RefCoords[2];
		bool rotateMatrix = calculateGrowthStrainsRotMat(v);
	}
	axpy_prod(trans(GrowthStrainsRotMat),StrainMat,tmpMat1);
	axpy_prod(tmpMat1,GrowthStrainsRotMat,StrainTissueMat);

}

void 	ShapeBase::updateTissueCoordPlasticStrain(){
	boost::numeric::ublas::matrix<double> tmpMat1(3,3);
	boost::numeric::ublas::matrix<double> StrainMat(3,3);
	tmpMat1 = boost::numeric::ublas::zero_matrix<double> (3,3);
	StrainMat = boost::numeric::ublas::zero_matrix<double> (3,3);
	StrainMat(0,0) = PlasticStrain(0);
	StrainMat(1,1) = PlasticStrain(1);
	StrainMat(2,2) = PlasticStrain(2);
	StrainMat(1,0) = PlasticStrain(3);
	StrainMat(1,2) = PlasticStrain(4);
	StrainMat(2,0) = PlasticStrain(5);
	StrainMat(0,1) = PlasticStrain(3);
	StrainMat(2,1) = PlasticStrain(4);
	StrainMat(0,2) = PlasticStrain(5);
	if(!GrowthStrainsRotMatUpToDate){
		//updating the rotation matrix:
		double* RefCoords;
		RefCoords = new double[9];
		calculateReferenceCoordSysAlignedToTissue(RefCoords);
		//Now this rotated vectors x, y, and (+)ve z defines the coordinate system you want the growth strains in. Get the rotation matrix, and calculate strains.
		//get (+)x vector in reference coordinates aligned to tissue in z.
		double* v = new double[3];
		v[0]=RefCoords[0];
		v[1]=RefCoords[1];
		v[2]=RefCoords[2];
		bool rotateMatrix = calculateGrowthStrainsRotMat(v);
	}
	axpy_prod(trans(GrowthStrainsRotMat),StrainMat,tmpMat1);
	axpy_prod(tmpMat1,GrowthStrainsRotMat,CurrPlasticStrainsInTissueCoordsMat);
}

void ShapeBase::getTissueCoordinaSystem(double* TissueCoords){
	if(!WorldToTissueRotMatUpToDate){
		//The element has not been grown yet, the tissue coordinates are not up-to-date:
		boost::numeric::ublas::matrix<double>  ReferenceCoordSysAlignedToTissue;
		boost::numeric::ublas::matrix<double>  ReferenceCoordSysAlignedOnCurrentShape;
		ReferenceCoordSysAlignedToTissue = boost::numeric::ublas::identity_matrix<double>(3,3);
		ReferenceCoordSysAlignedOnCurrentShape = boost::numeric::ublas::zero_matrix<double>(3,3);
		axpy_prod(trans(WorldToReferenceRotMat),ReferenceCoordSysAlignedToTissue,ReferenceCoordSysAlignedOnCurrentShape);
		//get the (+)ve z of the element on current coordinates
		double* v = new double[3];
		v[0]=ReferenceCoordSysAlignedOnCurrentShape(0,2);
		v[1]=ReferenceCoordSysAlignedOnCurrentShape(1,2);
		v[2]=ReferenceCoordSysAlignedOnCurrentShape(2,2);
		bool rotateMatrix = calculateWorldToTissueRotMat(v);
	}
	updateTissueCoordinateSystem(TissueCoords);
}

void 	ShapeBase::updateTissueCoordinateSystem(double* TissueCoords){
	//this is again only necessary for display
	boost::numeric::ublas::matrix<double> TissueCoordSystem;
	boost::numeric::ublas::matrix<double> tmpMat1;
	TissueCoordSystem =boost::numeric::ublas::identity_matrix<double>(3,3);
	tmpMat1 = boost::numeric::ublas::identity_matrix<double>(3,3);
	axpy_prod(WorldToTissueRotMat,tmpMat1,TissueCoordSystem);
	//0,1,2 will give the x coordinate, 3,4,5 will give the y coordinate etc.
	TissueCoords[0] = TissueCoordSystem(0,0);
	TissueCoords[1] = TissueCoordSystem(1,0);
	TissueCoords[2] = TissueCoordSystem(2,0);
	TissueCoords[3] = TissueCoordSystem(0,1);
	TissueCoords[4] = TissueCoordSystem(1,1);
	TissueCoords[5] = TissueCoordSystem(2,1);
	TissueCoords[6] = TissueCoordSystem(0,2);
	TissueCoords[7] = TissueCoordSystem(1,2);
	TissueCoords[8] = TissueCoordSystem(2,2);
}

bool 	ShapeBase::calculateWorldToTissueRotMat(double* v){
	WorldToTissueRotMatUpToDate = true;
	//define (+)z vector in world coordinates.
	double* u = new double[3];
	u[0] = 0;
	u[1] = 0;
	u[2] = 1;
	double c, s;
	calculateRotationAngleSinCos(u,v,c,s); //align u onto v
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
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
		delete[] rotAx;
		delete[] rotMat;
		delete[] u;
		//The rotation matrix to align coordinate systems is calculated. You will need to multiply the coordinate system
		return true;
	}
	else{
		WorldToTissueRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
		//The matrices are already almost aligned (less than 1 degree). Do not need to multiply the coordinate system
		delete[] u;
		return false;
	}
}

void 	ShapeBase::calculateReferenceCoordSysAlignedToTissue(double* RefCoordSys){
	boost::numeric::ublas::matrix<double>  IdentityCoordinateSystem;
	boost::numeric::ublas::matrix<double>  ReferenceCoordSysAlignedOnCurrentShape;
	IdentityCoordinateSystem  = boost::numeric::ublas::identity_matrix<double>(3,3);
	ReferenceCoordSysAlignedOnCurrentShape = boost::numeric::ublas::zero_matrix<double>(3,3);
	axpy_prod(trans(WorldToReferenceRotMat),IdentityCoordinateSystem,ReferenceCoordSysAlignedOnCurrentShape);
	//get the (+)ve z of the reference element on current coordinates
	double* v = new double[3];
	v[0]=ReferenceCoordSysAlignedOnCurrentShape(0,2);
	v[1]=ReferenceCoordSysAlignedOnCurrentShape(1,2);
	v[2]=ReferenceCoordSysAlignedOnCurrentShape(2,2);
	bool rotateMatrix = calculateWorldToTissueRotMat(v);
	boost::numeric::ublas::matrix<double>  ReferenceCoordSysAlignedToTissue;
	ReferenceCoordSysAlignedToTissue = boost::numeric::ublas::identity_matrix<double>(3,3);
	if(rotateMatrix){
		//rotate the (+)ve x &y of the reference aligned to current, with the rotation matrix R world to tissue
		axpy_prod(trans(WorldToTissueRotMat),ReferenceCoordSysAlignedOnCurrentShape,ReferenceCoordSysAlignedToTissue);
	}
	else{
		ReferenceCoordSysAlignedToTissue = ReferenceCoordSysAlignedOnCurrentShape;
	}
	RefCoordSys[0]= ReferenceCoordSysAlignedToTissue(0,0);
	RefCoordSys[1]= ReferenceCoordSysAlignedToTissue(1,0);
	RefCoordSys[2]= ReferenceCoordSysAlignedToTissue(2,0);
	RefCoordSys[3]= ReferenceCoordSysAlignedToTissue(0,1);
	RefCoordSys[4]= ReferenceCoordSysAlignedToTissue(1,1);
	RefCoordSys[5]= ReferenceCoordSysAlignedToTissue(2,1);
	RefCoordSys[6]= ReferenceCoordSysAlignedToTissue(0,2);
	RefCoordSys[7]= ReferenceCoordSysAlignedToTissue(1,2);
	RefCoordSys[8]= ReferenceCoordSysAlignedToTissue(2,2);
	delete[] v;

}

bool 	ShapeBase::calculateGrowthStrainsRotMat(double* v){
	GrowthStrainsRotMatUpToDate = true;
	GrowthStrainsRotMat = boost::numeric::ublas::identity_matrix<double>(3,3);
	//define (+)x vector in world coordinates.
	double* u = new double[3];
	u[0] = 1;
	u[1] = 0;
	u[2] = 0;
	double c, s;
	calculateRotationAngleSinCos(u,v,c,s);  //align u to v: align the growth in tissue to growth in local
	if (c<0.9998){
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
		constructRotationMatrix(c,s,rotAx,rotMat);
		GrowthStrainsRotMat(0,0)=rotMat[0];
		GrowthStrainsRotMat(0,1)=rotMat[1];
		GrowthStrainsRotMat(0,2)=rotMat[2];
		GrowthStrainsRotMat(1,0)=rotMat[3];
		GrowthStrainsRotMat(1,1)=rotMat[4];
		GrowthStrainsRotMat(1,2)=rotMat[5];
		GrowthStrainsRotMat(2,0)=rotMat[6];
		GrowthStrainsRotMat(2,1)=rotMat[7];
		GrowthStrainsRotMat(2,2)=rotMat[8];
		//You need to carry out a rotation:
		return true;
	}
	else{
		//there is no need to multiply the matrices
		return false;
	}
}

void 	ShapeBase::calculateGrowthInLocalCoordinates(double * strainsToAdd){
	//get the (+)ve z of the reference aligned onto current coordinate system (use transpose of rotation matrix you already have for alignment)
	double* RefCoords;
	RefCoords = new double[9];
	calculateReferenceCoordSysAlignedToTissue(RefCoords);
	//Now this rotated vectors x, y, and (+)ve z defines the coordinate system you want the growth strains in. Get the rotation matrix, and calculate strains.
	//get (+)x vector in reference coordinates aligned to tissue in z.
	double* v = new double[3];
	v[0]=RefCoords[0];
	v[1]=RefCoords[1];
	v[2]=RefCoords[2];
	bool rotateMatrix = calculateGrowthStrainsRotMat(v);
	boost::numeric::ublas::matrix<double>  CurrLocalGrowthToAdd;
	CurrLocalGrowthToAdd = boost::numeric::ublas::zero_matrix<double>(3,3);
	if (rotateMatrix){
		boost::numeric::ublas::matrix<double>  CurrGrowthToAddTissue;
		CurrGrowthToAddTissue = boost::numeric::ublas::zero_matrix<double>(3,3);
		boost::numeric::ublas::matrix<double>tmpMat1;
		tmpMat1 =  boost::numeric::ublas::zero_matrix<double>(3,3);
		//cout<<"constructing CurrGrowthToAddTissue"<<endl;
		CurrGrowthToAddTissue(0,0)= strainsToAdd[0];
		CurrGrowthToAddTissue(1,1)= strainsToAdd[1];
		CurrGrowthToAddTissue(2,2)= strainsToAdd[2];
		//boost::numeric::ublas::axpy_prod(GrowthStrainsRotMat,CurrGrowthToAddTissue,tmpMat1);
		//boost::numeric::ublas::axpy_prod(tmpMat1,trans(GrowthStrainsRotMat),CurrLocalGrowthToAdd);
		boost::numeric::ublas::matrix<double>R;
		R =  boost::numeric::ublas::identity_matrix<double>(3,3);
		R(2,2)=2.0;
		boost::numeric::ublas::matrix<double>RT;
		RT =  boost::numeric::ublas::identity_matrix<double>(3,3);
		RT(2,2)=0.5;
		boost::numeric::ublas::matrix<double>tmpMat2;
		tmpMat2 =  boost::numeric::ublas::zero_matrix<double>(3,3);
		boost::numeric::ublas::matrix<double>tmpMat3;
		tmpMat3 =  boost::numeric::ublas::zero_matrix<double>(3,3);
		boost::numeric::ublas::axpy_prod(R,GrowthStrainsRotMat,tmpMat1);
		boost::numeric::ublas::axpy_prod(tmpMat1,RT,tmpMat2);
		boost::numeric::ublas::axpy_prod(tmpMat2,CurrGrowthToAddTissue,tmpMat3);
		//boost::numeric::ublas::axpy_prod(R,trans(GrowthStrainsRotMat),tmpMat1);
		//boost::numeric::ublas::axpy_prod(tmpMat1,RT,tmpMat2);
		boost::numeric::ublas::axpy_prod(tmpMat3,trans(tmpMat2),CurrLocalGrowthToAdd);
	}
	else{
		CurrLocalGrowthToAdd(0,0)= strainsToAdd[0];
		CurrLocalGrowthToAdd(1,1)= strainsToAdd[1];
		CurrLocalGrowthToAdd(2,2)= strainsToAdd[2];
	}
	//writing as a upper triangular
	LocalGrowthStrainsMat(0,0) = ( (1.0 + LocalGrowthStrainsMat(0,0)) * (1.0 + CurrLocalGrowthToAdd(0,0)) ) - 1.0;
	LocalGrowthStrainsMat(1,1) = ( (1.0 + LocalGrowthStrainsMat(1,1)) * (1.0 + CurrLocalGrowthToAdd(1,1)) ) - 1.0;
	LocalGrowthStrainsMat(2,2) = ( (1.0 + LocalGrowthStrainsMat(2,2)) * (1.0 + CurrLocalGrowthToAdd(2,2)) ) - 1.0;
	LocalGrowthStrainsMat(0,1) = ( (1.0 + LocalGrowthStrainsMat(0,1)) * (1.0 + CurrLocalGrowthToAdd(0,1)) ) - 1.0;
	LocalGrowthStrainsMat(0,2) = ( (1.0 + LocalGrowthStrainsMat(0,2)) * (1.0 + CurrLocalGrowthToAdd(0,2)) ) - 1.0;
	LocalGrowthStrainsMat(1,2) = ( (1.0 + LocalGrowthStrainsMat(1,2)) * (1.0 + CurrLocalGrowthToAdd(1,2)) ) - 1.0;

	//cout<<"Element: "<<Id<<endl;
	//displayMatrix(LocalGrowthStrainsMat,"LocalGrowthStrainsMat");

	delete[] RefCoords;
}

void 	ShapeBase::updatePositionsAlignedToReferenceWithBuffers(){
	boost::numeric::ublas::matrix<double> PositionsMat(nDim,nNodes);
	boost::numeric::ublas::matrix<double> tmpMat1(nDim,nNodes);
	boost::numeric::ublas::matrix<double> PositionsOnReferenceMat(nDim,nNodes);
	PositionsOnReferenceMat = boost::numeric::ublas::zero_matrix<double> (nDim,nNodes);
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

void 	ShapeBase::bringPositionAlignedToReferenceToOrigin(double* refCentre){
	double centre[3] = {0.0,0.0,0.0};
	refCentre[0]=0.0;
	refCentre[1]=0.0;
	refCentre[2]=0.0;
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
	//threshold = 1000;
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

bool	ShapeBase::calculateDisplacementGradientRotationMatrix(double** RefNormalised, double* rotMat){
	using namespace boost::numeric::ublas;
	/*const int nMult = nNodes*nDim;
	boost::numeric::ublas::vector<double> displacement(nMult);
	boost::numeric::ublas::vector<double> refpos(nMult);

	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]- RefNormalised[i][j];
			refpos(counter) =  RefNormalised[i][j];
			counter++;
		}
	}*/
	boost::numeric::ublas::matrix<double>S(nDim,nDim);
	boost::numeric::ublas::matrix<double>P(nDim,nDim);
	boost::numeric::ublas::vector<double>dude(nDim*nDim);
	boost::numeric::ublas::vector<double>dXde(nDim*nDim);
	P = boost::numeric::ublas::zero_matrix<double>(nDim,nDim);
	S = boost::numeric::ublas::zero_matrix<double>(nDim,nDim);
	dXde = boost::numeric::ublas::zero_vector<double>(nDim*nDim);
	dude = boost::numeric::ublas::zero_vector<double>(nDim*nDim);
	//Bo is the shape functions derivative stack integrated over volume.
	//The deformation gradient matrix is
	calculatedudEdXde(RefNormalised, dude, dXde);
	//S(0,0) = e(0); S(1,1) = e(1); S(2,2) = e(2);
	//S(1,0) = e(3); S(2,0) = e(5); S(2,1) = e(4);
	//S(0,1) = e(3); S(0,2) = e(5); S(1,2) = e(4);
	S(0,0) = dude(0); S(0,1) = dude(1); S(0,2) = dude(2);
	S(1,0) = dude(3); S(1,1) = dude(4); S(1,2) = dude(5);
	S(2,0) = dude(6); S(2,1) = dude(7); S(2,2) = dude(8);
	P(0,0) = dXde(0); P(0,1) = dXde(1); P(0,2) = dXde(2);
	P(1,0) = dXde(3); P(1,1) = dXde(4); P(1,2) = dXde(5);
	P(2,0) = dXde(6); P(2,1) = dXde(7); P(2,2) = dXde(8);
	boost::numeric::ublas::matrix<double>InvP(nDim,nDim);
	InvP = boost::numeric::ublas::zero_matrix<double>(nDim,nDim);
	//double detP = 0.0;
	bool inverted = InvertMatrix(P, InvP);
	boost::numeric::ublas::matrix<double>F(nDim,nDim);
	F = boost::numeric::ublas::zero_matrix<double>(nDim,nDim);
	boost::numeric::ublas::axpy_prod(S,InvP,F);
	boost::numeric::ublas::matrix<double>FIdentity(nDim,nDim);
	FIdentity = boost::numeric::ublas::identity_matrix<double>(nDim,nDim);
	F = F + FIdentity;
	gsl_matrix * Sgsl = gsl_matrix_alloc (3, 3);
	gsl_matrix * V = gsl_matrix_alloc (3, 3);
	gsl_matrix * R = gsl_matrix_alloc (3, 3);
	gsl_vector * Sig = gsl_vector_alloc (3);
	gsl_vector * workspace = gsl_vector_alloc (3);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			gsl_matrix_set (Sgsl, i, j, F(i,j));
		}
	}
	//Singular Value Decomposition of covariance matrix S
	int a  =  gsl_linalg_SV_decomp (Sgsl, V, Sig, workspace);
	boost::numeric::ublas::matrix<double>Vublas(3,3);
	boost::numeric::ublas::matrix<double>UT(3,3);
	boost::numeric::ublas::matrix<double>CurrentRotMat(3,3);
	CurrentRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
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
	//cout<<"det: "<<det<<endl;
	//displayMatrix(CurrentRotMat,"CurrentRotMat");
	if (det<0){
		cout<<"Error! Flipped element, Id: "<<Id<<endl;
		//writing a matrix that is identity, except for the last diagonal element, which is det
		boost::numeric::ublas::matrix<double> weights(3,3);
		weights = boost::numeric::ublas::identity_matrix<double>(3,3);
		boost::numeric::ublas::matrix<double> tempMat(3,3);
		tempMat = boost::numeric::ublas::zero_matrix<double>(3,3);
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
	double threshold = 0.017; //this is sine 1 degrees
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
void 	ShapeBase::calculatedudEdXde(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde){
	if (ShapeType == 1 || ShapeType == 2 || ShapeType == 3){
		dudEdXde3D(RefNormalised, dude,dXde);
	}
	else if (ShapeType == 4){
		dudEdXde2D(RefNormalised, dude,dXde);
	}
}

void 	ShapeBase::dudEdXde3D(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde){
	const int nMult = nNodes*nDim;
	boost::numeric::ublas::vector<double> displacement(nMult);
	boost::numeric::ublas::vector<double> refpos(nMult);

	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]- RefNormalised[i][j];
			refpos(counter) =  RefNormalised[i][j];
			counter++;
		}
	}
	boost::numeric::ublas::axpy_prod(Bo,displacement,dude);
	boost::numeric::ublas::axpy_prod(Bo,refpos,dXde);
}

void 	ShapeBase::dudEdXde2D(double** RefNormalised, boost::numeric::ublas::vector<double>& dude, boost::numeric::ublas::vector<double>& dXde){
	int dim =nDim-1; //calculating  for 2D
	const int nMult = nNodes*dim;
	boost::numeric::ublas::vector<double> displacement(nMult);
	boost::numeric::ublas::vector<double> refpos(nMult);

	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<dim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]- RefNormalised[i][j];
			refpos(counter) =  RefNormalised[i][j];
			counter++;
		}
	}
	boost::numeric::ublas::vector<double>dude2D(dim*dim);
	boost::numeric::ublas::vector<double>dXde2D(dim*dim);
	dXde2D = boost::numeric::ublas::zero_vector<double>(dim*dim);
	dude2D = boost::numeric::ublas::zero_vector<double>(dim*dim);
	//displayMatrix(Bo,"Bo");
	//displayMatrix(displacement,"displacement");
	//displayMatrix(refpos,"refpos");
	boost::numeric::ublas::axpy_prod(Bo,displacement,dude2D);
	boost::numeric::ublas::axpy_prod(Bo,refpos,dXde2D);
	//displayMatrix(dXde2D,"dXde2D");
	//now transferring the information on 2D to a 3D matrix
	dude(0) = dude2D(0); dude(1) = dude2D(1);
	dude(3) = dude2D(2); dude(4) = dude2D(3);
	dude(8) = 1.0;
	dXde(0) = dXde2D(0); dXde(1) = dXde2D(1);
	dXde(3) = dXde2D(2); dXde(4) = dXde2D(3);
	dXde(8) = 1.0;
}
/*
bool 	ShapeBase::calculateAlignmentRotationMatrix(double** RefNormalised, double* rotMat){
	//Calculating the covarience matrix:
	boost::numeric::ublas::matrix<double>S(3,3);
	S = boost::numeric::ublas::zero_matrix<double>(3,3);
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
	CurrentRotMat = boost::numeric::ublas::zero_matrix<double>(3,3);
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
		tempMat = boost::numeric::ublas::zero_matrix<double>(3,3);
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
*/


bool 	ShapeBase::checkPackingToThisNodeViaState(int ColumnarLayerDiscretisationLayers, Node* NodePointer){
	if(IsAblated){
		//if the element is ablated, do not pack against it
		return false;
	}
	if(ColumnarLayerDiscretisationLayers>1){
		//If the columnar layer is discretised into multiple layers, the apical elements should be checked against apical nodes,
		// and basal nodes should be checked against basal elements. The midline elements should not have packing, BUT on  a single layer tissue, all is midline, therefore
		// this check would not be valid.
		if (tissuePlacement == 2){	//tissue placement of the element is midline in a multi-layered columnar layer, it should not pack to anything
			return false;
		}
		if (NodePointer->tissuePlacement == 1){
			//node is apical, should pack to apical elements of the columnar layer only - and all of the peripodial membrane
			if (tissueType == 0 && tissuePlacement == 1){ //tissue type of the element is columnar, tissue placement is basal
				return false;
			}
		}
		else if (NodePointer->tissuePlacement == 0){
			//node is basal, should pack to apical elements of the columnar layer only - and all of the peripodial membrane
			//BUT, all peripodial membrane nodes are put in as basal nodes, the node itself can be on the peripodial membrane,
			//in which case, it should pack to the element regardless (midline elemetns are already eliminated above)
			if (NodePointer->tissueType == 0){  //tissue type of the node is columnar
				if (tissueType == 0 && tissuePlacement == 0){ //tissue type of the element is columnar, tissue placement is apical
					return false;
				}
			}
		}
	}
	//The node and element are positioned correctly to be able to pack, then does the element belong to the node?
	bool pointBelongsToElement = DoesPointBelogToMe(NodePointer->Id);
	if (pointBelongsToElement){
		return false;
	}
	return true;
}

bool 	ShapeBase::DoesPointBelogToMe(int IdNode){
	for (int i = 0; i<nNodes; ++i){
		if (NodeIds[i] == IdNode){
			return true;
		}
	}
	return false;
}

double 	ShapeBase::determinant3by3Matrix(double* rotMat){
	double det =0.0;
	det  =  rotMat[0]*(rotMat[4]*rotMat[8]-rotMat[5]*rotMat[7]);
	det -= rotMat[1]*(rotMat[3]*rotMat[8]-rotMat[5]*rotMat[6]);
	det += rotMat[2]*(rotMat[3]*rotMat[7]-rotMat[4]*rotMat[6]);
	return det;
}
double 	ShapeBase::determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det =0.0;
	det  =  Mat(0,0)*(Mat(1,1)*Mat(2,2)-Mat(1,2)*Mat(2,1));
	det -= Mat(0,1)*(Mat(1,0)*Mat(2,2)-Mat(1,2)*Mat(2,0));
	det += Mat(0,2)*(Mat(1,0)*Mat(2,1)-Mat(1,1)*Mat(2,0));
	return det;
}
double 	ShapeBase::determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det = Mat(0,0) * Mat(1,1) - Mat(0,1) * Mat(1,0);
	return det;
}

void 	ShapeBase::updatePositionsAlignedToReferenceForRK(){
	//cout<<"Id: "<<Id<<" updating positions aligned to reference"<<endl;
	updatePositionsAlignedToReferenceWithBuffers();
	const int dim = nDim;
	double* refCentre = new double[dim];
	bringPositionAlignedToReferenceToOrigin(refCentre);
	for (int i=0; i<nNodes; ++i){
		PositionsAlignedToReference[i][0] += refCentre[0];
		PositionsAlignedToReference[i][1] += refCentre[1];
		PositionsAlignedToReference[i][2] += refCentre[2];
	}
	//cout<<"finalised update"<<endl;
}

void 	ShapeBase::alignElementOnReference(){
	//if (tissueType == 1){
	//	cout<<" Element : "<<Id<<endl;
	//	displayMatrix(WorldToReferenceRotMat,"WorldToReferenceRotMat_BeforeAlignment");
	//}
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
		//bool calculateRotation = calculateAlignmentRotationMatrix(RefNormalised, rotMat);
		/*if (Id == 12 ) {
			cout<<"rotMat from node alignment: "<<endl;
			for (int i=0; i<9 ;++i){
				cout<<rotMat[i]<<" ";
			}
			cout<<endl;
		}*/
		//SV decomposition alignment will work nicely on 3D elements, BUT
		//a 2D element will not be able correct for rotations in z plane.
		//Now I will manually correct for alignment in z plane for 2D elements, then move on the SV decomposition:
		if (ShapeType == 4){
			correctFor2DAlignment();
			//if (tissueType == 1){
			//	cout<<" Element : "<<Id<<endl;
			//	displayMatrix(WorldToReferenceRotMat,"WorldToReferenceRotMat_After2DAlignment");
			//}
		}
		//Now continuing on SV decomposition
		bool calculateRotation = calculateDisplacementGradientRotationMatrix(RefNormalised, rotMat);
		/*if (Id == 12 &&  calculateRotation ) {
			cout<<"rotMat from strain decomposition: "<<endl;
			for (int i=0; i<9 ;++i){
				cout<<rotMat[i]<<" ";
			}
			cout<<endl;
		}*/
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
			tmpMat = boost::numeric::ublas::zero_matrix<double>(3,3);
			boost::numeric::ublas::axpy_prod(CurrentRotMat,WorldToReferenceRotMat, tmpMat);
			WorldToReferenceRotMat = tmpMat;
			//if (tissueType == 1){
			//	cout<<" Element : "<<Id<<endl;
			//	displayMatrix(WorldToReferenceRotMat,"WorldToReferenceRotMat_after3DAlignment");
			//}
		}
		else{
			//alignement seems necessary, yet the rotation matrix was identity
			for (int i=0; i<nNodes; ++i){
				PositionsAlignedToReference[i][0] += refCentre[0];
				PositionsAlignedToReference[i][1] += refCentre[1];
				PositionsAlignedToReference[i][2] += refCentre[2];
			}
		}
	}
	else{
		//there have been no rotation, I want to centre the positions on reference now
		for (int i=0; i<nNodes; ++i){
			PositionsAlignedToReference[i][0] += refCentre[0];
			PositionsAlignedToReference[i][1] += refCentre[1];
			PositionsAlignedToReference[i][2] += refCentre[2];
		}
	}
	//if (tissueType == 1){
	//	cout<<" Element : "<<Id<<endl;
	//	displayMatrix(WorldToReferenceRotMat,"WorldToReferenceRotMat_afterAllAlignment");
	//}
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

void	ShapeBase::calculateRotationAxis(double* u, double* v,double* rotAx, double c){
	//aligning u onto v:
	if (c>-0.9998){
		crossProduct3D(u,v,rotAx);
		normaliseVector3D(rotAx);
	}
	else{
		//the angle is 180 degree, the standard rotation axis calculation will be wrong, I am rotating over x axis at all times;
		rotAx[0]= 1;rotAx[1]= 0;rotAx[2]= 0;
	}
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

void  ShapeBase::rotateReferenceElementByRotationMatrix(double* rotMat){
	//cout<<"rotating the reference matrix of element: "<<Id<<endl;
	for (int i=0; i<nNodes; ++i){
		double * u;
		u = new double[3];
		for (int j=0; j<nDim; ++j){
			u[j] = ReferenceShape->Positions[i][j];
		}
		rotateVectorByRotationMatrix(u,rotMat);
		for (int j=0; j<nDim; ++j){
			ReferenceShape->Positions[i][j] = u[j];
		}
	}
}

void	ShapeBase::calculateForces(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile){
	if (ShapeDim == 3){		//3D element
		calculateForces3D(RKId, SystemForces, Nodes, outputFile);
	}
	else if (ShapeDim == 2){	//2D element
		calculateForces2D(RKId, SystemForces, Nodes, outputFile);
	}
}

void	ShapeBase::calculateForces3D(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile){
	//cout<<"calculating forces"<<endl;
	const int nMult = nNodes*nDim;
	using namespace boost::numeric::ublas;
	boost::numeric::ublas::vector<double> displacement(nMult);
	if (RKId != 0 ){
		//If we are at the first RK step, the positions aligned to reference
		//will be up to date from alignment calculation.
		//In any other RK step (inside this if clause), the positions are updated.
		//I do not need to re-calculate alignment, the change will be negligable, but I still need to update the positions
		//using the same rotation matrices.
		//outputFile<<"  id: "<<Id<<" Updating positions aligned to reference"<<endl;
		updatePositionsAlignedToReferenceForRK();
	}
	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]-ReferenceShape->Positions[i][j];
			counter++;
		}
	}
	//outputFile<<"  id: "<<Id<<"   calculating strain"<<endl;
	Strain = zero_vector<double>(6);
	boost::numeric::ublas::axpy_prod(B,displacement,Strain);
	if (RKId == 0){
		//I need to keep record of this for saving and displaying purposes. This is the strain at the beginning of the time step.
		//This strain should be recorded as the strain of the current step.
		//Otherwise, the recorded strains will belong to artificial setup of RK step 4.
		//Necessary updates are done when saving or displaying is applicable
		RK1Strain = Strain;
	}
	//reading from the upper triangular:
	PlasticStrain(0)= LocalGrowthStrainsMat(0,0);
	PlasticStrain(1)= LocalGrowthStrainsMat(1,1);
	PlasticStrain(2)= LocalGrowthStrainsMat(2,2);
	PlasticStrain(3)= LocalGrowthStrainsMat(0,1);
	PlasticStrain(4)= LocalGrowthStrainsMat(0,2);
	PlasticStrain(5)= LocalGrowthStrainsMat(1,2);
	boost::numeric::ublas::vector<double> NetStrain;
	NetStrain= zero_vector<double>(6);
	NetStrain = Strain - PlasticStrain;
	Forces = zero_vector<double>(nMult);
	//outputFile<<"  id: "<<Id<<"   calculating forces"<<endl;
	boost::numeric::ublas::axpy_prod(BE,NetStrain,Forces);
	//Now I have the forces in tissue coordinate system, I need the forces in world coordinates:
	boost::numeric::ublas::matrix<double>forcesInReferenceCoordsMat(nDim,nNodes);
	forcesInReferenceCoordsMat = zero_matrix<double>(nDim,nNodes);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			forcesInReferenceCoordsMat(j,i)= Forces(counter);
			counter++;
		}
	}
	boost::numeric::ublas::matrix<double>forcesInWorldT(nDim,nNodes);
	forcesInWorldT = zero_matrix<double>(nDim,nNodes);
	boost::numeric::ublas::axpy_prod(trans(WorldToReferenceRotMat),forcesInReferenceCoordsMat,forcesInWorldT);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			Forces(counter) = forcesInWorldT(j,i);
			counter++;
		}
	}
	//cout<<"Element Id: "<<Id<<"Forces on node 351: "<<SystemForces[RKId][351][0]<<" "<<SystemForces[RKId][351][1]<<" "<<SystemForces[RKId][351][2]<<endl;
	//Now put the forces in world coordinates into system forces, in forces per volume format
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			if (!Nodes[NodeIds[i]]->FixedPos[j]){
				//SystemForces[RKId][NodeIds[i]][j] = SystemForces[RKId][NodeIds[i]][j] - (Forces(counter) / ReferenceShape->Volume);
				SystemForces[RKId][NodeIds[i]][j] = SystemForces[RKId][NodeIds[i]][j] - Forces(counter);
			}
			/*else{
				SystemForces[RKId][NodeIds[0]][j] = SystemForces[RKId][NodeIds[0]][j] + (Forces(counter) / ReferenceShape->Volume)/4;
				SystemForces[RKId][NodeIds[1]][j] = SystemForces[RKId][NodeIds[1]][j] + (Forces(counter) / ReferenceShape->Volume)/4;
				SystemForces[RKId][NodeIds[3]][j] = SystemForces[RKId][NodeIds[3]][j] + (Forces(counter) / ReferenceShape->Volume)/4;
				SystemForces[RKId][NodeIds[4]][j] = SystemForces[RKId][NodeIds[4]][j] + (Forces(counter) / ReferenceShape->Volume)/4;
			}*/
			counter++;
		}
	}
}

void	ShapeBase::updatePositionALignedToReferenceForDrawing(){
	alignElementOnReference();
}

void	ShapeBase::calculateForces2D(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile){
	//cout<<"calculating forces for 2D"<<endl;
	int dim = nDim-1;	//calculating forces and strains for 2D
	const int nMult = nNodes*dim;
	using namespace boost::numeric::ublas;
	boost::numeric::ublas::vector<double> displacement(nMult);
	//cout<<"updating positions aligned to ref"<<endl;
	if (RKId != 0 ){
		//If we are at the first RK step, the positions aligned to reference
		//will be up to date from alignment calculation.
		//In any other RK step (inside this if clause), the positions are updated.
		//I do not need to re-calculate alignment, the change will be negligable, but I still need to update the positions
		//using the same rotation matrices.
		//outputFile<<"  id: "<<Id<<" Updating positions aligned to reference"<<endl;
		//updatePositionsAlignedToReferenceForRK();
		//No correct alignemtn every turn:
		alignElementOnReference();
	}
	//cout<<"finalised positions aligned to ref"<<endl;
	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<dim; ++j){
			displacement(counter) = PositionsAlignedToReference[i][j]-ReferenceShape->Positions[i][j];
			counter++;
		}
	}
	//outputFile<<"  id: "<<Id<<"   calculating strain"<<endl;
	//cout<<"calculated displacement - RKId: "<<RKId<<" elementId : "<<Id<<endl;
	//displayMatrix(B,"B");
	//displayMatrix(displacement,"displacement");
	//outputFile<<"  id: "<<Id<<"   calculating strain"<<endl;
	Strain = zero_vector<double>(6);
	boost::numeric::ublas::vector<double> Strain2D(3);
	Strain2D = zero_vector<double>(3);
	boost::numeric::ublas::axpy_prod(B,displacement,Strain2D);
	//cout<<"calculated strain2D"<<endl;
	Strain(0)=Strain2D(0);
	Strain(1)=Strain2D(1);
	Strain(3)=Strain2D(2);
	if (RKId == 0){
		//I need to keep record of this for saving and displaying purposes. This is the strain at the beginning of the time step.
		//This strain should be recorded as the strain of the current step.
		//Otherwise, the recorded strains will belong to artificial setup of RK step 4.
		//Necessary updates are done when saving or displaying is applicable
		RK1Strain = Strain;
	}
	//cout<<"updated RK1Strain"<<endl;
	//reading from the upper triangular:
	PlasticStrain(0)= LocalGrowthStrainsMat(0,0);
	PlasticStrain(1)= LocalGrowthStrainsMat(1,1);
	PlasticStrain(2)= LocalGrowthStrainsMat(2,2);
	PlasticStrain(3)= LocalGrowthStrainsMat(0,1);
	PlasticStrain(4)= LocalGrowthStrainsMat(0,2);
	PlasticStrain(5)= LocalGrowthStrainsMat(1,2);
	//cout<<"updated PlasticStrain"<<endl;
	boost::numeric::ublas::vector<double> NetStrain;
	NetStrain= zero_vector<double>(3);
	NetStrain(0) = Strain(0) - PlasticStrain(0); //ex
	NetStrain(1) = Strain(1) - PlasticStrain(1); //ey
	NetStrain(2) = Strain(3) - PlasticStrain(3); //gxy  -- skipping z terms
	/*if (RKId == 0){
		cout<<"RKID: "<<RKId<<" Element: "<<Id<<" Local Strains: "<<Strain[0]<<" "<<Strain[1]<<" "<<Strain[2]<<" "<<Strain[3]<<" "<<Strain[4]<<" "<<Strain[5]<<endl;
		cout<<"RKID: "<<RKId<<" Element: "<<Id<<" Plastic Strains: "<<PlasticStrain[0]<<" "<<PlasticStrain[1]<<" "<<PlasticStrain[2]<<" "<<PlasticStrain[3]<<" "<<PlasticStrain[4]<<" "<<PlasticStrain[5]<<endl;
		cout<<"RKID: "<<RKId<<" Element: "<<Id<<" Net Strains: "<<NetStrain[0]<<" "<<NetStrain[1]<<" "<<NetStrain[2]<<endl;
	}*/
	//cout<<"calculated  NetStrain"<<endl;
	/*if(Id == 1166 || Id == 1164 || Id == 1162 || Id == 1168 || Id == 1170 || Id == 1095 || Id == 1093 || Id == 1091 || Id == 1097 || Id == 1099){
		cout<<"RK: "<<RKId<<" Element: "<<Id<<" Positions: "<<endl;
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				cout<<" 	"<<Positions[i][j]<<" ";
			}
			cout<<endl;
		}
		cout<<"RK: "<<RKId<<" Element: "<<Id<<" PositionsAlignedToReference: "<<endl;
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<dim; ++j){
				cout<<" 	"<<PositionsAlignedToReference[i][j]<<" ";
			}
			cout<<endl;
		}
		cout<<"RK: "<<RKId<<" Element: "<<Id<<" ReferenceShape->Positions: "<<endl;
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				cout<<" 	"<<ReferenceShape->Positions[i][j]<<" ";
			}
			cout<<endl;
		}
		double* v1;
		v1 = new double[3];
		double* v2;
		v2 = new double[3];
		double* vcross;
		vcross = new double[3];
		for (int i=0;i<3;i++){
			v1[i] =  PositionsAlignedToReference[1][i] -  PositionsAlignedToReference[0][i];
			v2[i] =  PositionsAlignedToReference[2][i] -  PositionsAlignedToReference[0][i];
		}
		crossProduct3D(v1,v2,vcross);
		cout<<"positions aligned to reference normal: "<<vcross[0]<<" "<<vcross[1]<<"  "<<vcross[2]<<endl;
		cout<<"Element "<<Id <<" ";
		displayMatrix(Strain2D,"Strain2D");
		cout<<"Element "<<Id <<" ";
		displayMatrix(Strain,"Strain");
		cout<<"Element "<<Id <<" ";
		displayMatrix(B,"B");
	}*/
	Forces = zero_vector<double>(nNodes*nDim);
	boost::numeric::ublas::vector<double> Forces2D;
	Forces2D = zero_vector<double>(nMult);
	//outputFile<<"  id: "<<Id<<"   calculating forces"<<endl;
	boost::numeric::ublas::axpy_prod(BE,NetStrain,Forces2D);
	//cout<<"calculated  Forces2D"<<endl;
	//boost::numeric::ublas::vector<double> ForcesFromk2D = zero_vector<double>(nMult);
	//boost::numeric::ublas::axpy_prod(k,displacement,ForcesFromk2D);
	//displayMatrix(Forces2D,"Forces2D______");
	//displayMatrix(ForcesFromk2D,"ForcesFromk2D_");
	//Now I have the forces in tissue coordinate system, I need the forces in world coordinates:
	boost::numeric::ublas::matrix<double>forcesInReferenceCoordsMat(nDim,nNodes);
	forcesInReferenceCoordsMat = zero_matrix<double>(nDim,nNodes);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<dim; ++j){
			forcesInReferenceCoordsMat(j,i)= Forces2D(counter);
			counter++;
		}
	}
	//cout<<"assigned forcesInReferenceCoordsMat"<<endl;
	boost::numeric::ublas::matrix<double>forcesInWorldT(nDim,nNodes);
	forcesInWorldT = zero_matrix<double>(nDim,nNodes);
	boost::numeric::ublas::axpy_prod(trans(WorldToReferenceRotMat),forcesInReferenceCoordsMat,forcesInWorldT);
	//cout<<"calculated forcesInReferenceCoordsMat"<<endl;
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			Forces(counter) = forcesInWorldT(j,i);
			counter++;
		}
	}
	//cout<<"calculated Forces"<<endl;
	//cout<<"Element Id: "<<Id<<"Forces on node 351: "<<SystemForces[RKId][351][0]<<" "<<SystemForces[RKId][351][1]<<" "<<SystemForces[RKId][351][2]<<endl;
	//Now put the forces in world coordinates into system forces, in forces per volume format
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			if (!Nodes[NodeIds[i]]->FixedPos[j]){
				SystemForces[RKId][NodeIds[i]][j] = SystemForces[RKId][NodeIds[i]][j] - Forces(counter);

			}
			counter++;
		}
	}
	//cout<<"updated SystemForces"<<endl;
	//displayMatrix(Strain,"Strain");
	//displayMatrix(PlasticStrain,"PlasticStrain");
	//displayMatrix(NetStrain,"NetStrain");
	//displayMatrix(Forces2D,"Forces2D");
	//displayMatrix(forcesInReferenceCoordsMat,"forcesInReferenceCoordsMat");
	//displayMatrix(forcesInWorldT,"forcesInWorldT");
	//displayMatrix(Forces,"Forces");
	//cout<<"element id: "<<Id<<" RK: "<<RKId<<"system Forces: "<<endl;
	//for (int i = 0; i<Nodes.size(); ++i){
	//	for (int j = 0; j<3; ++j){
	//		cout<<SystemForces[RKId][i][j]<<" ";
	//	}
	//	cout<<endl;
	//}
}

void	ShapeBase::fillNodeNeighbourhood(vector<Node*>& Nodes){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nNodes; ++j){
			if ( i !=j ){
				int n = Nodes[NodeIds[i]]->immediateNeigs.size();
				bool alreadyOnList = false;
				for (int k=0; k<n; ++k){
					if (NodeIds[j] == Nodes[NodeIds[i]]->immediateNeigs[k]){
						alreadyOnList = true;
						break;
					}
				}
				if (!alreadyOnList){
					Nodes[NodeIds[i]]->immediateNeigs.push_back(NodeIds[j]);
				}
			}
		}
	}
}


void	ShapeBase::updatePositions(int RKId, vector<Node*>& Nodes){
	if (RKId < 3){
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				Positions[i][j] = Nodes[NodeIds[i]]->RKPosition[j];
			}
		}
	}
	else{
		for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
			}
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
bool 	ShapeBase::InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse/*, double& det*/){
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

	/*det = 1.0;
	for(unsigned int i = 0; i < A.size1(); i++) {
		det *= A(i,i); // multiply by elements on diagonal
	    det = det * determinant_sign( pm );
	}*/
	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

/*int 	ShapeBase::determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}*/

void	ShapeBase::crossProduct3D(double* u, double* v, double* cross){
	cross[0] = u[1]*v[2] - u[2]*v[1];
	cross[1] = u[2]*v[0] - u[0]*v[2];
	cross[2] = u[0]*v[1] - u[1]*v[0];
}

double	ShapeBase::calculateMagnitudeVector3D(double* v){
	double mag = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	mag = pow(mag,0.5);
	return mag;
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

void 	ShapeBase:: assignVolumesToNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
		Nodes[NodeIds[i]]->mass +=ReferenceShape->Volume/nNodes;
	}
}

void 	ShapeBase:: assignElementToConnectedNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
		Nodes[NodeIds[i]]->connectedElementIds.push_back(Id);
		double weightfFraction = (ReferenceShape->Volume/nNodes)/Nodes[NodeIds[i]]->mass;
		Nodes[NodeIds[i]]->connectedElementWeights.push_back(weightfFraction);
	}
}

void 	ShapeBase:: removeMassFromNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
			Nodes[NodeIds[i]]->mass -=ReferenceShape->Volume/nNodes;
			//updating the weight fractions of the elements on the node due to elimination of the ablated element:
			int n = Nodes[NodeIds[i]]->connectedElementIds.size();
			double scaler = 1.0;
			for (int j=0;j<n;++j){
				if (Nodes[NodeIds[i]]->connectedElementIds[j]==Id){
					scaler = 1.0 - Nodes[NodeIds[i]]->connectedElementWeights[j];
					Nodes[NodeIds[i]]->connectedElementWeights[j]  = 0.0;
					break;
				}
			}
			for (int j=0;j<n;++j){
				Nodes[NodeIds[i]]->connectedElementWeights[j] /= scaler;
			}
			//All wiights are normlised as the sum will make 1.0. Now I do not want this element in the weighing,
			//it does not have a mass anymore, therefore I will multiply all the remaining weights with (1-w_ablated);
		}
}
/*
void ShapeBase::updateElementsNodePositions(int RKId, double ***SystemForces, vector <Node*>& Nodes, double dt){
	//Update Node positions:
	int n = nNodes;
	if (RKId < 3){
		//the first 3 RK steps, the velocity will be calculated, and the positions will be updated from normal positions to RKPositions, with half dt:
		double multiplier=0.0;
		if (RKId<2){
			multiplier =0.5;
		}
		else{
			multiplier =1.0;
		}
		for (int i=0;i<n;++i){
			for (int j=0; j<Nodes[NodeIds[i]]->nDim; ++j){
				Nodes[NodeIds[i]]->Velocity[RKId][j] = SystemForces[RKId][NodeIds[i]][j]/ Nodes[NodeIds[i]]->Viscosity ;
				Nodes[NodeIds[i]]->RKPosition[j] = Nodes[NodeIds[i]]->Position[j] + Nodes[NodeIds[i]]->Velocity[RKId][j]*multiplier*dt;
			}
		}
	}
	else{
		//this is the last RK step, I need to update the velocity only with RK, then I need to calculate the final positions
		//from 4 RK velocities:
		for (int i=0;i<n;++i){
		//	cout<<"Nodes "<<i<<" velocity: ";
			for (int j=0; j<Nodes[i]->nDim; ++j){
				Nodes[NodeIds[i]]->Velocity[RKId][j] = SystemForces[RKId][NodeIds[i]][j]/Nodes[NodeIds[i]]->Viscosity;
				//now I have 4 velocity data (corresponding to Runge-Kutta  k1, k2, k3, and k4)
				//writing  the velocity into v[0]
				//cout<<Nodes[i]->Velocity[0][j]<<" "<<Nodes[i]->Velocity[1][j]<<" "<<Nodes[i]->Velocity[2][j]<<" "<<Nodes[i]->Velocity[3][j]<<" ";
				Nodes[NodeIds[i]]->Velocity[0][j] = 1.0/6.0 * (Nodes[NodeIds[i]]->Velocity[0][j] + 2.0 * (Nodes[NodeIds[i]]->Velocity[1][j] + Nodes[NodeIds[i]]->Velocity[2][j]) + Nodes[NodeIds[i]]->Velocity[3][j]);
				Nodes[NodeIds[i]]->Position[j] += Nodes[NodeIds[i]]->Velocity[0][j]*dt;
			}
		//	cout<<endl;
		}
	}
	//Now I need to clear these forces, as they are already applied, no cumulative effect!
	for (int i=0;i<n;++i){
		for (int j=0; j<Nodes[i]->nDim; ++j){
			SystemForces[RKId][NodeIds[i]][j]=0.0;
		}
	}
};*/

