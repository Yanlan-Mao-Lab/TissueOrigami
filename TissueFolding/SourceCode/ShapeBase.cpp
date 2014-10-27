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
	//cout<<"calling calculate GrowthInLocalCoordinates"<<endl;
	calculateGrowthInLocalCoordinates(CurrGrowthStrainAddition);
	/*
	// The growth scale is the growth of this element should feel,
	// on the coordinates of the tissue, on Dorsal-Ventral, Anterior-Posterior, and Apical-Basal axes.
	// I need to get the orientation of this tissue coordinate system first.
	// This coordinate system, and all associated rotation matrices wiill change only when the
	// reference prism rotates, and it could have been updated in shape change functions up till this point.
	// There fore, I will re-calculate only if necessary
	if (!TissueCoordinateSystemUpToDate){
		calculateTissueCoordinateSystem();
		//Now I have the tissue coordinate system, I need to calculate the rotation matrix to rotate
		//the coordinate system of the reference prism to tissue coordinates
		calculateRotationMatrixReferenceToTissue();
	}
	//Now I have the rotation matrix and its transpose, I can rotate the strains:
	if (!CurrGrowthStrainsUpToDate){
		calculateCurrGrowthStrainsInTissueCoords();
	}
	updateGrowthStrainInTissueCoords(CurrGrowthStrainAddition);
	*/
}

void 	ShapeBase::changeShape(double shapechangescale, int axis){
	IsChangingShape = true;
	ChangedShapeInThePast = true;
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
	// I need to calculate the shape change in the coordinate system of the tissue,
	// This is why I need to have the rotation matrices between the
	// tissue coordinate system, local coordinate system and world coordinate system.
	// I will calculate these coordinate systems and corresponding rotation matrices]
	// ONLY if they are not up to date.
	// They will not be up-to date if the rotation if the reference element has changed
	// since last calculation.
	if (!TissueCoordinateSystemUpToDate){
		// I need to get the orientation of the tissue coordinate system first.
		calculateTissueCoordinateSystem();
		// Now I have the tissue coordinate system, I need to calculate the rotation matrix to rotate
		// the coordinate system of the reference prism to tissue coordinate system
		calculateRotationMatrixReferenceToTissue();
	}
	//Now I have the matrices necessary to calculate the shape changes, I need to rotate the shaoe change
	// that is already on the prism, to the tissue coordinate system. Since shape change function will be called
	// as many times on the element as there are active shape change inputs on its vicinity, I do not need to do this
	// rotation over and over again. It is different for growth changes, as I do all the update at once, and I do not
	// need such book-keeping.
	// CurrShapeChangeToAdd and the flag CurrShapeChangeStrainsUpToDate are reset at the beginnig of each time step, in function
	// resetCurrStepShapeChangeData
	if (!CurrShapeChangeStrainsUpToDate){
		calculateCurrShapeChangeStrainsInTissueCoords();
	}
	float effect =  (1.0 + CurrShapeChangeStrainsInTissueCoordsMat(axis,axis)) * shapechangescale;
	float compansation = -0.5 * effect;
	CurrShapeChangeToAdd[axis] += effect;
	CurrShapeChangeToAdd[compansatingaxes[0]] += compansation;
	CurrShapeChangeToAdd[compansatingaxes[1]] += compansation;
}

void 	ShapeBase::calculatePlasticStrain(){
	//If the element is growing or changing shape, I need to do that update first,
	//This also means the current strains on tissue coordinates will be up-to-date
	if (IsGrowing || IsChangingShape){
		if (IsGrowing){
			growShape();
		}
		if (IsChangingShape){
			//Here I need to update shape change
			updateShapeChangeStrainInTissueCoords();
		}
	}
	else {
		//The tissue is not changing its plastic deformation this step, but it may have rotate, and the
		// strains may not be correct any more
		if ( ( !CurrGrowthStrainsUpToDate && GrewInThePast ) || ( !CurrShapeChangeStrainsUpToDate && ChangedShapeInThePast ) ){
			if (!TissueCoordinateSystemUpToDate){
				calculateTissueCoordinateSystem();
				//Now I have the tissue coordinate system, I need to calculate the rotation matrix to rotate
				//the coordinate system of the reference prism to tissue coordinates
				calculateRotationMatrixReferenceToTissue();
			}
			if (!CurrGrowthStrainsUpToDate){
				calculateCurrGrowthStrainsInTissueCoords();
			}
			if (!CurrShapeChangeStrainsUpToDate){
				calculateCurrShapeChangeStrainsInTissueCoords();
			}
		}
	}
	//now sum up all the plastic strains
	updatePlasticStrainInTissueCoords();
	//Now I have the growth strain in growth axes, I need to convert it to world coordinates:
	//cout<<"entering updateGrowthStainInWorldAxes"<<endl;
	updatePlasticStainInWorldCoords();
	//cout<<"entering updateLocalPlasticStrain"<<endl;
	if (IsGrowing || IsChangingShape ){
		//The plastic strains on local coordinate system will only be affected if there is a plastic change on the element,
		//rotational changes will not affect the local coordinates, the local is on the reference shape frame in the first place.
		updateLocalPlasticStrains();
	}
	//cout<<"finished"<<endl;
}

void 	ShapeBase::calculateCurrGrowthStrainsInTissueCoords(){
/*	using namespace boost::numeric::ublas;
	matrix<double> tmpMat1(3,3);
	axpy_prod(RefToTissueRotMat,LocalGrowthStrainsMat,tmpMat1);
	axpy_prod(tmpMat1,RefToTissueRotMatT,CurrGrowthStrainsInTissueCoordsMat);
	CurrGrowthStrainsUpToDate = true;
*/
}

void 	ShapeBase::calculateCurrShapeChangeStrainsInTissueCoords(){
	using namespace boost::numeric::ublas;
	matrix<double> tmpMat1(3,3);
	axpy_prod(RefToTissueRotMat,LocalShapeChangeStrainsMat,tmpMat1);
	axpy_prod(tmpMat1,RefToTissueRotMatT,CurrShapeChangeStrainsInTissueCoordsMat);
	CurrShapeChangeStrainsUpToDate = true;
}

void 	ShapeBase::updateGrowthStrainInTissueCoords(double* strainsToAdd){
	//Now I have the current strains defining the size of the element in gorwth axis coordinates
	//I can grow the system, with respect to its current size now (10% growth will be 10% of current size, cumulative interest)
	//displayMatrix(CurrPlasticStrainsInTissueCoordsMat,"CurrPlasticStrainsInTissueCoordsMat");
	//cout<<"before addition"<<endl;
	//displayMatrix(CurrGrowthStrainsInTissueCoordsMat,"CurrGrowthStrainsInTissueCoordsMat");
	CurrGrowthStrainsInTissueCoordsMat(0,0) = ( (1.0 + CurrGrowthStrainsInTissueCoordsMat(0,0)) * (1.0 + strainsToAdd[0]) ) - 1.0;
	CurrGrowthStrainsInTissueCoordsMat(1,1) = ( (1.0 + CurrGrowthStrainsInTissueCoordsMat(1,1)) * (1.0 + strainsToAdd[1]) ) - 1.0;
	CurrGrowthStrainsInTissueCoordsMat(2,2) = ( (1.0 + CurrGrowthStrainsInTissueCoordsMat(2,2)) * (1.0 + strainsToAdd[2]) ) - 1.0;
	//cout<<"after addition"<<endl;
	//displayMatrix(CurrGrowthStrainsInTissueCoordsMat,"CurrGrowthStrainsInTissueCoordsMat");
}

void 	ShapeBase::updateShapeChangeStrainInTissueCoords(){
	//Now I have the current strains defining the size of the element in gorwth axis coordinates
	//I can grow the system, with respect to its current size now (10% growth will be 10% of current size, cumulative interest)
	//displayMatrix(CurrPlasticStrainsInTissueCoordsMat,"CurrPlasticStrainsInTissueCoordsMat");
	if (!CurrShapeChangeStrainsUpToDate){
		calculateCurrShapeChangeStrainsInTissueCoords();
	}
	CurrShapeChangeStrainsInTissueCoordsMat(0,0) += CurrShapeChangeToAdd[0];
	CurrShapeChangeStrainsInTissueCoordsMat(1,1) += CurrShapeChangeToAdd[1];
	CurrShapeChangeStrainsInTissueCoordsMat(2,2) += CurrShapeChangeToAdd[2];
	//displayMatrix(CurrShapeChangeStrainsInTissueCoordsMat,"CurrShapeChangeStrainsInTissueCoordsMat");
}

void 	ShapeBase::updatePlasticStrainInTissueCoords(){
	//This should be new growth strain in growth axes, not plastic strain to add
	CurrPlasticStrainsInTissueCoordsMat = CurrGrowthStrainsInTissueCoordsMat + CurrShapeChangeStrainsInTissueCoordsMat;
}

void 	ShapeBase::updateLocalPlasticStrains(){
/*	boost::numeric::ublas::matrix<double> tmpMat1(3,3);
	boost::numeric::ublas::matrix<double> tmpMat2(3,3);
	//	if I am in this function, at least on of IsGrowing or IsChangingShape is true, then I need to update
	//	total plastic strain regardless of details.
	axpy_prod(RefToTissueRotMatT,CurrPlasticStrainsInTissueCoordsMat,tmpMat1);
	axpy_prod(tmpMat1,RefToTissueRotMat,LocalPlasticStrainsMat);
	if (IsGrowing){
		axpy_prod(RefToTissueRotMatT,CurrGrowthStrainsInTissueCoordsMat,tmpMat1);
		axpy_prod(tmpMat1,RefToTissueRotMat,LocalGrowthStrainsMat);
	}
	if (IsChangingShape){
		axpy_prod(RefToTissueRotMatT,CurrShapeChangeStrainsInTissueCoordsMat,tmpMat1);
		axpy_prod(tmpMat1,RefToTissueRotMat,LocalShapeChangeStrainsMat);
	}
	//displayMatrix(LocalGrowthStrainsMat,"LocalGrowthStrainsMat");
*/
}

void 	ShapeBase::updateTissueCoordStrain(){
	boost::numeric::ublas::matrix<double> TissueToWorldRotMatT(3,3);
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
	TissueToWorldRotMatT = trans(TissueToWorldRotMat);
	axpy_prod(TissueToWorldRotMatT,StrainMat,tmpMat1);
	axpy_prod(tmpMat1,TissueToWorldRotMat,StrainTissueMat);

}

void 	ShapeBase::updatePlasticStainInWorldCoords(){
	boost::numeric::ublas::matrix<double> tmpMat1(3,3);
	boost::numeric::ublas::matrix<double> tmpMat2(3,3);
	axpy_prod(TissueToWorldRotMat,CurrPlasticStrainsInTissueCoordsMat,tmpMat1);
	boost::numeric::ublas::matrix<double> TissueToWorldRotMatT(3,3);
	TissueToWorldRotMatT = trans(TissueToWorldRotMat);
	axpy_prod(tmpMat1,TissueToWorldRotMatT,tmpMat2);
	PlasticStrain(0)= tmpMat2(0,0);
	PlasticStrain(1)= tmpMat2(1,1);
	PlasticStrain(2)= tmpMat2(2,2);
	PlasticStrain(3)= tmpMat2(1,0);
	PlasticStrain(4)= tmpMat2(1,2);
	PlasticStrain(5)= tmpMat2(2,0);
	//displayMatrix(PlasticStrain,"PlasticStrain");
}

void 	ShapeBase::calculateGrowthInLocalCoordinates(double * strainsToAdd){
	//get the (+)ve x of the reference aligned onto current coordinate system (use transpose of rotation matrix you already have for alignment)
	boost::numeric::ublas::matrix<double>  ReferenceCoordSystem;
	boost::numeric::ublas::matrix<double>  tmpMat1;
	ReferenceCoordSystem = boost::numeric::ublas::identity_matrix<double>(3,3);
	tmpMat1 = boost::numeric::ublas::zero_matrix<double>(3,3);
	//tmp1 will be reference coordinate system on current shape
	axpy_prod(trans(WorldToReferenceRotMat),ReferenceCoordSystem,tmpMat1);
	//if(Id == 0){
	//	displayMatrix(tmpMat1,"ReferenceCoordSystem - on current shape");
	//}
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
	calculateRotationAngleSinCos(u,v,c,s);
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
		//if(Id == 0){
		//	displayMatrix(WorldToTissueRotMat,"WorldToTissueRotMat");
		//}
		//rotate the (+)ve x &y of the reference aligned to current, with the rotation matrix R world to tissue
		axpy_prod(WorldToTissueRotMat,tmpMat1,ReferenceCoordSystem);
	}
	else{
		ReferenceCoordSystem = tmpMat1;
	}
	boost::numeric::ublas::matrix<double>  TissueCoordSystem;
	TissueCoordSystem = boost::numeric::ublas::identity_matrix<double>(3,3);
	axpy_prod(WorldToTissueRotMat,TissueCoordSystem,tmpMat1);
	TissueCoordSystem = tmpMat1;
	TissueCoordinateSystem[0] = TissueCoordSystem(0,0);
	TissueCoordinateSystem[1] = TissueCoordSystem(1,0);
	TissueCoordinateSystem[2] = TissueCoordSystem(2,0);
	TissueCoordinateSystem[3] = TissueCoordSystem(0,1);
	TissueCoordinateSystem[4] = TissueCoordSystem(1,1);
	TissueCoordinateSystem[5] = TissueCoordSystem(2,1);
	TissueCoordinateSystem[6] = TissueCoordSystem(0,2);
	TissueCoordinateSystem[7] = TissueCoordSystem(1,2);
	TissueCoordinateSystem[8] = TissueCoordSystem(2,2);

	//if(Id == 0){
	//	displayMatrix(ReferenceCoordSystem,"ReferenceCoordSystem - on tissue coordinates");
	//}
	//Now this rotated vectors x, y, and (+)ve z defines the coordinate system you want the growth strains in. Get the rotation matrix, and calculate strains.
	v[0] = 1;
	v[1] = 0;
	v[2] = 0;
	u[0]=ReferenceCoordSystem(0,0);
	u[1]=ReferenceCoordSystem(1,0);
	u[2]=ReferenceCoordSystem(2,0);
	calculateRotationAngleSinCos(u,v,c,s);
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
		//if(Id == 0){
		//	displayMatrix(CurrGrowthToAddTissue,"CurrGrowthToAddTissue");
		//}
	}
	else{
		CurrLocalGrowthToAdd(0,0)= strainsToAdd[0];
		CurrLocalGrowthToAdd(1,1)= strainsToAdd[1];
		CurrLocalGrowthToAdd(2,2)= strainsToAdd[2];
	}
	LocalGrowthStrainsMat(0,0) = ( (1.0 + LocalGrowthStrainsMat(0,0)) * (1.0 + CurrLocalGrowthToAdd(0,0)) ) - 1.0;
	LocalGrowthStrainsMat(1,1) = ( (1.0 + LocalGrowthStrainsMat(1,1)) * (1.0 + CurrLocalGrowthToAdd(1,1)) ) - 1.0;
	LocalGrowthStrainsMat(2,2) = ( (1.0 + LocalGrowthStrainsMat(2,2)) * (1.0 + CurrLocalGrowthToAdd(2,2)) ) - 1.0;
	//if(Id == 0){
	//	displayMatrix(LocalGrowthStrainsMat,"LocalGrowthStrainsMat");
	//}


	/*boost::numeric::ublas::matrix<double>  CurrGrowthToAddTissue;
	boost::numeric::ublas::matrix<double>  CurrLocalGrowthToAdd;
	boost::numeric::ublas::matrix<double>  tmpMat1;
	boost::numeric::ublas::matrix<double>  R;
	CurrGrowthToAddTissue = boost::numeric::ublas::zero_matrix<double>(3,3);
	CurrLocalGrowthToAdd= boost::numeric::ublas::zero_matrix<double>(3,3);
	tmpMat1 =  boost::numeric::ublas::zero_matrix<double>(3,3);
	R =  boost::numeric::ublas::zero_matrix<double>(3,3);
	cout<<"constructing CurrGrowthToAddTissue"<<endl;
	CurrGrowthToAddTissue(0,0)= strainsToAdd[0];
	CurrGrowthToAddTissue(1,1)= strainsToAdd[1];
	CurrGrowthToAddTissue(2,2)= strainsToAdd[2];

	boost::numeric::ublas::axpy_prod(trans(WorldToReferenceNormalRotMat),trans(WorldToReferenceXaxisRotMat),tmpMat1);
	boost::numeric::ublas::axpy_prod(tmpMat1,WorldToReferenceNormalRotMat,R);

	boost::numeric::ublas::axpy_prod(R,CurrGrowthToAddTissue,tmpMat1);
	boost::numeric::ublas::axpy_prod(tmpMat1,R,CurrLocalGrowthToAdd);
	LocalGrowthStrainsMat(0,0) = ( (1.0 + LocalGrowthStrainsMat(0,0)) * (1.0 + CurrLocalGrowthToAdd(0,0)) ) - 1.0;
	LocalGrowthStrainsMat(1,1) = ( (1.0 + LocalGrowthStrainsMat(1,1)) * (1.0 + CurrLocalGrowthToAdd(1,1)) ) - 1.0;
	LocalGrowthStrainsMat(2,2) = ( (1.0 + LocalGrowthStrainsMat(2,2)) * (1.0 + CurrLocalGrowthToAdd(2,2)) ) - 1.0;
	displayMatrix(LocalGrowthStrainsMat,"LocalGrowthStrainsMat");
	displayMatrix(CurrGrowthToAddTissue,"CurrGrowthToAddTissue");*/
}

void 	ShapeBase::calculateTissueCoordinateSystem(){
/*	//I assume the orientation of the tissue stays as the dorsal-ventral axis is aligned with x = {1,0,0} and
	//anterior-posterior axis is aligned with y = {0,1,0}. I need
	//to align the z axis with the current apicobasal axis.
	//The current apical-basal axis is the z axis of the reference element.
	double* v = new double[3];
	//calculateZVecForTissueCoordAlignment(v);
	double* u = new double[3];
	u[0]=0.0;
	u[1]=0.0;
	u[2]=1.0;

	double c, s;
	calculateRotationAngleSinCos(u,v,c,s);
	double *rotAx;
	rotAx = new double[3];
	double *rotMat;
	rotMat = new double[9]; //matrix is written in one row
	calculateRotationAxis(u,v,rotAx);
	constructRotationMatrix(c,s,rotAx,rotMat);
	updateRotationMatrixTissueToWorld(rotMat);
	rotateWorldsCoordinatesByRotationMatrix(TissueCoordinateSystem,rotMat);
	//cout<<"TissueCoordinateSystem:"<<endl;
	//cout<<TissueCoordinateSystem[0]<<" 	"<<TissueCoordinateSystem[1]<<" 	"<<TissueCoordinateSystem[2]<<endl;
	//cout<<TissueCoordinateSystem[3]<<" 	"<<TissueCoordinateSystem[4]<<" 	"<<TissueCoordinateSystem[5]<<endl;
	//cout<<TissueCoordinateSystem[6]<<" 	"<<TissueCoordinateSystem[7]<<" 	"<<TissueCoordinateSystem[8]<<endl;
	delete[] u;
	delete[] v;
	delete[] rotAx;
	delete[] rotMat;
	*/
}
void 	ShapeBase::calculateWorldToTissueRotationMatrix(){
	cout<<"inside calculateWorldToTissueRotationMatrix"<<endl;
	//I assume the orientation of the tissue stays as the dorsal-ventral axis is aligned with x = {1,0,0} and
	//anterior-posterior axis is aligned with y = {0,1,0}. I need
	//to align the z axis with the current apicobasal axis.
	//The current apical-basal axis is the z axis of the reference element.
	double* v = new double[3];
	calculateZVecForTissueCoordAlignment(v);
	double* u = new double[3];
	u[0]=0.0;
	u[1]=0.0;
	u[2]=1.0;

	double c, s;
	calculateRotationAngleSinCos(u,v,c,s);
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
	cout<<"finalised calculateWorldToTissueRotationMatrix"<<endl;
}


void 	ShapeBase::calculatePositionsOnTissueCoordinateSystem(){
	cout<<"inside calc. PositionsOnTissueCoordinateSystem"<<endl;
	calculateWorldToTissueRotationMatrix();
	boost::numeric::ublas::matrix<double> pos(nDim,nNodes);
	boost::numeric::ublas::matrix<double> tissuePosT(nDim,nNodes);
	for (int i=0; i<nNodes; i++){
		for (int j=0; j<nDim; j++){
			pos(j,i) = Positions[i][j];
		}
	}
	boost::numeric::ublas::axpy_prod(WorldToTissueRotMat,pos,tissuePosT);
	for (int i=0; i<nNodes; i++){
		for (int j=0; j<nDim; j++){
			PositionsInTissueCoord[i][j] = tissuePosT(j,i);
		}
	}
	cout<<"finalised calc. PositionsOnTissueCoordinateSystem"<<endl;
};

void 	ShapeBase::updateRotationMatrixTissueToWorld(double* rotMat){
		//this is the transpose of the rotation matrix I am using to rotate the world to tissue axes
		TissueToWorldRotMat(0,0)=rotMat[0];
		TissueToWorldRotMat(1,0)=rotMat[1];
		TissueToWorldRotMat(2,0)=rotMat[2];
		TissueToWorldRotMat(0,1)=rotMat[3];
		TissueToWorldRotMat(1,1)=rotMat[4];
		TissueToWorldRotMat(2,1)=rotMat[5];
		TissueToWorldRotMat(0,2)=rotMat[6];
		TissueToWorldRotMat(1,2)=rotMat[7];
		TissueToWorldRotMat(2,2)=rotMat[8];
}

void 	ShapeBase::calculateRotationMatrixReferenceToTissue(){
	//Now I need the rotation matrix to rotate the strains on growth coordinate system,
	//to strains in reference shape coordinate system;
	//The z axes are already aligned.
	//I will define the reference shape coordinate system as x-axis as the axis from node 0 to 1 for a prism
	//This selection is arbitrary, as long as it is consistent, I will be fine.
	double* u = new double[3];
	calculateXVecForTissueCoordAlignment(u);
	double* v = new double[3];
	v[0]=TissueCoordinateSystem[0];
	v[1]=TissueCoordinateSystem[1];
	v[2]=TissueCoordinateSystem[2];
	double c, s;
	calculateRotationAngleSinCos(u,v,c,s);
	double *rotAx;
	rotAx = new double[3];
	double *rotMat;
	rotMat = new double[9]; //matrix is written in one row
	calculateRotationAxis(u,v,rotAx);
	constructRotationMatrix(c,s,rotAx,rotMat);
	RefToTissueRotMat(0,0)=rotMat[0];
	RefToTissueRotMat(0,1)=rotMat[1];
	RefToTissueRotMat(0,2)=rotMat[2];
	RefToTissueRotMat(1,0)=rotMat[3];
	RefToTissueRotMat(1,1)=rotMat[4];
	RefToTissueRotMat(1,2)=rotMat[5];
	RefToTissueRotMat(2,0)=rotMat[6];
	RefToTissueRotMat(2,1)=rotMat[7];
	RefToTissueRotMat(2,2)=rotMat[8];
	RefToTissueRotMatT = trans(RefToTissueRotMat);
	delete[] u;
	delete[] v;
	delete[] rotAx;
	delete[] rotMat;

	//ReferenceShape->CurrAlignmentSide;
}

void 	ShapeBase::rotateWorldsCoordinatesByRotationMatrix(double* NewCoordinates, double *rotMat){
	double* u;
	u = new double[3];
	u[0]=1.0;
	u[1]=0.0;
	u[2]=0.0;
	rotateVectorByRotationMatrix(u,rotMat);
	NewCoordinates[0] = u[0];
	NewCoordinates[1] = u[1];
	NewCoordinates[2] = u[2];
	u[0]=0.0;
	u[1]=1.0;
	u[2]=0.0;
	rotateVectorByRotationMatrix(u,rotMat);
	NewCoordinates[3] = u[0];
	NewCoordinates[4] = u[1];
	NewCoordinates[5] = u[2];
	u[0]=0.0;
	u[1]=0.0;
	u[2]=1.0;
	rotateVectorByRotationMatrix(u,rotMat);
	NewCoordinates[6] = u[0];
	NewCoordinates[7] = u[1];
	NewCoordinates[8] = u[2];
	delete[] u;
}
void 	ShapeBase::updatePositionsAlignedToReferenceWithBuffers(){
	//cout<<"starting updatePositionsAlignedToReferenceWithBuffers"<<endl;
	boost::numeric::ublas::matrix<double> PositionsMat(nDim,nNodes);
	boost::numeric::ublas::matrix<double> tmpMat1(nDim,nNodes);
	boost::numeric::ublas::matrix<double> PositionsOnReferenceMat(nDim,nNodes);
	//cout<<"declared the matrices"<<endl;
	for (int i=0;i<nNodes; ++i){
		for (int j=0;j<nDim; ++j){
			PositionsMat(j,i)=Positions[i][j];
		}
	}
	//cout<<"Wrote positions on a matrix"<<endl;
	axpy_prod(WorldToReferenceRotMat,PositionsMat,PositionsOnReferenceMat);
	//cout<<"Calculated rotated positions matrix"<<endl;
	for (int i=0;i<nNodes; ++i){
		for (int j=0;j<nDim; ++j){
			PositionsAlignedToReference[i][j] = PositionsOnReferenceMat(j,i);
		}
	}
	//cout<<"finished updatePositionsAlignedToReferenceWithBuffers"<<endl;
}

void 	ShapeBase::normaliseShapePositions(double** RefNormalised, double* refCentre){
	/*cout<<"PositionsAlignedToReference before normalisation"<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			cout<<PositionsAlignedToReference[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"Reference before normalisation"<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			cout<<ReferenceShape ->Positions[i][j]<<" ";
		}
		cout<<endl;
	}*/
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
	normaliseShapePositions(RefNormalised,refCentre);
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

void 	ShapeBase::alignReference(){
	//cout<<"aligning reference:"<<endl;

	updateAlignmentTurn();
	updateReferenceShapeBaseFromBuffer();
	calculateNormals();
	double *v,*u;
	v = new double[3];
	u = new double[3];
	for (int i=0;i<nDim; ++i){
		v[i] = CurrentNormal[i];
		u[i] = ReferenceShape->CurrentNormal[i];
	}
	double c, s;
	calculateRotationAngleSinCos(u,v,c,s);
	double *rotAx;
	rotAx = new double[3];
	double *rotMat;
	rotMat = new double[9]; //matrix is written in one row
	if (c<0.9998){ //only rotate if the vectors are more than 1 degree apart
		updatedReference = true;
		TissueCoordinateSystemUpToDate = false;
		CurrShapeChangeStrainsUpToDate = false;
		CurrGrowthStrainsUpToDate = false;
		//cout<<"Aligning the plane, angle big enough - cos: "<<c<<" sin: "<<s<<endl;
		calculateRotationAxis(u,v,rotAx);
		constructRotationMatrix(c,s,rotAx,rotMat);
		//double centre[3];
		for (int i=0; i<nNodes; ++i){
			u[0] = ReferenceShape->Positions[i][0];
			u[1] = ReferenceShape->Positions[i][1];
			u[2] = ReferenceShape->Positions[i][2];
			rotateVectorByRotationMatrix(u,rotMat);
			ReferenceShape->Positions[i][0] = u[0];
			ReferenceShape->Positions[i][1] = u[1];
			ReferenceShape->Positions[i][2] = u[2];
		}
		//Aligned the planes, now align the sides:
	}
	//else{ cout<<"skipping plane alignment, degree is too small: cosine: "<<c<<" sine: "<<s<<endl;}
	double* RefSide = new double[3];
	double* ShapeSide;
	ShapeSide = new double[3];

	getCurrentAlignmentSides(RefSide, ShapeSide);
	calculateRotationAngleSinCos(RefSide,ShapeSide,c,s);
	if (c<0.9998){  //only rotate if the vectors are more than 1 degree apart
		updatedReference = true;
		TissueCoordinateSystemUpToDate = false;
		CurrShapeChangeStrainsUpToDate = false;
		CurrGrowthStrainsUpToDate = false;
		//cout<<"Aligning the side, angle big enough - cos: "<<c<<" sin: "<<s<<endl;
		calculateRotationAxis(RefSide,ShapeSide,rotAx);
		constructRotationMatrix(c,s,rotAx,rotMat);
		for (int i=0; i<nNodes; ++i){
			u[0] = ReferenceShape->Positions[i][0];
			u[1] = ReferenceShape->Positions[i][1];
			u[2] = ReferenceShape->Positions[i][2];
			rotateVectorByRotationMatrix(u,rotMat);
			ReferenceShape->Positions[i][0] = u[0];
			ReferenceShape->Positions[i][1] = u[1];
			ReferenceShape->Positions[i][2] = u[2];
		}
		//Now they are aligned, but I do not know if they face the same direction, so I need to check that:
		/*
		bool facecorrection = areSidesFacingSameDirection(RefSide, ShapeSide);
		if (!facecorrection){
			//face is in the wrong direction, I need to rotate the element a further 180 degrees
			constructRotationMatrix(-1.0,0.0,rotAx,rotMat);
			for (int i=0; i<nNodes; ++i){
				u[0] = ReferenceShape->Positions[i][0];
				u[1] = ReferenceShape->Positions[i][1];
				u[2] = ReferenceShape->Positions[i][2];
				rotateVectorByRotationMatrix(u,rotMat);
				ReferenceShape->Positions[i][0] = u[0];
				ReferenceShape->Positions[i][1] = u[1];
				ReferenceShape->Positions[i][2] = u[2];
			}
			cout<<"Reference positions after face correction: "<<endl;
			for (int i = 0; i<nNodes; ++i){
				for (int j = 0; j<nDim; ++j){
					cout<<ReferenceShape->Positions[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;

		}
		*/
	}
	//else{ cout<<"skipping side alignment, degree is too small: cosine: "<<c<<" sine: "<<s<<endl;}
	//aligning the position to the curr shape centre:
	double translocation[3] = {0.0,0.0,0.0};
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			translocation[j] +=  Positions[i][j] - ReferenceShape ->Positions[i][j];
		}
	}
	for (int j = 0; j<nDim; ++j){
		translocation[j] /= nNodes;
	}
	//translocation[0] -= 5.0;
	for (int i = 0; i<nNodes; ++i){
			for (int j = 0; j<nDim; ++j){
				ReferenceShape->Positions[i][j] +=translocation[j];
			}
	}
	delete[] u;
	delete[] v;
	delete[] rotAx;
	delete[] rotMat;
	delete[] RefSide;
	delete[] ShapeSide;
	//cout<<"finalised aligning reference, update reference? "<<updatedReference<<endl;
}

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

void	ShapeBase::calculateRotationAngleSinCos(double* u, double* v, double& c, double& s){
	//aligning u onto v:
	c = dotProduct3D(u,v);
	//if (Id <2){cout<<"u: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<" v: "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;}
	if (c > 1.0){
		//if (Id <2){cout<<"c was above 1: "<<c<<endl;}
		c = 1.0;
		s = 0.0;

	}
	else if( c<-1.0){
		//if (Id <2){cout<<"c was below -1: "<<c<<endl;}
		c = -1.0;
		s = 0.0;
	}
	else{
		double tet = acos(c);
		s = sin(tet);
		//if (Id <2){cout<<"c is ok: "<<c<<" tet "<<tet<<" s "<<s<<endl;}
	}
	//double s2c2 = c*c + s*s;
}

void	ShapeBase::calculateRotationAxis(double* u, double* v,double* rotAx){
	//aligning u onto v:
	crossProduct3D(u,v,rotAx);
	normaliseVector3D(rotAx);
}

void	ShapeBase::constructRotationMatrix(double c, double s, double* rotAx, double* rotMat){
	//cout<<" Rotation axis: "<<rotAx[0]<<" "<<rotAx[1]<<" "<<rotAx[2]<<endl;

	rotMat[0] = c + rotAx[0]*rotAx[0]*(1 - c);
	rotMat[1] = rotAx[0]*rotAx[1]*(1 - c) - rotAx[2]*s;
	rotMat[2] = rotAx[0]*rotAx[2]*(1 - c) + rotAx[1]*s;

	rotMat[3] = rotAx[1]*rotAx[0]*(1 - c) + rotAx[2]*s;
	rotMat[4] = c + rotAx[1]*rotAx[1]*(1 - c);
	rotMat[5] = rotAx[1]*rotAx[2]*(1 - c) - rotAx[0]*s;

	rotMat[6] = rotAx[2]*rotAx[0]*(1 - c) - rotAx[1]*s;
	rotMat[7] = rotAx[2]*rotAx[1]*(1 - c) + rotAx[0]*s;
	rotMat[8] = c + rotAx[2]*rotAx[2]*(1 - c);

	/*cout<<"Rotation Matrix: "<<endl;
	cout<<rotMat[0]<<" "<<rotMat[1]<<" "<<rotMat[2]<<endl;
	cout<<rotMat[3]<<" "<<rotMat[4]<<" "<<rotMat[5]<<endl;
	cout<<rotMat[6]<<" "<<rotMat[7]<<" "<<rotMat[8]<<endl;*/
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,double* rotMat){
	double x = rotMat[0]*u[0]+rotMat[1]*u[1]+rotMat[2]*u[2];
	double y = rotMat[3]*u[0]+rotMat[4]*u[1]+rotMat[5]*u[2];
	double z = rotMat[6]*u[0]+rotMat[7]*u[1]+rotMat[8]*u[2];
	//cout<<"before rotation: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
	u[0] = x;
	u[1] = y;
	u[2] = z;
	//cout<<"after rotation: "<<u[0]<<" "<<u[1]<<" "<<u[2]<<endl;
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
	//displayMatrix(k,"k");
	//displayMatrix(B,"B");
	//displayMatrix(BE,"BE");
	//cout<<"multiplying to get the contribution of the forces"<<endl;
	//displayMatrix(GrowthStrain,"GrowthStrain");
	//displayMatrix(ShapeChangeStrain,"ShapeChangeStrain");
	//displayMatrix(Strain,"Strain");
	//displayMatrix(LocalGrowthStrainsMat,"LocalGrowthStrainsMat");
	PlasticStrain(0)= LocalGrowthStrainsMat(0,0);
	PlasticStrain(1)= LocalGrowthStrainsMat(1,1);
	PlasticStrain(2)= LocalGrowthStrainsMat(2,2);
	PlasticStrain(3)= LocalGrowthStrainsMat(1,0);
	PlasticStrain(4)= LocalGrowthStrainsMat(1,2);
	PlasticStrain(5)= LocalGrowthStrainsMat(2,0);
	//displayMatrix(PlasticStrain,"PlasticStrain");
	boost::numeric::ublas::axpy_prod(BE,PlasticStrain,PlasticStrainForces);

	//cout<<"getting the net forces"<<endl;
	//cout<<"total forces:"<<endl;
	//displayMatrix(Forces,"Forces");
	Forces = Forces - PlasticStrainForces;
	//cout<<"forces after balancing growth contribution: "<<endl;
	/*if(Id < 2){
		//displayMatrix(PlasticStrain,"PlasticStrain");
		displayMatrix(displacement,"displacement");
		displayMatrix(Strain,"Strain");
		displayMatrix(Forces,"Forces");
		//displayMatrix(PlasticStrainForces,"PlasticStrainForces");
	}*/
	//Now I have the forces in tissue coordinate system, I need the forces in world coordinates:
	boost::numeric::ublas::matrix<double>forcesInReferenceCoordsMat(nDim,nNodes);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			forcesInReferenceCoordsMat(j,i)= Forces(counter);
			counter++;
		}
	}
	//boost::numeric::ublas::matrix<double>tempMat(nDim,nNodes);
	boost::numeric::ublas::matrix<double>forcesInWorldT(nDim,nNodes);
	//cout<<"calculating conversion"<<endl;
	//boost::numeric::ublas::axpy_prod(trans(WorldToReferenceSideRotMat),forcesInReferenceCoordsMat,tempMat);
	//boost::numeric::ublas::axpy_prod(trans(WorldToReferenceNormalRotMat),tempMat,forcesInWorldT);
	boost::numeric::ublas::axpy_prod(trans(WorldToReferenceRotMat),forcesInReferenceCoordsMat,forcesInWorldT);
	//cout<<"finalised conversion, writing on forces"<<endl;
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			Forces(counter) = forcesInWorldT(j,i);
			counter++;
		}
	}
	/*if(Id < 2){
		displayMatrix(WorldToReferenceRotMat,"WorldToReferenceRotMat");
		displayMatrix(Forces,"Forces in world coordinates");
	}*/
	//cout<<"finalised writing on forces"<<endl;
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
	/*
	cout<<"Element: "<<Id<<endl;
	cout<<"positions: "<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			cout<<Positions[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<"Reference positions: "<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			cout<<ReferenceShape->Positions[i][j]<<" ";
		}
		cout<<endl;
	}
	displayMatrix(displacement,"displacement");
	displayMatrix(Forces,"Forces");
	cout<<"SystemForces:"<<endl;
	int nSysNodes = Nodes.size();
	for (int i=0;i<nSysNodes;++i){
		for (int j=0;j<3;++j){
			cout<<SystemForces[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;*/
	//cout<<"finalised calculating forces"<<endl;
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
