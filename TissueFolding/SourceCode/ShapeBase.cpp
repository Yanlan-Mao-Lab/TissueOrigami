#include "ShapeBase.h"
#include "Node.h"
#include <sstream>

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>



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

void ShapeBase::getPos(gsl_matrix* Pos){
    for (int i=0; i<nNodes; ++i){
        for (int j =0; j<nDim; ++j){
            gsl_matrix_set (Pos, i, j, Positions[i][j]);
        }
    }
}

double 	ShapeBase::getYoungModulus(){
	return E;
}

double 	ShapeBase::getPoissonRatio(){
	return v;
}

double* ShapeBase::getGrowthRate(){
	//cout<<"Element "<<Id<<" Growth rate: "<<GrowthRate[0]<<" "<<GrowthRate[1]<<" "<<GrowthRate[2]<<endl;
	return GrowthRate;
}

gsl_matrix* ShapeBase::getFg(){
    gsl_matrix* tmpFg =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpFg,Fg);
    return tmpFg;
}

void ShapeBase::createMatrixCopy(gsl_matrix* dest, gsl_matrix* src){
    int m = src->size1;
    int n = src->size2;
    gsl_matrix_set_zero(dest);
    double tmp= 0.0;
    for (int i=0; i<m; ++i){
        for (int j=0 ; j<n; ++j){
            tmp = gsl_matrix_get(src,i,j);
            gsl_matrix_set(dest,i,j,tmp);
        }
    }
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

double ShapeBase::getPeripodialness(){
	return peripodialGrowthWeight;
}

double ShapeBase::getColumnarness(){
	return columnarGrowthWeight;
}

void ShapeBase::getRelativePosInColumnarBoundingBox(double* relativePos){
	relativePos[0] =  columnarRelativePosInBoundingBox[0];
	relativePos[1] =  columnarRelativePosInBoundingBox[1];
}

void ShapeBase::getRelativePosInPeripodialBoundingBox(double* relativePos){
	relativePos[0] =  peripodialRelativePosInBoundingBox[0];
	relativePos[1] =  peripodialRelativePosInBoundingBox[1];
}

void ShapeBase::convertRelativePosToGridIndex(double* relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY){
	//cout<<"relpos: "<<relpos[0]<<" "<<relpos[1]<<endl;
	relpos[0] *= (float) (nGridX-1);
	relpos[1] *= (float) (nGridY-1);
	indexX = floor(relpos[0]);
	fracX  = relpos[0] - indexX;
	indexY = floor(relpos[1]);
	fracY  = relpos[1] - indexY;
	//cout<<" indexX "<<indexX<<" fracX "<<fracX<<" indexY "<<indexY<<" fracY "<<fracY<<endl;
	if (indexX >= nGridX-1) { //this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexX = nGridX-2;
		fracX = 1.0;
		//cout<<" in if 1, indexX: "<<indexX<<" fracX: "<<fracX<<endl;
	}else if (indexX<0){
		indexX = 0;
		fracX = 0.0;
		//cout<<" in if 2, indexX: "<<indexX<<" fracX: "<<fracX<<endl;
	}
	if (indexY >= nGridY-1) {//this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexY = nGridY-2;
		fracY = 1.0;
		//cout<<" in if 3, indexY: "<<indexY<<" fracY: "<<fracY<<endl;
	}else if (indexY<0){
		indexY = 0;
		fracY = 0.0;
		//cout<<" in if 4, indexY: "<<indexY<<" fracY: "<<fracY<<endl;
	}
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
	for (int i = 0; i<nNodes; ++i){
		Positions[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
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
	bool hasPeripodialNode = false;
	bool hasLinkerNode = false;
	for (int i = 0; i<nNodes; ++i){
		if (Nodes[NodeIds[i]]->tissueType == 0){
			hasColumnarNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissueType == 1){
			hasPeripodialNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissueType == 2){
			hasLinkerNode = true;
		}
	}
	if (hasLinkerNode){
		tissueType = 2;
	}
	else if (hasPeripodialNode){
		//ASK LINKER ZONE BEFORE THIS, SOME LINKER ELEMENTS CAN HAVE LINKER NODES, NO COLUMNAR ELEMENT OR PERIPODIAL ELEMENT SHOULD HAVE A LINKER NODE
		tissueType = 1;
		setGrowthWeightsViaTissuePlacement( 1.0);//default is set to be columnar, I will not set this for linkers, as they are set in the initiation of peripodial membrane
	}
	else if (hasColumnarNode){
		//ASK PERIPODIAL MEMBRANE BEFORE THIS, SOME PERIPODIAL ELEMENTS CAN HAVE COLUMNAR NODES, AND SOME LINKER ELEMENTS CAN HAVE COLUMNAR NODES. NO COLUMNAR ELEMENT SHOULD HAVE A PERIPODIAL NODE
		tissueType = 0;
	}
	else {
		cerr<<"Element is not placed into tissue correctly, Id: "<<Id<<endl;
	}
	//cout<<"Element : "<<Id<<" hasColumnarNode: "<<hasColumnarNode<<" hasPeripodialmNode "<<hasPeripodialNode<<" tissueType: "<<tissueType<<endl;
}

void 	ShapeBase::setGrowthWeightsViaTissuePlacement (double periWeight){
	peripodialGrowthWeight = periWeight;
	columnarGrowthWeight = 1.0 - peripodialGrowthWeight;
	//cout<<" Element: "<<Id<<" peripodialness: "<<peripodialGrowthWeight<<" columnarness: "<<columnarGrowthWeight<<endl;
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

void ShapeBase::setFg(gsl_matrix* currFg){
    gsl_matrix_memcpy (Fg, currFg);
    gsl_matrix* tmpFgForInversion =gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion, Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    gsl_matrix_free(tmpFgForInversion);
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

bool 	ShapeBase::readNodeIdData(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;
		if (NodeIds[i] != savedId){
			return false;
		}
	}
	return true;
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

bool	ShapeBase::readReferencePositionData(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			if (ReferenceShape -> Positions[i][j] != savedPos){
				//the positions are not equal, it may be an issue of rounding, my satisfactory precision is 2%
				float percentError = (ReferenceShape -> Positions[i][j] - savedPos) / ReferenceShape -> Positions[i][j]*100.0;
				if (percentError>2.0 || percentError< -2.0){
					cout<<"ReferenceShape->Positions: "<<ReferenceShape -> Positions[i][j]<<" savedPos "<<savedPos<<" percent Error: "<<percentError<<endl;
					return false;
				}
			}
		}
	}
	return true;
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
	StrainMag = 0.0;
	if (type == 0){
		//this is the average strain
        for (int i=0; i<3; ++i){
           StrainMag += gsl_matrix_get(Strain,i,0) ;
        }
		StrainMag /= 3;
	}
	else if (type == 1){
        StrainMag = gsl_matrix_get(Strain,0,0);
	}
	else if (type == 2){
        StrainMag = gsl_matrix_get(Strain,1,0);
	}
	else if (type == 3){
        StrainMag = gsl_matrix_get(Strain,2,0);
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

void 	ShapeBase::getPysProp(int type, float &PysPropMag, double dt){
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
        double timescale = 60.0*60.0/dt; //converting groth rate from per dt to per hour
		for (int i =0 ; i< nDim ; ++i){
            PysPropMag += growth[i]*timescale;
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
	for (int i=0;i<6;++i){
			CurrGrowthStrainAddition[i]  += growthscale[i];
	}
}

void 	ShapeBase::changeShapeByFsc(double dt){
    gsl_matrix* FscIncrement = gsl_matrix_calloc(nDim,nDim); ///< The increment of shape change that will be induced this step
    if (rotatedGrowth){
    	double rTemp[3] = {0.0,0.0,0.0};
		for (int i = 0; i<3; ++i){
			rTemp[i] = gsl_matrix_get(GrowthStrainsRotMat,i,0)*ShapeChangeRate[0]+gsl_matrix_get(GrowthStrainsRotMat,i,1)*ShapeChangeRate[1]+gsl_matrix_get(GrowthStrainsRotMat,i,2)*ShapeChangeRate[2];
			if ( (ShapeChangeRate[i] <0 && rTemp[i] >0.0 ) || (ShapeChangeRate[i] >0 && rTemp[i] < 0.0) ){
				rTemp[i] *= -1.0;
			}
		}
		for (int i = 0; i<3; ++i){
			ShapeChangeRate[i] = rTemp[i];
		}
	}
    for (int i=0; i<3 ;++i){
    	gsl_matrix_set(FscIncrement,i,i, exp(ShapeChangeRate[i]*dt));
    }
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FscIncrement, Fsc, 0.0, temp1);
	gsl_matrix_memcpy(Fsc, temp1);
	gsl_matrix * tmpFscForInversion = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFscForInversion,Fsc);
	bool inverted = InvertMatrix(tmpFscForInversion, InvFsc);
	if (!inverted){
		cerr<<"Fsc not inverted!!"<<endl;
	}
	//double detFsc = determinant3by3Matrix(Fsc);
	//freeing matrices allocated in this function
	gsl_matrix_free(FscIncrement);
	gsl_matrix_free(temp1);
	gsl_matrix_free(tmpFscForInversion);
}

void 	ShapeBase::growShapeByFg(double dt){
    gsl_matrix* FgIncrement = gsl_matrix_calloc(nDim,nDim);
    //I need to rotate the growth rate vector, not the growth increment matrix!
    //When rate is zero, the growth increment vector should be identity.
    //Such that the total growth gradient (Fg) will stay the same.
    //If you rotate a zero growth increment, then it is not identity anymore!
    if (rotatedGrowth){
        double rTemp[3] = {0.0,0.0,0.0};
        for (int i = 0; i<3; ++i){
            rTemp[i] = gsl_matrix_get(GrowthStrainsRotMat,i,0)*GrowthRate[0]+gsl_matrix_get(GrowthStrainsRotMat,i,1)*GrowthRate[1]+gsl_matrix_get(GrowthStrainsRotMat,i,2)*GrowthRate[2];
            if (rTemp[i] <0.0){
                rTemp[i] *= -1.0;
            }
        }
        for (int i = 0; i<3; ++i){
            GrowthRate[i] = rTemp[i];
        }
    }
    for (int i=0; i<3 ;++i){
        gsl_matrix_set(FgIncrement,i,i, exp(GrowthRate[i]*dt));
    }
    //bibap
    //cout<<"inside Fg, Id : "<<Id<<endl;
    //gsl_matrix_set(FgIncrement,0,1, M_PI/12.0);
    //gsl_matrix_set(FgIncrement,1,0, M_PI/12.0);
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FgIncrement, Fg, 0.0, temp1);
    gsl_matrix_memcpy(Fg, temp1);

    gsl_matrix * tmpFgForInversion = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion,Fg);

    //gsl_matrix * tmpFgForInversion =  gsl_matrix_alloc(nDim, nDim);
    //gsl_matrix_memcpy (tmpFgForInversion, Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
        cerr<<"Fg not inverted!!"<<endl;
    }
    double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
    VolumePerNode = GrownVolume/nNodes;
    //freeing matrices allocated in this function
    gsl_matrix_free(FgIncrement);
    gsl_matrix_free(temp1);
    gsl_matrix_free(tmpFgForInversion);
}

void 	ShapeBase::CalculateGrowthRotationByF(){
    gsl_matrix* rotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(rotMat);
    rotatedGrowth = false;
    rotatedGrowth = calculate3DRotMatFromF(rotMat);
    if (rotatedGrowth){
        rotatedGrowth = disassembleRotationMatrixForZ(rotMat);
        if (rotatedGrowth){
            gsl_matrix_transpose(rotMat);
            gsl_matrix_memcpy(GrowthStrainsRotMat,rotMat);
        }
    }
    //freeing matrices allocated in this function
    gsl_matrix_free(rotMat);
}

bool 	ShapeBase::disassembleRotationMatrixForZ(gsl_matrix* rotMat){
    //For a rotation matrix R = [r11 r12 r13
    //                          r21 r22 r23
    //                          r31 r32 r33]
    //The rotations in X, Y and Z are calculated as follows:
    // tethaX = atan2(r32,r33);
    // tethaY = atan2(-r31, sqrt(r32 * r32 + r33 * r33)
    // tethaZ = atan2(r21,r11)
    // then Rx = [ 1              0              0
    //             0              cos(tethaX)   -sin(tethaX)
    //             0              sin(tethaX)    cos(tethaX)]
    //      Ry = [ cos(tethaY)    0              sin(tethaY)
    //             0              1              0
    //            -sin(tethaY)    0              cos(tethaY)]
    //      Rz = [ cos(tethaZ)   -sin(tethaZ)    0
    //             sin(tethaZ)    cos(tethaZ)    0
    //             0              0              1]
    double tethaZ = atan2(gsl_matrix_get(rotMat,1,0),gsl_matrix_get(rotMat,0,0));
    if (tethaZ > 0.017 || tethaZ < -0.017){ //0.017 rad is ~1 degrees
        //rotation is more than 1 degrees, element incremental growth should be rotated
        double c = cos(tethaZ);
        double s = sin(tethaZ);
        gsl_matrix_set(rotMat,0,0,  c );
        gsl_matrix_set(rotMat,0,1, -1.0*s);
        gsl_matrix_set(rotMat,0,2,  0.0);
        gsl_matrix_set(rotMat,1,0,  s);
        gsl_matrix_set(rotMat,1,1,  c);
        gsl_matrix_set(rotMat,1,2,  0.0);
        gsl_matrix_set(rotMat,2,0,  0.0);
        gsl_matrix_set(rotMat,2,1,  0.0);
        gsl_matrix_set(rotMat,2,2,  1.0);
        return true;
    }
    else{
        return false;   //rotation is less than 1 degrees;
    }
}

bool 	ShapeBase::calculate3DRotMatFromF(gsl_matrix* rotMat){
    //if (Id == 0) {displayMatrix(TriPointF,"TriPointF");}
    gsl_matrix * Sgsl = gsl_matrix_alloc (3, 3);
    gsl_matrix * V = gsl_matrix_alloc (3, 3);
    gsl_matrix * R = gsl_matrix_alloc (3, 3);
    gsl_vector * Sig = gsl_vector_alloc (3);
    gsl_vector * workspace = gsl_vector_alloc (3);
    createMatrixCopy (Sgsl,TriPointF);
    //Singular Value Decomposition of covariance matrix S
    int a  =  gsl_linalg_SV_decomp (Sgsl, V, Sig, workspace);

    gsl_matrix_transpose(Sgsl); //Sgsl ended up as U, I need U^T to calculate rotation matrix as : V*U^T
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, V, Sgsl,0.0, rotMat);
    double det = determinant3by3Matrix(rotMat);
    if (det<0){
        cout<<"Error! Flipped element, Id: "<<Id<<endl;
    }
    gsl_matrix_free (Sgsl);
    gsl_matrix_free (V);
    gsl_matrix_free (R);
    gsl_vector_free (Sig);
    gsl_vector_free (workspace);
    //Now I need to check if there is only numerical error accumulationg on rotMat, or there is an actual rotation (above 1 degrees):
    double threshold = 0.017; //this is sine 1 degrees
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            if(i != j){
                if (gsl_matrix_get(rotMat,i,j)>threshold || gsl_matrix_get(rotMat,i,j)< (-1.0*threshold)) {
                    return true;
                }
            }
        }
    }
    return false; //none of the off - diagonal terms of the matrix are above the threshold, the current rotation is only niumerical error.
}

void 	ShapeBase::calculateRelativePosInBoundingBox(double columnarBoundingBoxXMin, double columnarBoundingBoxYMin, double columnarBoundingBoxLength, double columnarBoundingBoxWidth, double peripodialBoundingBoxXMin, double peripodialBoundingBoxYMin, double peripodialBoundingBoxLength, double peripodialBoundingBoxWidth){
	columnarRelativePosInBoundingBox = getCentre();
	if (tissueType != 0){ //the tissue is not columnar, so there is peripodial membrane
		peripodialRelativePosInBoundingBox[0] = columnarRelativePosInBoundingBox[0];
		peripodialRelativePosInBoundingBox[1] = columnarRelativePosInBoundingBox[1];
		peripodialRelativePosInBoundingBox[0] = (peripodialRelativePosInBoundingBox[0] - peripodialBoundingBoxXMin) / peripodialBoundingBoxLength;
		peripodialRelativePosInBoundingBox[1] = (peripodialRelativePosInBoundingBox[1] - peripodialBoundingBoxYMin) / peripodialBoundingBoxWidth;
	}
	columnarRelativePosInBoundingBox[0] = (columnarRelativePosInBoundingBox[0] - columnarBoundingBoxXMin) / columnarBoundingBoxLength;
	columnarRelativePosInBoundingBox[1] = (columnarRelativePosInBoundingBox[1] - columnarBoundingBoxYMin) / columnarBoundingBoxWidth;
	//cout<<"Element: "<<Id<<" RelPos: "<<columnarRelativePosInBoundingBox[0]<<" "<<columnarRelativePosInBoundingBox[1]<<" "<<peripodialRelativePosInBoundingBox[0]<<" "<<peripodialRelativePosInBoundingBox[1]<<endl;
	//double* a = new double[3];
	//a = getRelativePosInBoundingBox();
	//cout<<" a: "<< a[0]<<" "<<a[1]<<endl;
	//delete[] a;
}

void 	ShapeBase::displayRelativePosInBoundingBox(){
	if (tissueType == 0){
		cout<<"Element: "<<Id<<"  rel Pos from element columnar: "<<columnarRelativePosInBoundingBox[0]<<" "<<columnarRelativePosInBoundingBox[1]<<endl;
	}
	else if (tissueType == 1){
		cout<<"Element: "<<Id<<"  rel Pos from element peripodial: "<<peripodialRelativePosInBoundingBox[0]<<" "<<peripodialRelativePosInBoundingBox[1]<<endl;
	}
	else{
		cout<<"Element: "<<Id<<"  rel Pos from element columnar: "<<columnarRelativePosInBoundingBox[0]<<" "<<columnarRelativePosInBoundingBox[1]<<" peripodial: "<<peripodialRelativePosInBoundingBox[0]<<" "<<peripodialRelativePosInBoundingBox[1]<<endl;
	}
}

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

double 	ShapeBase::determinant3by3Matrix(gsl_matrix* Mat){
    double det =0.0;
    gsl_matrix_get(Mat,0,0);
    det  =  gsl_matrix_get(Mat,0,0)*(gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,1));
    det -= gsl_matrix_get(Mat,0,1)*(gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,0));
    det += gsl_matrix_get(Mat,0,2)*(gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,1)-gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,0));
    return det;
}
double 	ShapeBase::determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det = Mat(0,0) * Mat(1,1) - Mat(0,1) * Mat(1,0);
	return det;
}

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

void	ShapeBase::calculateForces(double **SystemForces, vector <Node*>& Nodes, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile){
    if (ShapeDim == 3){		//3D element
        calculateForces3D(SystemForces, Nodes, recordForcesOnFixedNodes, FixedNodeForces, outputFile);
    }
}

void	ShapeBase::calculateForces3D(double **SystemForces, vector <Node*>& Nodes,  bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile){
    int dim = nDim;
    int n = nNodes;
    //calculating F and B in a 3 point gaussian:

    gsl_matrix* TriPointg  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix_set_zero(TriPointF);
    gsl_matrix* currg = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currF = gsl_matrix_calloc(dim,dim);
    //The point order is established in shape function derivative calculation!
    //Make sure the weights fir in with the order - eta zeta nu:
    //double points[3][3]={{1.0/6.0,1.0/6.0,0.0},{2.0/3.0,1.0/6.0,0.0},{1.0/6.0,2.0/3.0,0.0}};
    double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
    for (int iter =0; iter<3;++iter){
        //cout<<"Calculating gauss point: "<<eta<<" "<<nu<<" "<<zeta<<endl;
        calculateCurrNodalForces(currg, currF, iter);
        gsl_matrix_scale(currg,weights[iter]);
        gsl_matrix_add(TriPointg, currg);
        gsl_matrix_scale(currF,weights[iter]);
        gsl_matrix_add(TriPointF, currF);
        //displayMatrix(currg,"currg");
    }
    //displayMatrix(TriPointg,"TriPointg");
    //cout<<"Element Id: "<<Id<<"Forces on node 351: "<<SystemForces[351][0]<<" "<<SystemForces[351][1]<<" "<<SystemForces[351][2]<<endl;
    //Now put the forces in world coordinates into system forces, in forces per volume format
    int counter = 0;
    for (int i = 0; i<nNodes; ++i){
        for (int j = 0; j<nDim; ++j){
            if (!Nodes[NodeIds[i]]->FixedPos[j]){
                SystemForces[NodeIds[i]][j] = SystemForces[NodeIds[i]][j] - gsl_matrix_get(TriPointg,counter,0);
            }
            else if(recordForcesOnFixedNodes){
                FixedNodeForces[NodeIds[i]][j] = FixedNodeForces[NodeIds[i]][j] - gsl_matrix_get(TriPointg,counter,0);
            }
            counter++;
        }
    }
    /*cout<<"SystemForces " <<endl;
    for (int i = 0; i<nNodes; ++i){
        for (int j = 0; j<nDim; ++j){
            cout<<SystemForces[NodeIds[i]][j]<<" ";
        }
        cout<<endl;
    }*/
    //freeing matrices allocated in this function
    gsl_matrix_free(TriPointg);
    gsl_matrix_free(currg);
    gsl_matrix_free(currF);
    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(Fg, "Fg");
}

gsl_matrix* ShapeBase::calculateEForNodalForces(gsl_matrix* Fe, gsl_matrix* FeT){
    //calculating E (E = 1/2 *(Fe^T*Fe-I):
    gsl_matrix * E =  gsl_matrix_alloc(nDim, nDim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FeT, Fe,0.0, E);
    gsl_matrix * I = gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set_identity(I);
    gsl_matrix_sub(E,I);
    gsl_matrix_scale(E, 0.5);

    /*
    gsl_matrix * E =  gsl_matrix_calloc(nDim, nDim);
    gsl_matrix_add(E,Fe);
    gsl_matrix_add(E,FeT);
    gsl_matrix_scale(E, 0.5);
    gsl_matrix * I = gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set_identity(I);
    gsl_matrix_sub(E,I);
*/
    gsl_matrix_free(I);
    return E;
}

gsl_matrix* ShapeBase::calculateSForNodalForces(gsl_matrix* E){
    //calculating S: (S = D:E)
    gsl_matrix_set_zero(Strain);
    gsl_matrix* compactS = gsl_matrix_calloc(6,1);
    gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
    gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
    gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
    gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
    gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
    gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, D, Strain,0.0, compactS);

    gsl_matrix * S =  gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set(S,0,0,gsl_matrix_get(compactS,0,0));
    gsl_matrix_set(S,1,1,gsl_matrix_get(compactS,1,0));
    gsl_matrix_set(S,2,2,gsl_matrix_get(compactS,2,0));
    gsl_matrix_set(S,1,0,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,1,2,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,2,0,gsl_matrix_get(compactS,5,0));
    gsl_matrix_set(S,0,1,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,2,1,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,0,2,gsl_matrix_get(compactS,5,0));

    //displayMatrix(D,"D");
    //displayMatrix(E,"E");
    //displayMatrix(S, "S");
    //displayMatrix(Strain,"Strain");
    gsl_matrix_free(compactS);
/*
    double Idouble[3][3] = {{1.0,0.0,0.0} , {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    double D81[3][3][3][3];
    for (int I = 0; I<3; ++I){
        for (int J = 0; J<3; ++J){
            for (int K = 0; K<3; ++K){
                for (int L = 0; L<3; ++L){
                    D81[I][J][K][L]=lambda*Idouble[K][L]*Idouble[I][J] + mu * (Idouble[I][K]*Idouble[J][L] + Idouble[I][L]*Idouble[J][K] );
                 }
            }
        }
    }

    gsl_matrix * S2 =  gsl_matrix_calloc(nDim, nDim);
    for (int I = 0; I<3; ++I){
        for (int J = 0; J<3; ++J){
            for (int K = 0; K<3; ++K){
                for (int L = 0; L<3; ++L){
                    double value = gsl_matrix_get(S2,I,J);
                    value += D81[I][J][K][L]*gsl_matrix_get(E,K,L);
                    gsl_matrix_set(S2,I,J,value);
                }
            }
        }
    }
    displayMatrix(S2,"S2");*/
    return S;
}

gsl_matrix* ShapeBase::calculateCompactStressForNodalForces(gsl_matrix* Fe, gsl_matrix* S, gsl_matrix* FeT, gsl_matrix* Stress){
    //calculating stress (stress = detFe^-1 Fe S Fe^T):
    double detFe = determinant3by3Matrix(Fe);
    gsl_matrix * tmpMat1 =  gsl_matrix_calloc(nDim, nDim);
    //cout<<"detFe: "<<detFe<<endl;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Fe, S,0.0, tmpMat1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmpMat1, FeT,0.0, Stress);
    gsl_matrix_scale(Stress, 1.0/detFe);
    gsl_matrix * compactStress =  gsl_matrix_calloc(6,1);
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(Stress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(Stress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(Stress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(Stress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(Stress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(Stress,0,2));

    gsl_matrix_free(tmpMat1);
    return compactStress;
}

gsl_matrix* ShapeBase::calculateInvJShFuncDerSWithFe(gsl_matrix* currFe, gsl_matrix* InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithFe){
	//I want InvJe, normally J InvDXde = F, I can get Je from
	// Je InvDXde = Fe
	// but I can also get InvJe directly from:
	// InvJe Je InvdXde = InvJe Fe => I InvdXde = InvJe Fe => InvdXde InvFe = InvJe I => InvJe = InvdXde InvFe
	gsl_matrix * tmpFeforInversion =  gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvFe = gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvJe = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFeforInversion,currFe);
	bool inverted = InvertMatrix(tmpFeforInversion, InvFe);
	if (!inverted){
		cerr<<"Fe not inverted!!"<<endl;
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvDXde, InvFe,0.0, InvJe);

	int dim2 = nDim*nDim;
	//Generating the inverse Jacobian(elastic) stack:
	gsl_matrix * InvJacobianElasticStack =  gsl_matrix_calloc(dim2,dim2);
	for (int i =0; i<nDim; i++){
		for (int m=0; m<nDim; ++m){
			for (int n=0; n<3; ++n){
				gsl_matrix_set(InvJacobianElasticStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJe,n,m));
			}
		}
	}

	//I am calculating this for k calculation, in case there is growth. Under conditions that there is no growth, this function is not necessary,
	//the values of invJShFuncDerSWithF and  invJShFuncDerS will be equal
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianElasticStack, ShapeFuncDerStack,0.0, invJShFuncDerSWithFe);
	gsl_matrix_free(tmpFeforInversion);
	gsl_matrix_free(InvFe);
	gsl_matrix_free(InvJe);
	gsl_matrix_free(InvJacobianElasticStack);
}

gsl_matrix* ShapeBase::calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix* B, gsl_matrix *invJShFuncDerS){
    //calculating the transpose of B:
    //gsl_matrix* tmpMat2 = gsl_matrix_calloc(6,nDim*nDim);
    //calculating B:
    // gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, CoeffMat, InvJacobianStack,0.0, tmpMat2);
    // gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmpMat2, ShapeFuncDerStack,0.0, B);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianStack, ShapeFuncDerStack,0.0, invJShFuncDerS);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, CoeffMat, invJShFuncDerS,0.0, B);
    //displayMatrix(B,"B");
    //displayMatrix(InvJacobianStack,"InvJacobianStack");
    //displayMatrix(ShapeFuncDerStack,"ShapeFuncDerStack");
    //displayMatrix(invJShFuncDerS,"invJShFuncDerS");
    //displayMatrix(ShapeFuncDerivatives[0],"ShapeFuncDerivatives[0]");
    //displayMatrix(ShapeFuncDerivatives[1],"ShapeFuncDerivatives[1]");
    //displayMatrix(ShapeFuncDerivatives[2],"ShapeFuncDerivatives[2]");

    //generating B^T:
    gsl_matrix * BT = gsl_matrix_alloc(nNodes*nDim,6);
    gsl_matrix_transpose_memcpy(BT,B);
    //displayMatrix(BT,"BT");
    //gsl_matrix_free(tmpMat2);
    return BT;
}

gsl_matrix* ShapeBase::calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian){
    int dim2 = nDim*nDim;
    //invrting the Jacobian:
    gsl_matrix* tmpJacobianForInversion =  gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* InvJacobian = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpJacobianForInversion,Jacobian);
    bool inverted = InvertMatrix(tmpJacobianForInversion, InvJacobian);
    if (!inverted){
        cerr<<"Jacobian not inverted!!"<<endl;
    }
    //displayMatrix(Jacobian,"Jacobian");
    //Generating the inverse Jacobian stack:
    gsl_matrix * InvJacobianStack =  gsl_matrix_calloc(dim2,dim2);
    for (int i =0; i<nDim; i++){
        for (int m=0; m<nDim; ++m){
            for (int n=0; n<3; ++n){
                gsl_matrix_set(InvJacobianStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJacobian,n,m));
            }
        }
    }
    gsl_matrix_free(tmpJacobianForInversion);
    gsl_matrix_free(InvJacobian);
    return InvJacobianStack;
}

void ShapeBase::calculateImplicitKElastic(){
    //cout<<"calculating implicit K elastic for element: "<<Id<<endl;
    int dim = nDim;
    int n = nNodes;
    //calculating K in a 3 point gaussian:
    gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
    gsl_matrix_set_zero(TriPointKe);
    double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
    for (int iter =0; iter<3;++iter){
        gsl_matrix_set_zero(currK);
        //cout<<"Calculating gauss point: "<<iter<<endl;
        calculateElasticKIntegral1(currK,iter);
        calculateElasticKIntegral2(currK,iter);
        gsl_matrix_scale(currK,weights[iter]);
        gsl_matrix_add(TriPointKe, currK);
    }
    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(TriPointKe,"TriPointKe0");


/*
    gsl_matrix* disp = gsl_matrix_calloc(18,1);
    for(int ii =0; ii<6; ++ii){
        for (int jj=0; jj<3; jj++){
            gsl_matrix_set(disp,ii*3+jj,0,Positions[ii][jj]-ReferenceShape->Positions[ii][jj]);
        }
    }
    gsl_matrix* B = Bmatrices[0];
    gsl_matrix* Fe1E = gsl_matrix_calloc(6,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, B, disp,0.0, Fe1E);
    displayMatrix(Fe1E,"Fe1Ecurr");
    displayMatrix(FeMatrices[0],"Fecurr");
    gsl_matrix* C = CMatrices[0];
    gsl_matrix* sigmacurr = gsl_matrix_calloc(6,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, C, Fe1E,0.0, sigmacurr);
    displayMatrix(sigmacurr,"sigmacurr");

*/
    gsl_matrix_free(currK);
}

void ShapeBase::writeKelasticToMainKatrix(gsl_matrix* Ke){
    for (int a=0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            int NodeId1 = NodeIds[a];
            int NodeId2 = NodeIds[b];
            NodeId1 *= nDim;
            NodeId2 *= nDim;
            for (int i=0; i<nDim; ++i){
                for (int j=0; j<nDim; ++j){
                    double valueij = gsl_matrix_get(Ke,NodeId1+i,NodeId2+j) + gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                    gsl_matrix_set(Ke,NodeId1+i,NodeId2+j,valueij);
                }
            }
        }
    }
}

void	ShapeBase::calculateCMatrix(int pointNo){
    gsl_matrix* C = CMatrices[pointNo];
    gsl_matrix* Fe = FeMatrices[pointNo];
    double detFe = determinant3by3Matrix(Fe);
    double Idouble[3][3] = {{1.0,0.0,0.0} , {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    double D81[3][3][3][3];
    double C81[3][3][3][3];
    for (int I = 0; I<3; ++I){
        for (int J = 0; J<3; ++J){
            for (int K = 0; K<3; ++K){
                for (int L = 0; L<3; ++L){
                    D81[I][J][K][L] = lambda*Idouble[K][L]*Idouble[I][J] + mu * (Idouble[I][K]*Idouble[J][L] + Idouble[I][L]*Idouble[J][K] );
                    C81[I][J][K][L] = 0.0;
                }
            }
        }
    }
/*
    gsl_matrix* DD = gsl_matrix_calloc(6,6);
    gsl_matrix_set(DD,0,0, D81[0][0][0][0]);
    gsl_matrix_set(DD,0,1, D81[0][0][1][1]);
    gsl_matrix_set(DD,0,2, D81[0][0][2][2]);
    gsl_matrix_set(DD,0,3, 0.5*(D81[0][0][0][1] + D81[0][0][1][0]));
    gsl_matrix_set(DD,0,4, 0.5*(D81[0][0][0][2] + D81[0][0][2][0]));
    gsl_matrix_set(DD,0,5, 0.5*(D81[0][0][1][2] + D81[0][0][2][1]));

    gsl_matrix_set(DD,1,1, D81[1][1][1][1]);
    gsl_matrix_set(DD,1,2, D81[1][1][2][2]);
    gsl_matrix_set(DD,1,3, 0.5*(D81[1][1][0][1] + D81[1][1][1][0]));
    gsl_matrix_set(DD,1,4, 0.5*(D81[1][1][0][2] + D81[1][1][2][0]));
    gsl_matrix_set(DD,1,5, 0.5*(D81[1][1][1][2] + D81[1][1][2][1]));

    gsl_matrix_set(DD,2,2, D81[2][2][2][2]);
    gsl_matrix_set(DD,2,3, 0.5*(D81[2][2][0][1] + D81[2][2][1][0]));
    gsl_matrix_set(DD,2,4, 0.5*(D81[2][2][0][2] + D81[2][2][2][0]));
    gsl_matrix_set(DD,2,5, 0.5*(D81[2][2][1][2] + D81[2][2][2][1]));

    gsl_matrix_set(DD,3,3, 0.5*(D81[0][1][0][1] + D81[0][1][1][0]));
    gsl_matrix_set(DD,3,4, 0.5*(D81[0][1][0][2] + D81[0][1][2][0]));
    gsl_matrix_set(DD,3,5, 0.5*(D81[0][1][1][2] + D81[1][1][2][1]));

    gsl_matrix_set(DD,4,4, 0.5*(D81[0][2][0][2] + D81[0][2][2][0]));
    gsl_matrix_set(DD,4,5, 0.5*(D81[0][2][1][2] + D81[0][2][2][1]));

    gsl_matrix_set(DD,5,5, 0.5*(D81[1][2][1][2] + D81[1][2][2][1]));

    displayMatrix(DD,"DD");
    cout<<"lambda: "<<lambda<<" mu: "<<mu<<endl;*/



    for (int i = 0; i<3; ++i){
        for (int j = 0; j<3; ++j){
            for (int k = 0; k<3; ++k){
                for (int l = 0; l<3; ++l){
                    //sum over D(I,J,K,L)
                    for (int I = 0; I<3; ++I){
                        for (int J = 0; J<3; ++J){
                            for (int K = 0; K<3; ++K){
                                for (int L = 0; L<3; ++L){
                                    C81[i][j][k][l] += gsl_matrix_get(Fe,i,I)*gsl_matrix_get(Fe,j,J)*gsl_matrix_get(Fe,k,K)*gsl_matrix_get(Fe,l,L)*D81[I][J][K][L];
                                }
                            }
                        }
                    }
                    C81[i][j][k][l] /= detFe;
                    //cout<<"C["<<i<<"]["<<j<<"]["<<k<<"]["<<l<<"]["<<"]: "<<C81[i][j][k][l]<<endl;
                }
            }
        }
    }
    gsl_matrix_set_zero(C);
    gsl_matrix_set(C,0,0, C81[0][0][0][0]);
    gsl_matrix_set(C,0,1, C81[0][0][1][1]);
    gsl_matrix_set(C,0,2, C81[0][0][2][2]);
    gsl_matrix_set(C,0,3, 0.5*(C81[0][0][0][1] + C81[0][0][1][0]));
    gsl_matrix_set(C,0,4, 0.5*(C81[0][0][0][2] + C81[0][0][2][0]));
    gsl_matrix_set(C,0,5, 0.5*(C81[0][0][1][2] + C81[0][0][2][1]));

    gsl_matrix_set(C,1,1, C81[1][1][1][1]);
    gsl_matrix_set(C,1,2, C81[1][1][2][2]);
    gsl_matrix_set(C,1,3, 0.5*(C81[1][1][0][1] + C81[1][1][1][0]));
    gsl_matrix_set(C,1,4, 0.5*(C81[1][1][0][2] + C81[1][1][2][0]));
    gsl_matrix_set(C,1,5, 0.5*(C81[1][1][1][2] + C81[1][1][2][1]));

    gsl_matrix_set(C,2,2, C81[2][2][2][2]);
    gsl_matrix_set(C,2,3, 0.5*(C81[2][2][0][1] + C81[2][2][1][0]));
    gsl_matrix_set(C,2,4, 0.5*(C81[2][2][0][2] + C81[2][2][2][0]));
    gsl_matrix_set(C,2,5, 0.5*(C81[2][2][1][2] + C81[2][2][2][1]));

    gsl_matrix_set(C,3,3, 0.5*(C81[0][1][0][1] + C81[0][1][1][0]));
    gsl_matrix_set(C,3,4, 0.5*(C81[0][1][0][2] + C81[0][1][2][0]));
    gsl_matrix_set(C,3,5, 0.5*(C81[0][1][1][2] + C81[0][1][2][1]));

    gsl_matrix_set(C,4,4, 0.5*(C81[0][2][0][2] + C81[0][2][2][0]));
    gsl_matrix_set(C,4,5, 0.5*(C81[0][2][1][2] + C81[0][2][2][1]));

    gsl_matrix_set(C,5,5, 0.5*(C81[1][2][1][2] + C81[1][2][2][1]));

    //Swapping positions of 4 and 5:
    //record Rows 4 and 5 of C matirx into a smaller matrix
    double CCol45 [6][2];
    for (int i=0;i<6; ++i){
        for (int j=0;j<2; ++j){
            CCol45[i][j] = gsl_matrix_get(C,i,j+4);
        }
    }
    //Now swap
    for (int i=0;i<6; ++i){
        gsl_matrix_set(C,i,4, CCol45[i][1]);
        gsl_matrix_set(C,i,5, CCol45[i][0]);
    }
    //Now swap rows - :
    double CRow45 [2][6];
    for (int i=0;i<6; ++i){
       CRow45[0][i] = gsl_matrix_get(C,4,i);
       CRow45[1][i] = gsl_matrix_get(C,5,i);
    }
    //Now swap
    for (int i=0;i<6; ++i){
        gsl_matrix_set(C,4,i, CRow45[1][i]);
        gsl_matrix_set(C,5,i, CRow45[0][i]);
    }
    //filling the symmetry:
    for (int i=0; i<6; ++i){
        for (int j=i+1; j<6; ++j){
            gsl_matrix_set(C,j,i,gsl_matrix_get(C,i,j));
        }
    }
    //displayMatrix(C,"Cmatrix");
}

void	ShapeBase::calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix *ExternalNodalForces){
    gsl_matrix_set_zero(ExternalNodalForces);
    int nodeNo = 0;
    for (int i=0; i<nNodes; i++){
        if (NodeIds[i] == nodeId){
            nodeNo = i;
            break;
        }
    }
    for (int pointNo = 0; pointNo<3; pointNo++){
        gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
        gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
        gsl_matrix* B = Bmatrices[pointNo];
        consturctBaTBb(B, BaT,Bb,nodeNo,0);
        gsl_matrix* NodeForces = gsl_matrix_calloc(3,1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Externalstress,0.0, NodeForces);
        gsl_matrix_scale(NodeForces,1.0/3.0);
        gsl_matrix_scale(NodeForces,detFs[pointNo]);
        gsl_matrix_scale(NodeForces,detdXdes[pointNo]);
        gsl_matrix_add(ExternalNodalForces,NodeForces);
        gsl_matrix_free(BaT);
        gsl_matrix_free(Bb);
    }
    //displayMatrix(ExternalNodalForces,"ExternalNodalForces");
}
void	ShapeBase::calculateElasticKIntegral1(gsl_matrix* currK,int pointNo){
    gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix * invJShFuncDerSWithFe = invJShapeFuncDerStackwithFe[pointNo];

    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];
    gsl_matrix* Fe = FeMatrices[pointNo];
	double detFe = determinant3by3Matrix(Fe);
    //finished calculating 4th order tensor D
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            gsl_matrix* Keab = gsl_matrix_calloc(3,3);
            double DNa[3] = {0.0,0.0,0.0};
            double DNb[3] = {0.0,0.0,0.0};

            for (int i=0;i<nDim;++i){
                // original version: DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                // original version: DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            	DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            }
            //cout<<" DNb from Fe: "<<DNb[0]<<" "<<DNb[1]<<" "<<DNb[2]<<" DNb from F: "<<DNbold[0]<<" "<<DNbold[1]<<" "<<DNbold[2]<<endl;
            //writing Kab:
            for (int i = 0 ; i<nDim; ++i){
                for (int k=0; k<nDim; ++k){
                    double value = 0;
                    //the sum over j,l,I,J,K,L, to get Kab(i,k):
                    for (int j = 0; j<nDim; ++j){
                        for (int l=0; l<nDim; ++l){

                            for (int I=0; I<nDim; ++I){
                                for (int J=0; J<nDim; ++J){
                                    for (int K=0; K<nDim; ++K){
                                        for (int L=0; L<nDim; ++L){
                                            value += (gsl_matrix_get(Fe,i,I)*gsl_matrix_get(Fe,j,J)*gsl_matrix_get(Fe,k,K)*gsl_matrix_get(Fe,l,L)*D81[I][J][K][L]*DNb[l]*DNa[j]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    value *= detF*detdXde;
                    //value /= detF;
                    value /= detFe;
                    value += gsl_matrix_get(Keab,i,k);
                    gsl_matrix_set(Keab,i,k,value);
                }
            }
            //now I have Kab for current gauss point, I need to write in into currK:
            for (int i=0; i<nDim; ++i){
                for (int j=0; j<nDim; ++j){
                    double value = gsl_matrix_get(currK,a*nDim+i, b*nDim+j);
                    value += gsl_matrix_get(Keab,i, j);
                    gsl_matrix_set(currK,a*nDim+i, b*nDim+j,value);
                }
            }
            gsl_matrix_free(Keab);
        }
    }
}


void	ShapeBase::calculateElasticKIntegral2(gsl_matrix* currK,int pointNo){
    gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix * Stress = elasticStress[pointNo];
    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];

    gsl_matrix * DNaT = gsl_matrix_calloc(1,nDim);
    gsl_matrix * DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix * Keab2 = gsl_matrix_calloc(1,1);
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            for (int i=0;i<nDim;++i){
                gsl_matrix_set(DNaT,0,i,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix * tmp1 = gsl_matrix_calloc(1,nDim);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, DNaT, Stress,0.0, tmp1);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp1, DNb,0.0, Keab2);
            double value = gsl_matrix_get(Keab2,0,0)*detF*detdXde;
            for (int i=0; i<nDim; ++i){
                int index1 = a*nDim+i;
                int index2 = b*nDim+i;
                double addedValue = gsl_matrix_get(currK,index1,index2) + value;//adding the calculated value to current K matirx
                gsl_matrix_set(currK,index1,index2,addedValue);
            }
            gsl_matrix_free(tmp1);
        }
    }
    gsl_matrix_free(DNaT);
    gsl_matrix_free(DNb);
    gsl_matrix_free(Keab2);
}

void ShapeBase::consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b){
    for (int i=0; i<6; ++i){
        for (int j=0; j<nDim; ++j){
            //double value = gsl_matrix_get(B,i,a*dim+j);
            gsl_matrix_set(BaT,j,i,gsl_matrix_get(B,i,a*nDim+j)); //transpose of Ba
            //value  = gsl_matrix_get(B,i,b*dim+j);
            gsl_matrix_set(Bb,i,j,gsl_matrix_get(B,i,b*nDim+j)); //Bb
            //displayMatrix(Bb,"Bb");
        }
    }
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

void 	ShapeBase::cleanMyosinForce(){
	for (int i=0; i<nNodes; ++i){
		MyoForce[i][0] = 0;
		MyoForce[i][1] = 0;
		MyoForce[i][2] = 0;
	}
}

bool ShapeBase::calculateIfInsideActiveStripe(double initialPoint,double endPoint, double stripeSize1, double stripeSize2){
	//All nodes and the centre of the element should be inside the active zone:
	//starting from the centre:
	double* c;
	c = new double[3];
	c = getCentre();
	for (int i=0; i<nNodes+1; ++i){
		//getting the node position
		double x;
		double y;
		double bufferFrac;
		if (i == 0 ){
			x = c[0];
			y = c[1];
			bufferFrac = 0.0;
		}
		else{
			x = Positions[i-1][0];
			y = Positions[i-1][1];
			bufferFrac = 0.1;
		}
		//is the node inside the active region:
		bool xInActiveZone = false;
		bool yInActiveZone = false;
		bool loopComplete = false;
		double lowEnd = initialPoint - bufferFrac*stripeSize1;
		double highEnd = lowEnd + stripeSize1 + 2.0*bufferFrac*stripeSize1;
		while (!xInActiveZone && !loopComplete){
			if(stripeSize1 == 0 || highEnd>endPoint){
				highEnd = endPoint + bufferFrac*stripeSize1;
				loopComplete = true;
			}
			//cout<<" x: "<<x<<" low: "<<lowEnd<<" high: "<<highEnd<<endl;
			if (x>= lowEnd && x<=highEnd){
				xInActiveZone = true;
			}
			lowEnd += stripeSize1*2.0;
			highEnd = lowEnd + stripeSize1 + 2.0*bufferFrac*stripeSize1;
		}
		if (xInActiveZone){
			//if x is in active zone, I will move on to check y:
			loopComplete = false;
			lowEnd = initialPoint;
			highEnd = lowEnd + stripeSize2 + 2.0*bufferFrac*stripeSize2;
			while (!yInActiveZone && !loopComplete){
				if(stripeSize2 == 0 || highEnd>endPoint){
					highEnd = endPoint + bufferFrac*stripeSize2;
					loopComplete = true;
				}
				if (y>= lowEnd && y<=highEnd){
					yInActiveZone = true;
				}
				lowEnd += stripeSize2*2.0;
				highEnd = lowEnd + stripeSize2 + 2.0*bufferFrac*stripeSize2;
			}
		}
		//cout<<"Element: "<<Id<<" point ( 0 for centre, i-1 for node): "<<i<<" pos: "<<x<<" "<<y<<" is inside: "<<xInActiveZone<<" "<<yInActiveZone<<endl;
		if (!xInActiveZone || !yInActiveZone ){
			//if this node is not in the active zone, then the element is not in the active zone
			delete[] c;
			return false;
		}
	}
	//I did not return the function in any of the nodes, then all nodes must be inside the active zone:
	delete[] c;
	return true;
};

double ShapeBase::getCmyosinUniformForNode (int TissuePlacement){
	if(TissuePlacement == 1) {	//apical node
		return cMyoUniform[0];
	}
	if(TissuePlacement == 0) {	//basal node
		return cMyoUniform[1];
	}
	return 0.0;
}

double ShapeBase::getCmyosinUnipolarForNode (int TissuePlacement){
	if(TissuePlacement == 1) {
		return cMyoUnipolar[0];
	}
	if(TissuePlacement == 0) {
		return cMyoUnipolar[1];
	}
	return 0.0;
}

void ShapeBase::getMyosinLevels (double *cMyo){
	cMyo[0] = cMyoUniform[0];
	cMyo[1] = cMyoUniform[1];
	cMyo[2] = cMyoUnipolar[0];
	cMyo[3] = cMyoUnipolar[1];
}

void ShapeBase::getEquilibriumMyosinLevels (double* cMyoEq){
	cMyoEq[0] = cMyoUniformEq[0];
	cMyoEq[1] = cMyoUniformEq[1];
	cMyoEq[2] = cMyoUnipolarEq[0];
	cMyoEq[3] = cMyoUnipolarEq[1];
}

void ShapeBase::setMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1){
	cMyoUniform[0] = cUni0;
	cMyoUniform[1] = cUni1;
	cMyoUnipolar[0] = cPol0;
	cMyoUnipolar[1] = cPol1;
}

void ShapeBase::setEquilibriumMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1){
	cMyoUniformEq[0] = cUni0;
	cMyoUniformEq[1] = cUni1;
	cMyoUnipolarEq[0] = cPol0;
	cMyoUnipolarEq[1] = cPol1;
}

void	ShapeBase::updateUniformEquilibriumMyosinConcentration(bool isApical, double cEqUniform){
	if (isApical){
		cMyoUniformEq[0] = cEqUniform;
	}
	else{
		cMyoUniformEq[1] = cEqUniform;
	}
}

void	ShapeBase::updateUnipolarEquilibriumMyosinConcentration(bool isApical, double cEqUnipolar, double orientationX, double orientationY){
	int indice = 1;
	if (isApical){
		indice = 0;
	}
	cMyoUnipolarEq[indice] = cEqUnipolar;
	gsl_matrix_set(myoPolarityDir,indice,0,orientationX);
	gsl_matrix_set(myoPolarityDir,indice,1,orientationY);
	gsl_matrix_set(myoPolarityDir,indice,2,0.0);
}

void	ShapeBase::updateMyosinConcentration(double dt, double kMyo){
	//gsl_matrix_set(myoPolarityDir,0,0,1.0);
	//gsl_matrix_set(myoPolarityDir,0,1,0.0);
	//gsl_matrix_set(myoPolarityDir,0,2,0.0);
	double thresholdValue = 1E-8, thresholdFraction= 0.01;
	//the value of kMyo is taken form my thesis
	double currMyoDt[3] = {dt,dt*2.0,dt/2.0};
	double cFinal[3];
	//First value is the final value with the current time step,
	//second is with currTimeStep*2 and
	//third is with 0.5 currTimeStep;
	double cInitial, cEq;
	for (int myoIter =0; myoIter<4; myoIter++){
		if (myoIter == 0){
			cInitial = cMyoUniform[0];
			cEq = cMyoUniformEq[0];
		}
		else if (myoIter == 1){
			cInitial = cMyoUniform[1];
			cEq = cMyoUniformEq[1];
		}
		else if (myoIter == 2){
			cInitial = cMyoUnipolar[0];
			cEq = cMyoUnipolarEq[0];
		}
		else if (myoIter == 3){
			cInitial = cMyoUnipolar[1];
			cEq = cMyoUnipolarEq[1];
		}
		bool converged = false;
		while (!converged){
			int steps[3] = {dt/currMyoDt[0],dt/currMyoDt[1],dt/currMyoDt[2]};
			//cout<<"steps: "<<steps[0]<<" "<<steps[1]<<" "<<steps[2]<<" currMyoDt: "<<currMyoDt[0]<<" "<<currMyoDt[1]<<" "<<currMyoDt[2]<<endl;
			for (int j=0; j<3; ++j){
				cFinal[j] = cInitial;
				for (int i =0 ;i<steps[j]; ++i){
					cFinal[j] += (cEq - cFinal[j])*kMyo*currMyoDt[j];
				}
			}
			//check if the value of the current dt and half current dt are below the threshold:
			//cout<<"cFinal: "<<cFinal[0]<<" "<<cFinal[1]<<" "<<cFinal[2]<<endl;
			double diff = fabs((cFinal[0] - cFinal[2]));
			//cout<<" diff is : "<<diff<<" threshold is :"<<thresholdValue<<endl;
			if ( diff < thresholdValue ){
				converged = true;
				//cout<<"converged by difference"<<endl;
			}
			else if( fabs(diff / cFinal[2]) < thresholdFraction){
				converged = true;
				//cout<<"converged by fraction"<<endl;
			}
			else{
				currMyoDt[1] = currMyoDt[0];
				currMyoDt[0] = currMyoDt[2];
				currMyoDt[2] *= 0.5;
				//cout<<"updated currMyoDt: "<<currMyoDt[0]<<" "<<currMyoDt[1]<<" "<<currMyoDt[2]<<endl;
			}
			//need to implement increasing this
		}
		if (myoIter == 0){
			cMyoUniform[0] = cFinal[2];
		}
		else if (myoIter == 1){
			cMyoUniform[1] = cFinal[2];
		}
		else if (myoIter == 2){
			cMyoUnipolar[0] = cFinal[2];
		}
		else if (myoIter == 3){
			cMyoUnipolar[1] = cFinal[2];
		}
	}
	//cout<<"Element: "<<Id<<" EQ myosin levels: "<<cMyoUniformEq[0]<<" "<<cMyoUniformEq[1]<<" "<<cMyoUnipolarEq[0]<<" "<<cMyoUnipolarEq[1]<<endl;
	//cout<<"Element: "<<Id<<" myosin levels: "<<cMyoUniform[0]<<" "<<cMyoUniform[1]<<" "<<cMyoUnipolar[0]<<" "<<cMyoUnipolar[1]<<endl;
}

void 	ShapeBase::setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
	ShapeChangeRate[0] = x;
	ShapeChangeRate[1] = y;
	ShapeChangeRate[2] = z;
	ShapeChangeRate[3] = xy;
	ShapeChangeRate[4] = yz;
	ShapeChangeRate[5] = xz;
}

void 	ShapeBase::updateGrowthRate(double growthx, double growthy, double growthz){
    GrowthRate[0] += growthx;
    GrowthRate[1] += growthy;
    GrowthRate[2] += growthz;
}

void 	ShapeBase::updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
	ShapeChangeRate[0] += x;
	ShapeChangeRate[1] += y;
	ShapeChangeRate[2] += z;
	ShapeChangeRate[3] += xy;
	ShapeChangeRate[4] += yz;
	ShapeChangeRate[5] += xz;
}

bool 	ShapeBase::InvertMatrix(gsl_matrix* input, gsl_matrix* inverse){
    // Define the dimension n of the matrix
    // and the signum s (for LU decomposition)
    int s;

    // Define all the used matrices
    gsl_permutation * perm = gsl_permutation_alloc (input->size1);

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (input, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert (input, perm, inverse);

    return true;
}

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

void	ShapeBase::crossProduct3D(double* u, double* v, double* cross){
	cross[0] = u[1]*v[2] - u[2]*v[1];
	cross[1] = u[2]*v[0] - u[0]*v[2];
	cross[2] = u[0]*v[1] - u[1]*v[0];
}

void	ShapeBase::crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross){
	gsl_vector_set(cross,0, ( gsl_vector_get(u,1)*gsl_vector_get(v,2) - gsl_vector_get(u,2)*gsl_vector_get(v,1) ) );
	gsl_vector_set(cross,1, ( gsl_vector_get(u,2)*gsl_vector_get(v,0) - gsl_vector_get(u,0)*gsl_vector_get(v,2) ) );
	gsl_vector_set(cross,2, ( gsl_vector_get(u,0)*gsl_vector_get(v,1) - gsl_vector_get(u,1)*gsl_vector_get(v,0) ) );
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
void	ShapeBase::normaliseVector3D(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	double z = gsl_vector_get(v,2);
	double mag2 = x*x + y*y + z*z;
	if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
		double mag = pow(mag2,0.5);
		gsl_vector_scale(v,1.0/mag);
	}
}

double	ShapeBase::getNormVector3D(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	double z = gsl_vector_get(v,2);
	double mag2 = x*x + y*y + z*z;
	return pow(mag2,0.5);
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

void 	ShapeBase::displayMatrix(gsl_matrix* mat, string matname){
    int m = mat->size1;
    int n = mat->size2;
    cout<<matname<<": "<<endl;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            cout.precision(4);
            cout.width(6);
            cout<<gsl_matrix_get(mat,i,j)<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void 	ShapeBase::displayMatrix(gsl_vector* mat, string matname){
    int m = mat->size;
    cout<<matname<<": "<<endl;

    for (int i =0; i<m; i++){
        cout.precision(4);
        cout.width(6);
        cout<<gsl_vector_get(mat,i)<<endl;
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
        Nodes[NodeIds[i]]->mass += VolumePerNode;
	}
}
void 	ShapeBase:: assignSurfaceAreaToNodes(vector <Node*>& Nodes){
    double multiplier = 1.0;
    if (ShapeType ==1 ){ multiplier = 0.5;}
    for (int i=0; i<nNodes; i++){
        Nodes[NodeIds[i]]->surface +=ReferenceShape->BasalArea/(multiplier*nNodes);
	}
}

void 	ShapeBase::calculateZProjectedAreas(){
    double Threshold = 1E-5;
    int id0 = 0, id1 = 1, id2 = 2; // this is correct for basal side, I will change it for apical calculation
    for (int tissueSide = 0; tissueSide<2; tissueSide++){
        if ( tissueSide == 1){
            //I am calculating basal area,
            id0 = 3;
            id1 = 4;
            id2 = 5;
        }
        double sideVec1[2];
        double sideVec2[2];
        double Side1 = 0.0;
        double Side2 = 0.0;
        double costet = 0.0;
        double Area = 0.0;
        for (int i = 0; i<2; ++i){
            sideVec1[i]= Positions[id1][i] - Positions[id0][i];
            sideVec2[i]= Positions[id2][i] - Positions[id0][i];
            costet += sideVec1[i] * sideVec2[i];
            Side1  += sideVec1[i] * sideVec1[i];
            Side2  += sideVec2[i] * sideVec2[i];
        }
        if (Side1 > Threshold && Side2 > Threshold){
            Side1 = pow(Side1,0.5);
            Side2 = pow(Side2,0.5);
            costet /= (Side1*Side2);
            double sintet = pow((1-costet*costet),0.5);
            Area = Side1* Side2 * sintet / 2.0;
        }
        if(tissueSide == 0){
            ZProjectedBasalArea = Area;
        }
        else{
            ZProjectedApicalArea = Area;
        }
    }
}




void 	ShapeBase::assignZProjectedAreas(vector <Node*> Nodes){
    if (ShapeType == 1 ){ //only written for prisms
        for (int i=0; i<3; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedBasalArea/3.0;
        }
        for (int i=3; i<6; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedApicalArea/3.0;
        }
    }
}

void 	ShapeBase:: assignElementToConnectedNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
		Nodes[NodeIds[i]]->connectedElementIds.push_back(Id);
		double weightfFraction = (ReferenceShape->Volume/nNodes)/Nodes[NodeIds[i]]->mass;
		Nodes[NodeIds[i]]->connectedElementWeights.push_back(weightfFraction);
	}
}

void 	ShapeBase::removeMassFromNodes(vector <Node*>& Nodes){
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
			//All weights are normlised as the sum will make 1.0. Now I do not want this element in the weighing,
			//it does not have a mass anymore, therefore I will multiply all the remaining weights with (1-w_ablated);
		}
}

void 	ShapeBase::checkDisplayClipping(double xClip, double yClip, double zClip){
	IsClippedInDisplay=false;
	for (int j=0; j<nNodes; ++j){
		 if(Positions[j][0]>xClip){
			 IsClippedInDisplay = true;
			 return;
		 }
		 if(Positions[j][1]<yClip){
			 IsClippedInDisplay = true;
			 return;
		 }
		 if(Positions[j][2]>zClip){
			 IsClippedInDisplay = true;
			 return;
		 }
	 }
}
