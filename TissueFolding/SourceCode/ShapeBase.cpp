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

int		ShapeBase::getNodeId(int i){
	return NodeIds[i];
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
double 	ShapeBase::getInternalViscosity(){
	return internalViscosity;
}
void ShapeBase::updateInternalViscosityTest(){
	double d[2] = {0.0,0.0};
	for (int i = 0; i<nNodes; ++i ){
		d[0] += Positions[i][0];
		d[1] += Positions[i][1];
	}
	d[0] /= nNodes; d[1] /= nNodes;
	double dmag = d[0]*d[0] + d[1]*d[1];
	dmag = pow(dmag,0.5);
	internalViscosity = dmag*originalInternalViscosity;
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


void ShapeBase::getRelativePositionInTissueInGridIndex(int nGridX, int nGridY , int& IndexX, int& IndexY, double& FracX, double& FracY){
	//cout<<"inside getRelativePositionInTissueInGridIndex"<<endl;
	double* reletivePos = new double[2];
	getRelativePosInBoundingBox(reletivePos);
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
	delete[] reletivePos;
}

void ShapeBase::getInitialRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY){
	double* reletivePos = new double[2];
	getInitialRelativePosInBoundingBox(reletivePos);
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
	delete[] reletivePos;

}

bool ShapeBase::isGrowthRateApplicable( int sourceTissue, double& weight){
	//wight is the weight of hte current tissue in linker sites
	if (sourceTissue == 0){//columnar layer growth
		if (tissueType == 0){ //columnar
			weight = 1.0;
			return true;
		}
		else if(tissueType == 2){ //linker
			weight = columnarGrowthWeight;
			return true;
		}
	}
	else if (sourceTissue == 1){//peripodial membrane growth
		if (tissueType == 1){ //peripodial
			weight = 1.0;
			return  true;
		}
		else if ( tissueType == 2) { //linker
			weight = peripodialGrowthWeight;
			return true;
		}
	}
	return false;
}

void ShapeBase::calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue){
	double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue, tissueWeight);
	if (continueCalaculation){
		gsl_matrix_set(increment,0,0,exp(x*tissueWeight*dt));
		gsl_matrix_set(increment,1,1,exp(y*tissueWeight*dt));
		gsl_matrix_set(increment,2,2,exp(z*tissueWeight*dt));
		gsl_matrix* temp = gsl_matrix_calloc(3,3);
		gsl_matrix* rotMatT = gsl_matrix_calloc(3,3);
		gsl_matrix_transpose_memcpy(rotMatT,rotMat);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, increment);
		gsl_matrix_free(temp);
		gsl_matrix_free(rotMatT);
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}
/*
void ShapeBase::calculateFgFromGridCorners(double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue,  int IndexX, int IndexY, double FracX, double FracY){
	double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue,tissueWeight);
	if (continueCalaculation){
		gsl_matrix* corner0 = gsl_matrix_calloc(3,3);
		gsl_matrix* corner1 = gsl_matrix_calloc(3,3);
		gsl_matrix* corner2 = gsl_matrix_calloc(3,3);
		gsl_matrix* corner3 = gsl_matrix_calloc(3,3);
		double growth0[3], growth1[3], growth2[3], growth3[3];
		gsl_matrix* rotMat0;
		gsl_matrix* rotMat1;
		gsl_matrix* rotMat2;
		gsl_matrix* rotMat3;
		gsl_matrix* rotMatT = gsl_matrix_calloc(3,3);
		for (int axis =0; axis<3; axis++){
			growth0[axis] = currGF->getGrowthMatrixElement(IndexX,IndexY,axis)*tissueWeight*(1.0-FracX)*(1.0-FracY);
			gsl_matrix_set(corner0,axis,axis,exp(growth0[axis]*dt));
			growth1[axis] = currGF->getGrowthMatrixElement(IndexX+1,IndexY,axis)*tissueWeight*FracX*(1.0-FracY);
			gsl_matrix_set(corner1,axis,axis,exp(growth1[axis]*dt));
			growth2[axis] = currGF->getGrowthMatrixElement(IndexX,IndexY+1,axis)*tissueWeight*(1.0-FracX)*FracY;
			gsl_matrix_set(corner2,axis,axis,exp(growth2[axis]*dt));
			growth3[axis] = currGF->getGrowthMatrixElement(IndexX+1,IndexY+1,axis)*tissueWeight*FracX*FracY;
			gsl_matrix_set(corner3,axis,axis,exp(growth3[axis]*dt));
		}
		rotMat0 = currGF->getXyShearRotationsMatrixElement(IndexX,IndexY);
		rotMat1 = currGF->getXyShearRotationsMatrixElement(IndexX+1,IndexY);
		rotMat2 = currGF->getXyShearRotationsMatrixElement(IndexX,IndexY+1);
		rotMat3 = currGF->getXyShearRotationsMatrixElement(IndexX+1,IndexY+1);
		gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
		gsl_matrix_transpose_memcpy(rotMatT,rotMat0);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat0, corner0, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, corner0);
		gsl_matrix_transpose_memcpy(rotMatT,rotMat1);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat1, corner1, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, corner1);
		gsl_matrix_transpose_memcpy(rotMatT,rotMat2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat2, corner2, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, corner2);
		gsl_matrix_transpose_memcpy(rotMatT,rotMat3);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat3, corner3, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, corner3);
		//applying all the growhts to increment:
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, corner2, corner3, 0.0, increment);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, corner1, increment, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, corner0, temp, 0.0, increment);
		gsl_matrix_free(temp);
		gsl_matrix_free(corner0);
		gsl_matrix_free(corner1);
		gsl_matrix_free(corner2);
		gsl_matrix_free(corner3);
		gsl_matrix_free(rotMatT);
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}*/


void ShapeBase::calculateFgFromGridCorners(int gridGrowthsInterpolationType, double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue,  int IndexX, int IndexY, double FracX, double FracY){
	double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue,tissueWeight);
	if (continueCalaculation){
		//taking growth data around 4 grid points
		//
		// Grid shape:
		//    [point 2] ------------- [point 3]
		//       |                        |
		//       |<--fracX----> (o)       |
		//       |               |        |
		//       |               |        |
		//       |             fracY      |
		//       |               |        |
		//    [point 0] ------------- [point 1]
		//
		double *growth0, *growth1, *growth2, *growth3;
		growth0 = new double[3];
		growth1 = new double[3];
		growth2 = new double[3];
		growth3 = new double[3];
		double *angles;
		angles = new double[4];
		bool *angleEliminated;
		angleEliminated = new bool[4];
		currGF->getGrowthProfileAt4Corners(IndexX, IndexY, growth0, growth1, growth2, growth3, angles, angleEliminated);
		double growth[3];
		double angle;
		if (gridGrowthsInterpolationType == 0){
			//using the growth rate at the grid point:
			//if fraction is below 0,5, I will use the index availabe.
			//If it is above 0.5, it is within the range of next groid point, I will use index +1:
			if(FracX > 0.5){
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth3[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth1[i];
					}
				}
			}
			else{
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth2[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth0[i];
					}
				}
			}

		}
		else if (gridGrowthsInterpolationType == 1){
			//calculating the angle fraction eliminated, if any:
			double FracEliminated = 0.0;
			if (angleEliminated[0]){ FracEliminated += (1.0-FracX)*(1.0-FracY);	}
			if (angleEliminated[1]){ FracEliminated += FracX*(1.0-FracY);		}
			if (angleEliminated[2]){ FracEliminated += (1.0-FracX)*FracY;		}
			if (angleEliminated[3]){ FracEliminated += FracX*FracY;				}
			//taking the linear interpolation of 4 angles at 4 grid points:
			angle = angles[0]*(1.0-FracX)*(1.0-FracY)+angles[1]*FracX*(1.0-FracY)+angles[2]*(1.0-FracX)*FracY+angles[3]*FracX*FracY;
			if (FracEliminated>0){
				if (FracEliminated >= 0.9999999){
					angle = 0.0; //if all the angles should be eliminated because all corners have low aspect ratio, then angle is arbitrary, selected as zero
				}
				else{
					angle /= (1.0-FracEliminated); //normalising the sum to the eliminated averaging
				}
			}
			//taking the linear interpolation of 4 growth rates at 4 grid points
			for (int axis =0; axis<3; axis++){
				growth[axis]  = growth0[axis]*(1.0-FracX)*(1.0-FracY)+growth1[axis]*FracX*(1.0-FracY)+growth2[axis]*(1.0-FracX)*FracY+growth3[axis]*FracX*FracY;
				growth[axis] *= tissueWeight;
				//if (tissuePlacement == 0){
					//Prevented basal element z growth! Correct this later!!!
					//growth[2] = 0;
				//}
			}
		}
		//write the increment from obtained growth:
		for (int axis =0; axis<3; axis++){
			gsl_matrix_set(increment,axis,axis,exp(growth[axis]*dt));
		}
		//Rotate the growth if the angel is not zero:
		if (angle != 0.0){
			gsl_matrix* rotMat  = gsl_matrix_calloc(3,3);
			gsl_matrix* rotMatT = gsl_matrix_calloc(3,3);
			double c = cos(angle);
			double s = sin(angle);
			gsl_matrix_set(rotMat,0,0,  c );
			gsl_matrix_set(rotMat,0,1, -1.0*s);
			gsl_matrix_set(rotMat,0,2,  0.0);
			gsl_matrix_set(rotMat,1,0,  s);
			gsl_matrix_set(rotMat,1,1,  c);
			gsl_matrix_set(rotMat,1,2,  0.0);
			gsl_matrix_set(rotMat,2,0,  0.0);
			gsl_matrix_set(rotMat,2,1,  0.0);
			gsl_matrix_set(rotMat,2,2,  1.0);
			gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
			gsl_matrix_transpose_memcpy(rotMatT,rotMat);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, rotMatT, 0.0, increment);
			gsl_matrix_free(temp);
			gsl_matrix_free(rotMat);
			gsl_matrix_free(rotMatT);
		}
		/*if ( Id == 176 ) {//|| Id == 1 || Id == 195 || Id == 152 || Id == 48 || Id == 71){
			cout<<"Element: "<<Id<<endl;
			cout<<" growth0: "<<growth0[0]<<" "<<growth0[1]<<" "<<growth0[2]<<endl;
			cout<<" growth1: "<<growth1[0]<<" "<<growth1[1]<<" "<<growth1[2]<<endl;
			cout<<" growth2: "<<growth2[0]<<" "<<growth2[1]<<" "<<growth2[2]<<endl;
			cout<<" growth3: "<<growth3[0]<<" "<<growth3[1]<<" "<<growth3[2]<<endl;
			cout<<" angles: "<<angles[0]<<" "<<angles[1]<<" "<<angles[2]<<" "<<angles[3]<<endl;
			cout<<" angleEliminated: "<<angleEliminated[0]<<" "<<angleEliminated[1]<<" "<<angleEliminated[2]<<" "<<angleEliminated[3]<<endl;

			cout<<"	angle: "<<angle<<" in degrees: "<<angle*180.0/M_PI<<" growth: "<<growth[0]<<" "<<growth[1]<<" "<<growth[2]<<endl;
			cout<<"	pos base:  "<<Positions[0][0]<<" "<<Positions[0][1]<<" "<<Positions[0][2]<<endl;
			cout<<"	        :  "<<Positions[1][0]<<" "<<Positions[1][1]<<" "<<Positions[1][2]<<endl;
			cout<<"	        :  "<<Positions[2][0]<<" "<<Positions[2][1]<<" "<<Positions[2][2]<<endl;
			cout<<"	grid values: "<<IndexX<<" "<<IndexY<<" "<<FracX<<" "<<FracY<<endl;
			cout<<"	relative pos: "<<relativePosInBoundingBox[0]<<" "<<relativePosInBoundingBox[1]<<endl;
			cout<<" FracEliminated: "<<FracEliminated<<endl;
			cout<<" angles: ";
			for (int i=0; i<4; ++i){cout<<angles[i]<<" ";}
			cout<<endl<<" angleEliminated: ";
			for (int i=0; i<4; ++i){cout<<angleEliminated[i]<<" ";}
			cout<<endl;
			displayMatrix(increment,"currIncrement");
		}*/
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}

void ShapeBase::updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial ){
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	if (tissueType == 0){//columnar layer element, no peripodial application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 1){//peripodial layer element, no columnar application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 2){//linker between columnar and peripodial layer element, the growths are already weighted, need to apply both
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, temp, 0.0, growthIncrement);
	}
	for (int i=0; i<3; ++i){
		//this is used for display purposes of the simulation. As each new value is added to growthIncrement, I an update this directly, and the result is already cumulative of multiple growth functions
		GrowthRate[i] = gsl_matrix_get(growthIncrement,i,i);
	}
	gsl_matrix_free(temp);
	//if (Id == 38 || Id == 39){
	//	cout<<"Element: "<<Id<<endl;
	//	displayMatrix(growthIncrement,"growthIncrement");
	//}
}
/*
void ShapeBase::getRelativePosInColumnarBoundingBox(double* relativePos){
	relativePos[0] =  columnarRelativePosInBoundingBox[0];
	relativePos[1] =  columnarRelativePosInBoundingBox[1];
}

void ShapeBase::getRelativePosInPeripodialBoundingBox(double* relativePos){
	relativePos[0] =  peripodialRelativePosInBoundingBox[0];
	relativePos[1] =  peripodialRelativePosInBoundingBox[1];
}*/

void ShapeBase::getRelativePosInBoundingBox(double* relativePos){
	relativePos[0] =  relativePosInBoundingBox[0];
	relativePos[1] =  relativePosInBoundingBox[1];
}

void ShapeBase::setInitialRelativePosInBoundingBox(){
	initialRelativePosInBoundingBox[0] = relativePosInBoundingBox[0];
	initialRelativePosInBoundingBox[1] = relativePosInBoundingBox[1];
}

void ShapeBase::getInitialRelativePosInBoundingBox(double* relativePos){
	relativePos[0] =  initialRelativePosInBoundingBox[0];
	relativePos[1] =  initialRelativePosInBoundingBox[1];
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
	spansWholeTissue = false;
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
				spansWholeTissue = true;
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
		//ASK LINKER ZONE BEFORE THIS, SOME LINKER ELEMENTS CAN HAVE LINKER NODES AND OTHER TISSUE NODES, NO COLUMNAR ELEMENT OR PERIPODIAL ELEMENT SHOULD HAVE A LINKER NODE
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

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal, double viscosityMid){
	this -> internalViscosity = viscosityMid;
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
	}
	if (Id>754 && Id<769){
		cout<<"element: "<<Id<<" tissue placement: "<<tissuePlacement<<" internal viscosity: "<<internalViscosity<<endl;
	}
	this -> originalInternalViscosity = internalViscosity;
}

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal){
	this -> internalViscosity = 0.5*(viscosityApical+viscosityBasal);
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
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
		//DV
        StrainMag = gsl_matrix_get(Strain,0,0);
	}
	else if (type == 2){
		//AP
        StrainMag = gsl_matrix_get(Strain,1,0);
	}
	else if (type == 3){
		//AB
        StrainMag = gsl_matrix_get(Strain,2,0);
	}
	else if (type == 4){
		//xy
        StrainMag = gsl_matrix_get(Strain,3,0);
	}
	else if (type == 5){
		//yz
        StrainMag = gsl_matrix_get(Strain,4,0);
	}
	else if (type == 3){
		//xz
        StrainMag = gsl_matrix_get(Strain,5,0);
	}
	else{
		return;
	}
}

void 	ShapeBase::getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag){
	PysPropMag = 0.0;
	if (type == 0){
		PysPropMag = Nodes[NodeIds[NodeNo]] -> externalViscosity[2];
	}
}

void 	ShapeBase::getPysProp(int type, float &PysPropMag, double dt){
	if (type == 1){
		PysPropMag = getInternalViscosity();
	}
	else if (type == 2){
		PysPropMag = getYoungModulus();
	}
	else if (type == 3 ){
		PysPropMag = getPoissonRatio();
	}
	else if (type == 4){
		double* growth;
		growth = getGrowthRate();
        double timescale = 60.0*60.0; //reporting per hour
        for (int i =0 ; i< nDim-1 ; ++i){ //reporting only x & y
        //for (int i =0 ; i< nDim ; ++i){
			//growth is in form exp(r*dt), get r first, then adjust the time scale, and report the exponential form still:
			//And I want to rate of volumetric growth, that is x*y*z
			double value = exp(log(growth[i])/dt*timescale);
			PysPropMag *= value;
		}
        //converting to percentage increase per hour:
        PysPropMag -= 1.0;
        PysPropMag *= 100.0;
	}
	else if (type == 5){
		double* shapechange;
		shapechange = getShapeChangeRate();
		PysPropMag = shapechange[2];
	}
}

void 	ShapeBase::displayIdentifierColour(){
	cout <<" IdentifierColour:  "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}
/*
void 	ShapeBase::resetCurrStepShapeChangeData(){
	for (int i=0;i<3;++i){
		CurrShapeChangeToAdd[i] = 0.0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	IsChangingShape = false;
}
*/
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
    if (rotatedGrowth){
        gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
        gsl_matrix* GrowthStrainsRotMatT = gsl_matrix_calloc(nDim,nDim);
        gsl_matrix_transpose_memcpy(GrowthStrainsRotMatT,GrowthStrainsRotMat);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, GrowthStrainsRotMatT, growthIncrement, 0.0, temp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, GrowthStrainsRotMat, 0.0, growthIncrement);
    	//gsl_matrix_memcpy(growthIncrement, temp);
    	gsl_matrix_free(temp);
    	gsl_matrix_free(GrowthStrainsRotMatT);
    }
    //incrementing Fg with current growth rate:
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, growthIncrement, Fg, 0.0, temp1);
    gsl_matrix_memcpy(Fg, temp1);
    gsl_matrix * tmpFgForInversion = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion,Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
        cerr<<"Fg not inverted!!"<<endl;
    }
    double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
    VolumePerNode = GrownVolume/nNodes;
    //freeing matrices allocated in this function
    gsl_matrix_free(temp1);
    gsl_matrix_free(tmpFgForInversion);
	//if (Id == 38) {displayMatrix(Fg,"element38Fg");}
}

void	ShapeBase::calculatePlasticDeformation(bool volumeConserved, double rate){
	gsl_matrix* TriPointFe = gsl_matrix_calloc(3,3);
	gsl_matrix* FplasticIncrement = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(FplasticIncrement);
    //gsl_matrix_set_identity(TriPointFe);
	double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
	for (int iter =0; iter<3;++iter){
		gsl_matrix* currFe =  gsl_matrix_calloc(3,3);
		createMatrixCopy(currFe,FeMatrices[iter]);
		gsl_matrix_scale(currFe,weights[iter]);
		gsl_matrix_add(TriPointFe, currFe);
	}
	double p[3] = {gsl_matrix_get(TriPointFe,0,0),gsl_matrix_get(TriPointFe,1,1),gsl_matrix_get(TriPointFe,2,2)};
	for (int i=0;i<3;++i){
		p[i] -= 1.0;
		p[i] *= rate;
		p[i] += 1.0;
		gsl_matrix_set(FplasticIncrement,i,i,p[i]);
	}
	if (volumeConserved){
		double det = determinant3by3Matrix(FplasticIncrement);
		double scale = 1.0/pow (det,1.0/3.0);
		gsl_matrix_scale(FplasticIncrement,scale);
	}
	gsl_matrix* temp = gsl_matrix_calloc(3,3);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FplasticIncrement, Fplastic, 0.0, temp);
    gsl_matrix_memcpy(Fplastic, temp);
	gsl_matrix * tmpFplasticForInversion = gsl_matrix_calloc(3,3);
	createMatrixCopy(tmpFplasticForInversion,Fplastic);
	bool inverted = InvertMatrix(tmpFplasticForInversion, invFplastic);
	if (!inverted){
		cerr<<"Fplastic not inverted!!"<<endl;
	}
	gsl_matrix_free(FplasticIncrement);
    gsl_matrix_free(TriPointFe);
    gsl_matrix_free(temp);
    gsl_matrix_free(tmpFplasticForInversion);
}

void 	ShapeBase::CalculateGrowthRotationByF(){
    gsl_matrix* rotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(rotMat);
    //rotatedGrowth = false;
    //updating the F for the current shape positions
    //(not using leftovers from previous iteration)
    calculateTriPointFForRatation();
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

void 	ShapeBase::calculateTriPointFForRatation(){
	gsl_matrix_set_zero(TriPointF);
    gsl_matrix* currF = gsl_matrix_calloc(nDim,nDim);
    //The point order is established in shape function derivative calculation!
    //Make sure the weights fir in with the order - eta zeta nu:
    //double points[3][3]={{1.0/6.0,1.0/6.0,0.0},{2.0/3.0,1.0/6.0,0.0},{1.0/6.0,2.0/3.0,0.0}};
    double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
    for (int iter =0; iter<3;++iter){
        //cout<<"Calculating gauss point: "<<eta<<" "<<nu<<" "<<zeta<<endl;
        calculateCurrTriPointFForRotation(currF,iter);
        gsl_matrix_scale(currF,weights[iter]);
        gsl_matrix_add(TriPointF, currF);
    }
}

bool 	ShapeBase::disassembleRotationMatrixForZ(gsl_matrix* rotMat){
    //For a rotation matrix R = [r11 r12 r13
    //                          r21 r22 r23
    //                          r31 r32 r33]
    //The rotations in X, Y and Z are calculated as follows:
    // tethaX = atan2(r32,r33);
    // tethaY = atan2(-r31, sqrt(r32 * r32 + r33 * r33))
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
    //rotatedGrowth_tethaZ = tethaZ;
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
        isFlipped = true;
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
    return false; //none of the off - diagonal terms of the matrix are above the threshold, the current rotation is only numerical error.
}

void 	ShapeBase::calculateRelativePosInBoundingBox(double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth){
	relativePosInBoundingBox = getCentre();
	relativePosInBoundingBox[0] = (relativePosInBoundingBox[0] -boundingBoxXMin) / boundingBoxLength;
	relativePosInBoundingBox[1] = (relativePosInBoundingBox[1] - boundingBoxYMin) / boundingBoxWidth;
}
/*
void 	ShapeBase::calculateRelativePosInBoundingBox(double columnarBoundingBoxXMin, double columnarBoundingBoxYMin, double columnarBoundingBoxLength, double columnarBoundingBoxWidth, double peripodialBoundingBoxXMin, double peripodialBoundingBoxYMin, double peripodialBoundingBoxLength, double peripodialBoundingBoxWidth){
	columnarRelativePosInBoundingBox = getCentre();
	if (tissueType != 0){ //the tissue is not columnar, so there is peripodial membrane
		peripodialRelativePosInBoundingBox[0] = columnarRelativePosInBoundingBox[0];
		peripodialRelativePosInBoundingBox[1] = columnarRelativePosInBoundingBox[1];
		peripodialRelativePosInBoundingBox[0] = (peripodialRelativePosInBoundingBox[0] - peripodialBoundingBoxXMin) / peripodialBoundingBoxLength;
		peripodialRelativePosInBoundingBox[1] = (peripodialRelativePosInBoundingBox[1] - peripodialBoundingBoxYMin) / peripodialBoundingBoxWidth;
		columnarRelativePosInBoundingBox[0] = peripodialRelativePosInBoundingBox[0];
		columnarRelativePosInBoundingBox[1] = peripodialRelativePosInBoundingBox[1];
	}
	else{
		columnarRelativePosInBoundingBox[0] = (columnarRelativePosInBoundingBox[0] - columnarBoundingBoxXMin) / columnarBoundingBoxLength;
		columnarRelativePosInBoundingBox[1] = (columnarRelativePosInBoundingBox[1] - columnarBoundingBoxYMin) / columnarBoundingBoxWidth;
	}
	//cout<<"Element: "<<Id<<" RelPos: "<<columnarRelativePosInBoundingBox[0]<<" "<<columnarRelativePosInBoundingBox[1]<<" "<<peripodialRelativePosInBoundingBox[0]<<" "<<peripodialRelativePosInBoundingBox[1]<<endl;
	//double* a = new double[3];
	//a = getRelativePosInBoundingBox();
	//cout<<" a: "<< a[0]<<" "<<a[1]<<endl;
	//delete[] a;
}
*/
void 	ShapeBase::displayRelativePosInBoundingBox(){
		cout<<"Element: "<<Id<<"  relative position in the tissue bounding box: "<<relativePosInBoundingBox[0]<<" "<<relativePosInBoundingBox[1]<<endl;
}

bool 	ShapeBase::checkPackingToThisNodeViaState(int ColumnarLayerDiscretisationLayers, Node* NodePointer){
	if(IsAblated){
		//if the element is ablated, do not pack against it
		return false;
	}
	if (tissueType == 2){
		//the element is on the lateral section of the tissue, linking peripodial to columnar, no packing on this element
		return false;
	}
	if(ColumnarLayerDiscretisationLayers>1){
		//If the columnar layer is discretised into multiple layers, the apical elements should be checked against apical nodes,
		// and basal nodes should be checked against basal elements. The midline elements should not have packing, BUT on  a single layer tissue, all is midline, therefore
		// this check would not be valid.
		if (tissuePlacement == 2){	//tissue placement of the element is midline in a multi-layered columnar layer, it should not pack to anything
			return false;
		}
		if (NodePointer->tissuePlacement != tissuePlacement == 1){
			//apical nodes pack to apical elements and basal nodes pack to basal elements only
			return false;
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
	if (c>-0.99998){
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

void	ShapeBase::rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat){
	double x = gsl_matrix_get(rotMat,0,0)*u[0]+gsl_matrix_get(rotMat,0,1)*u[1]+gsl_matrix_get(rotMat,0,2)*u[2];
	double y = gsl_matrix_get(rotMat,1,0)*u[0]+gsl_matrix_get(rotMat,1,1)*u[1]+gsl_matrix_get(rotMat,1,2)*u[2];
	double z = gsl_matrix_get(rotMat,2,0)*u[0]+gsl_matrix_get(rotMat,2,1)*u[1]+gsl_matrix_get(rotMat,2,2)*u[2];
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

void	ShapeBase::calculateForces(vector <Node*>& Nodes, gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile){
    if (ShapeDim == 3){		//3D element
        calculateForces3D(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile);
    }
}

void ShapeBase::writeInternalForcesTogeAndgv(gsl_matrix* ge, gsl_matrix* gvInternal, double** SystemForces, vector <Node*>& Nodes){
    //now all the forces are written on SysyemForces
    //Now I will add the forces into ge, this step can be made faster by separating calculate forces function into two,
    //and filling up either ge or System forces depending on the solution method:
	for (int i = 0; i< nNodes; ++i){
        for ( int j=0; j<nDim; j++){
        	int indexI = nDim*NodeIds[i]+j;
        	double elementalvalue = gsl_matrix_get(ElementalElasticSystemForces,i,j);
        	double matrixValue = gsl_matrix_get(ge,indexI,0);
            gsl_matrix_set(ge, indexI,0,matrixValue + elementalvalue);
        	elementalvalue = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
        	matrixValue = gsl_matrix_get(gvInternal,indexI,0);
            gsl_matrix_set(gvInternal, indexI,0,matrixValue + elementalvalue);
        }
    }
    int counter = 0;
    for (int i = 0; i<nNodes; ++i){
        for (int j = 0; j<nDim; ++j){
            if (!Nodes[NodeIds[i]]->FixedPos[j]){
                SystemForces[NodeIds[i]][j] = SystemForces[NodeIds[i]][j] + gsl_matrix_get(ElementalElasticSystemForces,i,j) + gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
            }
            /*else if(recordForcesOnFixedNodes){
                FixedNodeForces[NodeIds[i]][j] = FixedNodeForces[NodeIds[i]][j] - gsl_matrix_get(TriPointg,counter,0);
            }*/
            counter++;
        }
    }
}



void	ShapeBase::calculateForces3D(vector <Node*>& Nodes,  gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile){
    int dim = nDim;
    int n = nNodes;
    //calculating F and B in a 3 point gaussian:

    gsl_matrix* TriPointge  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* TriPointgv  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix_set_zero(TriPointF);
    gsl_matrix_set_zero(ElementalElasticSystemForces);
    gsl_matrix_set_zero(ElementalInternalViscousSystemForces);
    gsl_matrix* currge = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currgv = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currF = gsl_matrix_calloc(dim,dim);
    //The point order is established in shape function derivative calculation!
    //Make sure the weights fir in with the order - eta zeta nu:
    //double points[3][3]={{1.0/6.0,1.0/6.0,0.0},{2.0/3.0,1.0/6.0,0.0},{1.0/6.0,2.0/3.0,0.0}};
    double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
    for (int iter =0; iter<3;++iter){
        //cout<<"Calculating Gauss point: "<<eta<<" "<<nu<<" "<<zeta<<endl;
    	calculateCurrNodalForces(currge, currgv, currF, displacementPerDt, iter);
        gsl_matrix_scale(currge,weights[iter]);
        gsl_matrix_add(TriPointge, currge);
        gsl_matrix_scale(currgv,weights[iter]);
        gsl_matrix_add(TriPointgv, currgv);
        gsl_matrix_scale(currF,weights[iter]);
        gsl_matrix_add(TriPointF, currF);
        //displayMatrix(currg,"currg");
        //if (Id == 0) {displayMatrix(currgv,"currgv");}
    }
    int counter = 0;
    for (int i = 0; i<nNodes; ++i){
            for (int j = 0; j<nDim; ++j){
            	if (!Nodes[NodeIds[i]]->FixedPos[j]){
            		double value = gsl_matrix_get(ElementalElasticSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointge,counter,0);
            		gsl_matrix_set(ElementalElasticSystemForces,i,j,value);
            		value = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointgv,counter,0);
            		gsl_matrix_set(ElementalInternalViscousSystemForces,i,j,value);
				}
				/*else if(recordForcesOnFixedNodes){
					FixedNodeForces[NodeIds[i]][j] = FixedNodeForces[NodeIds[i]][j] - gsl_matrix_get(TriPointg,counter,0);
				}*/
				counter++;
            }
    }
    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(ElementalElasticSystemForces,"ElementalElasticSystemForces");
/*
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
*/
    //freeing matrices allocated in this function
    gsl_matrix_free(TriPointge);
    gsl_matrix_free(TriPointgv);
    gsl_matrix_free(currge);
    gsl_matrix_free(currgv);
    gsl_matrix_free(currF);
    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(Fg, "Fg");
}


gsl_matrix* ShapeBase::calculateCauchyGreenDeformationTensor(gsl_matrix* Fe, gsl_matrix* FeT){
	//calculating C (C = (Fe^T*Fe):
	gsl_matrix * C =  gsl_matrix_alloc(nDim, nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FeT, Fe,0.0, C);
	return C;
}

gsl_matrix* ShapeBase::calculateEForNodalForcesKirshoff(gsl_matrix* C){
    //calculating E ( E = 1/2 *(Fe^T*Fe-I) ; E = 1/2 *(C-I):):
    gsl_matrix * E =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(E,C);
    gsl_matrix * I = gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set_identity(I);
    gsl_matrix_sub(E,I);
    gsl_matrix_scale(E, 0.5);
    gsl_matrix_free(I);
    return E;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesKirshoff(gsl_matrix* E){
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

    gsl_matrix_free(compactS);
    return S;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesNeoHookean(gsl_matrix* invC, double lnJ){
	//S = mu (I - C^-1) + lambda (lnJ) C^-1
	gsl_matrix * S =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(S,invC);
	//displayMatrix(S,"S1");
	gsl_matrix * I = gsl_matrix_alloc(nDim, nDim);
	gsl_matrix_set_identity(I);
    gsl_matrix_sub(I,invC);  //(I - C^-1)
    gsl_matrix_scale(I, mu); // mu (I - C^-1)
    gsl_matrix_scale(S, lambda*lnJ); //lambda (lnJ) C^-1
    gsl_matrix_add(S,I); // mu (I - C^-1) + lambda (lnJ) C^-1
	return S;
}

void ShapeBase::updateLagrangianElasticityTensorNeoHookean(gsl_matrix* invC, double lnJ, int pointNo){
	//calculating 4th order tensor C, for convenience the matrix D81 is used for both Kirshoff materials and neo-Hookean materials in the code.
	//The documentation lists  Lagrangian Elasticity Tensor with C for neo-Hookean, and with D for Kirshoff materials.
    //lambda is Lame s first parameter and mu is the shear modulus .
	double multiplier = 2*(mu - lambda*lnJ);
	for (int I = 0; I<nDim; ++I){
        for (int J = 0; J<nDim; ++J){
            for (int K = 0; K<nDim; ++K){
                for (int L = 0; L<nDim; ++L){
                	double Iijkl = 0.5* (gsl_matrix_get(invC,I,K)*gsl_matrix_get(invC,J,L) + gsl_matrix_get(invC,I,L)*gsl_matrix_get(invC,J,K));
                    D81[pointNo][I][J][K][L] = lambda*gsl_matrix_get(invC,I,J)*gsl_matrix_get(invC,K,L) + multiplier * Iijkl;
                	//D81[I][J][K][L] = lambda*gsl_matrix_get(invC,I,J)*gsl_matrix_get(invC,K,L) + multiplier * gsl_matrix_get(invC,I,K)*gsl_matrix_get(invC,J,L);
                }
            }
        }
    }
}

gsl_matrix* ShapeBase::calculateCompactStressForNodalForces(double detFe, gsl_matrix* Fe, gsl_matrix* S, gsl_matrix* FeT, gsl_matrix* Stress){
    //calculating stress (stress = detFe^-1 Fe S Fe^T):
    //double detFe = determinant3by3Matrix(Fe);
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

gsl_matrix* ShapeBase::calculateVelocityGradientTensor(gsl_matrix* B, gsl_matrix* displacementPerDt){
	/**
	 * Inputs:
	 * -# The elemental B matrix (6 , ShapeBase#nDim x ShapeBase#nNodes).
	 * -# The displacement of all nodes of the system, divided by the time
	 * step (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * Output:
	 * -# Velocity gradient tensor in Voigt notation (6, 1).
	 *
	 * This function calculates the velocity gradient tensor from elemental B matrix and elemental displacement.
	 * The elemental B matrix is composed of a stack of B matrices for each node of the element:
	 * 	\f{eqnarray*}{
        	\textbf{B}  &=& \left[ \left[ \textbf{B}_{0} \right] \left[ \textbf{B}_{1} \right] ... \left[ \textbf{B}_{nNode}\right] \right]\\
						&=& \left[
		\begin{bmatrix}
			\partial_x N_0 	& 0 				& 0	\\
			0 				& \partial_y N_0  	& 0	\\
			0 				& 0 				& \partial_z N_0 \\
			\partial_y N_0 	& \partial_x N_0 	& 0\\
			\partial_z N_0 	& 0 				& \partial_x N_0 \\
			0 				& \partial_z N_0 	& \partial_y N_0
		\end{bmatrix}

		\begin{bmatrix}
			\partial_x N_1 	& 0 				& 0	\\
			0 				& \partial_y N_1  	& 0	\\
			0 				& 0 				& \partial_z N_1 \\
			\partial_y N_1 	& \partial_x N_1 	& 0\\
			\partial_z N_1 	& 0 				& \partial_x N_1 \\
			0 				& \partial_z N_1 	& \partial_y N_1
		\end{bmatrix}

		...

		\begin{bmatrix}
			\partial_x N_{nNode} 	& 0 					& 0	\\
			0 						& \partial_y N_{nNode}  & 0	\\
			0 						& 0 					& \partial_z N_{nNode} \\
			\partial_y N_{nNode} 	& \partial_x N_{nNode} 	& 0\\
			\partial_z N_{nNode} 	& 0 					& \partial_x N_{nNode} \\
			0 						& \partial_z N_{nNode} 	& \partial_y N_{nNode}
		\end{bmatrix}
		\right]
		\f}
	 * The elemental displacement matrix is extracted from the system displacement matrix via
	 * the function ShapeBase#constructElementalDisplacementMatrix. The displacement is calculated
	 * as the displacement of a node from its position at the end of last time step, \f$ u_{n}\f$ to the position
	 * at the current Newton-Raphson iteration \f$ u_{k}\f$. With the velocities (displacement per time step),
	 * and the \f$\textbf{B}\f$ matrix, velocity gradient tensor can be calculated through:
	 * \f{eqnarray*}{
	 	 	 	 \boldsymbol{l} & = \boldsymbol{B}  \boldsymbol{v_{n+1}}\nonumber \\
								& = \boldsymbol{B} \frac{{u_{n+1}^{k} - u_{n}}} {\delta t}.
		\f}
	 *
	 * Procedure:
	 * - construct the ElementalDisplacementMatrix.
	 */
	gsl_matrix* elementalDisplacementMatrix = constructElementalDisplacementMatrix(displacementPerDt);
	/**
	 * - Allocate the velocity gradient tensor in Voigt notation.
	 */
	gsl_matrix* l =  gsl_matrix_calloc(6,1);
	/**
	 * - calculate velocity gradient tensor.
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, B, elementalDisplacementMatrix,0.0, l);
    /**
	 * - free allocated memory.
	 */
	gsl_matrix_free(elementalDisplacementMatrix);
	//displayMatrix(l,"l_forForceCalc");
    /**
	 * - return velocity gradient tensor.
	 */
	return l;
}

gsl_matrix* ShapeBase::calculateRateOfDeformationTensor(gsl_matrix* l){
	/**
	 * Inputs:
	 * -# The velocity gradient tensor given in Voigt notation (6 x 1).
	 *
	 * Output:
	 * -# Rate of deformation tensor (3 x 3).
	 *
	 * This function primarily rearranges the elements of the
	 * velocity gradient tensor given in Voigt notation, to from the rate of deformation tensor
	 *
	 * Procedure:
	 * - Allocate the memory for rate of deformation tensor
	 */
	gsl_matrix* d =  gsl_matrix_calloc(3,3);
	/**
	 * - Write the terms of velocity gradient tensor into rate of deformation tensor
	 */
	gsl_matrix_set(d,0,0,gsl_matrix_get(l,0,0) );
	gsl_matrix_set(d,1,1,gsl_matrix_get(l,1,0) );
	gsl_matrix_set(d,2,2,gsl_matrix_get(l,2,0) );
	gsl_matrix_set(d,0,1,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,2,1,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,0,2,0.5*gsl_matrix_get(l,5,0));
	gsl_matrix_set(d,1,0,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,1,2,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,2,0,0.5*gsl_matrix_get(l,5,0));
	/**
	 * - Return rate of deformation tensor
	 */
	return d;
}

void ShapeBase::calculateViscousStress(gsl_matrix* d, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# The rate of deformation matrix (ShapeBase#nDim x ShapeBase#nDim).
	 * -# the viscous stress of current Gauss point, the result will be written on this matrix
	 *
	 * This function will calculate the internal viscous stress of the element using rate of deformation matrix
	 * and ShapeBase#internalViscosity, \f$\eta\f$ , via:
	 * \f[\sigma^{v} = \eta  \textbf{d} \f]
	 *
	 * Procedure:
	 * - Copy rate of deformation tensor over to viscous stress tensor.
	 *
	 */
	createMatrixCopy(viscousStress, d);
	/**
	 *  - Scale with the internal viscosity to obtain viscous stress tensor
	 *
	 */
	gsl_matrix_scale(viscousStress,internalViscosity);
}

gsl_matrix* ShapeBase::constructElementalDisplacementMatrix(gsl_matrix* displacement){
	/**
	 * Inputs:
	 * -# The displacement matrix (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * This function calculates the elemental displacement matrix from the
	 * displacement matrix of the whole system, given as input. In
	 * current usage, under normal circumstances, the input matrix is displacement
	 * divided by time step. The displacement is calculated by Simulation#calculateDisplacementMatrix
	 * Both matrices, the displacement matrix of the whole system and the
	 * elemental displacement matrix are in vector form:
	 *
   	    \f$ displacement =
   	    \begin{bmatrix}
			\Delta x_{0}\\
			\Delta y_{0}\\
			\Delta z_{0}\\
			... ,\\
			\Delta x_{N}\\
			\Delta y_{N}\\
			\Delta z_{N}
			\end{bmatrix}
   	    \f$

	 */
	gsl_matrix* elementalDisplacementMatrix = gsl_matrix_calloc(nDim*nNodes,1);
	for (int i=0; i<nNodes; ++i){
		int index = NodeIds[i];
		for (int j=0; j<nDim; ++j){
			double value = gsl_matrix_get(displacement,index*nDim+j,0);
			gsl_matrix_set(elementalDisplacementMatrix,i*nDim+j,0,value);
		}
	}
	return elementalDisplacementMatrix;
}

void ShapeBase::calculateViscousForces(gsl_matrix*  gv, gsl_matrix*  BTdetFdetdXde, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# Elemental matrix for internal viscous forces (ShapeBase#nDim x ShapeBase#nNodes, 1)
	 * the resulting forces will be written on this matrix
	 * -# Transpose of elemental B matrix, multiplied by the determinant
	 * of the deformation gradient, \f$\textbf{F}\f$, and the determinant of \f$ \delta \textbf{X}/\delta \boldsymbol{\xi}\f$
	 * -#  The viscous stresses calculated in ShapeBase#calculateViscousStress
	 *
	 * This function will calculate the elemental viscous forces from viscous stress, via:
	 * 	\f{eqnarray*}{
        	\textbf{g}^v &=& \int_{V} \textbf{B}^{T} \sigma^{v} dV \\
          	  	  	  	 &=& det(\textbf{F})det\left( \frac{\delta \textbf{X} }{\delta \boldsymbol{\xi}} \right) \textbf{B}^{T} \sigma^{v}
		\f}
	 * Procedure:
	 * - Allocate the memory for stress in Voigt notation
	 * */
	gsl_matrix * compactStress =  gsl_matrix_calloc(6,1);
	/**
	 * - Write stress in Voigt notation
	 */
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(viscousStress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(viscousStress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(viscousStress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(viscousStress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(viscousStress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(viscousStress,0,2));
	/**
	 * - Calculate nodal forces
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BTdetFdetdXde, compactStress,0.0, gv);
	/**
	 * - Free memory
	 */
    gsl_matrix_free(compactStress);
}

void ShapeBase::calculateImplicitKElastic(){
    //cout<<"calculating implicit K elastic for element: "<<Id<<endl;
    int dim = nDim;
    int n = nNodes;
    if (IsAblated){
    	gsl_matrix_set_zero(TriPointKe);
    }
    else{
    	//calculating Kelastic in a 3 point gaussian:
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
		/*if (Id == 0){
			displayMatrix(currK,"currK");

		}*/
	    gsl_matrix_free(currK);
    }
}

void ShapeBase::calculateImplicitKViscous(gsl_matrix* displacementPerDt, double dt){
	//This is a function called over each element
	//First function is the function to calculate the sum of all integrals,
	//each integral calculation will be listed below, individually.
	//The inputs are:
	//1) the displacement per dt for all nodes (uk-un)/dt,
	//2) and is the time step the object has access to all
	//the necessary matrices to carry out the calculation.
	int dim = nDim;
	int n = nNodes;
	if (IsAblated || internalViscosity == 0){
	    //This is for efficiency, I do not calculate if the
	    //current element is laser ablated
		gsl_matrix_set_zero(TriPointKv);
	}
	else{
		//calculating Kviscous in a 3 point gaussian:
		//assign a temporary matix, 18 x 18 for a prism (6 nodes, 3D).
		gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
		//set the temporary matrix to zero.
		gsl_matrix_set_zero(TriPointKv);
		//the weights of the gauss points
		double weights[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
		//define the matrices for velocity gradient, and the
		//term in parentheses in calculation of the first integral.
		//These will be calculated once per Gauss point for the element.
		gsl_matrix* velocityGradient = gsl_matrix_calloc(dim,dim);
		gsl_matrix* paranthesisTermForKv1 =  gsl_matrix_calloc(dim,dim);
	    //loop over Gauss points
		for (int iter =0; iter<3;++iter){
			//set the temporary matrix to zero.
			gsl_matrix_set_zero(currK);
			//set the velocity gradient to zero.
			gsl_matrix_set_zero(velocityGradient);
			//set the parentheses term for first integral to zero.
			gsl_matrix_set_zero(paranthesisTermForKv1);
			//calculate the velocity gradient:
			calculateVelocityGradient(velocityGradient, displacementPerDt, iter);
			//calculate ( I / dt - velocityGradient) for first term:
			gsl_matrix* paranthesisTermForKv1 = gsl_matrix_alloc(nDim,nDim);
			//set the term to identity:
			gsl_matrix_set_identity(paranthesisTermForKv1);
			//divide the term by dt to obtain I/dt:
			gsl_matrix_scale(paranthesisTermForKv1,1.0/dt);
			//substract velocity gradient to obtain ( I / dt - velocityGradient):
			gsl_matrix_sub(paranthesisTermForKv1,velocityGradient);
			/*if (Id == 0){
				displayMatrix(velocityGradient,"velocityGradient");
				displayMatrix(paranthesisTermForKv1,"paranthesisTermForKv1");
				displayMatrix(viscousStress[iter],"viscousStress");
			}*/
			//calculate the first integral:
			calculateViscousKIntegral1(currK, paranthesisTermForKv1, iter);
			//calculate the second integral:
			calculateViscousKIntegral2(currK, iter);
		    //scaling the resulting temporary matrix with Gauss point weight.
			gsl_matrix_scale(currK,weights[iter]);
		    //Adding the temporary matrix to the elemental Kviscous.
			gsl_matrix_add(TriPointKv, currK);
		}
		//free the memory allocated in this function
		gsl_matrix_free(currK);
		gsl_matrix_free(velocityGradient);
		gsl_matrix_free(paranthesisTermForKv1);
	}
};

void ShapeBase::writeKelasticToMainKatrix(gsl_matrix* K){
    for (int a=0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            int NodeId1 = NodeIds[a];
            int NodeId2 = NodeIds[b];
            NodeId1 *= nDim;
            NodeId2 *= nDim;
            for (int i=0; i<nDim; ++i){
                for (int j=0; j<nDim; ++j){
                    double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
					valueij	+= gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                    gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
                }
            }
        }
    }
}

void ShapeBase::writeKviscousToMainKatrix(gsl_matrix* K){
	if (internalViscosity != 0){
		for (int a=0; a<nNodes; ++a){
			for (int b=0; b<nNodes; ++b){
				int NodeId1 = NodeIds[a];
				int NodeId2 = NodeIds[b];
				NodeId1 *= nDim;
				NodeId2 *= nDim;
				for (int i=0; i<nDim; ++i){
					for (int j=0; j<nDim; ++j){
						double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
						valueij	+= gsl_matrix_get(TriPointKv,a*nDim+i,b*nDim+j);
						gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
					}
				}
			}
		}
    }
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
void	ShapeBase::calculateElasticKIntegral1(gsl_matrix* currElementalK,int pointNo){
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
                                            value += (gsl_matrix_get(Fe,i,I)*gsl_matrix_get(Fe,j,J)*gsl_matrix_get(Fe,k,K)*gsl_matrix_get(Fe,l,L)*D81[pointNo][I][J][K][L]*DNb[l]*DNa[j]);
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
                    double value = gsl_matrix_get(currElementalK,a*nDim+i, b*nDim+j);
                    value += gsl_matrix_get(Keab,i, j);
                    gsl_matrix_set(currElementalK,a*nDim+i, b*nDim+j,value);
                }
            }
            gsl_matrix_free(Keab);
        }
    }
}


void	ShapeBase::calculateElasticKIntegral2(gsl_matrix* currElementalK,int pointNo){
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
                double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value;//adding the calculated value to current K matirx
                gsl_matrix_set(currElementalK,index1,index2,addedValue);
            }
            gsl_matrix_free(tmp1);
        }
    }
    gsl_matrix_free(DNaT);
    gsl_matrix_free(DNb);
    gsl_matrix_free(Keab2);
}

void ShapeBase::calculateVelocityGradient( gsl_matrix* velocityGradient, gsl_matrix* displacementPerDt, int pointNo){
	gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
	for (int c=0; c<nNodes; ++c){
		gsl_matrix* delVc = gsl_matrix_calloc(nDim,nDim);
		//get \DelNc^T
		gsl_matrix * DNc = gsl_matrix_calloc(nDim,1);
		for (int i=0;i<nDim;++i){
			gsl_matrix_set(DNc,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*c));
		}
		//calculate velocity of node c:
		gsl_matrix* vc = gsl_matrix_calloc(nDim,1);
		int id = NodeIds[c];
		gsl_matrix_set(vc,0,0,gsl_matrix_get(displacementPerDt,id*nDim,0));
		gsl_matrix_set(vc,1,0,gsl_matrix_get(displacementPerDt,id*nDim+1,0));
		gsl_matrix_set(vc,2,0,gsl_matrix_get(displacementPerDt,id*nDim+2,0));
		//calculate vc *DNc^T:
		calculateOuterProduct(vc,DNc,delVc);
		//add the nodal calculation to velocity gradient:
		gsl_matrix_add(velocityGradient, delVc);
		//free memory allocated in this loop:
		gsl_matrix_free(delVc);
		gsl_matrix_free(DNc);
		gsl_matrix_free(vc);
	}
	/*if (Id == 0 ){
		displayMatrix(velocityGradient,"DelV");
	}*/
}

void ShapeBase::calculateViscousKIntegral1(gsl_matrix* currElementalK, gsl_matrix* paranthesisTermForKv1, int pointNo){
//	gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
	gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
    gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
    gsl_matrix* BaTBb = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* KvabInt1 = gsl_matrix_calloc(nDim,nDim); //First integral in calculation of Kv for nodes a and b
    gsl_matrix* B = Bmatrices[pointNo];
    for (int a=0;a<nNodes;++a){
    	for (int b=0; b<nNodes; ++b){
    		consturctBaTBb(B, BaT,Bb,a,b);
    		//Bb matrix should have the last three rows as 0.5, this matrix stems from the
    		//"d" matrix, as opposed to viscous stresses, therefore the definition
    		//should have 0.5 on off diagonal terms, therefore the last three rows of Bb.
    		for (int i=3; i<6; ++i){
    			for (int j=0; j<nDim; ++j){
    				double value = gsl_matrix_get(Bb, i,j);
    				gsl_matrix_set(Bb,i,j,0.5*value);
    			}
    		}
    		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Bb,0.0, BaTBb);
    		//if (a == 3 && b ==3){
    		//	displayMatrix(BaTBb,"BaTBb_33");
    		//}
    		//the paranthesis term is: ( I / dt - velocityGradient)
    		//calculate all multiplication: BaT*Bb* (I/dt - \Del v):
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaTBb, paranthesisTermForKv1,0.0, KvabInt1);
			/*if (Id == 0 && a == 3 && b ==3){
				displayMatrix(BaTBb,"BaTBb_33");
				displayMatrix(KvabInt1,"KvabInt1_beforeVolumeIntegraiton_33");
				cout<<"detF: "<<detFs[pointNo]<<" detdXde: "<<detdXdes[pointNo]<<endl;
			}*/
			//volume integration:
    	    gsl_matrix_scale(KvabInt1,detFs[pointNo]);
    	    gsl_matrix_scale(KvabInt1,detdXdes[pointNo]);
    	    //scaling by viscosity:
    	    gsl_matrix_scale(KvabInt1,internalViscosity);
    	    /*if (Id == 0 && a == 3 && b ==3){
    	    	displayMatrix(KvabInt1,"KvabInt1_afterVolumeIntegraiton_33");
    	    }*/
    	    //Now KabvInt1 is a 3x3 (dim x dim) matrix, while currK is 18 x 18 (dim*n_node x dim*n_nodes)
    	    //currK is composed of 3x3 Kab matrices placed into the blocks (a,b), with indexing from zero,
    	    //for a = 1 and b=2, Kab will be placed in the 3x3 block covering indices (3,4,5) x (6,7,8) on K matrix.
			/*if (a == 3 && b == 3){
				cout<<"a: "<<a<<" b: "<<b<<endl;
				displayMatrix(KvabInt1,"KvabInt1");
			}*/
    	    for (int i=0; i<nDim; ++i){
        	    for (int j=0; j<nDim; ++j){
        	    	int index2 = a*nDim+i;
        	    	int index1 = b*nDim+j;
        	    	double value = gsl_matrix_get(KvabInt1,i,j);
        	    	double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
        	    	gsl_matrix_set(currElementalK,index1,index2,addedValue);
        	    }
			}
    	}
    }
    //if (Id ==0 ){displayMatrix(KvabInt1,"KvabInt1");}
    gsl_matrix_free(BaT);
    gsl_matrix_free(Bb);
    gsl_matrix_free(BaTBb);
    gsl_matrix_free(KvabInt1);
}

void ShapeBase::calculateViscousKIntegral2(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix * invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix * Stress = viscousStress[pointNo];
    gsl_matrix * DNa = gsl_matrix_calloc(nDim,1);
    gsl_matrix * DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix * KvabInt2 = gsl_matrix_calloc(nDim,nDim);
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            for (int i=0;i<nDim;++i){
                gsl_matrix_set(DNa,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix* DNaDNbOuterProduct = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix* DNbDNaOuterProduct = gsl_matrix_calloc(nDim, nDim);
            calculateOuterProduct(DNa, DNb, DNaDNbOuterProduct);
            calculateOuterProduct(DNb, DNa, DNbDNaOuterProduct);
            gsl_matrix* paranthesisTerm = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix_memcpy(paranthesisTerm,DNaDNbOuterProduct);
            gsl_matrix_sub(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_add(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_scale(paranthesisTerm, 0.5);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Stress, paranthesisTerm,0.0, KvabInt2);
			/*if (Id == 0){
				//displayMatrix(paranthesisTerm,"paranthesisTerm");
				displayMatrix(Stress,"Stress");
				//displayMatrix(KvabInt2,"KvabInt2-beforeVolumeIntegration");
			}*/
			gsl_matrix_scale(KvabInt2, detFs[pointNo]);
			gsl_matrix_scale(KvabInt2, detdXdes[pointNo]);
			/*if (a == 3 && b == 3){
				cout<<"a: "<<a<<" b: "<<b<<endl;
				displayMatrix(DNaDNbOuterProduct,"DNaDNbOuterProduct");
				displayMatrix(DNbDNaOuterProduct,"DNbDNaOuterProduct");
				displayMatrix(paranthesisTerm,"paranthesisTerm");
				displayMatrix(Stress,"Stress");
				cout<<"detFs[pointNo]: "<<detFs[pointNo]<<endl;
				cout<<"detdXdes[pointNo]: "<<detdXdes[pointNo]<<endl;
				displayMatrix(KvabInt2,"KvabInt2");
			}*/
			for (int i=0; i<nDim; ++i){
				for (int j=0; j<nDim; ++j){
					int index2 = a*nDim+i;
					int index1 = b*nDim+j;
					double value = gsl_matrix_get(KvabInt2,i,j);
					double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
					gsl_matrix_set(currElementalK,index1,index2,addedValue);
				}
			}
            gsl_matrix_free(DNaDNbOuterProduct);
            gsl_matrix_free(DNbDNaOuterProduct);
        }
    }
    gsl_matrix_free(DNa);
    gsl_matrix_free(DNb);
    gsl_matrix_free(KvabInt2);
	/*if (Id == 0){
		displayMatrix(DNa,"DNa");
		displayMatrix(DNaT,"DNaT");
		displayMatrix(DNb,"DNb");
		displayMatrix(DNbT,"DNbT");
		displayMatrix(DNaDNbT,"DNaDNbT");
		displayMatrix(DNbDNaT,"DNbDNaT");
		//displayMatrix(paranthesisTerm,"paranthesisTerm");
		displayMatrix(Stress,"Stress");
		displayMatrix(KvabInt2,"KvabInt2-beforeVolumeIntegration");
	}*/
}

void	ShapeBase::calculateOuterProduct(gsl_matrix* a, gsl_matrix* b, gsl_matrix* outerProduct){
	//cout<<"inside outer product"<<endl;
	int size1 = a->size1;
	int size2 = a->size2;
	if ((int) b->size2 != size2){
		cerr<<"matrix dimension mismatch in outer product calculation"<<endl;
	}
	gsl_matrix * bT = gsl_matrix_calloc(size2,size1);
	gsl_matrix_transpose_memcpy (bT,b);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, a, bT,0.0, outerProduct);
	gsl_matrix_free(bT);
	//cout<<"finalised outer product"<<endl;
}

gsl_matrix*	ShapeBase::calculateSymmetricisedTensorProduct(gsl_matrix* a, gsl_matrix* b){
	int size1 = a->size1;
	int size2 = a->size2;
	if ((int) b->size1 != size1){
		cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<endl;
	}
	if ((int) b->size2 != size2){
		cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<endl;
	}
	//generating the transposes:
	//gsl_matrix * aT = gsl_matrix_calloc(size2,size1);
	//gsl_matrix * bT = gsl_matrix_calloc(size2,size1);
	//gsl_matrix_transpose_memcpy (aT,a);
	//gsl_matrix_transpose_memcpy (bT,b);
	//calculating individual outer products a x b = a bT
	gsl_matrix * abOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(a,b,abOuterProduct);
	//gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, a, bT,0.0, abOuterProduct);
	gsl_matrix * baOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(b,a,abOuterProduct);
	//gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, b, aT,0.0, baOuterProduct);
	//calculating the averaged contraction term:
	gsl_matrix * averagedContraction = gsl_matrix_calloc(size1,size1);
	gsl_matrix_add(averagedContraction,abOuterProduct);
	gsl_matrix_add(averagedContraction,baOuterProduct);
	gsl_matrix_scale(averagedContraction,0.5);
	gsl_matrix_free(abOuterProduct);
	gsl_matrix_free(baOuterProduct);
	return averagedContraction;
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

void	ShapeBase::updateReferencePositionsToCurentShape(){
	for (int i=0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
			ReferenceShape ->Positions[i][j] = Positions[i][j];
		}
	}
}

void 	ShapeBase::setGrowthRate(double dt, double rx, double ry, double rz){
	GrowthRate[0] = exp(rx*dt);
	GrowthRate[1] = exp(ry*dt);
	GrowthRate[2] = exp(rz*dt);
	//if (Id ==0) {displayMatrix(growthIncrement, "Element0growthIncrement_initialSetting");}
}

void 	ShapeBase::setGrowthRateExpFromInput(double x, double y, double z){
	GrowthRate[0] = x;
	GrowthRate[1] = y;
	GrowthRate[2] = z;
}

void 	ShapeBase::updateGrowthIncrementFromRate(){
	gsl_matrix_set_identity(growthIncrement);
	gsl_matrix_set(growthIncrement,0,0,GrowthRate[0]);
	gsl_matrix_set(growthIncrement,1,1,GrowthRate[1]);
	gsl_matrix_set(growthIncrement,2,2,GrowthRate[2]);
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
	//cout<<"Element: "<<Id<<" eq upon signal: "<<cMyoUniformEq[0]<<" "<<cMyoUniformEq[1]<<endl;
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

void 	ShapeBase::calculatePrincipalStrainAxesOnXYPlane(double& e1, double &e2, double& tet){
	//principal strains:
	//e1,e2 = (exx + eyy) /2  +- sqrt( ( (exx - eyy)/2 ) ^2 + exy ^2)
	//extension is taken to be positive, therefore the most extended axis will be e1.
	//angle of the strains (direction of e1):
	// tan (2*tetha) = (2 exy ) / ( exx - eyy )
	double exx = gsl_matrix_get(Strain,0,0);
	double eyy = gsl_matrix_get(Strain,1,0);
	double exy = gsl_matrix_get(Strain,3,0)/2.0;
	double difference = (exx - eyy)/2.0;
	//double sumTerm = (exx + eyy) /2.0 ;
	//double sqrootTerm = pow ( difference*difference +  exy*exy, 0.5);
	//e1 = sumTerm + sqrootTerm;
	//e2 = sumTerm - sqrootTerm;
	double tan2Tet = exy/difference;
	tet = atan2(exy,difference)/2.0;
	//tet = atan(tan2Tet)/2;
	//I can calculate e1 and e2, but I will not know which direction they
	//correspond to, I should instead, calculate the rotation from a quaternian,
	//and therefore obtain the main strain direction:
	double c = cos(tet);
	double s = sin(tet);
	gsl_matrix* Rot = gsl_matrix_calloc(2,2);
	gsl_matrix* RotT = gsl_matrix_calloc(2,2);
	gsl_matrix* tmp = gsl_matrix_calloc(2,2);
	gsl_matrix* newStrain = gsl_matrix_calloc(2,2);
	gsl_matrix_set(Rot,0,0,c);
	gsl_matrix_set(Rot,0,1,s);
	gsl_matrix_set(Rot,1,0,(-1.0)*s);
	gsl_matrix_set(Rot,1,1,c);
	gsl_matrix_set(RotT,0,0,c);
	gsl_matrix_set(RotT,1,0,s);
	gsl_matrix_set(RotT,0,1,(-1.0)*s);
	gsl_matrix_set(RotT,1,1,c);
	gsl_matrix_set(newStrain,0,0,exx);
	gsl_matrix_set(newStrain,0,1,exy);
	gsl_matrix_set(newStrain,1,0,exy);
	gsl_matrix_set(newStrain,1,1,eyy);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Rot, newStrain, 0.0, tmp);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp, RotT, 0.0, newStrain);
	e1 = gsl_matrix_get(newStrain,0,0);
	e2 = gsl_matrix_get(newStrain,1,1);
	//cout<<"Id: "<<Id<<" tan2Tet "<<tan2Tet<<" 2Tet "<<atan2(exy,difference)<<" tet: "<<tet;
	if (e2>e1){
		//the main strain is in the direction perpendicular to the calculated angle.
		double tmp = e1;
		e1 = e2;
		e2 = tmp;
		tet = tet + M_PI/2.0;
	}
	//cout<<"after correction tet: "<<tet<<endl;
	//displayMatrix(Strain, "strain");
	//displayMatrix(newStrain, "newStrain");
	gsl_matrix_free(Rot);
	gsl_matrix_free(RotT);
	gsl_matrix_free(newStrain);
	gsl_matrix_free(tmp);

}

bool	ShapeBase::checkIfXYPlaneStrainAboveThreshold(double thres){
	double exx = gsl_matrix_get(Strain,0,0);
	double eyy = gsl_matrix_get(Strain,1,0);
	double exy = gsl_matrix_get(Strain,3,0)/2.0;
	//cout<<"Element: "<<Id<<" exx, eyy, exy: "<<exx<<" "<<eyy<<" "<<exy<<endl;
	double thresSquared = thres*thres;
	double diff= exx - eyy;
	if (diff*diff > thresSquared){
		return true;
	}
	if (exy*exy > thresSquared){
		return true;
	}
	return false;
}


void	ShapeBase::updateEquilibriumMyoWithFeedbackFromZero(double MyosinFeedbackCap){
	if (tissueType == 0 ){//|| tissueType == 1){
		//the feedback is applied only in peripodial membrane or the disc proper.
		//Linker zones are not affected.
		if (spansWholeTissue || tissuePlacement == 1){
			//The feedback is only on the apical surface of the tissue:
			//This feedback does not focus on the aspect ratio of the tissue.
			//Only stretch in one direction can result in myosin feedback.
			//If the tissue is under compression, with a high aspect ratio, there will be no feedback.
			double e1 = 0.0, e2 = 0.0, tet = 0.0;
			bool calculatePrincipalStrainDirection = checkIfXYPlaneStrainAboveThreshold(1E-5);
			if (calculatePrincipalStrainDirection){
				calculatePrincipalStrainAxesOnXYPlane(e1, e2, tet);
				double lowThres = 0.05;
				double upThres = 0.2;
				if(e1<lowThres){
					cMyoUnipolarEq[0] = 0.0;
				}else {
					//calculate concentration:
					if (e1>upThres){
						cMyoUnipolarEq[0] = MyosinFeedbackCap;
					}
					else{
						cMyoUnipolarEq[0] = e1*1000.0;
					}
					//give the current direction:
					gsl_matrix_set(myoPolarityDir,0,0,cos(tet));
					gsl_matrix_set(myoPolarityDir,0,1,sin(tet));
					gsl_matrix_set(myoPolarityDir,0,2,0.0);
					//cout<<"Element: "<<Id<<" e1: "<<e1<<" e2 "<<e2 <<" tet: "<<tet<<" cmyo: "<<cMyoUnipolarEq[0]<<endl;
					//displayMatrix(Strain, "strain");
					//displayMatrix(myoPolarityDir, "myoPolarityDir");
				}
			}
		}
		/*double exx = gsl_matrix_get(Strain,0,0);
		double eyy = gsl_matrix_get(Strain,1,0);
		double feedbackStrain = 0.0;
		int directionIndex = 0; //by default the myosin response is in x
		if (eyy > exx){
			//if the stretch is higher in yy, then change the direction
			//and assign the strain to be used in feedback calculation.
			directionIndex = 1;
			feedbackStrain = eyy;
		}
		else{
			feedbackStrain = exx;
		}

		if(feedbackStrain<lowThres){
			cMyoUnipolarEq[0] = 0.0;
		}else if (feedbackStrain>upThres){
			cMyoUnipolarEq[0] = MyosinFeedbackCap;
			//reset direction:
			gsl_matrix_set(myoPolarityDir,0,0,0);
			gsl_matrix_set(myoPolarityDir,0,1,0.0);
			gsl_matrix_set(myoPolarityDir,0,2,0.0);
			//give the current direction:
			gsl_matrix_set(myoPolarityDir,0,directionIndex,1.0);
		}
		else{
			cMyoUnipolarEq[0] = feedbackStrain*1000.0;
			//reset direction:
			gsl_matrix_set(myoPolarityDir,0,0,0);
			gsl_matrix_set(myoPolarityDir,0,1,0.0);
			gsl_matrix_set(myoPolarityDir,0,2,0.0);
			//give the current direction:
			gsl_matrix_set(myoPolarityDir,0,directionIndex,1.0);
		}*/

	}
}

void	ShapeBase::updateEquilibriumMyoWithFeedbackFromFixedTotal(double totalMyosinLevel){
	if (tissueType == 0 ){//|| tissueType == 1){
		//the feedback is applied only in the disc proper.
		//Linker zones are not affected.
		if (spansWholeTissue || tissuePlacement == 1){
			//The feedback is only on the apical surface of the tissue:
			//This feedback does not focus on the aspect ratio of the tissue.
			//Only stretch in one direction can result in myosin feedback.
			//If the tissue is under compression, with a high aspect ratio, there will be no feedback.
			double e1 = 0.0, e2 = 0.0, tet = 0.0;
			bool calculatePrincipalStrainDirection = checkIfXYPlaneStrainAboveThreshold(1E-5);
			if (calculatePrincipalStrainDirection){
				calculatePrincipalStrainAxesOnXYPlane(e1, e2, tet);
				//The experimental curve between strain and the ratio of polarised to non-polar myosin:
				//Y = 0.6164* X + 1.038
				//X = strain ( 50% strain is 0.5, defined in engineering strain terms),
				//Y = ratio of polar/uniform myo
				//The strain in the model is Green strain, I need to do the conversion:
				//Exx = ((1+exx)^2 - 1);
				//Then the formulation is:
				// Y = 0.6164 * ( (2*Exx +1)^0.5 - 1 ) + 1.038
				//In the experiments, the strain is applied by the stretcher, in one direction only.
				//In the model, this should be the residual strain, as is the difference of the strains
				//of the major and minor axis, as long as, e1 is positive (there is stretch)
				if (e1>0){
					//there is stretch
					double eEffectiveGreen = e1 - e2; //eEffective (eEffective = e1 - e2), is the difference between strains of major and minor axes.
					double eEffectiveEngineering = pow(2*eEffectiveGreen +1 , 0.5) - 1;
					double ratioOfPolarToUniformMyo = 0.6164*eEffectiveEngineering + 1.038;
					cMyoUniformEq[0] = totalMyosinLevel / (1+ratioOfPolarToUniformMyo);
					cMyoUnipolarEq[0] = totalMyosinLevel - cMyoUniformEq[0];
					//give the current direction:
					gsl_matrix_set(myoPolarityDir,0,0,cos(tet));
					gsl_matrix_set(myoPolarityDir,0,1,sin(tet));
					gsl_matrix_set(myoPolarityDir,0,2,0.0);
					//cout<<"Element: "<<Id<<" e1: "<<e1<<" e2 "<<e2 <<" tet: "<<tet<<" cmyo: "<<cMyoUnipolarEq[0]<<endl;
					//displayMatrix(Strain, "strain");
					//displayMatrix(myoPolarityDir, "myoPolarityDir");
				}
			}
			else{
				cMyoUnipolarEq[0] = 0.0;
				cMyoUniformEq[0] = totalMyosinLevel;
			}
			//cout<<" element: "<<Id<< " e1, e2: "<<e1<<" "<<e2<<" uniform: "<<cMyoUniformEq[0]<<" polar: "<<cMyoUnipolarEq[0]<<endl;
		}

	}
}

void	ShapeBase::updateMyosinConcentration(double dt, double kMyo, bool thereIsMyosinFeedback, double MyosinFeedbackCap){
	double thresholdValue = 1E-8, thresholdFraction= 0.01;
	//the value of kMyo is taken form my thesis
	double currMyoDt[3] = {dt,dt*2.0,dt/2.0};
	double cFinal[3];
	//First value is the final value with the current time step,
	//second is with currTimeStep*2 and
	//third is with 0.5 currTimeStep;
	double cInitial, cEq;
	//0 for any polarity below
	if (thereIsMyosinFeedback){
		//updateEquilibriumMyoWithFeedbackFromZero(MyosinFeedbackCap); // this function adds in polarised myosin to the tissue, independent of the uniform myosin levels
		updateEquilibriumMyoWithFeedbackFromFixedTotal(MyosinFeedbackCap); //this function redistributes a total pool of myosin in between uniform and polarised myosin levels
	}
	for (int myoIter =0; myoIter<4; myoIter++){
		if (myoIter == 0){ //apical uniform
			cInitial = cMyoUniform[0];
			cEq = cMyoUniformEq[0];
		}
		else if (myoIter == 1){//basal uniform
			cInitial = cMyoUniform[1];
			cEq = cMyoUniformEq[1];
		}
		else if (myoIter == 2){//apical polar
			cInitial = cMyoUnipolar[0];
			cEq = cMyoUnipolarEq[0];
		}
		else if (myoIter == 3){//basal polar
			cInitial = cMyoUnipolar[1];
			cEq = cMyoUnipolarEq[1];
		}
		bool converged = false;
		while (!converged){
			int steps[3] = {dt/currMyoDt[0],dt/currMyoDt[1],dt/currMyoDt[2]};
			for (int j=0; j<3; ++j){
				cFinal[j] = cInitial;
				for (int i =0 ;i<steps[j]; ++i){
					cFinal[j] += (cEq - cFinal[j])*kMyo*currMyoDt[j];
				}
			}
			//check if the value of the current dt and half current dt are below the threshold:
			double diff = fabs((cFinal[0] - cFinal[2]));
			if ( diff < thresholdValue ){
				converged = true;
			}
			else if( fabs(diff / cFinal[2]) < thresholdFraction){
				converged = true;
			}
			else{
				currMyoDt[1] = currMyoDt[0];
				currMyoDt[0] = currMyoDt[2];
				currMyoDt[2] *= 0.5;
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
	//if (cMyoUnipolar[0] > 0 || cMyoUnipolar[1]>0){
		//cout<<" Element: "<<Id<<" unipolar myo:  "<<cMyoUnipolar[0]<<" "<<cMyoUnipolar[1]<<" eq: "<<cMyoUnipolarEq[0]<<" "<<cMyoUnipolarEq[1]<<" polarity dir: "<<endl;
		//displayMatrix(myoPolarityDir,"myoPolarityDir");
		//cout<<" Element: "<<Id<<" uniform myo:  "<<cMyoUniform[0]<<" "<<cMyoUniform[1]<<" eq: "<<cMyoUniformEq[0]<<" "<<cMyoUniformEq[1]<<endl;
	//}
}

void ShapeBase::adjustCMyosinFromSave(){
	//if (cMyoUnipolarEq[0] < 1E-5){cMyoUnipolar[0] = 0.0;}
	//if (cMyoUnipolarEq[1] < 1E-5){cMyoUnipolar[1] = 0.0;}
	cMyoUnipolar[0] = 0.0;
	cMyoUnipolar[1] = 0.0;
}

void 	ShapeBase::setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
	ShapeChangeRate[0] = x;
	ShapeChangeRate[1] = y;
	ShapeChangeRate[2] = z;
	ShapeChangeRate[3] = xy;
	ShapeChangeRate[4] = yz;
	ShapeChangeRate[5] = xz;
}

//R Fginc R^T version:
void ShapeBase::calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix){
	gsl_matrix_set_identity(resultingGrowthIncrement);
	gsl_matrix_set(resultingGrowthIncrement,0,0,exp(growthx*dt));
	gsl_matrix_set(resultingGrowthIncrement,1,1,exp(growthy*dt));
	gsl_matrix_set(resultingGrowthIncrement,2,2,exp(growthz*dt));
	gsl_matrix* RotT = gsl_matrix_calloc(3,3);
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);;
	gsl_matrix_transpose_memcpy(RotT,ShearAngleRotationMatrix);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShearAngleRotationMatrix, resultingGrowthIncrement, 0.0, temp);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, RotT, 0.0, resultingGrowthIncrement);
	gsl_matrix_free(RotT);
	gsl_matrix_free(temp);

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

bool 	ShapeBase::InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse){
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
		if(tissueType == 2){
			Nodes[NodeIds[i]]->hasLateralElementOwner = true;
		}
	}
}

void 	ShapeBase::removeMassFromNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
			Nodes[NodeIds[i]]->mass -= VolumePerNode;
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
	IsXSymmetricClippedInDisplay=false;
	IsYSymmetricClippedInDisplay=false;
	for (int j=0; j<nNodes; ++j){
		 if((-1.0)*Positions[j][0]>xClip){
			 IsXSymmetricClippedInDisplay = true;
		 }
		 if((-1.0)*Positions[j][1]<yClip){
			 IsYSymmetricClippedInDisplay = true;
		 }
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
