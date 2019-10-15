#include "Prism.h"
#include "ReferenceShapeBase.h"
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

Prism::Prism(int* tmpNodeIds, vector<Node*>& Nodes, int CurrId, bool thereIsPlasticDeformation){
	//cout<<"constructing prism"<<endl;
	nNodes = 6;
	nDim = 3;
	Id = CurrId;
	ShapeDim = 3;	//3D shape
	NodeIds = new int[6];
	IdentifierColour = new int[3];
	E = 10.0;
	v = 0.3;
    internalViscosity = 0;
    lambda = E*v /(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);
    stiffnessMultiplier = 1.0;
    minimumValueOfStiffnessMultiplier = 0.00001;
    maximumValueOfStiffnessMultiplier = 100000;

    D = gsl_matrix_calloc(6,6);
    MyoForce = new double*[6];
	GrowthRate = new double[3];
	ShapeChangeRate  = new double[6];
	//CurrGrowthStrainAddition = new double[6];

    ApicalNormalForPacking =  new double[3];
	BasalNormalForPacking =  new double[3];
	relativePosInBoundingBox = new double[3];
	initialRelativePosInBoundingBox = new double[3];
	initialRelativePositionInZ = 1.0;
	//columnarRelativePosInBoundingBox = new double[3];
	//peripodialRelativePosInBoundingBox = new double[3];

	for (int i=0; i<nNodes; ++i){
		MyoForce[i] = new double[3];
		for (int j=0; j<3; j++){
			MyoForce[i][j] = 0.0;
		}
	}
	cMyoUniform[0] = 0.0;
	cMyoUniform[1] = 0.0;
	cMyoUnipolar[0] = 0.0;
	cMyoUnipolar[1] = 0.0;
	cMyoUniformEq[0] = 0.0;
	cMyoUniformEq[1] = 0.0;
	cMyoUnipolarEq[0] = 0.0;
	cMyoUnipolarEq[1] = 0.0;

	apicalNormalCurrentShape = new double[nDim];
	myoPolarityDir = gsl_matrix_calloc(2,3);
	for (int i=0; i<3; ++i){
		GrowthRate[i] = 0;
		ShapeChangeRate[i] =0;
		ApicalNormalForPacking[i] = 0;
		BasalNormalForPacking[i] = 0;
		relativePosInBoundingBox[i] = 0;
		initialRelativePosInBoundingBox[i] = 0;
		apicalNormalCurrentShape[i] = 0;
		//columnarRelativePosInBoundingBox[i] =0;
		//peripodialRelativePosInBoundingBox[i] =0;
	}
	columnarGrowthWeight = 1.0;
	peripodialGrowthWeight = 0.0;
	//CurrShapeChangeStrainsUpToDate = false;
	//CurrGrowthStrainsUpToDate = false;
	//IsGrowing = false;
	isFlipped = false;
	IsChangingShape = false;
	//GrewInThePast = false;
	//ChangedShapeInThePast = false;
	IsAblated = false;
	atSymetricityBoundary = false;
	IsClippedInDisplay = false;
	IsXSymmetricClippedInDisplay = false;
	IsYSymmetricClippedInDisplay = false;
	capElement = false;
    rotatedGrowth = false;
    //rotatedGrowth_tethaZ = 0;
	setIdentificationColour();
	setShapeType("Prism");
	ReferenceShape = new ReferenceShapeBase("Prism",Id);

	readNodeIds(tmpNodeIds);
	setPositionMatrix(Nodes);
	setReferencePositionMatrix();
    setInitialEdgeLenghts();

	//setGrowthTemplateMatrix();
	setCoeffMat();
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);

	numberOfGaussPoints =6;
	const int nGauss = numberOfGaussPoints;
	gaussWeights = new double[nGauss];
	gaussPoints = new double*[nGauss];
    for (int i=0; i<numberOfGaussPoints; ++i){
    	gaussPoints[i] = new double [nDim];
    }
    if (numberOfGaussPoints == 3){
    	gaussPoints[0][0]= 1.0/6.0;
		gaussPoints[0][1]= 1.0/6.0;
		gaussPoints[0][2]= 0;
		gaussPoints[1][0]= 2.0/3.0;
		gaussPoints[1][1]= 1.0/6.0;
		gaussPoints[1][2]= 0;
		gaussPoints[2][0]= 1.0/6.0;
		gaussPoints[2][1]= 2.0/3.0;
		gaussPoints[2][2]= 0;
		gaussWeights[0] = 1.0/3.0;
		gaussWeights[1] = 1.0/3.0;
		gaussWeights[2] = 1.0/3.0;
    }
    if (numberOfGaussPoints == 6){
		double oneOverSqrtThree = 1.0/ pow(3,0.5);
		gaussPoints[0][0]= 1.0/6.0;
		gaussPoints[0][1]= 1.0/6.0;
		gaussPoints[0][2]= oneOverSqrtThree;
		gaussPoints[1][0]= 2.0/3.0;
		gaussPoints[1][1]= 1.0/6.0;
		gaussPoints[1][2]= oneOverSqrtThree;
		gaussPoints[2][0]= 1.0/6.0;
		gaussPoints[2][1]= 2.0/3.0;
		gaussPoints[2][2]= oneOverSqrtThree;
		gaussPoints[3][0]= 1.0/6.0;
		gaussPoints[3][1]= 1.0/6.0;
		gaussPoints[3][2]= -oneOverSqrtThree;
		gaussPoints[4][0]= 2.0/3.0;
		gaussPoints[4][1]= 1.0/6.0;
		gaussPoints[4][2]= -oneOverSqrtThree;
		gaussPoints[5][0]= 1.0/6.0;
		gaussPoints[5][1]= 2.0/3.0;
		gaussPoints[5][2]= -oneOverSqrtThree;
		gaussWeights[0] = 1.0/6.0;
		gaussWeights[1] = 1.0/6.0;
		gaussWeights[2] = 1.0/6.0;
		gaussWeights[3] = 1.0/6.0;
		gaussWeights[4] = 1.0/6.0;
		gaussWeights[5] = 1.0/6.0;
    }
    ShapeFuncDerivatives = new gsl_matrix*[nGauss];
    ShapeFuncDerStacks = new gsl_matrix*[nGauss];
    InvdXdes = new gsl_matrix*[nGauss];
    detdXdes = new double[nGauss];
    Bmatrices = new gsl_matrix*[nGauss];
    FeMatrices = new gsl_matrix*[nGauss];
    detFs = new double[nGauss];
    invJShapeFuncDerStack = new gsl_matrix*[nGauss];
    invJShapeFuncDerStackwithFe  = new gsl_matrix*[nGauss];
    elasticStress = new gsl_matrix*[nGauss];
    viscousStress = new gsl_matrix*[nGauss];
    D81 = new double****[nGauss];
    for (int i=0; i<numberOfGaussPoints; ++i){
        ShapeFuncDerivatives[i] = gsl_matrix_calloc(nDim, nNodes);
        ShapeFuncDerStacks[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        InvdXdes[i] = gsl_matrix_calloc(nDim, nDim);
        detdXdes[i] = 0.0;
        Bmatrices[i] = gsl_matrix_calloc(nNodes,nDim*nNodes);
        FeMatrices[i] = gsl_matrix_calloc(3,3);
        invJShapeFuncDerStack[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        invJShapeFuncDerStackwithFe[i] = gsl_matrix_calloc(nDim*nDim, nDim*nNodes);
        elasticStress[i] = gsl_matrix_calloc(3,3);
        viscousStress[i] = gsl_matrix_calloc(3,3);
        //need to reach:  D81[gauss points][4][4][4][4];
        D81[i] = new double***[nDim];
        for (int j=0; j<nDim; ++j){
        	 D81[i][j] = new double**[nDim];
        	 for (int k=0; k<nDim; ++k){
        		 D81[i][j][k] = new double*[nDim];
        		 for (int l=0; l<nDim; ++l){
        			 D81[i][j][k][l] = new double[nDim];
        		 }
			}
        }
    }
    Strain = gsl_matrix_calloc(6,1);
    GrowthStrainsRotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(GrowthStrainsRotMat);
    Fg = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fg);
    InvFg = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFg);
    Fsc = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(Fsc);
    InvFsc = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(InvFsc);



    growthIncrement = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(growthIncrement);
    plasticDeformationIncrement= gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(plasticDeformationIncrement);
    shapeChangeIncrement = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(shapeChangeIncrement);
    zRemodellingSoFar = 1.0;
    TriPointF = gsl_matrix_calloc(3,3);
	gsl_matrix_set_identity(TriPointF);
    TriPointKe = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
    TriPointKv = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
    ElementalElasticSystemForces = gsl_matrix_calloc(nNodes,nDim);
    ElementalInternalViscousSystemForces = gsl_matrix_calloc(nNodes,nDim);
	RotatedElement = false;    

	CurrShapeChangeToAdd[0] = 0;
	CurrShapeChangeToAdd[1] = 0;
	CurrShapeChangeToAdd[2] = 0;

	VolumePerNode = 0;
    ZProjectedBasalArea=0.0;
    ZProjectedApicalArea=0.0;
    BasalArea=0.0;
    ApicalArea=0.0;
    exposedLateralAreaApicalSide = 0;
    exposedLateralAreaApicalSide = 0;
    elementHasExposedApicalSurface = false;
    elementHasExposedBasalSurface = false;
    elementHasExposedLateralApicalSurface = false;
    elementHasExposedLateralBasalSurface = false;
    exposedApicalSurfaceNodeIds[0] = 0;
    exposedApicalSurfaceNodeIds[1] = 0;
    exposedApicalSurfaceNodeIds[2] = 0;
    exposedBasalSurfaceNodeIds[0] = 0;
    exposedBasalSurfaceNodeIds[1] = 0;
    exposedBasalSurfaceNodeIds[2] = 0;
    exposedLateralAreaApicalSideNodeIds[0] = 0;
    exposedLateralAreaApicalSideNodeIds[1] = 0;
    exposedLateralAreaApicalSideNodeIds[2] = 0;
    exposedLateralAreaApicalSideNodeIds[3] = 0;
    exposedLateralAreaBasalSideNodeIds[0] = 0;
    exposedLateralAreaBasalSideNodeIds[1] = 0;
    exposedLateralAreaBasalSideNodeIds[2] = 0;
    exposedLateralAreaBasalSideNodeIds[3] = 0;
    nLateralSurfaceAreaNodeNumber = 4;
    nSurfaceAreaNodeNumber = 3;

    cellsMigrating = false;
    remodellingPlaneRotationMatrix = gsl_matrix_calloc(3,3);
    gsl_matrix_set_identity(remodellingPlaneRotationMatrix);
    isECMMimicing = false;
    isECMMimimcingAtCircumference = false;
    atBasalBorderOfECM = false;
    isActinMimicing = false;
    atApicalBorderOfActin = false;
    insideEllipseBand = false;
    coveringEllipseBandId = -1;
    stiffnessPerturbationRateInSec = 0;
    ECMThicknessPlaneRotationalMatrix = gsl_matrix_calloc(3,3);
    basalNeigElementId = -1;

    isMutated = false;
    mutationGrowthRatePerSec = 0.0;

    thereIsGrowthRedistribution = false;
    growthRedistributionShrinksElement = false;
    growthRedistributionScale = 0.5;

    plasticDeformationHalfLifeMultiplier = 1.0;
}

Prism::~Prism(){
	for (int i=0; i<nNodes; ++i){
		delete[] Positions[i];
		delete[] MyoForce[i];
	}
	delete[] Positions;
	delete[] relativePosInBoundingBox;
	delete[] initialRelativePosInBoundingBox;
	delete[] NodeIds;
	delete[] IdentifierColour;
	delete[] MyoForce;
	delete[] GrowthRate;
	delete[] ShapeChangeRate;
	delete[] apicalNormalCurrentShape;
	delete	 ReferenceShape;

    //freeing matrices allocated
	gsl_matrix_free(growthIncrement);
	gsl_matrix_free(plasticDeformationIncrement);
	gsl_matrix_free(D);
	gsl_matrix_free(CoeffMat);
	gsl_matrix_free(Fg);
	gsl_matrix_free(InvFg);
	gsl_matrix_free(Fsc);
	gsl_matrix_free(InvFsc);
	gsl_matrix_free(TriPointF);
	gsl_matrix_free(Strain);
	gsl_matrix_free(TriPointKe);
	gsl_matrix_free(TriPointKv);
	gsl_matrix_free(GrowthStrainsRotMat);
    for (int i=0; i<numberOfGaussPoints; ++i){
    	delete[] gaussPoints[i];
        gsl_matrix_free (ShapeFuncDerivatives[i]);
        gsl_matrix_free (ShapeFuncDerStacks[i]);
        gsl_matrix_free (InvdXdes[i]);
        gsl_matrix_free (Bmatrices[i]);
        gsl_matrix_free (FeMatrices[i]);
        gsl_matrix_free (invJShapeFuncDerStack[i]);
        gsl_matrix_free (invJShapeFuncDerStackwithFe[i]);
        gsl_matrix_free (elasticStress[i]);
        gsl_matrix_free (viscousStress[i]);
        for (int j=0; j<nDim; ++j){
        	 for (int k=0; k<nDim; ++k){
        		 for (int l=0; l<nDim; ++l){
        			 delete[] D81[i][j][k][l];
        		 }
			}
        }
        for (int j=0; j<nDim; ++j){
        	 for (int k=0; k<nDim; ++k){
        		 delete[] D81[i][j][k];
			}
        }
        for (int j=0; j<nDim; ++j){
        	delete[] D81[i][j];
        }
        delete[] D81[i];
    }
    delete[] ShapeFuncDerivatives;
    delete[] ShapeFuncDerStacks;
    delete[] InvdXdes;
    delete[] detdXdes;
    delete[] Bmatrices;
    delete[] FeMatrices;
    delete[] detFs;
    delete[] invJShapeFuncDerStack;
    delete[] invJShapeFuncDerStackwithFe;
    delete[] elasticStress;
    delete[] viscousStress;
    delete[] D81;
    delete[] gaussPoints;
    delete[] gaussWeights;
    delete[] elementsIdsOnSameColumn;
    gsl_matrix_free(ECMThicknessPlaneRotationalMatrix);
}

void Prism::setCoeffMat(){
   /** The coefficient matrix collating the calculated stresses in the form of the vector
     * form of elemental strss and strain vectors:
     *
     * \f{eqnarray*}{
        \boldsymbol{c} =
        \begin{bmatrix}
        1	& 0	& 0	& 0	& 0	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 1	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 0	& 0	& 0	& 0	& 1	\\
        0	& 1	& 0	& 1	& 0	& 0	& 0	& 0	& 0	\\
        0	& 0	& 0	& 0	& 0	& 1	& 0	& 1	& 0	\\
        0	& 0	& 1	& 0	& 0	& 0	& 1	& 0	& 0	\\
        \end{bmatrix}
        \sim
        \begin{bmatrix}
        \epsilon_{x} \\
        \epsilon_{y} \\
        \epsilon_{z}\\
        \gamma_{x,y} \\
        \gamma_{x,z} \\
        \gamma_{y,z}
        \end{bmatrix}
        \sim
        \begin{bmatrix}
        \sigma_{x} \\
        \sigma_{y} \\
        \sigma_{z}\\
        \sigma_{x,y} \\
        \sigma_{x,z} \\
        \sigma_{y,z}
        \end{bmatrix}.
        \f}
      */
    CoeffMat = gsl_matrix_calloc(6, nDim*nDim);
    gsl_matrix_set(CoeffMat,0,0,1);
    gsl_matrix_set(CoeffMat,1,4,1);
    gsl_matrix_set(CoeffMat,2,8,1);
    gsl_matrix_set(CoeffMat,3,1,1);
    gsl_matrix_set(CoeffMat,3,3,1);
    gsl_matrix_set(CoeffMat,4,5,1);
    gsl_matrix_set(CoeffMat,4,7,1);
    gsl_matrix_set(CoeffMat,5,2,1);
    gsl_matrix_set(CoeffMat,5,6,1);
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
	delete[] vec1;
	delete[] vec2;
	delete[] view;
	delete[] normal;
}

void  Prism::calculateApicalNormalCurrentShape(){
	double * u = new double[3];
	double * v = new double[3];
	for (int i=0; i<nDim; ++i){
		u[i] = Positions[4][i] - Positions[3][i];
		v[i] = Positions[5][i] - Positions[3][i];
		apicalNormalCurrentShape[i] = 0.0;
	}
	crossProduct3D(u,v,apicalNormalCurrentShape);
	double dummy = normaliseVector3D(apicalNormalCurrentShape);
	for (int i=0; i<nDim; ++i){
		u[i] = Positions[0][i] - Positions[3][i];
	}
	double  dot = dotProduct3D(u,apicalNormalCurrentShape);
	if (dot<0){
		for (int i=0; i<nDim; ++i){
			apicalNormalCurrentShape[i] *=(-1.0);
		}
	}
	//cerr<<"		apical normal after direction correction: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
	delete[] v;
	delete[] u;
}

void  Prism::calculateBasalNormal(double * normal){
	double * u = new double[3];
	double * v = new double[3];
	for (int i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[1][i] - ReferenceShape->Positions[0][i];
		v[i] = ReferenceShape->Positions[2][i] - ReferenceShape->Positions[0][i];
		normal[i] = 0.0;
	}
	crossProduct3D(u,v,normal);
	double dummy = normaliseVector3D(normal);
	for (int i=0; i<nDim; ++i){
		u[i] = ReferenceShape->Positions[3][i] - ReferenceShape->Positions[0][i];
	}
	double  dot = dotProduct3D(u,normal);
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
	double * normal = new double[3];
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


void  Prism::setElasticProperties(double EApical, double EBasal, double EMid, double EECM, double v){
	/** This function takes the Young's moduli for apical and basal surfaces, the tissue midline, 
	 *  and the ECM. It also accepts the Poisson's ratio. The choice of the element being 
	 *  a peripodial or columnar layer is made before calling this function. /n
	 * 
	 * Initially, the element's Young's modulus is set to the tissue mid layer. First check is for ECM,
	 * if the element is representing ECM (ShapeBase#isECMMimicing), then the stiffness is set to 
	 * \f$ E_{ECM} \f$. If it is not an ECM element, the tissue placement is checked. 
	*/
	this -> E = EMid;
	if (isECMMimicing){
		this -> E = EECM;
	}
	else{
		/**  
		 * If the element is on the basal surface, or at the basal border of ECM, its stiffness is set
		 * to \f$ E_{basal} \f$. The double check is necessary, as in the case of a mesh with ECM, the 
		 * ECM elements will take up the basal flag (ShapeBase#tissuePlacement = 0). The actual cellular basal layer will be the first 
		 * layer of elements that are immediately above the ECM, that are flagged as being at 
		 * basal border of ECM (ShapeBase#atBasalBorderOfECM).
		 */
		if (tissuePlacement == 0 || atBasalBorderOfECM){
			this -> E = EBasal;
		}
		else if(tissuePlacement == 1 ){
		/**
		* If the element has an apical tissue placement (ShapeBase#tissuePlacement = 1), then 
		* its stiffness will be set to \f$ E_{apical} \f$. \n
		*/
			this -> E = EApical;
		}
	}
	/**
	* The poisson ratio is set to the user input, followed by a sanity check, such taht if the value is 
	* bounded to the interval [0.0, 0.5].
	*/
	this -> v = v; //poisson ratio
	if (this -> v>0.5){this -> v= 0.5;}
	else if (this -> v<0.0){this -> v = 0.0;}
	/**
	* Once the Young's modulus and the Poisson's ratio of the element are set, the Lame's parameters, 
	* ShapeBase#lambda and ShapeBase#mu, are calculated and the mechanical model is simulated over these parameters.
	* This detail highlights the importance of structured parameter updates. At any instance the 
	* the Young's modulus or Poisson's ratio of an element is altered during the simulation, Prism#updateElasticProperties
	* function must be called to update the Lame parameter, the shear modulus, then the Lagrangian elasticity tensor (ShapeBase#D81)
	* (the tensor for elastic force calculation, depending on the material model). 
	* Otherwise the changes will not be reflected in the mechanical model.
	*/
    lambda = E*v/(1+v)/(1-2.0*v);
    mu = E/2.0/(1+v);
    calculateDVector();
    calculateD81Tensor();
}


void  Prism::updateElasticProperties(){
    /**
     * This function updates the Lame's parameters and the related tensors for elastic stress calculation
     * with the current Young's moduli and Poisson's ratio. Of note, for ease of tracking the extent of change,
     * during modifications of stiffness, the initial Young's modulus is not altered, nad the modifications are 
     * stored in ShapeBase#stiffnessMultiplier.
     */
    lambda = stiffnessMultiplier*E*v/(1+v)/(1-2.0*v);
    mu = stiffnessMultiplier*E/2.0/(1+v);
    /** The two tensor updates are not necessary for a neo-Hookean material, as the calculation for
    * ShapeBase#D81 tensor is carried out at each Newton-Raphson iteration. At this point, the material model is
    * not known,  therefore the values will be updated. 
    */
    calculateDVector();
    calculateD81Tensor();
}

void Prism::fillLateralNeighbours(vector<Node*>& Nodes, vector<int>& lateralNeigbours){
	// to do: do i use this?
	//generate a list of element that are connected to this element
	vector <int> connectedElementList;
	vector <int> connectedElementEncounterCounter;
	for (int i=0; i<nNodes; ++i){
		int nConnectedElement = Nodes[NodeIds[i]]->connectedElementIds.size();
		for (int j=0; j<nConnectedElement; ++j){
			int currElementId = Nodes[NodeIds[i]]->connectedElementIds[j];
			if (currElementId != this->Id){
				//the connected element of the node is not this element itself
				int currListSize = connectedElementList.size();
				bool alreadyRecorded = false;
				for (int k=0; k<currListSize; ++k){
					if (connectedElementList[k] == currElementId){
						connectedElementEncounterCounter[k]++;
						alreadyRecorded = true;
						break;
					}
				}
				if (!alreadyRecorded){
					connectedElementList.push_back(currElementId);
					connectedElementEncounterCounter.push_back(1);
				}
			}
		}
	}
	//If the element id have been encoiuntered 4 times, then it is a lateral neighbour:
	int currListSize = connectedElementList.size();
	for (int k=0; k<currListSize; ++k){
		if (connectedElementEncounterCounter[k] == 4){
			lateralNeigbours.push_back(connectedElementList[k]);
		}
	}
}


void Prism::calculateDVector(){
	double multiplier = E/((1+v)*(1-2*v));
	gsl_matrix_set(D,0,0,  multiplier*(1-v));
	gsl_matrix_set(D,0,1,  multiplier*v);
	gsl_matrix_set(D,0,2,  multiplier*v);
	gsl_matrix_set(D,1,0,  multiplier*v);
	gsl_matrix_set(D,1,1,  multiplier*(1-v));
	gsl_matrix_set(D,1,2,  multiplier*v);
	gsl_matrix_set(D,2,0,  multiplier*v);
	gsl_matrix_set(D,2,1,  multiplier*v);
	gsl_matrix_set(D,2,2,  multiplier*(1-v));
	gsl_matrix_set(D,3,3,  multiplier*(1-2*v)/2);
	gsl_matrix_set(D,4,4,  multiplier*(1-2*v)/2);
	gsl_matrix_set(D,5,5,  multiplier*(1-2*v)/2);
}

void Prism::calculateD81Tensor(){
 /**
    * This function calculates the Lagrangian elasticity tensor for a neo-hookean material.
    * The Lagrangian elasticity tensor depends on the deformation gradient through the
    * Cauchy-Green Deformation tensor, and the physical material properties, namely Lame's second
    * parameter \f$ \lambda \f$ and the sheer modulus \f$ \mu \f$. \n
    * \f$\boldsymbol {\mathcal D} \f$ can be obtained in the Voigt notation:
    * \f{eqnarray*}{
            \mathcal{D}_{ijkl}  = \lambda \boldsymbol{C}^{-1}_{ij} \boldsymbol{C}^{-1}_{kl} + 2 \left( \mu - \lambda ln(J)\right) \mathcal{I}_{ijkl}
      \f}
    * where \f$ \mathcal{I}_{ijkl} \f$ is defined as:
      \f{eqnarray*}{
            \mathcal{I}_{ijkl}  = \frac{1}{2} \left[ \boldsymbol{C}^{-1}_{ik} \boldsymbol{C}^{-1}_{jl} + \boldsymbol{C}^{-1}_{il} \boldsymbol{C}^{-1}_{jk}\right]
      \f}
    * where the indices $ i,j,k,l $ go through 3 dimensions.\n
    *
    */
    	// lambda is Lame s first parameter and mu is the shear modulus .
	double Idouble[3][3] = {{1.0,0.0,0.0} , {0.0,1.0,0.0}, {0.0,0.0,1.0}};
	for (int pointNo = 0; pointNo<3; pointNo++){
		for (int I = 0; I<nDim; ++I){
			for (int J = 0; J<nDim; ++J){
				for (int K = 0; K<nDim; ++K){
					for (int L = 0; L<nDim; ++L){
						D81[pointNo][I][J][K][L] = lambda*Idouble[K][L]*Idouble[I][J] + mu * ( Idouble[I][K]*Idouble[J][L] + Idouble[I][L]*Idouble[J][K] );
					}
				}
			}
		}
	}
}

void Prism::getCurrRelaxedShape(gsl_matrix* CurrRelaxedShape){
    /** This function writes the current relaxed shape of the element in the gsl_matrix
     * pointed by the provided input gsl_matrix pointer. The matrix must be allocated before
     * the call to the function.
     *
     */
	for (int i =0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
            gsl_matrix_set(CurrRelaxedShape,i,j,ReferenceShape->Positions[i][j]);
		}
	}
}

void Prism::setShapeFunctionDerivatives(gsl_matrix* ShapeFuncDer, double eta, double zeta, double nu){
   /** This function calculates and writes the shape function derivatives on
     * the gsl_matrix pointed by the provided input gsl_matrix pointer,
     * for the provided barycentric coordinates. The matrix must be allocated prior to
     * calling the function. \n
     * For the prism, the shape functions and their derivatives are given as:
     *
     * \f{tabular}{{|c|c|c|c|c|}
            \multicolumn{5}{|c|}{ $\lambda = 1 - \xi - \eta$ }            \\
            \multicolumn{5}{|c|}{ $\alpha = \left( 1-\zeta \right) /2$  } \\
            \multicolumn{5}{|c|}{ $b = \left(1+\zeta\right)/2$ }          \\
            Node & Shape function.    & Shape function         	 & Shape function 	  	   & Shape function.         \\
                 & 			        & derivative wrt $\xi$      	 &  derivative wrt $\zeta$      & derivative wrt $\eta$ \\
                 & 		$N$ 	        & $\frac{\partial{N}}{ \partial \xi}$     & $\frac{\partial{N}}{\partial \zeta}$  &  $\frac{\partial{N}}{\partial \eta}$ \\
            1       & $\lambda \alpha$  & $- \alpha$         		 & $- \alpha$	  	   	    & $ \frac{-\lambda} {2}$        \\
            2	 & $\xi \alpha$           & $\alpha$     		 	& 0  	  	   			    & $\frac{-\xi} {2} $        \\
            3  	& $\eta \alpha$  	& 0         		 		& $\alpha$ 	  	   	    & $\frac{-\nu} {2} $        \\
            4    	& $\lambda b$  		&  $-b$      			& $-b$ 	  	   		    & $\frac{\lambda} {2} $         \\
            5	& $\xi b$  			& $b$         		 	&0 	  	   			    & $\frac{\xi} {2} $         \\
            6 	 & $\eta b$  		& 0    		 		& $b$   	  	   		    & $\frac{\eta} {2} $        \\
         \f}
     *
     * The shape function derivatives matrix has the form:
     * 	\f{eqnarray*}{
        \begin{bmatrix}
            \frac{\partial N_0}{ \partial \xi}      & ... 	& \frac{\partial N_n}{ \partial \xi}    \\
            \frac{\partial N_0}{ \partial \zeta}    & ...   & \frac{\partial N_n}{ \partial \zeta}  \\
            \frac{\partial N_0}{ \partial \eta}     & ...   & \frac{\partial N_n}{ \partial \eta}
        \end{bmatrix}
        \f}
     */
	double alpha  = (1 - zeta)/2;
	double beta = (1 + zeta)/2;
	double lambda = 1-eta-nu;

    gsl_matrix_set(ShapeFuncDer,0,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,0,1,  alpha);
    gsl_matrix_set(ShapeFuncDer,0,2,  0);
    gsl_matrix_set(ShapeFuncDer,0,3, -beta);
    gsl_matrix_set(ShapeFuncDer,0,4,  beta);
    gsl_matrix_set(ShapeFuncDer,0,5, 0);

    gsl_matrix_set(ShapeFuncDer,1,0, -alpha);
    gsl_matrix_set(ShapeFuncDer,1,1,  0);
    gsl_matrix_set(ShapeFuncDer,1,2,  alpha);
    gsl_matrix_set(ShapeFuncDer,1,3, -beta);
    gsl_matrix_set(ShapeFuncDer,1,4,  0);
    gsl_matrix_set(ShapeFuncDer,1,5,  beta);

    gsl_matrix_set(ShapeFuncDer,2,0, -lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,1, -eta/2);
    gsl_matrix_set(ShapeFuncDer,2,2, -nu/2);
    gsl_matrix_set(ShapeFuncDer,2,3,  lambda/2);
    gsl_matrix_set(ShapeFuncDer,2,4,  eta/2);
    gsl_matrix_set(ShapeFuncDer,2,5,  nu/2);
}

void Prism::setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack){
    /** This function writes the input gsl_matrix ShapeFuncDer into a stack form, in the
     * matrix ShapeFuncDerStack. Both matrices must be allocated prior to calling the function.
     * This stack is for direct multiplication with B matrix during calculation of stresses.
     * with the notation  \f$ N_{i, \xi} \f$ being the short hand for
     * \f$ \partial \boldsymbol{x_{i}} / \partial {\xi} \f$
     * for node i, the stack for of the shape function derivatives in open form is:
     *
     *
     *\f{eqnarray*}{
        \begin{bmatrix}
            N_{0,\xi} 		& 0 			& 0			& N_{1,\xi} 	& \cdots &N_{n-1,\xi} 		& 0 			& 0			\\
            N_{0,\zeta} 	& 0 			& 0 		& N_{1,\zeta}  	& \cdots &N_{n-1,\zeta}       & 0 			& 0			\\
            N_{0,\nu} 		& 0			& 0 			& N_{1,\nu}  	& \cdots &N_{n-1,\nu} 		& 0 			& 0			\\
            0 			& N_{0,\xi} 	& 0 			& 0 			& \cdots &0                 & N_{n-1,\xi} 	& 0			\\
            0 			& N_{0,\zeta} 	& 0 			& 0 	 		& \cdots &0                 & N_{n-1,\zeta} 	& 0			\\
            0 			& N_{0,\nu} 	& 0 			& 0  			& \cdots &0                 & N_{n-1,\nu} 	& 0			\\
            0 			& 0			& N_{0,\xi}         & 0  			& \cdots &0                 &0			& N_{n-1,\xi} 	\\
            0 			& 0			& N_{0,\zeta}       & 0 	 		& \cdots &0                 &0			& N_{n-1,\zeta}	\\
            0 			& 0			& N_{0,\nu}         & 0   			& \cdots &0                 &0			& N_{n-1,\nu} 	.
        \end{bmatrix}
        \f}
     *
     */	
     	int n = nNodes;
	int dim = nDim;
	for (int i=0; i<n;++i){
        for (int k=0; k<dim; ++k){
            gsl_matrix_set(ShapeFuncDerStack,k,i*dim, gsl_matrix_get(ShapeFuncDer,k,i));
        }
        //subrange(ShapeFuncDerStack, 0,dim,i*dim,i*dim+1) = subrange(ShapeFuncDer,0,dim,i,i+1);
	}
    for (int k =dim; k<2*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
    for (int k =2*dim; k<3*dim; ++k){
        for (int m =1; m<dim*n; ++m){
            gsl_matrix_set(ShapeFuncDerStack,k,m,gsl_matrix_get(ShapeFuncDerStack,k-dim,m-1));
        }
    }
}

void Prism::calculateElementShapeFunctionDerivatives(){
    /** This function calculates the Shape Function Derivatives for all Gauss points,
     * and then calculates the derivative of reference coordiantes with respect to barycentric coordinates
     * \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$. The determnant of htis derivative is recorded for next steps of calculation.\n
     * First the matrices for the reference shape and the \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$
     * for one gauss point are allocated.
     */
    gsl_matrix* CurrRelaxedShape = gsl_matrix_calloc(nNodes, nDim);
    gsl_matrix* dXde  = gsl_matrix_calloc(nDim, nDim);
    /** The reference shape cordinates are obtained through function ShapeBase#getCurrRelaxedShape.
    */
    getCurrRelaxedShape(CurrRelaxedShape);
    for (size_t iter =0; iter<numberOfGaussPoints;++iter){
	/**
	 * Looping over all the gauss points as listed in constructor Prism#Prism,
	 * the shape function derivatives are obtained in stack form via ShapeBase#setShapeFunctionDerivatives and
	 * ShapeBase#setShapeFunctionDerivativeStack functions. Then \f$ \frac{\partial \boldsymbol{X}}{\partial \boldsymbol{\zeta}} \f$,
	 *its determinant recorded in the array  ShapeBase#detdXdes.
	 */
        double eta  = gaussPoints[iter][0];
        double nu   = gaussPoints[iter][1];
        double zeta = gaussPoints[iter][2];
        //calculating shape function derivatives:
        setShapeFunctionDerivatives(ShapeFuncDerivatives[iter],eta,zeta,nu);
        setShapeFunctionDerivativeStack(ShapeFuncDerivatives[iter],ShapeFuncDerStacks[iter]);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDerivatives[iter], CurrRelaxedShape, 0.0, dXde);
        gsl_matrix_transpose(dXde);
        detdXdes[iter] = determinant3by3Matrix(dXde);
        bool inverted = InvertMatrix(dXde, InvdXdes[iter]);
        if (!inverted){
            cerr<<"dXde not inverted at point "<<iter<<"!!"<<endl;
        }
    }
    gsl_matrix_free(CurrRelaxedShape);
    gsl_matrix_free(dXde);
}

void Prism::calculateCurrTriPointFForRotation(gsl_matrix *currF,int pointNo){
    /**
     * This function calcultes the defrmation gradient \f$ \boldsymbol{F} \f$ for the input
     * gauss point number pointNo. The calculated deformation gradient is recorded in the in input
     * gsl_matrix pointer currF, the matrix must be allocated prior to calling the function. \n
     * The deformation gradient is the derivative of current nodal positions of an element ( \f$ \boldsymbol{x} \f$ )
     * with respect to the reference shape positions ( \f$ \boldsymbol{X} \f$ ).
     *
     * \f{eqnarray*}{
        \frac{\partial \boldsymbol{x}} {\partial \boldsymbol{X}}
        \f}
     *
     */
     	const int n = nNodes;
	const int dim = nDim;
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    getPos(CurrShape);
    gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
	gsl_matrix* InvdXde = InvdXdes[pointNo];
	gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
	gsl_matrix_transpose(Jacobian);
	//calculating F:
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);
    gsl_matrix_free(CurrShape);
    gsl_matrix_free(Jacobian);
}




void Prism::calculateCurrNodalForces(gsl_matrix *currge, gsl_matrix *currgv, gsl_matrix *currF, gsl_matrix* displacementPerDt, int pointNo){
    const int n = nNodes;
    const int dim = nDim;
    gsl_matrix* currFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    //Getting the current shape positions matrix:
    getPos(CurrShape);

    //calculating dx/de (Jacobian) and reading out dX/de, shape function derivaties:
    gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
    gsl_matrix* ShapeFuncDerStack = ShapeFuncDerStacks[pointNo];
    gsl_matrix* InvdXde = InvdXdes[pointNo];
    gsl_matrix* B = Bmatrices[pointNo]; //I will update and use this value, it is not a constant.
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* invJShFuncDerSWithFe =invJShapeFuncDerStackwithFe[pointNo];
    gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
    gsl_matrix_transpose(Jacobian);
    //calculating F:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);

    //calculating Fe:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currF, InvFg, 0.0, currFe);	///< Removing growth

    createMatrixCopy(FeMatrices[pointNo], currFe); // storing Fe for use in implicit elastic K calculation.

    //setting material type:
    bool KirshoffMaterial = false;
    bool neoHookeanMaterial = !KirshoffMaterial;

    //calculating Cauchy-Green Deformation Tensor
    gsl_matrix* E;
    gsl_matrix* S;
    gsl_matrix* C = calculateCauchyGreenDeformationTensor(currFe);

    //calculating S:
    double detFe = determinant3by3Matrix(currFe);

    double lnJ = log(detFe);
    if(std::isnan(lnJ)){
    	cout<<"element: "<<Id<<" lnJ is nan, detFe: "<<detFe<<endl;
    	cout<<" Element positions: "<<endl;
    	displayPositions();
    }
    if (KirshoffMaterial){
    	//calculating E (E = 1/2 *(Fe^T*Fe-I):
    	E = calculateEForNodalForcesKirshoff(C);
    	//calculating S: (S = D:E)
    	S = calculateSForNodalForcesKirshoff(E);
    }else if (neoHookeanMaterial){
    	gsl_matrix* tmpCforInversion =  gsl_matrix_calloc(nDim,nDim);
		gsl_matrix* InvC = gsl_matrix_calloc(nDim,nDim);
		createMatrixCopy(tmpCforInversion,C);
		bool inverted = InvertMatrix(tmpCforInversion, InvC);
		if (!inverted){
			cerr<<"C not inverted!!"<<endl;
		}
    	S = calculateSForNodalForcesNeoHookean(InvC,lnJ);
    	updateLagrangianElasticityTensorNeoHookean(InvC,lnJ,pointNo);

    	//I would like to keep a record of strains, therefore I am repeating this calculation here,
    	//it does not contribute to force calculation
    	E = calculateEForNodalForcesKirshoff(C);
    	gsl_matrix_set_zero(Strain);

		gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
		gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
		gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
		gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
		gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
		gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
		gsl_matrix_free(tmpCforInversion);
		gsl_matrix_free(InvC);
    }
    //calculating elastic stress (elastic stress = detFe^-1 Fe S Fe^T):
    gsl_matrix_set_zero(elasticStress[pointNo]);
    gsl_matrix* compactStress  = calculateCompactStressForNodalForces(detFe, currFe,S, elasticStress[pointNo]);

    //Now from elastic stress, I will calculate nodal forces via B.
    //Calculating the inverse Jacobian stack matrix:
    gsl_matrix* InvJacobianStack = calculateInverseJacobianStackForNodalForces(Jacobian);

    //Calculating currB^T:
    detFs[pointNo] = determinant3by3Matrix(currF);
    gsl_matrix* currBT = calculateBTforNodalForces(InvJacobianStack,ShapeFuncDerStack, B, invJShFuncDerS);


    //Calculate invJShapeFuncDerStackwithFe for K calculation (using F in inverse jacobian calculation rather than Fe):
    calculateInvJShFuncDerSWithFe(currFe, InvdXde, ShapeFuncDerStack, invJShFuncDerSWithFe);

    //calculating nodal elastic forces as B^T compactStress detF
    gsl_matrix_scale(currBT,detFs[pointNo]);
    gsl_matrix_scale(currBT,detdXdes[pointNo]);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currBT, compactStress,0.0, currge);


    //calculating viscous stress:
    gsl_matrix_set_zero(viscousStress[pointNo]);
    if (internalViscosity != 0){
    	gsl_matrix* l = calculateVelocityGradientTensor(B, displacementPerDt);
    	gsl_matrix* d = calculateRateOfDeformationTensor(l);
    	calculateViscousStress(d,viscousStress[pointNo]);
    	//The Bt I am giving into this function has already been scaled to include volume integration.
    	calculateViscousForces(currgv, currBT,viscousStress[pointNo]);
        gsl_matrix_free(l);
        gsl_matrix_free(d);
    }
    //freeing the matrices allocated in this function
    gsl_matrix_free(C);
    gsl_matrix_free(E);
    gsl_matrix_free(S);
    gsl_matrix_free(compactStress);
    gsl_matrix_free(InvJacobianStack);
    gsl_matrix_free(currBT);
    gsl_matrix_free(currFe);
    gsl_matrix_free(CurrShape);
    gsl_matrix_free(Jacobian);

    //gsl_matrix_free(viscousStress);
}



void Prism::calculateReferenceVolume(){
	int nTriangularFaces = 8;
	int** triangularFaces = new int*[8];
	for (int i =0;i<nTriangularFaces; ++i){
		triangularFaces[i] = new int[3];
	}
	//two bases:
	triangularFaces[0][0] = 0; triangularFaces[0][1] = 1; triangularFaces[0][2] = 2;
	triangularFaces[1][0] = 3; triangularFaces[1][1] = 4; triangularFaces[1][2] = 5;
	//divide the rectangular sides into trianges:
	triangularFaces[2][0] = 0; triangularFaces[2][1] = 1; triangularFaces[2][2] = 3;
	triangularFaces[3][0] = 1; triangularFaces[3][1] = 3; triangularFaces[3][2] = 4;
	triangularFaces[4][0] = 1; triangularFaces[4][1] = 2; triangularFaces[4][2] = 4;
	triangularFaces[5][0] = 2; triangularFaces[5][1] = 4; triangularFaces[5][2] = 5;
	triangularFaces[6][0] = 0; triangularFaces[6][1] = 2; triangularFaces[6][2] = 3;
	triangularFaces[7][0] = 2; triangularFaces[7][1] = 3; triangularFaces[7][2] = 5;
	//calculating the centre point of the reference shape so that I will have the
	//mid point
	double* centre = new double[3];
	centre[0] = 0.0;
	centre[1] = 0.0;
	centre[2] = 0.0;
	for (int i = 0; i < nNodes; ++i){
		for (int j = 0 ;j < nDim; ++j){
			centre[j] += ReferenceShape->Positions[i][j];
		}
	}
	for (int j = 0 ;j < nDim; ++j){
		centre[j] /= nNodes;
	}
	ReferenceShape->Volume = calculateVolumeForInputShapeStructure(ReferenceShape->Positions, nTriangularFaces, triangularFaces,centre);
	GrownVolume = ReferenceShape->Volume;
	VolumePerNode = GrownVolume/nNodes;

	//calculating the basal area:
	double* vec1 = new double [(const int) nDim];
	double* vec2 = new double [(const int) nDim];
	double* baseVec = new double [(const int) nDim];
	for (int j=0 ;j < nDim; ++j){
			vec1 [j] = ReferenceShape->Positions[1][j] - ReferenceShape->Positions[0][j];
			vec2 [j] = ReferenceShape->Positions[2][j] - ReferenceShape->Positions[0][j];
	}
	crossProduct3D(vec1,vec2,baseVec);
	double normBaseVec= calculateMagnitudeVector3D (baseVec);
	double baseArea= normBaseVec/2;
	ReferenceShape->BasalArea = baseArea;
	delete [] vec1;
	delete [] vec2;
	delete [] baseVec;
	delete [] centre;
	for (int i =0;i<nTriangularFaces; ++i){
		delete [] triangularFaces[i];
	}
	delete [] triangularFaces;
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
        vec[i] = Positions[id1][i] - Positions[id0][i];
	}
}

void 	Prism::setInitialEdgeLenghts(){
	int node0 = 0, node1 = 1;
	for (int i=0;i<3; ++i){
		if (i ==1 ){
			node1 = 2;
		}
		if (i ==2 ){
			node0 = 1;
		}
		double dx = Positions[node0][0]-Positions[node1][0];
		double dy = Positions[node0][1]-Positions[node1][1];
		double dz = Positions[node0][2]-Positions[node1][2];
		initialBasalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
	}
	node0 = 3;
	node1 = 4;
	for (int i=0;i<3; ++i){
		if (i ==1 ){
			node1 = 5;
		}
		if (i ==2 ){
			node0 = 4;
		}
		double dx = Positions[node0][0]-Positions[node1][0];
		double dy = Positions[node0][1]-Positions[node1][1];
		double dz = Positions[node0][2]-Positions[node1][2];
		initialApilcalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
	}
}

bool Prism::areNodesDirectlyConnected(int node0, int node1){
	int * node0Placement;
	node0Placement = std::find (NodeIds, NodeIds+6, node0);
	int node1Placement=0;
	if (node0Placement != NodeIds+6){
		//node is on id list:
		ptrdiff_t pos = distance(NodeIds, node0Placement);
		if(pos >=3){
			node1Placement = pos-3;
		}
		else{
			node1Placement = pos+3;
		}
		return NodeIds[node1Placement] == node1;
	}
	return false;
}

void Prism::checkEdgeLenghtsForBinding(vector <int>&masterIds, vector <int>&slaveIds){
	/**
	* This function returns the nodes to collapse together if any of the element surfaces
	* are shrinking below the set threshold of 10% of the original side length (set in thresholdFraction variable).\n
	*/
	double thresholdFraction = 0.1; //fraction of original length the edge shrunk to
	double thresholdAngle = M_PI*8.0/9.0;// M_PI*8.0/9.0; //160 degrees, this was 100000 for most of simulations, I did not collapse based on angle
	double currentBasalEdgeLengthsSq[3] = {0,0,0};
	double currentApicalEdgeLengthsSq[3] = {0,0,0};
	int node0 = 0, node1 = 1;
	/**
	* The top side of the element will only be checked if the element top surface is at the exposed top surface of the tissue (apical
	* surface, ShapeBase#tissuePlacement = 0 or the element spans the whole tissue, ShapeBase#spansWholeTissue = true).
	*/
	if (tissuePlacement == 0 || spansWholeTissue == true){
		//check angles first, and fix if necessary:
		double* vec01 = new double[3];
		double* vec10 = new double[3];
		double* vec02 = new double[3];
		double* vec12 = new double[3];
		vec01[0] = Positions[1][0]-Positions[0][0];
		vec01[1] = Positions[1][1]-Positions[0][1];
		vec01[2] = Positions[1][2]-Positions[0][2];
		vec02[0] = Positions[2][0]-Positions[0][0];
		vec02[1] = Positions[2][1]-Positions[0][1];
		vec02[2] = Positions[2][2]-Positions[0][2];
		vec12[0] = Positions[2][0]-Positions[1][0];
		vec12[1] = Positions[2][1]-Positions[1][1];
		vec12[2] = Positions[2][2]-Positions[1][2];
		double L01 = normaliseVector3D(vec01);
		double L02 = normaliseVector3D(vec02);
		double L12 = normaliseVector3D(vec12);
		vec10[0] = vec01[0]*-1;
		vec10[1] = vec01[1]*-1;
		vec10[2] = vec01[2]*-1;
		double dot0 =dotProduct3D(vec01,vec02);
		double dot1 =dotProduct3D(vec12,vec10);
		double tet0 = acos(dot0);
		double tet1 = acos(dot1);
		if (tet0 > thresholdAngle){
			//snap to closest
			if (L01 < L02){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[1]);
			}
			else{
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[2]);
			}
		}
		if (tet1 > thresholdAngle){
			if (L01 < L12){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[1]);
			}
			else{
				masterIds.push_back(NodeIds[1]);
				slaveIds.push_back(NodeIds[2]);
			}
			//checkBasalLengths = false;
		}
		if ( M_PI - tet0 - tet1 > thresholdAngle){
			if (L02 < L12){
				masterIds.push_back(NodeIds[0]);
				slaveIds.push_back(NodeIds[2]);
			}
			else{
				masterIds.push_back(NodeIds[1]);
				slaveIds.push_back(NodeIds[2]);
			}
			//checkBasalLengths = false;
		}
		delete[] vec01;
		delete[] vec10;
		delete[] vec02;
		delete[] vec12;

		//check basal side only for the most basal section, which will belong to a basal element or an element that spans the whole tissue.
		for (int i=0;i<3; ++i){
			if (i ==1 ){
				node1 = 2;
			}
			if (i ==2 ){
				node0 = 1;
			}
			double dx = Positions[node0][0]-Positions[node1][0];
			double dy = Positions[node0][1]-Positions[node1][1];
			double dz = Positions[node0][2]-Positions[node1][2];
			currentBasalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
		}
	}
	/**
	* The bottom side of the element will only be checked if the element bottom surface is at the exposed bottom surface of the tissue (basal
	* surface, ShapeBase#tissuePlacement = 1 or the element spans the whole tissue, ShapeBase#spansWholeTissue = true).
	*/
	if (tissuePlacement == 1 || spansWholeTissue == true){
		//check angles first, and fix if necessary:
		double* vec34 = new double[3];
		double* vec43 = new double[3];
		double* vec35 = new double[3];
		double* vec45 = new double[3];
		vec34[0] = Positions[4][0]-Positions[3][0];
		vec34[1] = Positions[4][1]-Positions[3][1];
		vec34[2] = Positions[4][2]-Positions[3][2];
		vec35[0] = Positions[5][0]-Positions[3][0];
		vec35[1] = Positions[5][1]-Positions[3][1];
		vec35[2] = Positions[5][2]-Positions[3][2];
		vec45[0] = Positions[5][0]-Positions[4][0];
		vec45[1] = Positions[5][1]-Positions[4][1];
		vec45[2] = Positions[5][2]-Positions[4][2];
		double L34 = normaliseVector3D(vec34);
		double L35 = normaliseVector3D(vec35);
		double L45 = normaliseVector3D(vec45);
		vec43[0] = vec34[0]*-1;
		vec43[1] = vec34[1]*-1;
		vec43[2] = vec34[2]*-1;
		double dot3 =dotProduct3D(vec34,vec35);
		double dot4 =dotProduct3D(vec45,vec43);
		double tet3 = acos(dot3);
		double tet4 = acos(dot4);

		if (tet3 > thresholdAngle){
			//snap to closest
			if (L34 < L35){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[4]);
			}
			else{
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[5]);
			}
		}
		if (tet4 > thresholdAngle){
			if (L34 < L45){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[4]);
			}
			else{
				masterIds.push_back(NodeIds[4]);
				slaveIds.push_back(NodeIds[5]);
			}
			//checkBasalLengths = false;
		}
		if ( M_PI - tet3 - tet4 > thresholdAngle){
			if (L35 < L45){
				masterIds.push_back(NodeIds[3]);
				slaveIds.push_back(NodeIds[5]);
			}
			else{
				masterIds.push_back(NodeIds[4]);
				slaveIds.push_back(NodeIds[5]);
			}
			//checkBasalLengths = false;
		}
		delete[] vec34;
		delete[] vec43;
		delete[] vec35;
		delete[] vec45;

		node0 = 3;
		node1 = 4;
		for (int i=0;i<3; ++i){
			if (i ==1 ){
				node1 = 5;
			}
			if (i ==2 ){
				node0 = 4;
			}
			double dx = Positions[node0][0]-Positions[node1][0];
			double dy = Positions[node0][1]-Positions[node1][1];
			double dz = Positions[node0][2]-Positions[node1][2];
			currentApicalEdgeLengthsSq[i] = dx*dx+dy*dy+dz*dz;
		}
	}

	bool node1isSlave = false;
	bool node2isSlave = false;
	bool node4isSlave = false;
	bool node5isSlave = false;
	if (tissuePlacement == 0 || spansWholeTissue == true){
		if (currentBasalEdgeLengthsSq[0]<thresholdFraction*initialBasalEdgeLengthsSq[0]){
			//node 0 is too close to node 1, bind
			int master = NodeIds[0];
			int slave  = NodeIds[1];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node1isSlave = true;
		}
		if (currentBasalEdgeLengthsSq[1]<thresholdFraction*initialBasalEdgeLengthsSq[1]){
			//node 0 is too close to node 2, bind
			int master = NodeIds[0];
			int slave  = NodeIds[2];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node2isSlave = true;
		}
		if (currentBasalEdgeLengthsSq[2]<thresholdFraction*initialBasalEdgeLengthsSq[2]){
			//node 1 is too close to node 2, bind
			if (!node1isSlave || !node2isSlave){
				int master = NodeIds[1];
				if (node1isSlave){
					master = NodeIds[0];
				}
				int slave  = NodeIds[2];
				masterIds.push_back(master);
				slaveIds.push_back(slave);
			}
		}
	}
	if (tissuePlacement == 1 || spansWholeTissue == true){
		if (currentApicalEdgeLengthsSq[0]<thresholdFraction*initialApilcalEdgeLengthsSq[0]){
			//node 3 is too close to node 4, bind
			int master = NodeIds[3];
			int slave  = NodeIds[4];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node4isSlave = true;
		}
		if (currentApicalEdgeLengthsSq[1]<thresholdFraction*initialApilcalEdgeLengthsSq[1]){
			//node 3 is too close to node 5, bind
			int master = NodeIds[3];
			int slave  = NodeIds[5];
			masterIds.push_back(master);
			slaveIds.push_back(slave);
			node5isSlave = true;
		}
		if (currentApicalEdgeLengthsSq[2]<thresholdFraction*initialApilcalEdgeLengthsSq[2]){
			if (!node4isSlave || !node5isSlave){
				int master = NodeIds[4];
				if (node1isSlave){
					master = NodeIds[3];
				}
				int slave  = NodeIds[5];
				masterIds.push_back(master);
				slaveIds.push_back(slave);
			}
		}
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
	int counter = 0;
	int pairs[3][2] = {{3,4},{3,5},{4,5}};
	if (tissueType == 1){
		//peripodial
		pairs[0][0] = 0;
		pairs[0][1] = 1;
		pairs[1][0] = 0;
		pairs[1][1] = 2;
		pairs[2][0] = 1;
		pairs[2][1] = 2;
	}
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		double dsum2 = (dx*dx + dy*dy+dz*dz);
		if(dsum2 > 0){
			//the average side length can be zero for elements with collapsed nodes.
			//I do not want to add the zero length in calculation of average length.
			//I will handle the case where this function returns zero outside, in averaging all elements.
			dsum += pow((dx*dx + dy*dy+dz*dz),0.5);
			counter++;
		}
	}
	if (counter>0){
		dsum /= counter;
	}
	return dsum;
}

double Prism::getBasalSideLengthAverage(){
	double dx,dy,dz;
	double dsum =0.0;
	int counter = 0;
	int pairs[3][2] = {{0,1},{0,2},{1,2}};
	if (tissueType == 1){
		//peripodial
		pairs[0][0] = 3;
		pairs[0][1] = 4;
		pairs[1][0] = 3;
		pairs[1][1] = 5;
		pairs[2][0] = 4;
		pairs[2][1] = 5;
	}
	for (int i=0; i<3; ++i){
		dx = Positions[pairs[i][0]][0] - Positions[pairs[i][1]][0];
		dy = Positions[pairs[i][0]][1] - Positions[pairs[i][1]][1];
		dz = Positions[pairs[i][0]][2] - Positions[pairs[i][1]][2];
		double dsum2 = (dx*dx + dy*dy+dz*dz);
		if(dsum2 > 0){
			//the average side length can be zero for elements with collapsed nodes.
			//I do not want to add the zero length in calculation of average length.
			//I will handle the case where this function returns zero outside, in averaging all elements.
			dsum += pow((dx*dx + dy*dy+dz*dz),0.5);
			counter++;
		}
	}
	if (counter>0){
		dsum /= counter;
	}
	return dsum;
}

void Prism::getApicalTriangles(vector <int> &ApicalTriangles){
	if (tissueType == 0){
		ApicalTriangles.push_back(NodeIds[3]);
		ApicalTriangles.push_back(NodeIds[4]);
		ApicalTriangles.push_back(NodeIds[5]);
	}
	else{
		ApicalTriangles.push_back(NodeIds[0]);
		ApicalTriangles.push_back(NodeIds[1]);
		ApicalTriangles.push_back(NodeIds[2]);
	}
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


bool Prism::IsThisNodeMyApical(int currNodeId){
	if (NodeIds[3] == currNodeId || NodeIds[4] == currNodeId || NodeIds[5] == currNodeId ){
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


void Prism::AddPackingToSurface(int tissueplacementOfPackingNode, double Fx, double Fy,double Fz, double **PackingForces, vector<Node*> &Nodes, bool& allCornersFixedX, bool& allCornersFixedY, bool& allCornersFixedZ){
	int Id0 = -1, Id1= -1, Id2 = -1;
	if (tissueplacementOfPackingNode == 0){//basal node, packing to basal surface:
		if (tissueType == 0){ //element is columnar, basal surface is bottom
			Id0 = 0;
			Id1 = 1;
			Id2 = 2;
		}
		else {//element is peripodial, basal surface is top
			Id0 = 3;
			Id1 = 4;
			Id2 = 5;
		}
	}
	if (tissueplacementOfPackingNode == 1){//apical node, packing to apical surface:
		if (tissueType == 0){ //element is columnar, apical surface is top
			Id0 = 3;
			Id1 = 4;
			Id2 = 5;
		}
		else {//element is peripodial, apical surface is bottom
			Id0 = 0;
			Id1 = 1;
			Id2 = 2;
		}
	}
	double F[3];
	//F[0] = Fx / 3.0;
	//F[1] = Fy / 3.0;
	//F[2] = Fz / 3.0;
	//checking how many corners are fixed
	int numberOfNonFixedCorners[3] = {0,0,0};
	for(int j=0; j<nDim; ++j){
		if (!Nodes[NodeIds[Id0]]->FixedPos[j]){
			numberOfNonFixedCorners[j]++;
		}
		if (!Nodes[NodeIds[Id1]]->FixedPos[j]){
			numberOfNonFixedCorners[j]++;
		}
		if (!Nodes[NodeIds[Id2]]->FixedPos[j]){
			numberOfNonFixedCorners[j]++;
		}
	}
	if (numberOfNonFixedCorners[0] == 0){
		allCornersFixedX = true;
		F[0] = 0.0;
	}
	else{
		F[0] = Fx / numberOfNonFixedCorners[0];
	}
	if (numberOfNonFixedCorners[1] == 0){
		allCornersFixedY = true;
		F[1] = 0.0;
	}
	else{
		F[1] = Fy / numberOfNonFixedCorners[1];
	}
	if (numberOfNonFixedCorners[2] == 0){
		allCornersFixedZ = true;
		F[2] = 0.0;
	}
	else{
		F[2] = Fz / numberOfNonFixedCorners[2];
	}

	for(int j=0; j<nDim; ++j){
		if (!Nodes[NodeIds[Id0]]->FixedPos[j]){
            PackingForces[NodeIds[Id0]][j] -= F[j];
		}
		if (!Nodes[NodeIds[Id1]]->FixedPos[j]){
			PackingForces[NodeIds[Id1]][j] -= F[j];
		}
		if (!Nodes[NodeIds[Id2]]->FixedPos[j]){
			PackingForces[NodeIds[Id2]][j] -= F[j];
		}
	}
}

bool Prism::IsPointCloseEnoughForPacking(double* Pos,  float threshold, int TissuePlacementOfPackingNode){
	int initial =0, final = 6;
	if (TissuePlacementOfPackingNode == 0){
		//tissue placement of the node is basal, should check it against basal surface nodes only
		if(tissueType == 0){ //element is columnar, the basal nodes are nodes 0-2:
			final = 3;
		}
		else{
			//element is peripodial, the basal nodes are 3-5
			initial = 3;
		}
	}
	else if (TissuePlacementOfPackingNode == 1){
		//tissue placement of the node is apical, should check it against apical surface nodes  only
		if(tissueType == 0){ //element is columnar, the apical nodes are nodes 3-5:
			initial = 3;
		}
		else{
			//element is peripodial, the  the apical nodes are nodes 0-2:
			final = 3;
		}
	}
	//cout<<" initial : "<<initial<<" final: "<<final<<endl;
	float dmin = threshold;
	float dminNeg = (-1.0)*threshold;
	for (int i=initial; i<final; ++i){
		//cout<<"	checking against node: "<<NodeIds[i]<<endl;
		float dx =100.0, dy = 100.0, dz = 100.0;
		dx = Pos[0]-Positions[i][0];
		//cout<<" dx: "<<dx<<endl;
		if ((dx >=0 && dx < dmin) || (dx <=0 && dx >dminNeg)){
			dy = Pos[1]-Positions[i][1];
			//cout<<" dy: "<<dy<<endl;
			if ((dy >=0 && dy < dmin) || (dy <=0 && dy >dminNeg)){
				dz = Pos[2]-Positions[i][2];
				//cout<<" dz: "<<dz<<endl;
				if ((dz >=0 && dz < dmin) || (dz <=0 && dz >dminNeg)){
					return true;
				}
			}
		}
	}
	return false;
}

void  Prism::getApicalNodePos(double* posCorner){
	if (tissueType == 0){//columnar element, apical surface is the top nodes
		posCorner[0] = Positions[3][0];
		posCorner[1] = Positions[3][1];
		posCorner[2] = Positions[3][2];
	}
	else{
		//peripodial element, apical surface is the bottom nodes
		//the linker elements never reach this question
		posCorner[0] = Positions[0][0];
		posCorner[1] = Positions[0][1];
		posCorner[2] = Positions[0][2];

	}
}

void  Prism::getBasalNodePos(double* posCorner){
	if (tissueType == 0){//columnar element, basal surface is the bottom nodes
		posCorner[0] = Positions[0][0];
		posCorner[1] = Positions[0][1];
		posCorner[2] = Positions[0][2];
	}
	else{
		//peripodial element, basal surface is the top nodes
		//the linker elements never reach this question
		posCorner[0] = Positions[3][0];
		posCorner[1] = Positions[3][1];
		posCorner[2] = Positions[3][2];

	}
}

void Prism::getApicalNodeIds(vector <int> &nodeIds){
	if (tissueType == 0){//columnar element, apical surface is the top nodes
		//I iwll rotate the nodes as I want the normal to point into the lumen:
		nodeIds.push_back(NodeIds[3]);
		nodeIds.push_back(NodeIds[5]);
		nodeIds.push_back(NodeIds[4]);
	}
	else{
		//peripodial element, apical surface is the bottom nodes
		//the linker elements never reach this question

		nodeIds.push_back(NodeIds[0]);
		nodeIds.push_back(NodeIds[1]);
		nodeIds.push_back(NodeIds[2]);
	}
}

void Prism::getApicalNodeIndicesOnElement(vector <int> &apicalNodeIndices){
	if (tissueType == 0){//columnar element, apical surface is the top nodes
		//I iwll rotate the nodes as I want the normal to point into the lumen:
		apicalNodeIndices.push_back(3);
		apicalNodeIndices.push_back(5);
		apicalNodeIndices.push_back(4);
	}
	else{
		//peripodial element, apical surface is the bottom nodes
		//the linker elements never reach this question
		apicalNodeIndices.push_back(0);
		apicalNodeIndices.push_back(1);
		apicalNodeIndices.push_back(2);
	}
}

void Prism::getBasalNodeIds(vector <int> &nodeIds){
	if (tissueType == 0){//columnar element, apical surface is the top nodes
		//I iwll rotate the nodes as I want the normal to point into the lumen on the apical side
		//and I want these to be consistent:
		nodeIds.push_back(NodeIds[0]);
		nodeIds.push_back(NodeIds[2]);
		nodeIds.push_back(NodeIds[1]);
	}
	else{
		//peripodial element, apical surface is the bottom nodes
		//the linker elements never reach this question
		nodeIds.push_back(NodeIds[3]);
		nodeIds.push_back(NodeIds[4]);
		nodeIds.push_back(NodeIds[5]);
	}
}

void Prism::getBasalCentre(double* centre){
	centre[0] = 0;
	centre[1] = 0;
	centre[2] = 0;
	if (tissueType == 0){
		for (int j=0; j<3; ++j){
			for (int i=0; i<3; ++i){ //columnar element, basal surface is the bottom nodes, (first 3 nodes)
				centre[j] += Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
	else{
		//
		//the linker elements never reach this question
		for (int j=0; j<3; ++j){
			for (int i=3; i<6; ++i){//peripodial element, basal surface is the top nodes, (last 3 nodes)
				centre[j] += Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
}

void Prism::getApicalCentre(double* centre){
	centre[0] = 0;
	centre[1] = 0;
	centre[2] = 0;
	if (tissueType == 0){
		for (int j=0; j<3; ++j){
			for (int i=3; i<6; ++i){ //columnar element, apical surface is the top nodes, (last 3 nodes)
				centre[j] += Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
	else{
		//
		//the linker elements never reach this question
		for (int j=0; j<3; ++j){
			for (int i=0; i<3; ++i){//peripodial element, apical surface is the bottom nodes, (first 3 nodes)
				centre[j] += Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
}

void Prism::getReferenceBasalCentre(double* centre){
	centre[0] = 0;
	centre[1] = 0;
	centre[2] = 0;
	if (tissueType == 0){
		for (int j=0; j<3; ++j){
			for (int i=0; i<3; ++i){ //columnar element, basal surface is the bottom nodes, (first 3 nodes)
				centre[j] += ReferenceShape->Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
	else{
		//
		//the linker elements never reach this question
		for (int j=0; j<3; ++j){
			for (int i=3; i<6; ++i){//peripodial element, basal surface is the top nodes, (last 3 nodes)
				centre[j] += ReferenceShape->Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
}

void Prism::getReferenceApicalCentre(double* centre){
	centre[0] = 0;
	centre[1] = 0;
	centre[2] = 0;
	if (tissueType == 0){
		for (int j=0; j<3; ++j){
			for (int i=3; i<6; ++i){ //columnar element, apical surface is the top nodes, (last 3 nodes)
				centre[j] += ReferenceShape->Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
	else{
		//
		//the linker elements never reach this question
		for (int j=0; j<3; ++j){
			for (int i=0; i<3; ++i){//peripodial element, apical surface is the bottom nodes, (first 3 nodes)
				centre[j] += ReferenceShape->Positions[i][j];
			}
			centre[j] /= 3.0;
		}
	}
}

double* Prism::getBasalMinViscosity(vector<Node*> Nodes){
	double* ExtVisc = new double[3];
	ExtVisc[0] = 10000000.0;
	ExtVisc[1] = 10000000.0;
	ExtVisc[2] = 10000000.0;

	if (tissueType == 0){
		for (int i=0; i<3; ++i){
			for (int j= 0; j<nDim; ++j){
				if (ExtVisc[j] > Nodes[NodeIds[i]]->externalViscosity[j]){
					ExtVisc[j] = Nodes[NodeIds[i]]->externalViscosity[j];
				}
			}
		}
	}
	else{
		for (int i=3; i<6; ++i){//peripodial element, apical surface is the bottom nodes, (first 3 nodes)
			for (int j= 0; j<nDim; ++j){
				if (ExtVisc[j] > Nodes[NodeIds[i]]->externalViscosity[j]){
					ExtVisc[j] = Nodes[NodeIds[i]]->externalViscosity[j];
				}
			}
		}
	}
	return ExtVisc;
}

double* Prism::getApicalMinViscosity(vector<Node*> Nodes){
	double* ExtVisc = new double[3];
	ExtVisc[0] = 10000000.0;
	ExtVisc[1] = 10000000.0;
	ExtVisc[2] = 10000000.0;
	if (tissueType == 0){
		for (int i=3; i<6; ++i){
			for (int j= 0; j<nDim; ++j){
				if (ExtVisc[j] > Nodes[NodeIds[i]]->externalViscosity[j]){
					ExtVisc[j] = Nodes[NodeIds[i]]->externalViscosity[j];
				}
			}
		}
	}
	else{
		for (int i=0; i<3; ++i){//peripodial element, apical surface is the bottom nodes, (first 3 nodes)
			for (int j= 0; j<nDim; ++j){
				if (ExtVisc[j] > Nodes[NodeIds[i]]->externalViscosity[j]){
					ExtVisc[j] = Nodes[NodeIds[i]]->externalViscosity[j];
				}
			}
		}
	}
	return ExtVisc;
}

void  	Prism::calculateMyosinForcesAreaBased(double forcePerMyoMolecule){
	if (cMyoUniform[0]>0){
		calculateApicalArea();
		distributeMyosinForcesAreaBased(true,true, forcePerMyoMolecule); //isIsotropic, isApical, force per myosin molecule to calculte the actual total force
		//there is unfiorm myosin activity on apical surface
	}
	if (cMyoUniform[1]>0){
		calculateBasalArea();
		distributeMyosinForcesAreaBased(true,false, forcePerMyoMolecule);
		//there is unfiorm myosin activity on basal surface
	}
	if (cMyoUnipolar[0]>0){
		calculateApicalArea();
		distributeMyosinForcesAreaBased(false,true, forcePerMyoMolecule);
		//there is unipolar myosin activity on apical surface
	}
	if (cMyoUnipolar[1]>0){
		calculateBasalArea();
		distributeMyosinForcesAreaBased(false,false, forcePerMyoMolecule);
		//there is unipolar myosin activity on basal surface
	}
}

void  	Prism::calculateMyosinForcesTotalSizeBased(double forcePerMyoMolecule){
	if (cMyoUniform[0]>0){
		distributeMyosinForcesTotalSizeBased(true,true, forcePerMyoMolecule); //isIsotropic, isApical, force per myosin molecule to calculte the actual total force
		//there is unfiorm myosin activity on apical surface
	}
	if (cMyoUniform[1]>0){
		distributeMyosinForcesTotalSizeBased(true,false, forcePerMyoMolecule);
		//there is unfiorm myosin activity on basal surface
	}
	if (cMyoUnipolar[0]>0){
		distributeMyosinForcesTotalSizeBased(false,true, forcePerMyoMolecule);
		//there is unipolar myosin activity on apical surface
	}
	if (cMyoUnipolar[1]>0){
		distributeMyosinForcesTotalSizeBased(false,false, forcePerMyoMolecule);
		//there is unipolar myosin activity on basal surface
	}
}

void 	Prism::distributeMyosinForcesAreaBased(bool isIsotropic, bool apical, double forcePerMyoMolecule){
	int id0, id1, id2;
	double forcemag;
	int currSurface = 0;
	if (apical){
		id0 = 3;
		id1 = 4;
		id2 = 5;
		currSurface = 0;
		if (isIsotropic){
			forcemag = cMyoUniform[0]*ApicalArea;
		}
		else{
			forcemag = cMyoUnipolar[0]*ApicalArea;
		}
	}
	else{
		id0 = 0;
		id1 = 1;
		id2 = 2;
		currSurface = 1;
		if (isIsotropic){
			forcemag = cMyoUniform[1]*BasalArea;
		}
		else{
			forcemag = cMyoUnipolar[1]*BasalArea;
		}
	}
	forcemag *= forcePerMyoMolecule;
	double centre[3] = {0.0,0.0,0.0};
	gsl_vector* vec0 = gsl_vector_calloc(3);
	gsl_vector* vec1 = gsl_vector_calloc(3);
	gsl_vector* vec2 = gsl_vector_calloc(3);
	for (int i=0; i<3; ++i){
		centre[i] = (Positions[id0][i] + Positions[id1][i] +  Positions[id2][i])/3.0;
		//axisSum[i] = -2.0 * centre[i];
		gsl_vector_set(vec0,i,centre[i] - Positions[id0][i]);
		gsl_vector_set(vec1,i,centre[i] - Positions[id1][i]);
		gsl_vector_set(vec2,i,centre[i] - Positions[id2][i]);
	}
	if (!isIsotropic){
		//If I am doing a non-isotropic calculation, then I will project the vectors towards the centre
		// on top of the axis of interest:
		//First I need to align the myosin polarity direction on to the surface, the forces are on the surface!:
		//get the normal:
		gsl_vector* cross = gsl_vector_calloc(3);
		crossProduct3D(vec0,vec1,cross);
		normaliseVector3D(cross);
		//project vector on normal:
		//cout<<"myo polarity dir: "<<gsl_matrix_get(myoPolarityDir,currSurface,0)<<" "<<gsl_matrix_get(myoPolarityDir,currSurface,1)<<endl;
		double dotp = 0.0;
		for (int i =0; i<3; ++i){
			dotp += gsl_vector_get(cross,i)*gsl_matrix_get(myoPolarityDir,currSurface,i);
		}
		gsl_vector* alignedmyoPolarityDir = gsl_vector_calloc(3);
		for (int i =0; i<3; ++i){
			gsl_vector_set(alignedmyoPolarityDir,i,gsl_matrix_get(myoPolarityDir,currSurface,i) - dotp*gsl_vector_get(cross,i));
		}
		double proj[3] = {0,0,0};
		for (int i =0; i<3; ++i){
			proj[0] += gsl_vector_get(vec0,i)*gsl_vector_get(alignedmyoPolarityDir,i);
			proj[1] += gsl_vector_get(vec1,i)*gsl_vector_get(alignedmyoPolarityDir,i);
			proj[2] += gsl_vector_get(vec2,i)*gsl_vector_get(alignedmyoPolarityDir,i);
		}
		for (int i =0; i<3; ++i){
			gsl_vector_set(vec0,i,proj[0]*gsl_vector_get(alignedmyoPolarityDir,i));
			gsl_vector_set(vec1,i,proj[1]*gsl_vector_get(alignedmyoPolarityDir,i));
			gsl_vector_set(vec2,i,proj[2]*gsl_vector_get(alignedmyoPolarityDir,i));
		}
		gsl_vector_free(alignedmyoPolarityDir);
		gsl_vector_free(cross);
	}
	//now I will distribute the force proportional to the magnitude of each vector:
	double sumMag = 0.0;
	double x = gsl_vector_get(vec0,0);
	double y = gsl_vector_get(vec0,1);
	double z = gsl_vector_get(vec0,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	x = gsl_vector_get(vec1,0);
	y = gsl_vector_get(vec1,1);
	z = gsl_vector_get(vec1,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	x = gsl_vector_get(vec2,0);
	y = gsl_vector_get(vec2,1);
	z = gsl_vector_get(vec2,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	double scaleFactor = forcemag/sumMag;
	gsl_vector_scale(vec0,scaleFactor);
	gsl_vector_scale(vec1,scaleFactor);
	gsl_vector_scale(vec2,scaleFactor);

	for (int i=0; i<3; ++i){
		MyoForce[id0][i] += gsl_vector_get(vec0,i);
		MyoForce[id1][i] += gsl_vector_get(vec1,i);
		MyoForce[id2][i] += gsl_vector_get(vec2,i);
	}
	gsl_vector_free(vec0);
	gsl_vector_free(vec1);
	gsl_vector_free(vec2);
	/*cout<<" Element: "<<Id<<" forcemag: "<<forcemag<<" Myoforces: "<<endl;
	//double sum[3] = {0,0,0};
	for (int i=0; i<6; i++){
		for (int j=0; j<3; j++){
			cout<<MyoForce[i][j]<<" ";
			//sum[j] += MyoForce[i][j];
		}
		cout<<endl;
	}
	cout<<endl;
	displayMatrix(myoPolarityDir,"myoPolarityDir");
	*/
	//cout<<"sum: "<<sum[0]<<" "<<sum[1]<<" "<<sum[2]<<endl;


}

void 	Prism::distributeMyosinForcesTotalSizeBased(bool isIsotropic, bool apical, double forcePerMyoMolecule){
	int id0, id1, id2;
	double forcemag;
	int currSurface = 0;
	double AverageSurfaceArea = GrownVolume/(gsl_matrix_get(Fg,2,2)*ReferenceShape->height);
	if (apical){
		id0 = 3;
		id1 = 4;
		id2 = 5;
		currSurface = 0;
		if (isIsotropic){
			forcemag = cMyoUniform[0]*AverageSurfaceArea;
		}
		else{
			forcemag = cMyoUnipolar[0]*AverageSurfaceArea;
		}
	}
	else{
		id0 = 0;
		id1 = 1;
		id2 = 2;
		currSurface = 1;
		if (isIsotropic){
			forcemag = cMyoUniform[1]*AverageSurfaceArea;
		}
		else{
			forcemag = cMyoUnipolar[1]*AverageSurfaceArea;
		}
	}
	double centre[3] = {0.0,0.0,0.0};
	gsl_vector* vec0 = gsl_vector_calloc(3);
	gsl_vector* vec1 = gsl_vector_calloc(3);
	gsl_vector* vec2 = gsl_vector_calloc(3);
	for (int i=0; i<3; ++i){
		centre[i] = (Positions[id0][i] + Positions[id1][i] +  Positions[id2][i])/3.0;
		//axisSum[i] = -2.0 * centre[i];
		gsl_vector_set(vec0,i,centre[i] - Positions[id0][i]);
		gsl_vector_set(vec1,i,centre[i] - Positions[id1][i]);
		gsl_vector_set(vec2,i,centre[i] - Positions[id2][i]);
	}
	if (!isIsotropic){
		//If I am doing a non-isotropic calculation, then I will project the vectors towards the centre
		// on top of the axis of interest:
		//First I need to align the myosin polarity direction on to the surface, the forces are on the surface!:
		//get the normal:
		gsl_vector* cross = gsl_vector_calloc(3);
		crossProduct3D(vec0,vec1,cross);
		normaliseVector3D(cross);
		//project vector on normal:
		//cout<<"myo polarity dir: "<<gsl_matrix_get(myoPolarityDir,currSurface,0)<<" "<<gsl_matrix_get(myoPolarityDir,currSurface,1)<<endl;
		double dotp = 0.0;
		for (int i =0; i<3; ++i){
			dotp += gsl_vector_get(cross,i)*gsl_matrix_get(myoPolarityDir,currSurface,i);
		}
		gsl_vector* alignedmyoPolarityDir = gsl_vector_calloc(3);
		for (int i =0; i<3; ++i){
			gsl_vector_set(alignedmyoPolarityDir,i,gsl_matrix_get(myoPolarityDir,currSurface,i) - dotp*gsl_vector_get(cross,i));
		}
		double proj[3] = {0,0,0};
		for (int i =0; i<3; ++i){
			proj[0] += gsl_vector_get(vec0,i)*gsl_vector_get(alignedmyoPolarityDir,i);
			proj[1] += gsl_vector_get(vec1,i)*gsl_vector_get(alignedmyoPolarityDir,i);
			proj[2] += gsl_vector_get(vec2,i)*gsl_vector_get(alignedmyoPolarityDir,i);
		}
		for (int i =0; i<3; ++i){
			gsl_vector_set(vec0,i,proj[0]*gsl_vector_get(alignedmyoPolarityDir,i));
			gsl_vector_set(vec1,i,proj[1]*gsl_vector_get(alignedmyoPolarityDir,i));
			gsl_vector_set(vec2,i,proj[2]*gsl_vector_get(alignedmyoPolarityDir,i));
		}
		gsl_vector_free(alignedmyoPolarityDir);
		gsl_vector_free(cross);
	}
	//now I will distribute the force proportional to the magnitude of each vector:
	double sumMag = 0.0;
	double x = gsl_vector_get(vec0,0);
	double y = gsl_vector_get(vec0,1);
	double z = gsl_vector_get(vec0,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	x = gsl_vector_get(vec1,0);
	y = gsl_vector_get(vec1,1);
	z = gsl_vector_get(vec1,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	x = gsl_vector_get(vec2,0);
	y = gsl_vector_get(vec2,1);
	z = gsl_vector_get(vec2,2);
	sumMag += pow(x*x+y*y+z*z,0.5);
	double scaleFactor = forcemag/sumMag;
	gsl_vector_scale(vec0,scaleFactor);
	gsl_vector_scale(vec1,scaleFactor);
	gsl_vector_scale(vec2,scaleFactor);

	for (int i=0; i<3; ++i){
		MyoForce[id0][i] += gsl_vector_get(vec0,i);
		MyoForce[id1][i] += gsl_vector_get(vec1,i);
		MyoForce[id2][i] += gsl_vector_get(vec2,i);
	}
	gsl_vector_free(vec0);
	gsl_vector_free(vec1);
	gsl_vector_free(vec2);
}

void 	Prism::calculateBasalArea(){
    double Threshold = 1E-5;
	int id0 = 0, id1 = 1, id2 = 2; // this is correct for basal side, I will change it for apical calculation
	if (tissueType == 1){ //peroipodial element, the basal surface is actually on top
		id0 = 3;
		id1 = 4;
		id2 = 5;
	}
	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;

	for (int i = 0; i<3; ++i){
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
		double sinTetSq  = 1-costet*costet;
		double sintet = 0;
		if(sinTetSq>Threshold){
			sintet = pow(sinTetSq,0.5);
		}
		Area = Side1* Side2 * sintet / 2.0;
	}
	BasalArea = Area;
	//cout<<" Element "<<Id<<" basal area: "<<BasalArea<<endl;
	//ApicalArea = Area;
}

void 	Prism::calculateApicalArea(){
    double Threshold = 1E-5;
	int id0 = 3, id1 = 4, id2 = 5; // this is correct for basal side, I will change it for apical calculation
	if (tissueType == 1){ //peroipodial element, the basal surface is actually at the bottom
		id0 = 0;
		id1 = 1;
		id2 = 2;
	}
	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;
	for (int i = 0; i<3; ++i){
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
		double sinTetSq  = 1-costet*costet;
		double sintet = 0;
		if(sinTetSq>Threshold){
			sintet = pow(sinTetSq,0.5);
		}
		Area = Side1* Side2 * sintet / 2.0;
	}
	ApicalArea = Area;
}

void Prism::setBasalNeigElementId(vector<ShapeBase*>& elementsList){
    /**
     * This function sets the basal neighbours of all elemetns qualifying for a basal neighbour element.
     * The check is only carried out on the apical elements of the  coumnar tissue (ShapeBase#tissuePlacement = 0, ShapeBase#TissueType = 0).
     * The columnar elements bordering the basal extracellualar matrix (ShapeBase#atBasalBorderOfECM = true),
     * and the basal extracellular matrix itself (ShapeBase#isECMMimicing) are not
     * labelled, as these elements do not hold a basal neighbour that is a tissue piece.
     */
     	if (tissueType== 0 //columnar layer element
			&& !tissuePlacement == 0 //not checking for basal elements, to save time
			&& !atBasalBorderOfECM  //not checking for elements bordering ECM, same as basal
			&& !isECMMimicing		 //not checking for ECM mimicking elements
			){
		//I am setting this only for columnar elements.
		//I do not need to check for spanning the whole tissue, as those elements
		//will not have a basal neighbour of relevance.
		/** Basal nodes of a prism are hte  nodes recorded in ShapeBase#NodeIds array indices 0-2.
		*/
		int currBasalNodes[3] = {NodeIds[0], NodeIds[1],NodeIds[2]};
		for( vector<ShapeBase*>::iterator itElement=elementsList.begin(); itElement<elementsList.end(); ++itElement){
			bool isApicalOwner  = (*itElement)->IsThisNodeMyApical(currBasalNodes[0]);
			if (isApicalOwner){
				//node 0 is owned by this element, is node 1?
				isApicalOwner  = (*itElement)->IsThisNodeMyApical(currBasalNodes[1]);
				if (isApicalOwner){
					//node 0 & 1 are owned by this element, is node 2?
					isApicalOwner  = (*itElement)->IsThisNodeMyApical(currBasalNodes[2]);
					if (isApicalOwner){
						/** When all three nodes are owned as apical nodes by another element, this is the
						* neighbour I am looking for. If the element is NOT an explicit ECM element (ShapeBase#isECMMimicing = false),
						* then the new Id is recorded into the ShapeBase#basalNeigElementId.
						*/
						if (!(*itElement)-> isECMMimicing){
							basalNeigElementId = (*itElement)->getId();
						}
						/** Once the element is found, break the loop regardless of the ECM status. If the loop has
						*  reached an ECM element then this element does not have a basal bordering tissue. It should not
						* have this case anyway, as elemetns atBasalBorderOfECM are skipped.
						*/
						break;
					}
				}
			}
		}
	}
}

void Prism::constructElementStackList(const int discretisationLayers, vector<ShapeBase*>& elementsList){
	//This is a basal element. I will take my basal nodes, and find elements that have them as apical nodes.
	elementsIdsOnSameColumn	= new int[discretisationLayers];
	int currApicalNodes[3] = {NodeIds[3], NodeIds[4],NodeIds[5]};
	int currElementId = Id;
	//Starting the list with the basal element
	elementsIdsOnSameColumn[0] = currElementId;
	for (size_t i=1; i<  discretisationLayers; ++i ){
		//cout<<" starting to fill the list item: "<<i <<" of "<< discretisationLayers-1<<endl;
		for( vector<ShapeBase*>::iterator itElement=elementsList.begin(); itElement<elementsList.end(); ++itElement){
			int checkedElementId = (*itElement)->getId();
			if (checkedElementId != currElementId){ //the iterator is not this element
				//I will find the elemetns that have this elements apical nodes as basal. All three ndes should match;
				bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currApicalNodes[0]);
				if (IsBasalOwner){
					IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currApicalNodes[1]);
				}
				if (IsBasalOwner){
					IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currApicalNodes[2]);
				}
				//If IsBasalOwner is true at this point, then I found the matching element.
				if (IsBasalOwner){
					//I found the element:
					elementsIdsOnSameColumn[i] = checkedElementId;
					for (int j=0; j<3; ++j){
						currApicalNodes[j] = (*itElement)->getCorrecpondingApical(currApicalNodes[j]);
					}
					currElementId = checkedElementId;
					break;
				}
			}
		}
	}
	//Now I have the list for the basal element.
	//I should construct it for at midline and apical elements:
	for (size_t i=1; i<  discretisationLayers; ++i ){
		elementsList[elementsIdsOnSameColumn[i]]->elementsIdsOnSameColumn = new int[discretisationLayers];
		for (size_t j=0; j<  discretisationLayers; ++j ){
			elementsList[elementsIdsOnSameColumn[i]]->elementsIdsOnSameColumn[j] = elementsIdsOnSameColumn[j];
		}
	}
};


void Prism::assignExposedSurfaceAreaIndices(vector <Node*>& Nodes){
	if (tissueType == 0 || tissueType == 1){//columnar or peripodial Element
		//these are written for columnar tissue, should be swapped for peorpodial tissue:
		int apicalIds[3] = {3,4,5};
		int basalIds[3] = {0,1,2};
		if (tissueType == 1){//peripodial Element
			apicalIds[0]= 0;
			apicalIds[1]= 1;
			apicalIds[2]= 2;
			basalIds[0]= 3;
			basalIds[1]= 4;
			basalIds[2]= 5;
		}
		if (tissuePlacement == 1){ //apical element
			elementHasExposedApicalSurface = true;
			//If element is apical, I will calculate the apical surface
			exposedApicalSurfaceNodeIds[0] = apicalIds[0];
			exposedApicalSurfaceNodeIds[1] = apicalIds[1];
			exposedApicalSurfaceNodeIds[2] = apicalIds[2];
		}
		else if(tissuePlacement == 0){ //basal element
			elementHasExposedBasalSurface = true;
			//If element is basal, I will calculate the basal surface
			exposedBasalSurfaceNodeIds[0] = basalIds[0];
			exposedBasalSurfaceNodeIds[1] = basalIds[1];
			exposedBasalSurfaceNodeIds[2] = basalIds[2];
		}
		else if (tissuePlacement == 2 && spansWholeTissue) {//midline element, but spans whole tissue
			elementHasExposedApicalSurface = true;
			elementHasExposedBasalSurface = true;
			exposedApicalSurfaceNodeIds[0] = apicalIds[0];
			exposedApicalSurfaceNodeIds[1] = apicalIds[1];
			exposedApicalSurfaceNodeIds[2] = apicalIds[2];
			exposedBasalSurfaceNodeIds[0] = basalIds[0];
			exposedBasalSurfaceNodeIds[1] = basalIds[1];
			exposedBasalSurfaceNodeIds[2] = basalIds[2];
		}
	}
	else if (tissueType == 2){//linker element
		int numberOfApicalNodes = 0;
		int numberOfBasalNodes = 0;
		for (int i=0; i<nNodes; ++i){
			if (Nodes[NodeIds[i]]->tissuePlacement == 0){ //basal Node
				exposedLateralAreaBasalSideNodeIds[numberOfBasalNodes] = i;
				numberOfBasalNodes++;
			}
			if (Nodes[NodeIds[i]]->tissuePlacement == 1){ //apical Node
				exposedLateralAreaApicalSideNodeIds[numberOfApicalNodes] = i;
				numberOfApicalNodes++;
			}
		}
		if (numberOfApicalNodes>3){
			elementHasExposedLateralApicalSurface = true;
		}
		if (numberOfBasalNodes>3){
			elementHasExposedLateralBasalSurface = true;
		}
	}
}

