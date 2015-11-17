#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
//this was the working version in linux. It should be working with correct addition of the path to INCLUDEPATH in .pro file
//#include </usr/include/gsl/gsl_matrix.h>
//#include </usr/include/gsl/gsl_linalg.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "Node.h"
#include "ReferenceShapeBase.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"

using namespace std;

class ShapeBase{
private:
	void ParentErrorMessage(string functionName);					 	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string inputs.
	bool ParentErrorMessage(string functionName, bool returnValue); 	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and bool inputs, returning bool.
	double ParentErrorMessage(string functionName, double returnValue);	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and double inputs, returning double.
	int ParentErrorMessage(string functionName, int returnValue);		///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and int inputs, returning int.
protected:
	int 	ShapeType;				///< The integer defining the type of the shape, Prisms shape type = 1;
	int 	nNodes;					///< The number of nodes of the element, it is based on ShapeBase#ShapeType
	int 	nDim;					///< The number of dimensions for the positions of each of the nodes of the element
	int* 	IdentifierColour;		///< The unique identifier colour of the element, this is used for "picking" in the visual interface.
	double* GrowthRate;				///< Growth rate recording for display purposes only. The recorded growth rate in x, y, and z  coordinates, does not record shear deformation induced in growth. Recorded in exponential form through time step, converted to rate per hour for display within the visual interface
	gsl_matrix* growthIncrement;	///< The matrix (3,3) representing the incremental growth in current time step. Reset to identity at the beginning of each time step, updated in growth functions, and utilised to update Fg.
	double  columnarGrowthWeight;	///< The fraction defining how close to the columnar layer the element is. 1.0 for columnar layer, 0.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
	double  peripodialGrowthWeight;	///< The fraction defining how close to the peripodial membrane the element is. 0.0 for columnar layer, 1.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
	double* ShapeChangeRate;		///< Shape change rate of the elements, only orthagonal shape changes are allowed (x, y, z). Shape changes will be scaled to conserve volume, thus three values will not be independent.
    bool    rotatedGrowth;			///< The boolean stating if the element has rotated from the growth axis, hence the calculated growth requires further rotation to follow tissue axes.
	double* columnarRelativePosInBoundingBox;	///< The relative position on x-y plane, within the bounding box of the columnar layer (x,y).
	double* peripodialRelativePosInBoundingBox; ///< The relative position on x-y plane, within the bounding box of the peripodial membrane layer (x,y).
    gsl_matrix **ShapeFuncDerivatives;
    gsl_matrix **ShapeFuncDerStacks;
    gsl_matrix **InvdXdes;
    double* detdXdes;
    gsl_matrix **Bmatrices;
    gsl_matrix **CMatrices;
    gsl_matrix **FeMatrices;
    gsl_matrix **invJShapeFuncDerStack;
    gsl_matrix **invJShapeFuncDerStackwithFe;
    gsl_matrix **elasticStress;
	gsl_matrix* TriPointF;
    double* detFs;

    double ZProjectedBasalArea,ZProjectedApicalArea;
    double BasalArea, ApicalArea;
    double cMyoUniform[2]; //apical, basal
    double cMyoUnipolar[2]; //apical basal
    double cMyoUniformEq[2]; //apical, basal
    double cMyoUnipolarEq[2]; //apical basal
    gsl_matrix* myoPolarityDir;



    void 	setShapeType(string TypeName);
	void 	readNodeIds(int* tmpNodeIds);
	void 	setPositionMatrix(vector<Node*>& Nodes);
	void 	setTissuePlacement(vector<Node*>& Nodes);
	void 	setTissueType(vector<Node*>& Nodes);
	void 	setReferencePositionMatrix();
	void 	setIdentificationColour();
	void 	rotateReferenceElementByRotationMatrix(double* rotMat);
	bool 	InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);
    bool 	InvertMatrix(gsl_matrix* input, gsl_matrix* inverse);

    int 	determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm);

	double  dotProduct3D(double* u, double* v);
	void 	updateNodeIdsFromSave(ifstream& file);
	void 	updateReferencePositionMatrixFromSave(ifstream& file);
	virtual void calculateReferenceVolume(){ParentErrorMessage("calculateReferenceVolume");};
	bool 	calculateGrowthStrainsRotMat(double* v);
	void	calculateForces3D(double **SystemForces, vector <Node*>& Nodes,  bool recordForcesOnFixedNodes, double **FixedNodeForces, ofstream& outputFile);
	gsl_matrix* calculateEForNodalForces(gsl_matrix* F, gsl_matrix* Fe);
    gsl_matrix* calculateSForNodalForces(gsl_matrix* E);
    gsl_matrix* calculateCompactStressForNodalForces(gsl_matrix* Fe, gsl_matrix* S, gsl_matrix* FeT, gsl_matrix *Stress);
    gsl_matrix* calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian);
    gsl_matrix* calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix *B, gsl_matrix* invJShFuncDerS);
    gsl_matrix* calculateInvJShFuncDerSWithFe(gsl_matrix * currFe, gsl_matrix * InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithF);

    void	calculateCMatrix(int pointNo);
    void    consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b);
    void    calculateElasticKIntegral1(gsl_matrix* currK,int pointNo);
    void	calculateElasticKIntegral2(gsl_matrix* currK,int pointNo);

    bool 	disassembleRotationMatrixForZ(gsl_matrix* rotMat);
    bool 	calculate3DRotMatFromF(gsl_matrix* rotMat);

    gsl_matrix* D;
    gsl_matrix* CoeffMat;
    double D81[4][4][4][4];

    //boost::numeric::ublas::vector<double> Forces;

	double E, v;
    double lambda, mu;
    gsl_matrix* Fg;			///< Growth matrix
    gsl_matrix* InvFg;		///< Inverse of growth matrix
    gsl_matrix* Fsc;		///< Shape change matrix
    gsl_matrix* InvFsc;		///< Inverse of shape change matrix

    gsl_matrix* TriPointKe;

public:
     int 	Id;
	int		ShapeDim;
	int* 	NodeIds;
	virtual ~ShapeBase(){
			//while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
			};
	double** Positions;
	ReferenceShapeBase* ReferenceShape;
    gsl_matrix* Strain;

    //bool 	IsGrowing;
	bool 	IsChangingShape;
	bool	ApicalNormalForPackingUpToDate;
	bool	BasalNormalForPackingUpToDate;
	int 	tissuePlacement; //1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
	int 	tissueType;	///< The tissue type is 0 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
	bool	IsAblated;
	bool	IsClippedInDisplay;
	double 	CurrShapeChangeToAdd[3];
	double* ApicalNormalForPacking;
	double* BasalNormalForPacking;
    double  GrownVolume;
	double  VolumePerNode;
	bool 	capElement;
	double** MyoForce;

	int 	getId();
	string 	getName();
	int 	getShapeType();
	int 	getNodeNumber();
	int* 	getNodeIds();
	int 	getDim();
	int* 	getIdentifierColour();
	double* getCentre();
	double	getPeripodialness();
	double	getColumnarness();
	void	getRelativePositionInTissueInGridIndex(int nGridX, int nGridY, double*columnarReletivePos, double* peripodialRelativePos, int& columnarIndexX, int& peripodialIndexX, int& columnarIndexY, int& peripodialIndexY, double& columnarFracX, double& peripodialFracX, double& columnarFracY, double& peripodialFracY);
	bool 	isGrowthRateApplicable(int sourceTissue, double& weight);
	void 	calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue);
	void 	calculateFgFromGridCorners(double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue, int IndexX, int IndexY, double FracX, double dFracY);
	void 	updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial);
	void	calculateRelativePosInBoundingBox(double columnarBoundingBoxXMin, double columnarBoundingBoxYMin, double columnarBoundingBoxLength, double columnarBoundingBoxWidth, double peipodialBoundingBoxXMin, double peipodialBoundingBoxYMin, double peipodialBoundingBoxLength, double peipodialBoundingBoxWidth);
	void	displayRelativePosInBoundingBox();
	void	getRelativePosInColumnarBoundingBox(double* relativePos);
	void	getRelativePosInPeripodialBoundingBox(double* relativePos);
	void 	convertRelativePosToGridIndex(double* relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY);
	void 	getStrain(int type, float &StrainMag);
	void 	getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag);
    void 	getPysProp(int type, float &PysPropMag, double dt);
	double 	getYoungModulus();
	double 	getPoissonRatio();
	double* getGrowthRate();
	double* getShapeChangeRate();
	double** getReferencePos();
    void    getPos(gsl_matrix* Pos);
    gsl_matrix* getFg();
	void 	displayName();
	void	displayNodeIds();
	void 	displayPositions();
	void 	displayReferencePositions();
	void 	displayIdentifierColour();
    void    setFg(gsl_matrix* currFg);
	void 	setGrowthWeightsViaTissuePlacement (double periWeight);
    virtual void setElasticProperties(double EApical,double EBasal, double EMid,double v){ParentErrorMessage("setElasticProperties");};
	virtual void calculateBasalNormal(double * normal){ParentErrorMessage("calculateBasalNormal");};
	virtual void AlignReferenceBaseNormalToZ(){ParentErrorMessage("AlignReferenceBaseNormalToZ");};
	void 	calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix);
	void 	updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
	virtual void calculateReferenceStiffnessMatrix(){ParentErrorMessage("calculateReferenceStiffnessMatrix");};
    virtual void calculateElementShapeFunctionDerivatives(){ParentErrorMessage("calculateReferenceStiffnessMatrix");};
    virtual void calculateCurrNodalForces(gsl_matrix *gslcurrg, gsl_matrix *gslcurrF, int pointNo){ParentErrorMessage("calculateCurrNodalForces");};
    virtual void calculateCurrTriPointFForRotation(gsl_matrix *currF,int pointNo){ParentErrorMessage("calculateCurrTriPointFForRotation");};

    virtual void calculateApicalArea(){ParentErrorMessage("calculateApicalArea");};
    virtual void calculateBasalArea(){ParentErrorMessage("calculateBasalArea");};

    void 	calculateForces(double **SystemForces,  vector <Node*>& Nodes, bool recordForcesOnFixedNodes, double** FixedNodeForces,  ofstream& outputFile);
    void 	updatePositions(vector<Node*>& Nodes);
	void	updateReferencePositionsToCurentShape();
    void 	setGrowthRate(double dt, double x, double y, double z);
	void 	cleanMyosinForce();
	void	updateUniformEquilibriumMyosinConcentration(bool isApical, double cEqUniform);
	void	updateUnipolarEquilibriumMyosinConcentration(bool isApical, double cEqUnipolar, double orientationX, double orientationY);
	void	updateMyosinConcentration(double dt, double kMyo);
	bool 	calculateIfInsideActiveStripe(double initialPoint,double endPoint, double stripeSize1, double stripeSize2);
	double	getCmyosinUniformForNode (int TissuePlacement);
	double	getCmyosinUnipolarForNode (int TissuePlacement);
	void 	getMyosinLevels (double *cMyo);
	void 	getEquilibriumMyosinLevels (double *cMyoEq);
	void 	setMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1);
	void 	setEquilibriumMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1);
	virtual void calculateMyosinForces(double forcePerMyoMolecule){ParentErrorMessage("calculateMyosinForces");};
	virtual void distributeMyosinForce(bool isIsotropic, bool apical){ParentErrorMessage("distributeMyosin");};
	void 	setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
	void 	updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes);
	bool 	readNodeIdData(ifstream& file);
	bool	readReferencePositionData(ifstream& file);
	void 	convertPlasticStrainToGrowthStrain();
	virtual void  checkHealth(){ParentErrorMessage("checkHealth");};
	void 	resetCurrStepShapeChangeData();
    void    writeKelasticToMainKatrix(gsl_matrix* Ke);
    void    calculateImplicitKElastic();
    void	calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix* ExternalNodalForces);


	void 	updateShapeFromSave(ifstream& file);
	void 	displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
    void 	displayMatrix(gsl_matrix* mat, string matname);
    void 	displayMatrix(gsl_vector* mat, string matname);
    void createMatrixCopy(gsl_matrix *dest, gsl_matrix* src);
	double	calculateMagnitudeVector3D(double* v);
	void	normaliseVector3D(double* v);
	void	normaliseVector3D(gsl_vector* v);
	double	getNormVector3D(gsl_vector* v);
	double 	determinant3by3Matrix(double* rotMat);
	double 	determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat);
    double 	determinant3by3Matrix(gsl_matrix* Mat);
	double 	determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat);
	void	calculateRotationAngleSinCos(double* u, double* v, double& c, double& s);
	void	calculateRotationAxis(double* u, double* v,double* rotAx, double c);
	void	constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat);

	void    CalculateGrowthRotationByF();
	void 	calculateTriPointFForRatation();
    void 	growShapeByFg(double dt);
    void 	changeShapeByFsc(double dt);

	virtual double getApicalSideLengthAverage(){return ParentErrorMessage("getApicalSideLengthAverage",0.0);};
	virtual void getApicalTriangles(vector <int> &ApicalTriangles){ParentErrorMessage("getApicalTriangles");};
	virtual int getCorrecpondingApical(int currNodeId){return ParentErrorMessage("getCorrecpondingApical", -100);};
	virtual bool IsThisNodeMyBasal(int currNodeId){return ParentErrorMessage("IsThisNodeMyBasal", false);};
	virtual double getElementHeight(){return ParentErrorMessage("getElementHeight", 0.0);};
	virtual bool IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode){return ParentErrorMessage("IsPointCloseEnoughForPacking", false);};
	virtual void calculateNormalForPacking(int tissuePlacement){ParentErrorMessage("calculateNormalForPacking");};
	virtual void AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, double **SystemForces, double **PackingForces, vector<Node*> &Nodes){ParentErrorMessage("AddPackingToApicalSurface");};
	virtual void getApicalNodePos(double* posCorner){ParentErrorMessage("getApicalNodePos");};
	virtual void getBasalNodePos(double* posCorner){ParentErrorMessage("getBasalNodePos");};
	virtual bool IspointInsideTriangle(int tissueplacement, double x, double y,double z){return ParentErrorMessage("IspointInsideTriangle",false );};
	bool checkPackingToThisNodeViaState(int ColumnarLayerDiscretisationLAyers, Node* NodePointer);
	bool DoesPointBelogToMe(int IdNode);
	void growShape();
	void assignVolumesToNodes(vector <Node*>& Nodes);
	void assignSurfaceAreaToNodes(vector <Node*>& Nodes);
    void calculateZProjectedAreas();
    void assignZProjectedAreas(vector <Node*> Nodes);
	void assignElementToConnectedNodes(vector <Node*>& Nodes);
	void removeMassFromNodes(vector <Node*>& Nodes);


	void 	convertLocalStrainToTissueStrain(double* strainsToAdd);

	bool RotatedElement;
    gsl_matrix* GrowthStrainsRotMat;


	bool 	calculateAlignmentScore(double** RefNormalised);
	void 	bringShapePositionsToOrigin(double** RefNormalised, double* refCentre);
	void 	updateElementsNodePositions(int RKId, double ***SystemForces, vector <Node*>& Nodes, double dt);
	void 	updateReferencePositionMatrixFromMeshInput(ifstream& file);
	void	fillNodeNeighbourhood(vector<Node*>& Nodes);
	void 	checkDisplayClipping(double xClip, double yClip, double zClip);
	void	crossProduct3D(double* u, double* v, double* cross);
	void	crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross);
	void	alignGrowthCalculationOnReference();
	void	readNewGrowthRate(double* NewGrowth, double& ex, double&ey, double& ez, double& exy, double& exz, double& eyz);
	void	updateUniformOrRingGrowthRate(double* NewGrowth, int GrowthId);
	void	updateGridBasedGrowthRate(double* NewGrowth, int GrowthId, int i, int j);
};

#endif
