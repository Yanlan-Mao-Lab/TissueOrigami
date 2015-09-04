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
	void ParentErrorMessage(string functionName);
	bool ParentErrorMessage(string functionName, bool returnValue);
	double ParentErrorMessage(string functionName, double returnValue);
	int ParentErrorMessage(string functionName, int returnValue);
protected:
	int 	ShapeType;

	int 	nNodes;

	int 	nDim;
	int* 	IdentifierColour;
	double* GrowthRate;
	double  columnarGrowthWeight;
	double  peripodialGrowthWeight;
	double* ShapeChangeRate;
    bool    rotatedGrowth;
	double* CurrGrowthStrainAddition;
	bool 	CurrShapeChangeStrainsUpToDate;
	bool 	CurrGrowthStrainsUpToDate;
	bool 	GrewInThePast;
	bool 	ChangedShapeInThePast;
	double* RelativePosInBoundingBox;
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
    double* detFs;

    double ZProjectedBasalArea,ZProjectedApicalArea;
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
	void	calculateForces3D(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile);
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
    gsl_matrix* TriPointF;
    gsl_matrix* TriPointKe;

public:
    //boost::numeric::ublas::matrix<double> B;
    //boost::numeric::ublas::matrix<double> BE;
    int 	Id;
	int		ShapeDim;
    //boost::numeric::ublas::matrix<double> LocalGrowthStrainsMat;
	int* 	NodeIds;
	virtual ~ShapeBase(){
			//while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
			};
	double** Positions;
	ReferenceShapeBase* ReferenceShape;
    gsl_matrix* Strain;
    gsl_matrix* RK1Strain;

    bool 	IsGrowing;
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
    double GrownVolume;
	double VolumePerNode;
	bool capElement;

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
	void	calculateRelativePosInBoundingBox(double BoindingBoxXMin, double BoundingBoxYMin, double BoundingBoxLength, double BoundingBoxWidth);
	void	displayRelativePosInBoundingBox();
	void	getRelativePosInBoundingBox(double* relativePos);
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
	void 	updateGrowthRate(double scalex, double scaley, double scalez);
	void 	updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
	virtual void calculateReferenceStiffnessMatrix(){ParentErrorMessage("calculateReferenceStiffnessMatrix");};
    virtual void calculateElementShapeFunctionDerivatives(){ParentErrorMessage("calculateReferenceStiffnessMatrix");};
    virtual void calculateCurrNodalForces(gsl_matrix *gslcurrg, gsl_matrix *gslcurrF, int pointNo){ParentErrorMessage("gslcalculateCurrNodalForces");};


    void 	calculateForces(int RKId, double ***SystemForces, vector <Node*>& Nodes, ofstream& outputFile);
    void 	updatePositions(int RKId, vector<Node*>& Nodes);
	void 	setGrowthRate(double x, double y, double z);
	void 	setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
	void 	updateGrowthToAdd(double* growthscale);
	void 	updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes);
	bool 	readNodeIdData(ifstream& file);
	bool	readReferencePositionData(ifstream& file);
	void 	convertPlasticStrainToGrowthStrain();
	virtual void  checkHealth(){ParentErrorMessage("checkHealth");};
	void 	resetCurrStepGrowthData();
	void 	resetCurrStepShapeChangeData();
    void    writeKelasticToMainKatrix(gsl_matrix* Ke);
    void    calculateImplicitKElastic();
    void	calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix* ExternalNodalForces);


	//void 	calculatePlasticStrain();
	//void 	changeShape(double shapechangescale, int axis);
	void 	updateShapeFromSave(ifstream& file);
	void 	displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
    void 	displayMatrix(gsl_matrix* mat, string matname);
    void 	displayMatrix(gsl_vector* mat, string matname);
    void createMatrixCopy(gsl_matrix *dest, gsl_matrix* src);
	double	calculateMagnitudeVector3D(double* v);
	void	normaliseVector3D(double* v);
	double 	determinant3by3Matrix(double* rotMat);
	double 	determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat);
    double 	determinant3by3Matrix(gsl_matrix* Mat);
	double 	determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat);
	void	calculateRotationAngleSinCos(double* u, double* v, double& c, double& s);
	void	calculateRotationAxis(double* u, double* v,double* rotAx, double c);
	void	constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,double* rotMat);

	void    CalculateGrowthRotationByF();
    void 	growShapeByFg(double dt);
    void 	changeShapeByFsc(double dt);

	virtual double getApicalSideLengthAverage(){return ParentErrorMessage("getApicalSideLengthAverage",0.0);};
	virtual void getApicalTriangles(vector <int> &ApicalTriangles){ParentErrorMessage("getApicalTriangles");};
	virtual int getCorrecpondingApical(int currNodeId){return ParentErrorMessage("getCorrecpondingApical", -100);};
	virtual bool IsThisNodeMyBasal(int currNodeId){return ParentErrorMessage("IsThisNodeMyBasal", false);};
	virtual double getElementHeight(){return ParentErrorMessage("getElementHeight", 0.0);};
	virtual bool IsPointCloseEnoughForPacking(double* Pos,  float Peripodialthreshold, float Columnarthreshold, int TissuePlacementOfPackingNode, int TissueTypeOfPackingNode){return ParentErrorMessage("IsPointCloseEnoughForPacking", false);};
	virtual void calculateNormalForPacking(int tissuePlacement){ParentErrorMessage("calculateNormalForPacking");};
	virtual void AddPackingToSurface(int tissueplacement, double Fx, double Fy,double Fz, int RKId,  double ***SystemForces, double ***PackingForces, vector<Node*> &Nodes){ParentErrorMessage("AddPackingToApicalSurface");};
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


	//bool 	calculateAlignmentRotationMatrix(double** RefNormalised, double* rotMat);
	bool 	calculateAlignmentScore(double** RefNormalised);
	void 	bringShapePositionsToOrigin(double** RefNormalised, double* refCentre);
	void 	updateElementsNodePositions(int RKId, double ***SystemForces, vector <Node*>& Nodes, double dt);
	void 	updateReferencePositionMatrixFromMeshInput(ifstream& file);
	void	fillNodeNeighbourhood(vector<Node*>& Nodes);
	void 	checkDisplayClipping(double xClip, double yClip, double zClip);
	void	crossProduct3D(double* u, double* v, double* cross);
	void	alignGrowthCalculationOnReference();
	void	readNewGrowthRate(double* NewGrowth, double& ex, double&ey, double& ez, double& exy, double& exz, double& eyz);
	void	updateUniformOrRingGrowthRate(double* NewGrowth, int GrowthId);
	void	updateGridBasedGrowthRate(double* NewGrowth, int GrowthId, int i, int j);
};

#endif
