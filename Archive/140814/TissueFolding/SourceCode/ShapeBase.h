#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Node.h"
#include "ReferenceShapeBase.h"

using namespace std;

class ShapeBase{
private:
	void ParentErrorMessage();
protected:
	int 	ShapeType;
	int 	Id;
	int 	nNodes;
	int* 	NodeIds;
	int 	nDim;
	ReferenceShapeBase* ReferenceShape;
	int* 	IdentifierColour;
	double* GrowthRate;

	void 	setShapeType(string TypeName);
	void 	readNodeIds(int* tmpNodeIds);
	void 	setPositionMatrix(vector<Node*>& Nodes);
	void 	setReferencePositionMatrix();
	virtual void setRefShapePosBuffers(){ParentErrorMessage();};
	void 	setIdentificationColour();
	bool 	InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse, double& det);
	int 	determinant_sign(boost::numeric::ublas::permutation_matrix<std::size_t>& pm);
	void	crossProduct3D(double* u, double* v, double* cross);
	double  dotProduct3D(double* u, double* v);
	void	normaliseVector3D(double* v);
	void	calculateRotationQuaternian(double* u, double* v, double* Qrot);
	void	NormaliseQuaternian(double* Qrot);
	void	RotateVectorByQuaternian(double* u, double* Qrot);
	virtual void calculateNormals(){ParentErrorMessage();};
	virtual void setNormals(){ParentErrorMessage();};
	virtual void getCurrentAlignmentSides(double*, double*){ParentErrorMessage();};
	virtual void updateAlignmentTurn(){ParentErrorMessage();};
	virtual void updateReferenceShapeBaseFromBuffer(){ParentErrorMessage();};
	void 	updateNodeIdsFromSave(ifstream& file);
	void 	updateReferencePositionMatrixFromSave(ifstream& file);
	virtual void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace){ParentErrorMessage();};
	void	calculateRotationAngleSinCos(double* u, double* v, double& c, double& s);
	void	calculateRotationAxis(double* u, double* v,double* rotAx);
	void	constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,double* rotMat);
	bool	areSidesFacingSameDirection(double* RefSide, double* ShapeSide);
	boost::numeric::ublas::matrix<double> D;
	boost::numeric::ublas::matrix<int> CoeffMat;
	boost::numeric::ublas::matrix<double> k;
	boost::numeric::ublas::matrix<double> B;
	boost::numeric::ublas::vector<double> Forces;
	double E, v;

public:
	virtual ~ShapeBase(){
			//while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
			};

	double** Positions;
	boost::numeric::ublas::vector<double> Strain;
	int alingmentTurn;
	double* CurrentNormal;
	bool updatedReference;
	int getId();
	string getName();
	int getShapeType();
	int getNodeNumber();
	int* getNodeIds();
	int getDim();
	int* getIdentifierColour();
	double* getCentre();
	void getStrainColour(int type, float* StrainColour);
	void getStrain(int type, float &StrainMag);
	void getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag);
	void getPysProp(int type, float &PysPropMag);
	double getYoungModulus();
	double getPoissonRatio();
	double* getGrowthRate();
	double** getReferencePos();
	double*	 getReferenceNormal();
	void displayName();
	void displayPositions();
	void displayIdentifierColour();
	virtual void setViscosity(double ApicalVisc,double BasalVisc, vector <Node*>& Nodes){ParentErrorMessage();};
	virtual void setElasticProperties(double E,double v){ParentErrorMessage();};
	void updateGrowthRate(double scalex, double scaley, double scalez);
	virtual void calculateReferenceStiffnessMatrix(){ParentErrorMessage();};
	virtual void setStiffnessMatrixBuffers(){ParentErrorMessage();};
	void alignReference();
	void calculateForces(int RKid, double **SystemForces, vector <Node*>& Nodes);
	void updatePositions(vector<Node*>& Nodes);
	void setGrowthRate(double x, double y, double z);
	void growShape(double scalex, double scaley, double scalez);
	void updateShapeFromSave(ifstream& file);
	void displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
};

#endif
