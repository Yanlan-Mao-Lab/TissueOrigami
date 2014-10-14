#include "ShapeBase.h"
#include <sstream>

void	ShapeBase::setShapeType(string TypeName){
	if (TypeName == "Prism"){
		this->ShapeType = 1;
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

int 	ShapeBase::getDim(){
	return nDim;
}

string 	ShapeBase::getName(){
	string name;
	if (ShapeType == 1){
		name = "Prism";
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

void 	ShapeBase::readNodeIds(int* inpNodeIds){
	for (int i=0; i<nNodes; ++i){
		this->NodeIds[i] = inpNodeIds[i];
	}
}

void 	ShapeBase::displayName(){
	cout<<"Type: "<<this->ShapeType<<" Id: "<<this->Id<<endl;
}

void 	ShapeBase::setPositionMatrix(vector<double*>& Nodes){
	const int n = nNodes;
	const int dim = nDim;

	Positions = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		Positions[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]][j];
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

void 	ShapeBase::displayIdentifierColour(){
	cout <<" IdentifierColour:  "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}

void 	ShapeBase::alignReference(){
	//double CM[3] = {0.0,0.0,0.0};
	//double CMRef[3] = {0.0,0.0,0.0};
	double translation[3] = {0.0,0.0,0.0};
	cout<<"translation: "<<translation[0]<<" "<<translation[1]<<" "<<translation[1]<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			//CMRef[j] += Positions[i][j];
			//CMRef[j] += ReferenceShape -> Positions[i][j];
			translation[j] += ReferenceShape -> Positions[i][j] - Positions[i][j];
		}
	}
	for (int j = 0; j<nDim; ++j){
		translation[j] = translation[j]/nNodes;
	}
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			ReferenceShape -> Positions[i][j] -= translation[j];
		}
	}
}

void	ShapeBase::calculateForces(int RKid, double **SystemForces){
	const int nMult = nNodes*nDim;

	using namespace boost::numeric::ublas;
	boost::numeric::ublas::vector<double> displacement(nMult);
	int counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			displacement(counter) = Positions[i][j]-ReferenceShape->Positions[i][j];
			counter++;
		}
	}

	Forces = zero_vector<double>(nMult);
	boost::numeric::ublas::axpy_prod(k,displacement,Forces);
	counter = 0;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			if (i>2 || j!= 2 || FixedLowerZ == false){
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
	for (int i=0;i<8;++i){
		for (int j=0;j<3;++j){
			cout<<SystemForces[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;*/
}

void	ShapeBase::updatePositions(vector<double*>& Nodes){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
				Positions[i][j] = Nodes[NodeIds[i]][j];
		}
	}
}

void 	ShapeBase::growShape(float scale){
	//cout<<"growing shape: "<<Id<<endl;
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			ReferenceShape -> Positions[i][j] = ReferenceShape -> Positions[i][j] * scale;
		}
	}
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
	for(int i = 0; i < A.size1(); i++) {
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
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
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
