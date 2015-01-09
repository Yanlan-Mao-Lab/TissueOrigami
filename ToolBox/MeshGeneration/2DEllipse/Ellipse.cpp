using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include "math.h"
#include <vector>

class EllipseLayoutGenerator{
public:
	EllipseLayoutGenerator(double r1, double r2, double sideLen, double condensation);
	double pi;
	double r1;
	double r2;
	double r1Init;
	double r2Init;
	double condensation;
	double sideLen;
	double sideLenInit; 
	double curvedEndQuarterCircum;
	vector <double> posx,posy,posz;
	vector <double> InvardNormalsX, InvardNormalsY;
	vector <double> zVec0,zVec1,zVec2;
	vector <int> tissueType; //0: columnar layer, 1: peripodium
	vector <int> atBorder;	//1 id the node is at the circumference of the setup, 0 otherwise
	double tetha;
	double dtet;
	double Circumference;
	int Vborder;
	bool tethaAtZero;
	bool r1failed;	
	vector <int> Links0,Links1;	
	vector <int*>triangles;
	double dForNodeGeneration;
	
	void calculateCircumference();
	void calculateCurrentBorderNumber();
	void calculatedtet();
	void initiateTetha();
	void addRing();
	bool updateRadia();
	void addMidLine();
	void minimise(int n_iter);
	void minimiseWithTesselation2D();
	void minimiseWithTesselation3D();
	void Tesselate2D();
	void Tesselate3D();
	void readInTesselation2D();
	void readInTesselation3D();
	void minimisewitSpringMovement(int n_iter);
	void writeVectors2D(ofstream &vectorsForGnuplot);
	void writeNodes2D(ofstream &nodesForGnuplot);
	void writeMeshFileForSimulation(double zHeight, int zLayers);
	void writeTriangleMeshFileForSimulation(double zHeight);
	void addEquidistantRing();
	void addEquidistantRingMin();
	void addRectangle();
	void calculateAverageSideLength();
	void addZVecs(double curvedEndRadia, int nCurveLayers);
	void addInitialPosZ(double curvedEndRadia, int nCurveLayers);
	int calculateLargeRadiaAndLayers(double curvedEndRadia);
	void addZPosToCurve(int nLast, int counter, double dZ);
	void correctNormalsForAddedPoints();
};

EllipseLayoutGenerator::EllipseLayoutGenerator(double r1, double r2, double sideLen, double condensation){
	pi = 3.14;
	this->r1 = r1;
	this->r2 = r2;
	this->r1Init = r1;
	this->r2Init = r2;
	this->condensation = condensation;
	this->sideLen = sideLen*condensation;
	this->sideLenInit =this->sideLen;
	tethaAtZero=true;
	r1failed = false;
	curvedEndQuarterCircum = 0.0;
}
void EllipseLayoutGenerator::calculateCircumference(){
	double h = ((r1-r2)*(r1-r2))/((r1+r2)*(r1+r2));
	Circumference = pi * (r1 + r2) * (1.0 + (3*h / (10.0 + pow(4.0 - 3*h,0.5) )));  
}
void EllipseLayoutGenerator::calculateCurrentBorderNumber(){
	Vborder = Circumference / (sideLen*1.8);
}
void EllipseLayoutGenerator::calculatedtet(){
	int n  = Vborder/2;	
	dtet = pi/n;
}
void EllipseLayoutGenerator::initiateTetha(){
	if(tethaAtZero){
		tetha =0.0;
		tethaAtZero=false;
	}
	else{
		tetha =dtet/2.0;
		tethaAtZero=true;	
	}
}

void EllipseLayoutGenerator::addEquidistantRingMin(){	
	double CQ = Circumference/4.0;
	int n= CQ / sideLen;
	//n -=1;
	dtet= pi/2.0/n;
	double d = CQ / n;
	dForNodeGeneration  = d;
	//cerr<<"dForNodeGeneration: "<<dForNodeGeneration<<endl;
	double x0 = r1;
	double y0 = 0;
	vector <double> CurrPosX, CurrPosY;
	CurrPosX.push_back(x0);
	CurrPosY.push_back(y0);
	tetha = dtet;
	for (int i= 0; i<n-1; ++i){
		double x = r1*cos(tetha);
		double y = r2*sin(tetha);
		CurrPosX.push_back(x);
		CurrPosY.push_back(y);
		tetha += dtet;
	}
	CurrPosX.push_back(0);
	CurrPosY.push_back(r2);
	double ksp = 0.2;
	int niteration = 100000;
	for (int a = 0; a< niteration; ++a){
		if (a < 100) {ksp = 0.2;}
		else if (a<500) {ksp = 0.1;}
		else if (a<5000)  {ksp =0.05;}
		else if (a<10000)  {ksp =0.01;}
		else if (a<20000)  {ksp =0.005;}
		else  {ksp =0.001;}
		for (int i= 1; i<n; ++i){
			double vec[2]={0.0,0.0};
			double v[2]={0.0,0.0};
			//pre neig:
			vec[0] = CurrPosX[i-1] - CurrPosX[i];
			vec[1] = CurrPosY[i-1] - CurrPosY[i];	
			double dsq = vec[0]*vec[0] + vec[1]*vec[1];
			double dmag = pow(dsq,0.5);
			double F =  (dmag-d)*ksp;
			v[0] = F*vec[0]/dmag;
			v[1] = F*vec[1]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;
			//post neig:
			vec[0] = CurrPosX[i+1] - CurrPosX[i];
			vec[1] = CurrPosY[i+1] - CurrPosY[i];	
			dsq = vec[0]*vec[0] + vec[1]*vec[1];
			dmag = pow(dsq,0.5);
			//cerr<<"dmag of node : "<<i<<" "<<dmag<<endl;	
			F =  (dmag-d)*ksp;
			v[0] += F*vec[0]/dmag;
			v[1] += F*vec[1]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;
			//add energy for deviating from ellipse:
			dsq = CurrPosX[i]*CurrPosX[i] + CurrPosY[i]*CurrPosY[i];
			dmag = pow(dsq,0.5);
			double ellipseSum = (CurrPosX[i]/r1) * (CurrPosX[i]/r1) +  (CurrPosY[i]/r2) * (CurrPosY[i]/r2);
			if (ellipseSum > 1){
				ellipseSum = (-1.0)*pow(( ellipseSum-1),0.5);		
			}
			else{
				ellipseSum = pow(( 1.0 - ellipseSum),0.5);	
			}
			F = ellipseSum*ksp;
			v[0] += F*CurrPosX[i]/dmag;
			v[1] += F*CurrPosY[i]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;	
			CurrPosX[i] += v[0];
			CurrPosY[i] += v[1];
		}
	}
	//now cheking iof any of the non-neighbouring nodes are too close
	// I will remove them from the list, therefore I am nt checking first or last ones:
	double threshold = d*0.7;
	threshold *= threshold;	
	n = CurrPosX.size();
	int i=1;	
	while (i < n-1){
		bool deleteNode = false;
		for (int j=0; j< n; ++j){
			if (CurrPosX[i]<0 || CurrPosY[i]<0){deleteNode = true;}
			if (i != j && deleteNode ==false){
				double dx = CurrPosX[i] - CurrPosX[j];
				double dy = CurrPosY[i] - CurrPosY[j];
				double d2 = dx*dx + dy*dy;
				//cerr<<"i: "<<i<<" j "<<j<<" d: "<<d2 <<" threshold : "<<threshold <<" pos: "<<CurrPosX[i]<<" "<<CurrPosY[i]<<endl;
				if (d2 < threshold){deleteNode = true;}
			}
			if(deleteNode){
				vector<double>::iterator it;
  				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i = 0;		
				break;				
			}			
		}
		i++;	
	}
	//now doing the same check with already added points:
	n = CurrPosX.size();	
	i=1;	
	while (i < n-1){
		for (int j=0; j< posx.size(); ++j){
			double dx = CurrPosX[i] - posx[j];
			double dy = CurrPosY[i] - posy[j];
			double d2 = dx*dx + dy*dy;
			//cerr<<"i: "<<i<<" j "<<j<<" d: "<<d2 <<endl;
			if (d2 < threshold){
				vector<double>::iterator it;
				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i =0;			
				break;				
			}		
		}
		i++;	
	}
		
	//now rotating to add the remaining points:
	n = CurrPosX.size();
	for (int i = n-2; i>-1; --i){
		CurrPosX.push_back((-1.0)*CurrPosX[i]);
		CurrPosY.push_back(CurrPosY[i]);	
	}
	n = CurrPosX.size();
	for (int i = n-2; i>0; --i){
		CurrPosX.push_back(CurrPosX[i]);
		CurrPosY.push_back((-1.0)*CurrPosY[i]);		
	}
	n = CurrPosX.size();
	//calculating the normals pointing towards the inside of ellipse
	for (int i = 0; i<n; ++i){
		int index0 = i-1;
		int index1 = i+1;
		if (index0 == -1) { index0=n-1; }
		if (index1 ==  n) { index1 =0; }
		double vec1[2] = {CurrPosX[index0] - CurrPosX[i], CurrPosY[index0] - CurrPosY[i] };
		double vec2[2] = {CurrPosX[index1] - CurrPosX[i], CurrPosY[index1] - CurrPosY[i] };
		double mag = vec1[0] *vec1[0] + vec1[1]*vec1[1];
		if (mag> 1E-14){mag = pow(mag,0.5); vec1[0] /= mag; vec1[1] /= mag;}
		mag = vec2[0] *vec2[0] + vec2[1]*vec2[1];
		if (mag> 1E-14){mag = pow(mag,0.5); vec2[0] /= mag; vec2[1] /= mag;}
		vec1[0] += vec2[0]; vec1[1] += vec2[1];
		mag = vec1[0]*vec1[0] + vec1[1]*vec1[1];
		if (mag> 1E-14){mag = pow(mag,0.5); vec1[0] /= mag; vec1[1] /= mag;}
		InvardNormalsX.push_back(vec1[0]);
		InvardNormalsY.push_back(vec1[1]);
		//cout<<" point: "<<InvardNormalsX.size()-1<<" position: "<<CurrPosX[i]<<" "<<CurrPosY[i]<<"neigs: "<<CurrPosX[index0]<<" "<<CurrPosY[index0]<<"  "<<CurrPosX[index1]<<" "<<CurrPosY[index1]<<" normal: "<<InvardNormalsX[i]<<" "<<InvardNormalsY[i]<<endl;
	}
	n = CurrPosX.size();
	for (int i = 0; i<n; ++i){
		posx.push_back(CurrPosX[i]);
		posy.push_back(CurrPosY[i]);
	}
	//for (int i = 0; i<posx.size(); ++i){
	//	cout<<" point: "<<i<<" position: "<<posx[i]<<" "<<posy[i]<<" normal: "<<InvardNormalsX[i]<<" "<<InvardNormalsY[i]<<endl;
	//}
}

void EllipseLayoutGenerator::addRectangle(){	
	int n = r2 / sideLen;
	double d = r2 / n;
	dForNodeGeneration  = d;
	double x0 = r1;
	double y0 = 0;
	vector <double> CurrPosX, CurrPosY;
	CurrPosX.push_back(x0);
	CurrPosY.push_back(y0);
	tetha = dtet;
	for (int i= 0; i<n-1; ++i){
		double x = r1;
		double y = (i+1)*d;
		CurrPosX.push_back(x);
		CurrPosY.push_back(y);
	}
	CurrPosX.push_back(r1);
	CurrPosY.push_back(r2);
	n = r1 / sideLen;
	d = r1 / n;
	for (int i= 0; i<n-1; ++i){
		double x = r1 - (i+1)*d;
		double y = r2;
		CurrPosX.push_back(x);
		CurrPosY.push_back(y);
	}
	CurrPosX.push_back(0);
	CurrPosY.push_back(r2);
	
	//now rotating to add the remaining points:
	n = CurrPosX.size();
	for (int i = n-2; i>-1; --i){
		CurrPosX.push_back((-1.0)*CurrPosX[i]);
		CurrPosY.push_back(CurrPosY[i]);	
	}
	n = CurrPosX.size();
	for (int i = n-2; i>0; --i){
		CurrPosX.push_back(CurrPosX[i]);
		CurrPosY.push_back((-1.0)*CurrPosY[i]);		
	}
	n = CurrPosX.size();
	for (int i = 0; i<n; ++i){
		posx.push_back(CurrPosX[i]);
		posy.push_back(CurrPosY[i]);
	}
}



void EllipseLayoutGenerator::addRing(){
	double initialY = r2*sin(tetha);
	int firsindex = posx.size();
	while (tetha <= pi - 0.75*dtet){
	//while (tetha <= pi ){
		cerr<<"tetha: "<<tetha<<" dtet: "<<dtet<<" Vborder: "<<Vborder<<endl;
		double x = r1*cos(tetha);
		double y = r2*sin(tetha);
		posx.push_back(x);
		posy.push_back(y);
		tetha += dtet;
	}
	//adding the corner
	posx.push_back(-r1);
	posy.push_back(initialY);
	//taking mirror image:
	int n = posx.size();
	for (int i=firsindex+1; i<n-1;++i){
		posx.push_back(posx[i]);
		posy.push_back(-1.0*posy[i]);	
	}	
}

bool EllipseLayoutGenerator::updateRadia(){	
	if (r1>sideLen*2.0 && r2> sideLen*2.0){
		r1 -=sideLen;
		r2 -=sideLen;
		return true; 	
	}
	else{
		if (r1<=sideLen*1.5){
			r1failed = true;
		}
		return false;	
	}
}
void EllipseLayoutGenerator::addMidLine(){
	vector <double> CurrPosX, CurrPosY;
	CurrPosX.push_back(0.0);	
	CurrPosY.push_back(0.0);
	if (r1failed){
		//r1 was the shorter axe, I will add the midline along r2 (vertical)	:
		double y = sideLen;
		while (y<(r2-0.75*sideLen)){
			CurrPosX.push_back(0.0);
			CurrPosY.push_back(y);	
			CurrPosX.push_back(0.0);
			CurrPosY.push_back((-1.0)*y);
			y +=sideLen;	
		}
	}
	else{
		double x = sideLen;
		while (x<(r1-0.75*sideLen)){
			CurrPosX.push_back(x);
			CurrPosY.push_back(0.0);
			CurrPosX.push_back((-1.0)*x);
			CurrPosY.push_back(0.0);
			x +=sideLen;	
		}
	}
	//now eliminating any node that may be too close to already existing nodes:
	double threshold = dForNodeGeneration*0.7;
	threshold *= threshold;	
	int n = CurrPosX.size();	
	int i=1;	
	while (i < n-1){
		for (int j=0; j< posx.size(); ++j){
			double dx = CurrPosX[i] - posx[j];
			double dy = CurrPosY[i] - posy[j];
			double d2 = dx*dx + dy*dy;
			if (d2 < threshold){
				vector<double>::iterator it;
				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i--;			
				break;				
			}		
		}
		i++;	
	}
	n = CurrPosX.size();
	for (int i=0;i<n; ++i){
		posx.push_back(CurrPosX[i]);
		posy.push_back(CurrPosY[i]);
		InvardNormalsX.push_back(0.0);
		InvardNormalsY.push_back(0.0);
	}


}

void EllipseLayoutGenerator::minimise(int n_iter){
	double dt = 0.1;
	int n = posx.size();
	r1 = r1Init;
	r2 = r2Init;
	sideLen = sideLen/condensation;
	double r1sq = r1*r1;
	double r2sq = r2*r2;
	double sidesq = sideLen*sideLen;
	for (int i=0;i<n_iter; ++i){
		cerr<<"iteration: "<<i<<endl;		
		for (int j=0; j<n; ++j){
			double v[2] = {0.0,0.0};
			double d1 = posx[j]*posx[j] / r1sq  + posy[j]*posy[j]/r2sq;
			for (int k=0; k<n; ++k){
				if (j != k){
					double dx = posx[j] - posx[k];
					double dy = posy[j] - posy[k];
					double dsq = dx*dx + dy*dy;
					if (dsq < sidesq){
						double d = pow(dsq,0.5);
						double F = ( (sideLen - d)/sideLen ) / r2sq;
						v[0] += dx/d*F;	
						v[1] += dy/d*F;					
					}				
				}			
			}
			//calculated th pushing effect on the node, as its velocity v;
			double pos[2] = {posx[j]+v[0]*dt, posy[j]+v[1]*dt};
			double d2 = pos[0]*pos[0]/r1sq + pos[1]*pos[1] / r2sq;
			bool accept = false;
			if (d1>1){
				//the point was outside the ellipse
				if(d2<d1){
					//the new position is closer to ellipse:
					accept = true;				
				}			
			}
			else if (d2<1){
				//the point was inside the ellipse, the new position is also inside
				accept = true;			
			}	
			if (accept){
				posx[j] = pos[0];
				posy[j] = pos[1];			
			}
		}	
	}
}

void EllipseLayoutGenerator::minimiseWithTesselation2D(){
	Tesselate2D();
	readInTesselation2D();
	minimisewitSpringMovement(100);
}

void EllipseLayoutGenerator::minimiseWithTesselation3D(){
	//Tesselate3D();
	readInTesselation3D();
	//minimisewitSpringMovement(10);
}

void EllipseLayoutGenerator::Tesselate3D(){
	ofstream pointsForTesselation;
	pointsForTesselation.open("./Points.in",ofstream::trunc);
	for (int i =0; i< posx.size(); ++i){
		pointsForTesselation<<posx[i]<<" "<<posy[i]<<" 0.0 0.0 0.0 1.0"<<endl;
	}
	system("../PoissonRecon/Bin/PoissonRecon --in ./Points.in --out ./Points-out ");
	system("../PoissonRecon/plytomesh.pl Points-out.ply > Points.out");
}

void EllipseLayoutGenerator::Tesselate2D(){
	ofstream vectorsForGnuplot,nodesForGnuplot;;
	nodesForGnuplot.open("./Nodes1.out",ofstream::trunc);	
	writeNodes2D(nodesForGnuplot);
	ofstream pointsForTesselation;
	pointsForTesselation.open("./Points.node",ofstream::trunc);
	pointsForTesselation<<posx.size()<<" 2 0 1"<<endl; //dim, attribute number , border markers on or of
	for (int i =0; i< posx.size(); ++i){
		pointsForTesselation<<i<<" "<<posx[i]<<" "<<posy[i]<<endl;
	}
	double maxArea = (sideLen/condensation) * (sideLen/condensation)*pow(3.0,0.5)/2.0 * 1/2.0;
	maxArea = maxArea*1.5;
	cerr<<"Max Area: "<<maxArea<<endl;


	ostringstream Convert;
	Convert << maxArea; // Use some manipulators
	string maxAreaStr = Convert.str(); // Give the result to the string
	string sysCommand = "/home/melda/Desktop/MeshGenerate/triangle/triangle -q33a"+maxAreaStr+" ./Points.node  ";
	cerr<<"Running triangulation with: "<<sysCommand<<endl;
	system(sysCommand.c_str());
}
void EllipseLayoutGenerator::readInTesselation2D(){
	ifstream MeshFile;	
	MeshFile.open("./Points.1.ele", ifstream::in);
	//readHeader:
	int ntri;
	MeshFile>>ntri;
	int nodesPerTri;
	MeshFile>>nodesPerTri;
	int nAttributes;
	MeshFile>>nAttributes;
	for (int i=0; i<ntri; ++i){
		int triId;
		MeshFile>>triId;
		int* pnts;
		pnts = new int[3];
		for (int j=0;j<3;++j){
			MeshFile >> pnts[j];
		}
		triangles.push_back(pnts);
		Links0.push_back(pnts[0]);
		Links1.push_back(pnts[1]);
		Links0.push_back(pnts[1]);
		Links1.push_back(pnts[2]);
		Links0.push_back(pnts[2]);
		Links1.push_back(pnts[0]);
	}
	MeshFile.close();
	ifstream NodeFile;	
	NodeFile.open("./Points.1.node", ifstream::in);
	int nNode;
	int bordersMarked;
	NodeFile>>nNode;
	NodeFile>>nodesPerTri;	//dimensions
	NodeFile>>nAttributes;	//attribute number
	NodeFile>>bordersMarked;	//border markers on or off
	for (int i=0; i<nNode; ++i){
		int nodeId;
		NodeFile>>nodeId;
		//cerr<<"Reading node : "<<i<<"("<<nodeId<<") of "<<nNode<<endl;
		double pnts[2];
		for (int j=0;j<2;++j){
			NodeFile >> pnts[j];	
		}
		//reading the flag for border nodes, and generating the vector for it
		NodeFile>>bordersMarked;
		atBorder.push_back(bordersMarked);
		//cerr<<"		read pos: "<< pnts[0]<<" "<<pnts[1]<<" borders? "<<bordersMarked<<endl;
		if (nodeId<posx.size()){
			//cerr<<"		overwriting existing placement on vector"<<endl;
			posx[nodeId]=pnts[0];	
			posy[nodeId]=pnts[1];		
		}
		else{
			//cerr<<"		adding points to system"<<endl;
			posx.push_back(pnts[0]);	
			posy.push_back(pnts[1]);
			//will correct these normals later on
			InvardNormalsX.push_back(-100.0);
			InvardNormalsY.push_back(-100.0);		
		}
	}
	NodeFile.close();
	//writing for gnuplot vectors:
	ofstream vectorsForGnuplot,nodesForGnuplot;;
	vectorsForGnuplot.open("./Vectors1.out",ofstream::trunc);
	nodesForGnuplot.open("./Nodes2.out",ofstream::trunc);
	writeVectors2D(vectorsForGnuplot);	
	writeNodes2D(nodesForGnuplot);
	//cout<<"wrote nodes for gnuplot"<<endl;
}

void EllipseLayoutGenerator::writeVectors2D(ofstream &vectorsForGnuplot){
	for (int i =0; i< Links1.size(); ++i){
		double x = posx[Links0[i]];
		double y = posy[Links0[i]];
		double dx =  posx[Links1[i]] - x;
		double dy =  posy[Links1[i]] - y;
		vectorsForGnuplot<<x<<"  "<<y<<" "<<dx<<" "<<dy<<endl;
	}	
}
void EllipseLayoutGenerator::writeNodes2D(ofstream &nodesForGnuplot){
	for (int i =0; i< posx.size(); ++i){
		double x = posx[i];
		double y = posy[i];
		nodesForGnuplot<<x<<"  "<<y<<endl;
	}	
}

void EllipseLayoutGenerator::readInTesselation3D(){
	cout<<"reading file"<<endl;
	ifstream PLYMeshFile;
	PLYMeshFile.open("./Points.out", ifstream::in);
	if (!(PLYMeshFile.good() && PLYMeshFile.is_open())){
		cerr<<"Cannot open the save file to display:  ./Points.out"<<endl;
	}
	//readHeader:
	string currline;
	//skipping the header:
	getline(PLYMeshFile,currline);
	while(currline[0] == '#'){
		getline(PLYMeshFile,currline);
	}
	//Now I hold the first vertex line:
	istringstream currSStrem(currline);
	string InpType;		//is the input a vertex or a face
	currSStrem >> InpType;
	int tmp_int;
	int counter =0;
	while(InpType == "Vertex" && !PLYMeshFile.eof()){
		istringstream currSStrem(currline);
		currSStrem >> InpType;
		currSStrem >> tmp_int;
		double x,y,z;
		currSStrem  >> x;
		currSStrem  >> y;
		currSStrem  >> z;
		posx[counter] = x;
		posy[counter] = y;
		counter++;
		getline(PLYMeshFile,currline);
		istringstream currSStrem2(currline);
		currSStrem2 >> InpType;
	}
	//now the InpType is not Vertex, and I have the first line to read the first face from:
	while(InpType == "Face" && !PLYMeshFile.eof()){
		istringstream currSStrem(currline);
		currSStrem >> InpType;
		currSStrem >> tmp_int;
		int pnts[3];
		for (int i=0;i<3;++i){
			currSStrem  >> pnts[i];
		}
		Links0.push_back(pnts[0]);
		Links1.push_back(pnts[1]);
		Links0.push_back(pnts[1]);
		Links1.push_back(pnts[2]);
		Links0.push_back(pnts[2]);
		Links1.push_back(pnts[0]);
		getline(PLYMeshFile,currline);
		istringstream currSStrem2(currline);
		currSStrem2 >> InpType;
	}
	cout<<"finished reading file"<<endl;
}

void EllipseLayoutGenerator::minimisewitSpringMovement(int n_iter){
	double dt = 0.001;
	double ksp =10.0;
	const int n = posx.size();
	double F[n][2];
	for (int iter = 0; iter < n_iter; ++iter){
		for (int i=0; i<n; ++i){
			F[i][0]= 0.0; 
			F[i][1]= 0.0;
		}
		int nLink = Links0.size();
		for (int i=0; i<nLink;++i){
			double dx = posx[Links0[i]] - posx[Links1[i]];
			double dy = posy[Links0[i]] - posy[Links1[i]];
			double dsq = dx*dx + dy*dy;
			double d = 0; 		
			if (dsq > 1e-5){
				d= pow(dx*dx + dy*dy ,0.5);
				dx /=d;
				dy /=d;
			}
			else {
				dx = 0.5;
				dy = 0.5;			
			}
			double currForce = ksp * ( sideLen - d );
			double Fx = currForce*dx;
			double Fy = currForce*dy;
			F[Links0[i]][0] += Fx;
			F[Links0[i]][1] += Fy;
			F[Links1[i]][0] -= Fx;
			F[Links1[i]][1] -= Fy;
		}
		for (int i=0; i<n; ++i){
			cerr<<"Force on node : "<<i<<" "<< F[i][0]<<" "<< F[i][1]<<endl;	
		}
		for (int i=0; i<n; ++i){
			posx[i] += F[i][0]*dt;
			posy[i] += F[i][1]*dt;
		}
	}
	ofstream vectorsForGnuplot;
	vectorsForGnuplot.open("./Vectors2.out",ofstream::trunc);
	writeVectors2D(vectorsForGnuplot);	
}

void EllipseLayoutGenerator::calculateAverageSideLength(){
	int nTri = triangles.size();
	double SumSide =0.0;	
	for (int i =0; i<nTri; ++i){
		int IDcorner0, IDcorner1;
		int node0, node1;
		double dx, dy, side;		
		for (int j=0; j<3; ++j){		
			if (j==0){IDcorner0 = 0, IDcorner1 = 1;}
			else if (j==1){IDcorner0 = 1, IDcorner1 = 2;}
			else if (j==2){IDcorner0 = 2, IDcorner1 = 0;}
			node0 = triangles[i][IDcorner0];
			node1 = triangles[i][IDcorner1];
			dx = posx[node0] - posx[node1];
			dy = posy[node0] - posy[node1];
			side = pow(dx*dx + dy*dy , 0.5);
			SumSide +=side;
		}
	}
	SumSide /= (3.0 * nTri);
	cerr<<"Average side length of a triangle: "<<SumSide<<" -- desired length was: "<<sideLen/condensation<<endl;
}

/*void EllipseLayoutGenerator::writeMeshFileForSimulation(double zHeight, int zLayers){
	ofstream MeshFile;
	MeshFile.open("./MeshFile.out",ofstream::trunc);	
	MeshFile<<posx.size()*(zLayers+1)<<endl;
	cerr<<"posx.size(): "<<posx.size()<<" "<<posx.size()*(zLayers+1)<<endl;
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	0    0"<<endl;	//x-coord   y-coord   z-coord   basal-node identifier (0);
	}
	double dzHeight = zHeight/zLayers;
	double currzHeight = dzHeight;
	for (int layers=0; layers<zLayers; ++layers){
		int nodePositionIdentifier = 2; //0=basal. 1=apical, 2=midline
		if (layers == zLayers-1){nodePositionIdentifier=1;}
		for (int i=0;i<n;++i){
			MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<currzHeight<<"    "<<nodePositionIdentifier<<endl;	
		}
		currzHeight += dzHeight;
	}
	int nTri = triangles.size();
	MeshFile<<nTri*zLayers<<endl;
	int currOffset =0;
	for (int layers=0; layers<zLayers; ++layers){
		for (int i =0; i<nTri; ++i){
			double refpos[6][3];
			int nodes[6];
			for (int j=0;j<3;++j){			
				nodes[j]=triangles[i][j]+currOffset;
			}
			for (int j=3;j<6;++j){			
				nodes[j]=triangles[i][j-3]+currOffset+n;
			}	
			//writing Shapetype
			MeshFile<<"1	";
			//writing nodes of prism
			for (int j=0;j<6;++j){
				MeshFile<<nodes[j]<<"     ";
			}
			//writing positions for reference prism:
			//basal nodes:
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<dzHeight*layers<<"    ";		
			}
			//apical nodes:
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<dzHeight*(layers+1)<<"   ";		
			}
			MeshFile<<endl;
		}
		currOffset += n;
	}
}*/

void EllipseLayoutGenerator::writeMeshFileForSimulation(double zHeight, int zLayers){
	ofstream MeshFile;
	MeshFile.open("./MeshFile.out",ofstream::trunc);	
	MeshFile<<posx.size()*(zLayers+1)<<endl;
	cerr<<"posx.size(): "<<posx.size()<<" "<<posx.size()*(zLayers+1)<<endl;
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<posz[i]<<"    0	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	//x-coord   y-coord   z-coord   basal-node identifier(0) tissueType(columnar-0, peripodium-2) flag for border nodes
		
	}
	double dzHeight = zHeight/zLayers;
	double currzHeight = dzHeight;
	for (int layers=0; layers<zLayers; ++layers){
		int nodePositionIdentifier = 2; //0=basal. 1=apical, 2=midline
		if (layers == zLayers-1){nodePositionIdentifier=1;}
		for (int i=0;i<n;++i){
			MeshFile<<posx[i]+zVec0[i]*currzHeight<<"	"<<posy[i]+zVec1[i]*currzHeight<<"	"<<posz[i]+zVec2[i]*currzHeight<<"    "<<nodePositionIdentifier<<"	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	
		}
		currzHeight += dzHeight;
	}
	int nTri = triangles.size();
	MeshFile<<nTri*zLayers<<endl;
	int currOffset =0;
	for (int layers=0; layers<zLayers; ++layers){
		for (int i =0; i<nTri; ++i){
			double refpos[6][3];
			int nodes[6];
			for (int j=0;j<3;++j){			
				nodes[j]=triangles[i][j]+currOffset;
			}
			for (int j=3;j<6;++j){			
				nodes[j]=triangles[i][j-3]+currOffset+n;
			}	
			//writing Shapetype
			MeshFile<<"1	";
			//writing nodes of prism
			for (int j=0;j<6;++j){
				MeshFile<<nodes[j]<<"     ";
			}
			//writing positions for reference prism:
			//basal nodes:
			double currzHeight = dzHeight*layers;
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]+zVec0[triangles[i][j]]*currzHeight<<"   "<<posy[triangles[i][j]]+zVec1[triangles[i][j]]*currzHeight<<"   "<<posz[triangles[i][j]]+zVec2[triangles[i][j]]*currzHeight<<"    ";		
			}
			currzHeight = dzHeight*(layers+1);
			//apical nodes:
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]+zVec0[triangles[i][j]]*currzHeight<<"   "<<posy[triangles[i][j]]+zVec1[triangles[i][j]]*currzHeight<<"   "<<posz[triangles[i][j]]+zVec2[triangles[i][j]]*currzHeight<<"   ";		
			}
			MeshFile<<endl;
		}
		currOffset += n;
	}
}

void EllipseLayoutGenerator::writeTriangleMeshFileForSimulation(double zHeight){
	ofstream MeshFile;
	MeshFile.open("./TriangleFile.out",ofstream::trunc);	
	MeshFile<<posx.size()<<endl;
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	0	    2	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	//x-coord   y-coord   z-coord   mid-node identifier(0) tissueType(columnar-0,peripodium-2) flag for border nodes
		
	}
	int nTri = triangles.size();
	MeshFile<<nTri<<endl;
	for (int i =0; i<nTri; ++i){
		double refpos[3][3];
		int nodes[3];
		for (int j=0;j<3;++j){			
			nodes[j]=triangles[i][j];
		}
			
		//writing Shapetype as truiangle (4)
		MeshFile<<"4	";
		//writing nodes of triangle
		for (int j=0;j<3;++j){
			MeshFile<<nodes[j]<<"     ";
		}
		//writing the slab height for the triangle
		MeshFile<<zHeight<<"     ";
		//writing positions for reference triangle:
		for (int j=0;j<3;++j){
			MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<posz[triangles[i][j]]<<"    ";		
		}
		MeshFile<<endl;		
	}
}

void EllipseLayoutGenerator::addInitialPosZ(double curvedEndRadia, int nCurveLayers){
	int n = posx.size();
	for (int i=0; i<n; ++i){
		cout<<"node "<<i<<" "<<posx[i]<<" "<<posy[i]<<endl;	
	}
	double r12 = (r1Init + sideLen*0.5) * (r1Init + sideLen*0.5);
	double r22 = (r2Init + sideLen*0.5) * (r2Init + sideLen*0.5);
	double dZ =0;
	if (nCurveLayers>0){
		dZ = curvedEndRadia / ((double)  nCurveLayers);
	}
	for (int i=0;i<n;++i){
		//cout<<" assigning node: "<<i<<endl;		
		double distance = posx[i]*posx[i]/r12 + posy[i]*posy[i]/r22;
		if (distance<1){
			//cout<<"distance id: "<<distance<<"size of tissueType: "<<tissueType.size()<<endl;
			tissueType.push_back(0); //0: columnar layer; 2: peripodium
			posz.push_back(0.0);
			if (InvardNormalsX[i] == -100) {
				//this is a newly added point, and it does require setting the normals 
				//but it also is inside the columnar layer, where the normals will not be used
				//just as in the midline, I am setting the nomral to 0, and will not bother with it later on;
				InvardNormalsX[i]=0.0;
				InvardNormalsY[i]=0.0;
				
			}
		}	
	}
}

void EllipseLayoutGenerator::correctNormalsForAddedPoints(){
	int n = posx.size();
	double threshold = sideLen*1.25;	
	double t2 = threshold*threshold;
	for (int i=0; i<n; ++i){
		double summation = 0.0, sumX=0.0, sumY = 0.0;
		int counter = 0;
		if (InvardNormalsX[i] == -100) {
			for (int j=0; j<n; ++j){
				double dx = posx[i] - posx[j];
				double dy = posy[i] - posy[j];	
				double dz = posy[i] - posz[j];	
				double d2 = dx*dx + dy*dy + dz*dz;
				if (d2<t2 && InvardNormalsX[j] !=-100 && InvardNormalsY[j]!=-100 && !(InvardNormalsX[j]==0 && InvardNormalsY[j]==0)){
					//the node is in the viicinity, and it has a normal set, I will use it in interpolation
					d2 = pow(d2,0.5);
					summation += d2;
					sumX += d2*InvardNormalsX[j];
					sumY += d2*InvardNormalsY[j];
				}
			}
			InvardNormalsX[i] = sumX/summation;
			InvardNormalsY[i] = sumY/summation;
		}
	}
}

void EllipseLayoutGenerator::addZVecs(double curvedEndRadia, int nCurveLayers){
	//adding the zVectors of the columnar layer, 
	//the vectors are going from basal to apical layer (0,0,1);	
	int n = posx.size();
	double r12 = (r1Init + sideLen*0.5) * (r1Init + sideLen*0.5);
	double r22 = (r2Init + sideLen*0.5) * (r2Init + sideLen*0.5);
	double dZ =0;	
	if (nCurveLayers>0){
		dZ = curvedEndRadia / ((double)  nCurveLayers);
	}
	for (int i=0;i<n;++i){
		//cout<<"i: "<<i<<"n: "<<n<<" Xsize: "<<posx.size()<<" Ysize: "<<posy.size()<<" Zsize: "<<posz.size()<<" Xnorm: "<<InvardNormalsX.size()<<" Ynorm: "<<InvardNormalsY.size()<<endl;	
		double distance = posx[i]*posx[i]/r12 + posy[i]*posy[i]/r22;
		if (distance<1){		
			zVec0.push_back(0.0);
			zVec1.push_back(0.0);
			zVec2.push_back(1.0);
		}
		else{
			double rotAx [2] = {InvardNormalsY[i], -InvardNormalsX[i]}; //cross product of normal with (0,0,1);
			double sintet = (curvedEndRadia-posz[i])/curvedEndRadia;
			double costet = pow((2.0*curvedEndRadia-posz[i])*posz[i],0.5)/curvedEndRadia;
			double rotMat[3][3];
			rotMat[0][0] = costet + rotAx[0]*rotAx[0]*(1.0 - costet);
			rotMat[0][1] = rotAx[0]*rotAx[1]*(1.0 - costet);
			rotMat[0][2] = rotAx[1]*sintet;
			rotMat[1][0] = rotAx[0]*rotAx[1]*(1.0 - costet);
			rotMat[1][1] = costet + rotAx[1]*rotAx[1]*(1.0 - costet);
			rotMat[1][2] = (-1.0)*rotAx[0]*sintet;
			rotMat[2][0] = (-1.0)*rotAx[1]*sintet;
			rotMat[2][1] = rotAx[0]*sintet;
			rotMat[2][2] = costet;
			zVec0.push_back(rotMat[0][0]*InvardNormalsX[i] + rotMat[0][1]*InvardNormalsY[i]);
			zVec1.push_back(rotMat[1][0]*InvardNormalsX[i] + rotMat[1][1]*InvardNormalsY[i]);
			zVec2.push_back(rotMat[2][0]*InvardNormalsX[i] + rotMat[2][1]*InvardNormalsY[i]);
			//cout<<"node: "<<i<<" pos: "<<posx[i]<<" "<<posy[i]<<" "<<posz[i]<<" normal: "<<InvardNormalsX[i]<<" "<<InvardNormalsY[i]<<" zvec: "<<zVec0[i]<<" "<<zVec1[i]<<" "<<zVec2[i]<<endl;		
		}
	}

}

int EllipseLayoutGenerator::calculateLargeRadiaAndLayers(double curvedEndRadia){
	double QuarterCircumference = pi*curvedEndRadia/2.0;
	r1 = r1Init +  QuarterCircumference;
	r2 = r2Init +  QuarterCircumference;
	double tempNLayers = QuarterCircumference/sideLen;
	int NLayers = floor(tempNLayers);
	sideLen = QuarterCircumference/((double) NLayers);
	//cout<<"curvedEndRadia: "<<curvedEndRadia<<" r1: "<<r1<<" r2: "<<r2<<"  tempNLayers: "<<tempNLayers<< "NLayers: "<<NLayers<<" sideLen: "<<sideLen<<endl;
	return NLayers;
}

int main(int argc, char **argv)
{	
	double 	DVRadius = 20.0;
	double 	APRadius = 10.0;
	double 	ABHeight = 5.0;
	double 	sideLength = 5.0;
	int    	ABLayers =1;
	bool 	addCurvedEnd = false;
	double	curvedEndRadia = ABHeight*1.5; // starting from the basal side


	EllipseLayoutGenerator Lay01(DVRadius,APRadius,sideLength,0.9);
	bool calculateNextRing = true;
	int counter =0;
	while (calculateNextRing) {
		Lay01.calculateCircumference();
		Lay01.calculateCurrentBorderNumber();
		Lay01.calculatedtet();
		Lay01.addEquidistantRingMin();	
		//Lay01.addRectangle();
		calculateNextRing = Lay01.updateRadia();
	}
	//outside the ring calculation, now I will add points in the middle, in horizontal or vertical,
	//depending on which axis failed first, r1 or r2:
	Lay01.addMidLine();
	int nCurveLayers = 0;
	if(addCurvedEnd){
		nCurveLayers = Lay01.calculateLargeRadiaAndLayers(curvedEndRadia);
		int counter = 0;
		while (counter<nCurveLayers){
			int nLast = Lay01.posx.size();
			Lay01.calculateCircumference();
			Lay01.calculateCurrentBorderNumber();
			Lay01.calculatedtet();
			Lay01.addEquidistantRingMin();
			Lay01.updateRadia();
			counter++;
		}
	}
	//cout<<"before Tesselate2D, n =  "<<	Lay01.posx.size()<<endl;
	Lay01.Tesselate2D();
	Lay01.readInTesselation2D();
	Lay01.calculateAverageSideLength();
	Lay01.addInitialPosZ(curvedEndRadia, nCurveLayers);
	Lay01.correctNormalsForAddedPoints();
	Lay01.addZVecs(curvedEndRadia, nCurveLayers);
	Lay01.writeMeshFileForSimulation(ABHeight,ABLayers);	
	Lay01.writeTriangleMeshFileForSimulation(ABHeight);
	//output the points for plotting:
	int n=Lay01.posx.size();	
	cerr<<"r1: "<<Lay01.r1<<" r2: "<<Lay01.r2<<endl;
	Lay01.calculateAverageSideLength();
}

