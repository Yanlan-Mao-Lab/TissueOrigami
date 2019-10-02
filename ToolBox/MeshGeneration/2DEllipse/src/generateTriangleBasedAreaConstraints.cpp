using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>


int main(int argc, char **argv)
{	
	if (argc < 2){
		cerr<<"Please provide hte base input (e.g. Points.5)"<<endl;
		return 0;
	}
	if (argc < 3){
		cerr<<"Please provide the lower limit of triangles in x (relative pos, 0,1)"<<endl;
		return 0;
	}
	if (argc < 4){
		cerr<<"Please provide the higher limit of triangles in x (relative pos, 0,1)"<<endl;
		return 0;
	}	
	if (argc < 5){
		cerr<<"Please provide the higher limit of triangles in y (relative pos, 0.5,1)"<<endl;
		return 0;
	}
	if (argc < 6){
		cerr<<"Please provide the area limit for refined triangles"<<endl;
		return 0;
	}
	if (argc < 7){
		cerr<<"Please provide the buffer range outside refined triangles for double are limit (relative pos, 0,1)"<<endl;
		return 0;
	}
	
	const char* inpchar = argv[1];
	string inpstring = string(inpchar);
	string inpNodesNameString = inpstring +".node";
	string inpEleNameString = inpstring +".ele";
	string outAreaNameString = inpstring +".area";
	const char* inpchar2 = argv[2];
	double xLowLimit = atof(inpchar2);
	const char* inpchar3 = argv[3];
	double xHighLimit = atof(inpchar3);
	const char* inpchar4 = argv[4];
	double yHighLimit = atof(inpchar4);
	const char* inpchar5 = argv[5];
	double areaLimit = atof(inpchar5);
	const char* inpchar6 = argv[6];
	double buffer = atof(inpchar6);
	
	const char* inpNodesName = inpNodesNameString.c_str();
	const char* inpEleName = inpEleNameString.c_str();
	const char* outAreaName = outAreaNameString.c_str();
	cout<<"inpEleName: "<<inpEleName<<" inpNodesName: "<<inpNodesName<<" outAreaName: "<<outAreaName<<endl;
	ifstream inputMeshNodes;
	ifstream inputMeshEle;
	ofstream outputMeshArea;
	inputMeshNodes.open(inpNodesName,ifstream::in);
	inputMeshEle.open(inpEleName,ifstream::in);
	outputMeshArea.open(outAreaName,ifstream::out);
	if (!inputMeshNodes.good() || !inputMeshNodes.is_open()){
		cerr<<"input nodes cannot be opened: "<<inpNodesName<<endl;
		return 0;
	}
	if (!inputMeshEle.good() || !inputMeshEle.is_open()){
		cerr<<"input elements cannot be opened: "<<inpEleName<<endl;
		return 0;
	}
	if (!outputMeshArea.good() || !outputMeshArea.is_open()){
		cerr<<"output area cannot be opened: "<<outputMeshArea<<endl;
		return 0;
	}
	//Read x positions and generate bounding box in x:
	double minX = 100000.0, maxX = -100000.0, minY = 100000.0, maxY = -100000.0;
	int dummyInt;
	double dummyDouble;
	int nNodes;
	inputMeshNodes>>nNodes;
	cout<<"nNodes: "<<nNodes<<endl;
	const int n = nNodes;
	bool nodesWithinRange[n];
	bool nodesWithinBufferRange[n];
	double xValues[n],yValues[n];

	inputMeshNodes>>dummyInt;
	inputMeshNodes>>dummyInt;
	inputMeshNodes>>dummyInt;
	for (int i=0; i<nNodes; ++i){
		int id;
		double x,y,dummyDouble;
		inputMeshNodes>>id;
		inputMeshNodes>>x;
		inputMeshNodes>>y;
		inputMeshNodes>>dummyInt;
		//set bounding box in x
		if (x<minX){
			minX = x;
		}
		if (x>maxX){
			maxX = x;
		}
		if (y<minY){
			minY = y;
		}
		if (y>maxY){
			maxY = y;
		}
		//cout<<"i: "<<i<<" id: "<<id<<" x: "<<x<<" min-max X: "<<minX<<" "<<maxX<<endl;
		xValues[i] = x;
		yValues[i] = y;
	}
	double boundingBoxXSize = maxX-minX;
	double boundingBoxYSize = maxY-minY;
	for (int i=0; i<nNodes; ++i){
		cout<<" xValues["<<i<<"]: "<<xValues[i]<<" min-max X: "<<minX<<" "<<maxX<<" after refinement ";
		xValues[i] -=minX;
		xValues[i] /=boundingBoxXSize;
		yValues[i] -=minY;
		yValues[i] /=boundingBoxYSize;
		if (xValues[i]>=xLowLimit && xValues[i] <= xHighLimit
		    && yValues[i] <=yHighLimit){
			nodesWithinRange[i] = true;
			nodesWithinBufferRange[i] = false;
		}
		else {
			nodesWithinRange[i] = false;
			//check if they fall within buffer:
			if (xValues[i]>=xLowLimit-buffer && xValues[i] <= xLowLimit && yValues[i] <=yHighLimit){
				nodesWithinBufferRange[i] = true;
			}
			else if (xValues[i] >= xHighLimit && xValues[i] <= xHighLimit+buffer && yValues[i] <=yHighLimit){
				nodesWithinBufferRange[i] = true;
			}
			else{
				nodesWithinBufferRange[i] = false;	
			}
		}
		
		cout<<" x,y["<<i<<"]: "<<xValues[i]<<" "<<yValues[i]<<" WithinRange: "<<nodesWithinRange[i]<<endl;
	}
	//Now read in the triangles, and write area file:
	int nEle;
	const int constNEle = nEle;
	inputMeshEle>>nEle;
	inputMeshEle>>dummyInt;
	inputMeshEle>>dummyInt;
	outputMeshArea<<nEle<<endl;
	for (int i=0; i<nEle; ++i){
		//cout<<"i: "<<i<<endl;
		int id,node0,node1,node2;
		inputMeshEle>>id;	
		inputMeshEle>>node0;	
		inputMeshEle>>node1;	
		inputMeshEle>>node2;
		//check if within range:
		if (nodesWithinRange[node0] || nodesWithinRange[node1] || nodesWithinRange[node2]){
			outputMeshArea<<id<<" "<<areaLimit<<endl;
		}
		else if (nodesWithinBufferRange[node0] || nodesWithinBufferRange[node1] || nodesWithinBufferRange[node2]){
			outputMeshArea<<id<<" "<<areaLimit*2.0<<endl;
		}
		else{
			outputMeshArea<<id<<" -1"<<endl;
		}
	}
	//cout<<"outside loop"<<endl;
	inputMeshNodes.close();
	inputMeshEle.close();
	outputMeshArea.close();
}
