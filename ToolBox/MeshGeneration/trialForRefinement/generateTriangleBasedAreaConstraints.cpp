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
		cerr<<"Please provide the area limit for refined triangles"<<endl;
		return 0;
	}
	if (argc < 6){
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
	double areaLimit = atof(inpchar4);
	const char* inpchar5 = argv[5];
	double buffer = atof(inpchar5);
	
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
	double minX = 100000.0, maxX = -100000.0;
	int dummyInt;
	double dummyDouble;
	int nNodes;
	inputMeshNodes>>nNodes;
	cout<<"nNodes: "<<nNodes<<endl;
	const int n = nNodes;
	bool nodesWithinRange[n];
	bool nodesWithinBufferRange[n];
	double xValues[n];

	inputMeshNodes>>dummyInt;
	inputMeshNodes>>dummyInt;
	inputMeshNodes>>dummyInt;
	for (int i=0; i<nNodes; ++i){
		int id;
		double x,dummyDouble;
		inputMeshNodes>>id;
		inputMeshNodes>>x;
		inputMeshNodes>>dummyDouble;
		inputMeshNodes>>dummyInt;
		//set bounding box in x
		if (x<minX){
			minX = x;
		}
		if (x>maxX){
			maxX = x;
		}
		//cout<<"i: "<<i<<" id: "<<id<<" x: "<<x<<" min-max X: "<<minX<<" "<<maxX<<endl;
		xValues[i] = x;
	}
	double boundingBoxXSize = maxX-minX;
	for (int i=0; i<nNodes; ++i){
		cout<<" xValues["<<i<<"]: "<<xValues[i]<<" min-max X: "<<minX<<" "<<maxX<<" after refinement ";
		xValues[i] -=minX;
		xValues[i] /=boundingBoxXSize;
		if (xValues[i]>=xLowLimit && xValues[i] <= xHighLimit){
			nodesWithinRange[i] = true;
			nodesWithinBufferRange[i] = false;
		}
		else {
			nodesWithinRange[i] = false;
			//check if they fall within buffer:
			if (xValues[i]>=xLowLimit-buffer && xValues[i] <= xLowLimit){
				nodesWithinBufferRange[i] = true;
			}
			else if (xValues[i] >= xHighLimit && xValues[i] <= xHighLimit+buffer){
				nodesWithinBufferRange[i] = true;
			}
			else{
				nodesWithinBufferRange[i] = false;	
			}
		}
		
		cout<<" xValues["<<i<<"]: "<<xValues[i]<<" WithinRange: "<<nodesWithinRange[i]<<endl;
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
