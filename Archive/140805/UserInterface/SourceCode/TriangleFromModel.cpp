/*
 * TriangleFromModel.cpp
 *
 *  Created on: 20 Mar 2014
 *      Author: melda
 */

#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>

#include "TriangleFromModel.h"

using namespace std;

TriangleFromModel::TriangleFromModel(){
};

TriangleFromModel::~TriangleFromModel(){
};

void TriangleFromModel::setPoints(float **coord){
	for (int i=0;i<3;i++) {
		for (int j=0; j<3; j++){
			Coords[i][j]=coord[i][j];
		}
	}
};

void TriangleFromModel::setNo(int no){
	No = no;
	setIdentifyColour();
	setName();
}

void TriangleFromModel::setIdentifyColour(){
	//Triangles will use red:
	IdentifierColor[0] = 1+No;
	IdentifierColor[1] = 0;
	IdentifierColor[2] = 0;
	IdentifierColor[3] = 255;
	cerr<<"Color for triangle: "<<No<<" "<<IdentifierColor[0]<<" "<<IdentifierColor[1]<<" "<<IdentifierColor[2]<<" "<<IdentifierColor[3]<<endl;
}

void TriangleFromModel::setName(){
	stringstream convert;    // stream used for the conversion
	convert.fill('0');
	convert.width(4);
	convert <<No;      // insert the textual representation of 'Number' in the characters in the stream
	string NoString;          // string which will contain the result
	NoString = convert.str(); // set 'Result' to the contents of the stream
	name = "Triangle"+NoString;
	cerr<<"Name for triangle: "<<No<<" "<<name<<endl;
}
string TriangleFromModel::getName(){
	return name;
}


NodeFromModel::NodeFromModel(){
};

NodeFromModel::~NodeFromModel(){
};

void NodeFromModel::setPoints(float coord[3]){
	for (int i=0;i<3;i++) {
		Coords[i]=coord[i];
	}
};
void NodeFromModel::setNo(int no){
	No = no;
	setIdentifyColour();
	setName();
}

void NodeFromModel::setIdentifyColour(){
	//Nodes will use green:
	IdentifierColor[0] = 0;
	IdentifierColor[1] = 1+No;
	IdentifierColor[2] = 0;
	IdentifierColor[1] = 255;
}

void NodeFromModel::setName(){
	stringstream convert;    // stream used for the conversion
	convert.fill('0');
	convert.width(4);
	convert <<No;      // insert the textual representation of 'Number' in the characters in the stream
	string NoString;          // string which will contain the result
	NoString = convert.str(); // set 'Result' to the contents of the stream
	name = "Node"+NoString;
	cerr<<"Name for Node: "<<No<<" "<<name<<endl;
}
string NodeFromModel::getName(){
	return name;
}
