/*
 * TriangleFromModel.h
 *
 *  Created on: 20 Mar 2014
 *      Author: melda
 */

#include <string.h>
using namespace std;

class TriangleFromModel{
public:
	TriangleFromModel();
	~TriangleFromModel();
	float Coords[3][3];
	float DrawColour[4];
	void setPoints(float **coord);
	void setNo(int no);
	void setIdentifyColour();
	void setName();
	string getName();
	int No;
	int IdentifierColor[4];
	string name;
};

class NodeFromModel{
public:
	NodeFromModel();
	~NodeFromModel();
	float Coords[3];
	float DrawColour[4];
	void setPoints(float coord[3]);
	void setNo(int no);
	void setIdentifyColour();
	void setName();
	string getName();
	int No;
	int IdentifierColor[4];
	string name;
};

