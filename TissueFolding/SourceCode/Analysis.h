/*
 * Analysis.h
 *
 *  Created on: 25 Jul 2016
 *      Author: melda
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <vector>
#include "Node.h"
#include "ShapeBase.h"
using namespace std;

class Analysis{

protected:
   //double getApicalSideLengthAverage();
	int Dim ;		//dimensions of the system, should be 3D.
	int nNodes;		//number of nodes of the system to be analysed
	int nElements; //number of elements of the system to be analysed

	ofstream saveFileVolumeMaps;
	ofstream saveFileLenghtMeasurements;
	ofstream saveFileKinkPositions;
	//vector <int> ApicalContourLineAPNodeIds;
	vector <int> apicalContourLineDVNodeIds;
	vector <int> basalContourLineDVNodeIds;
	/*vector <int> BasalContourLineAPNodeIds;
	vector <int> MidlineContourLineAPNodeIds;
	vector <int> MidlineContourLineDVNodeIds;
*/

public:
	Analysis(int dim ,int nEle, string saveDirectoryToDisplayString,vector<Node*> &nodes);
	~Analysis();
	void calculateBoundingBoxSizeAndAspectRatio(int timeInSec,double boundingBoxLength, double boundingBoxWidth);
	void sortPositionMinToMax(vector<Node*> &nodes, int axisToSortWith, vector <int> &linkToArrayToSort );

	void setUpContourLinesDV(vector<Node*> &nodes);
	void calculateContourLineLengthsDV(vector<Node*> &nodes);
	void findApicalKinkPointsDV(int timeInSec, double boundingBoxXMin,  double boundingBoxLength, vector<Node*> &nodes);
	void setUpContourLinesAP(vector<Node*> &nodes);
	void calculateContourLineLengthsAP(vector<Node*> &nodes);

	void calculateTissueVolumeMap(vector<ShapeBase*> &elements, int timeInSec, double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth);
	};

#endif /* ANALYSIS_H_ */
