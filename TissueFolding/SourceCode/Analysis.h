/*
 * Analysis.h
 *
 *  Created on: 25 Jul 2016
 *      Author: melda
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <vector>
#include <memory>
#include "Node.h"
#include "ShapeBase.h"

using namespace std;

class Analysis{

protected:
   //double getApicalSideLengthAverage();
	size_t Dim ;		//dimensions of the system, should be 3D.
	size_t nNodes;		//number of nodes of the system to be analysed
	size_t nElements; //number of elements of the system to be analysed

	ofstream saveFileVolumeMaps;
	ofstream saveFileLenghtMeasurements;
	ofstream saveFileKinkPositions;
	ofstream saveFileFoldMarkers;
	ofstream saveFileCircumference;
	string saveDirectoryString;
	vector <int> apicalContourLineDVNodeIds;
	vector <int> basalContourLineDVNodeIds;
	bool** connectivityMap;
	void setUpConnectivityMap(const std::vector<std::unique_ptr<Node>>& nodes);
public:
	double yPosForSideDVLine;
	double relativeYPosSideDVLine;
	vector<double> apicalContourLineDVSelectedYPositionsX;
	vector<double> apicalContourLineDVSelectedYPositionsZ;
	Analysis(int dim, string saveDirectoryToDisplayString, const std::vector<std::unique_ptr<Node>>& nodes, double boundingBoxWidth);
	~Analysis();

	void calculateBoundingBoxSizeAndAspectRatio(int timeInSec,double boundingBoxLength, double boundingBoxWidth);
	void sortPositionMinToMax(const std::vector<std::unique_ptr<Node>>& nodes, int axisToSortWith, vector <int> &linkToArrayToSort );
	void sortPointsMinToMaxBasedFirstArray(vector<double> &x, vector<double> &z, vector<int> &baseNodeId0, vector <int> &baseNodeId1);
	void sortPointsMinToMaxBasedOnInitialPos(vector<double> &base, vector<double> &x, vector<double> &z);
	void setUpContourLinesDV(const std::vector<std::unique_ptr<Node>>& nodes, double boundingBoxWidth);
	void calculateContourLineLengthsDV(const std::vector<std::unique_ptr<Node>>& nodes);
	void setUpSideApicalDVContour(const std::vector<std::unique_ptr<Node>>& nodes, double boundingBoxWidth);
	void updateSideContourPosition(double boundingBoxSizeY);
	void findApicalKinkPointsDV(int timeInSec, double boundingBoxXMin,  double boundingBoxLength, double boundingBoxWidth, const std::vector<std::unique_ptr<Node>>& nodes);
	void setUpContourLinesAP(const std::vector<std::unique_ptr<Node>>& nodes);
	void calculateContourLineLengthsAP(const std::vector<std::unique_ptr<Node>>& nodes);
	void calculateTissueVolumeMap(vector<ShapeBase*> &elements, int timeInSec, double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth);
	void saveNodesOnFold(int timeInSec, const std::vector<std::unique_ptr<Node>>& Nodes);
	void saveApicalCircumferencePosition(int timeInSec, const std::vector<std::unique_ptr<Node>>& nodes);

};

#endif /* ANALYSIS_H_ */
