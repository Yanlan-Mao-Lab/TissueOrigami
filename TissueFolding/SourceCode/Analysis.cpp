/*
 * Analysis.cpp
 *
 *  Created on: 25 Jul 2016
 *      Author: melda
 */



#include "Analysis.h"

Analysis::Analysis(int dim ,int nEle, string saveDirectoryToDisplayString, vector<Node*> &nodes){
	Dim = dim;
	nElements = nEle;
	nNodes = nodes.size();

	string saveFileString = saveDirectoryToDisplayString +"/Save_AnalysisVolumeMaps";
	const char* name_saveFileVolumeMaps = saveFileString.c_str();;
	saveFileVolumeMaps.open(name_saveFileVolumeMaps, ofstream::out);
	if (!(saveFileVolumeMaps.good() && saveFileVolumeMaps.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileVolumeMaps<<endl;
	}
	string saveFileLenghtMeasurementsString = saveDirectoryToDisplayString +"/Save_AnalysisLengthMeasurements";
	const char* name_saveFileLenghtMeasurements = saveFileLenghtMeasurementsString.c_str();;
	saveFileLenghtMeasurements.open(name_saveFileLenghtMeasurements, ofstream::out);
	if (!(saveFileLenghtMeasurements.good() && saveFileLenghtMeasurements.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileLenghtMeasurements<<endl;
	}
	string saveFileKinkPositionsString = saveDirectoryToDisplayString +"/Save_AnalysisKinkPositions";
	const char* name_saveFileKinkPositions = saveFileKinkPositionsString.c_str();;
	saveFileKinkPositions.open(name_saveFileKinkPositions, ofstream::out);
	if (!(saveFileKinkPositions.good() && saveFileKinkPositions.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileKinkPositions<<endl;
	}
	setUpContourLinesDV(nodes);
}




void Analysis::calculateBoundingBoxSizeAndAspectRatio(int timeInSec, double boundingBoxLength, double boundingBoxWidth){
	// format is [time in sec] [DV length in microns] [AP length in microns] [aspect ratio DV/AP] [...] then contour lengths will follow in coming functions.
	saveFileLenghtMeasurements<<timeInSec<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxLength<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxWidth<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxLength/boundingBoxWidth<<" ";

}

void Analysis::setUpContourLinesDV(vector<Node*> &nodes){

	double posThreshold = 10E-2;
	double negThreshold = (-1.0)*posThreshold;
	for (int i=0 ; i<nNodes; ++i){
		if (nodes[i]->tissuePlacement == 0 ){ //tissue placement 0 is for basal nodes
			if (nodes[i]->Position[1] > negThreshold && nodes[i]->Position[1] < posThreshold){
				basalContourLineDVNodeIds.push_back(i);
			}
		}
		if (nodes[i]->tissuePlacement == 1 ){ //tissue placement 1 is for apical nodes
			if (nodes[i]->Position[1] > negThreshold && nodes[i]->Position[1] < posThreshold){
				apicalContourLineDVNodeIds.push_back(i);
			}
		}
	}
	sortPositionMinToMax(nodes, 0, basalContourLineDVNodeIds); //sort basalContourLineDVNodeIds with x position going from minimum to maximum
	sortPositionMinToMax(nodes, 0, apicalContourLineDVNodeIds); //sort apicalContourLineDVNodeIds with x position going from minimum to maximum
}


void Analysis::calculateContourLineLengthsDV(vector<Node*> &nodes){

	double currContourLength = 0;
	//calculating contour length for basal DV axis:
	int n = basalContourLineDVNodeIds.size();
	int id0 = basalContourLineDVNodeIds[0];
	for (int i=1;i<n;++i){
		int id1 = basalContourLineDVNodeIds[i];
		double dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
		double dy = nodes[id0]->Position[1] - nodes[id1]->Position[1];
		double dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
		double d = dx*dx + dy*dy + dz*dz;
		d = pow(d,0.5);
		currContourLength += d;
		id0 = id1;
	}
	//writing the basal contour length:
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<currContourLength<<" ";
	//calculating apical contour length:
	currContourLength = 0;
	n = apicalContourLineDVNodeIds.size();
	id0 = apicalContourLineDVNodeIds[0];
	for (int i=1;i<n;++i){
		int id1 = apicalContourLineDVNodeIds[i];
		double dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
		double dy = nodes[id0]->Position[1] - nodes[id1]->Position[1];
		double dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
		double d = dx*dx + dy*dy + dz*dz;
		d = pow(d,0.5);
		currContourLength += d;
		id0 = id1;
	}
	//writing the apical contour length:
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<currContourLength<<endl;

}

void Analysis::findApicalKinkPointsDV(int timeInSec,  double boundingBoxXMin, double boundingBoxLength, vector<Node*> &nodes){
	bool withRadiusOfCurvature = true;
	if (withRadiusOfCurvature){
		int n = apicalContourLineDVNodeIds.size();
		if (n<3){
			return;
		}
		double thresholdRCurvature = 2.0;
		int id0 = apicalContourLineDVNodeIds[0];
		int id1 = apicalContourLineDVNodeIds[1];
		for (int i=0;i<n-2;++i){
			int id2 = apicalContourLineDVNodeIds[i+2];
			//calculating the first slope between node id 0 and node id1:
			double dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
			double dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
			double m01 = dz/dx;
			//the slope for node id1 and node id2
			dx = nodes[id1]->Position[0] - nodes[id2]->Position[0];
			dz = nodes[id1]->Position[2] - nodes[id2]->Position[2];
			double m12 = dz/dx;
			double meanSlope = ( m01 + m12 ) * 0.5;
			double midPointX01 = (nodes[id0]->Position[0] + nodes[id1]->Position[0])*0.5;
			double midPointX12 = (nodes[id1]->Position[0] + nodes[id2]->Position[0])*0.5;
			double secondDerivative = (m01 - m12)/(midPointX01 - midPointX12);
			double radiusofCurvature = pow((1 + meanSlope*meanSlope ),1.5)/secondDerivative;
			double relativeX = (nodes[id1]->Position[0] - boundingBoxXMin)/boundingBoxLength;
			//cout<<timeInSec<<" "<<id1<<" "<<relativeX<<" "<<radiusofCurvature<<" "<<endl;
			id0 = id1;
			id1 = id2;
			if ((radiusofCurvature>0 && radiusofCurvature<thresholdRCurvature) || (radiusofCurvature<0 && radiusofCurvature>(-1.0)*thresholdRCurvature)){
				saveFileKinkPositions<<timeInSec<<" ";
				saveFileKinkPositions<<id1<<" ";
				saveFileKinkPositions.width(6);
				saveFileKinkPositions.precision(3);
				if (radiusofCurvature > 0){
					saveFileKinkPositions<<relativeX<<" -1.0"<<endl; // a positive normal means the vectors are forming a kink on the apical surface, writing kink first, peak second]
				}
				else{
					saveFileKinkPositions<<"-1.0 "<<relativeX<<endl;// a negative normal means the vectors are forming a peak on the apical surface, writing kink first, peak second]
				}
			}
		}
	}
	else{
		int n = apicalContourLineDVNodeIds.size();
		if (n>5){
			double thresholdRadian = 120.0 * (M_PI/180.0);
			double thresholdCosine =  cos(thresholdRadian);
			//I am taking 5 node regions, and calculating the angle that the vectors -2 , +2 neigs, and -1, +1 neigs make.
			// If the angle is sharp enough (cos above threshold), then I will look at the direction.
			// then I will check the direction of the curvature.
			//distance in microns, the threshold I would accept to qualify as a downward/upeard pattern between two nodes
			//as I am looking into 5 node clusters, this will require a kink to be 2*threshold deep before it is recorded.
			int id0 = apicalContourLineDVNodeIds[0];
			int id1 = apicalContourLineDVNodeIds[1];
			int id2 = apicalContourLineDVNodeIds[2];
			int id3 = apicalContourLineDVNodeIds[3];
			for (int i=1;i<n-1;++i){
				int id4 = apicalContourLineDVNodeIds[i+2];
				double vec0[3],vec1[3],vec3[3],vec4[3];
				for ( int a=0; a<3; ++a){
					 vec0[a] = nodes[id0]->Position[a] - nodes[id2]->Position[a];
					 vec1[a] = nodes[id1]->Position[a] - nodes[id2]->Position[a];
					 vec3[a] = nodes[id3]->Position[a] - nodes[id2]->Position[a];
					 vec4[a] = nodes[id4]->Position[a] - nodes[id2]->Position[a];
				}
				double dotP = vec0[0]*vec4[0]+vec0[1]*vec4[1]+vec0[2]*vec4[2];
				double mag0 = pow(vec0[0]*vec0[0]+vec0[1]*vec0[1]+vec0[2]*vec0[2],0.5);
				double mag4 = pow(vec4[0]*vec4[0]+vec4[1]*vec4[1]+vec4[2]*vec4[2],0.5);
				double cosTet = dotP / mag0 / mag4;
				if (cosTet > thresholdCosine){
					//the angle is small enough, check direction:
					//I am checking which direction the Y of hte normal vector points to, a kink should point ot
					// positive y.
					double crossPY = vec0[2]*vec4[0]-vec0[0]*vec4[2];
					int signCrossPY04 = 1;
					if (crossPY < 0) {
						signCrossPY04 = -1;
					}

					//Now I can move on to the next connection:
					dotP = vec1[0]*vec3[0]+vec1[1]*vec3[1]+vec1[2]*vec3[2];
					double mag1 = pow(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2],0.5);
					double mag3 = pow(vec3[0]*vec3[0]+vec3[1]*vec3[1]+vec3[2]*vec3[2],0.5);
					double cosTet = dotP / mag1 / mag3;
					if (cosTet > thresholdCosine){
						double crossPY = vec1[2]*vec3[0]-vec1[0]*vec3[2];
						int signCrossPY13 = 1;
						if (crossPY < 0) {
							signCrossPY13 = -1;
						}

						//now I have the sine and cosines of both comparisons, and they are both eligable as the angles are below the threshold.
						//I need to make sure they are facing the same direction, if yes, I can record the information:
						if (signCrossPY04*signCrossPY13 > 0){
							double relativeX = (nodes[id3]->Position[0] - boundingBoxXMin)/boundingBoxLength;
							saveFileKinkPositions<<timeInSec<<" ";
							saveFileKinkPositions<<id3<<" ";
							saveFileKinkPositions.width(6);
							saveFileKinkPositions.precision(3);
							if (signCrossPY04 > 0){
								saveFileKinkPositions<<relativeX<<" -1.0"<<endl; // a positive normal means the vectors are forming a kink on the apical surface, writing kink first, peak second]
							}
							else{
								saveFileKinkPositions<<"-1.0 "<<relativeX<<endl;// a negative normal means the vectors are forming a peak on the apical surface, writing kink first, peak second]
							}
						}
					}
				}
				id0 = id1;
				id1 = id2;
				id2 = id3;
				id3 = id4;
			}
		}
	}
	cout<<" finished findApicalKinkPointsDV "<<endl;
}

void Analysis::sortPositionMinToMax(vector<Node*> &nodes, int axisToSortWith, vector <int> &linkToArrayToSort ){
	int n = linkToArrayToSort.size();
	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			double pos1 = nodes[linkToArrayToSort[i]]->Position[axisToSortWith];
			double pos0 = nodes[linkToArrayToSort[i-1]]->Position[axisToSortWith];
			if(pos1<pos0){
				int temp=linkToArrayToSort[i-1];
				linkToArrayToSort[i-1]=linkToArrayToSort[i];
				linkToArrayToSort[i]=temp;
				swapped = true;
			}
		}
	}
}


void Analysis::calculateTissueVolumeMap	(vector<ShapeBase*> &elements, int timeInSec, double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth){
	//Format:
	//[time] [element id] [relative pos in X] [relative pos in Y] [reference shape volume] [current ideal volume] [current emergent volume]
	double totTissueIdealVolume = 0;
	double totTissueEmergentVolume = 0;
	vector<ShapeBase*>::iterator itEle;
	for (itEle=elements.begin(); itEle<elements.end(); ++itEle){
		//double* c = new double[3];
		//c = (*itEle)->getCentre();
    	(*itEle)->calculateRelativePosInBoundingBox(boundingBoxXMin, boundingBoxYMin,boundingBoxLength, boundingBoxWidth);
    	double* ReletivePos = new double[2];
		(*itEle)->getRelativePosInBoundingBox(ReletivePos);
    	double currEmergentVolume = (*itEle) -> calculateCurrentGrownAndEmergentVolumes();
		double currIdealVolume = (*itEle)->GrownVolume;
		totTissueIdealVolume  += currIdealVolume;
		totTissueEmergentVolume  += currEmergentVolume;
		saveFileVolumeMaps<<timeInSec<<" "<<(*itEle) ->Id<<" "<<ReletivePos[0]<<" "<<ReletivePos[1]<<" "<<(*itEle)->ReferenceShape->Volume<<" "<<currIdealVolume<<" "<<currEmergentVolume<<endl;
		delete[] ReletivePos;
	}
	cout<<" Time: "<<timeInSec<<" total tissue ideal volume: "<<totTissueIdealVolume<<" total tissue emergent volume:" <<totTissueEmergentVolume<<endl;
}
