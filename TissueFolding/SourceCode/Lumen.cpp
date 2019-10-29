/*
 * Lumen.cpp
 *
 *  Created on: 22 Jun 2018
 *      Author: melda
 */

#include "Lumen.h"
#include <algorithm>

using namespace std;
Lumen::Lumen(const std::vector <std::unique_ptr<ShapeBase>>& Elements, const std::vector <std::unique_ptr<Node>>& Nodes, double lumenBulkModulus, double lumenGrowthFold){
	Dim = 3;
	bulkModulus = lumenBulkModulus;
	initialIdealVolume = 0;
	currentIdealVolume = 0;
	currentVolume = 0;
	growthRate = 0;

	//lumenGrowthFold is the growth in 24 hours, the growth rate is in seconds:
    if (lumenGrowthFold>0){
        //calculate the rate only if the rate is non-zero
        growthRate = log(lumenGrowthFold)/3600/24;
    }
	rV = 0;
	nTriangleSize= 0;
    for (const auto & iterEle : Elements){
		if (iterEle->isECMMimimcingAtCircumference == false){
			if (iterEle->tissuePlacement != 3){ //not a lateral element

				if(iterEle->tissuePlacement == 1 ||(iterEle->tissuePlacement == 2 && iterEle->spansWholeTissue)){
					encapsulatingElementIds.push_back(iterEle->Id);
					nTriangleSize++;
					vector <int> apicalNodeIds;
					iterEle->getApicalNodeIds(apicalNodeIds);
					for (int nodeIdIterator =0;nodeIdIterator<3; ++nodeIdIterator){
						int currId = apicalNodeIds[nodeIdIterator];
						if (find(nodeIdsList.begin(), nodeIdsList.end(),currId)==nodeIdsList.end()){
							nodeIdsList.push_back(currId);
							Nodes[currId]->facingLumen = true;
						}
					}

				}

			}
		}
	}
	//set up the vector and matrix arrays:
	for (int i=0; i<nTriangleSize; ++i){
		//gsl_matrix* tmpMat1 = gsl_matrix_calloc(Dim, Dim);
		//gsl_matrix* tmpMat2 = gsl_matrix_calloc(Dim, Dim);
		//gsl_matrix* tmpMat3 = gsl_matrix_calloc(Dim, Dim);
		xCap1.push_back(gsl_matrix_calloc(Dim, Dim));
		xCap2.push_back(gsl_matrix_calloc(Dim, Dim));
		xCap3.push_back(gsl_matrix_calloc(Dim, Dim));
		x1.push_back(gsl_matrix_calloc(Dim, 1));
		x2.push_back(gsl_matrix_calloc(Dim, 1));
		x3.push_back(gsl_matrix_calloc(Dim, 1));
		g1.push_back(gsl_matrix_calloc(Dim, 1));
		g2.push_back(gsl_matrix_calloc(Dim, 1));
		g3.push_back(gsl_matrix_calloc(Dim, 1));
	}
	int nNode = nodeIdsList.size();
	Kv = gsl_matrix_calloc(Dim*nNode,Dim*nNode);
	KvNumerical = gsl_matrix_calloc(Dim*nNode,Dim*nNode);
	updateMatrices(Nodes, Elements);
	calculateCurrentVolume();
	currentIdealVolume = currentVolume*1.0;
	initialIdealVolume = currentIdealVolume;
	cout<<" in lumen constructor, currentIdealVolume "<<currentIdealVolume<<endl;

}


Lumen::~Lumen(){
	for (int i=0; i<nTriangleSize; ++i){
		gsl_matrix_free (xCap1[i]);
		gsl_matrix_free (xCap2[i]);
		gsl_matrix_free (xCap3[i]);
		gsl_matrix_free (x1[i]);
		gsl_matrix_free (x2[i]);
		gsl_matrix_free (x3[i]);
		gsl_matrix_free (g1[i]);
		gsl_matrix_free (g2[i]);
		gsl_matrix_free (g3[i]);
	}
	gsl_matrix_free (Kv);
	gsl_matrix_free (KvNumerical);
}

void	Lumen::updateMatrices(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	//cout<<"in update matrices"<<endl;
	for (int eleIndex=0; eleIndex<nTriangleSize; ++eleIndex){
		vector <int> apicalNodeIds;
		int elementId = encapsulatingElementIds[eleIndex];
		Elements[elementId]->getApicalNodeIds(apicalNodeIds);
		for (int k=0; k<3; ++k){//going over trienagle corners
			gsl_matrix* tmpPos; //pointer to the vector for position of corner k
			gsl_matrix* tmpPosCap; //pointer to the matrix of position of corner k
			if 		(k==0)	{tmpPos =x1[eleIndex]; tmpPosCap = xCap1[eleIndex];}
			else if (k==1)	{tmpPos =x2[eleIndex]; tmpPosCap = xCap2[eleIndex];}
			else if (k==2)	{tmpPos =x3[eleIndex]; tmpPosCap = xCap3[eleIndex];}

			gsl_matrix_set(tmpPos,0,0,Nodes[apicalNodeIds[k]]->Position[0]);
			gsl_matrix_set(tmpPos,1,0,Nodes[apicalNodeIds[k]]->Position[1]);
			gsl_matrix_set(tmpPos,2,0,Nodes[apicalNodeIds[k]]->Position[2]);

			gsl_matrix_set(tmpPosCap,0,1, (-1.0)*Nodes[apicalNodeIds[k]]->Position[2]);
			gsl_matrix_set(tmpPosCap,0,2, Nodes[apicalNodeIds[k]]->Position[1]);
			gsl_matrix_set(tmpPosCap,1,0, Nodes[apicalNodeIds[k]]->Position[2]);
			gsl_matrix_set(tmpPosCap,1,2, (-1.0)*Nodes[apicalNodeIds[k]]->Position[0]);
			gsl_matrix_set(tmpPosCap,2,0, (-1.0)*Nodes[apicalNodeIds[k]]->Position[1]);
			gsl_matrix_set(tmpPosCap,2,1, Nodes[apicalNodeIds[k]]->Position[0]);
		}
	}
}


void	Lumen::calculateCurrentVolume(){
	//cout<<"in calculate matrices"<<endl;
	double currentElementalVolume[(const int) nTriangleSize];
	#ifndef DO_NOT_USE_OMP
	/** If DO_NOT_USE_OMP is not defined,I will be using omp. This
	 * is necessary as omp is not set up on mac
	 */
	#pragma omp parallel for
	#endif
	for (int eleIndex=0; eleIndex<nTriangleSize; ++eleIndex){
		gsl_matrix* tmp = gsl_matrix_calloc(Dim, 1);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, xCap3[eleIndex], x1[eleIndex],0.0, tmp);
		//displayMatrix(xCap3[eleIndex],"xCap3");
		//displayMatrix(x1[eleIndex],"x1");
		//displayMatrix(x2[eleIndex],"x2");
		double dotP = 0;
		for (int k=0;k<3;++k){
			dotP += gsl_matrix_get(tmp,k,0)*gsl_matrix_get(x2[eleIndex],k,0);
		}
		currentElementalVolume[eleIndex] = 1.0/6.0*dotP;
	}
	currentVolume=0;
	for (int eleIndex=0; eleIndex<nTriangleSize; ++eleIndex){
		currentVolume +=currentElementalVolume[eleIndex];
	}
	//currentIdealVolume = currentVolume;
	rV =  (currentVolume - currentIdealVolume)/currentIdealVolume;
	//cout<<"current Volume of Lumen: "<<currentVolume<<" ideal volume: "<<currentIdealVolume <<" rV "<<rV<<endl;
}

void	Lumen::calculateResiduals(const std::vector <std::unique_ptr<Node>>& Nodes, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	//cout<<"in calculate residuals for Lumen"<<endl;
	double rVover6V0 =  rV /6.0 / currentIdealVolume;
	#ifndef DO_NOT_USE_OMP
	/** If DO_NOT_USE_OMP is not defined,I will be using omp. This
	 * is necessary as omp is not set up on mac
	 */
	#pragma omp parallel for
	#endif
	for (int eleIndex=0; eleIndex<nTriangleSize; ++eleIndex){
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, xCap2[eleIndex], x3[eleIndex],0.0, g1[eleIndex]);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, xCap3[eleIndex], x1[eleIndex],0.0, g2[eleIndex]);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, xCap1[eleIndex], x2[eleIndex],0.0, g3[eleIndex]);
		//now write to elastic g:
		vector <int> apicalNodeIndices;
		int eleId = encapsulatingElementIds[eleIndex];
		Elements[eleId]->getApicalNodeIndicesOnElement(apicalNodeIndices);
		for (int i=0; i<3; ++i){//triangle corners for apical nodes
			int indexOnElement = apicalNodeIndices[i];
			gsl_matrix* currentg;
			if (i == 0){currentg = g1[eleIndex];}
			else if (i == 1){currentg = g2[eleIndex];}
			else if (i == 2){currentg = g3[eleIndex];}
			for (int j=0; j<Dim; ++j){
				if (!Nodes[Elements[eleId]->NodeIds[indexOnElement]]->FixedPos[j]){
					//if (encapsulatingElements[eleIndex]->NodeIds[indexOnElement] == 176){
					//	cout<<" from element: "<<encapsulatingElements[eleIndex]->Id<<" node "<<encapsulatingElements[eleIndex]->NodeIds[indexOnElement]<<" g"<<i+1<<"["<<j<<"]: "<<gsl_matrix_get(currentg,j,0)<<endl;
					//}
					double valueToAdd= -1.0*bulkModulus*rVover6V0*gsl_matrix_get(currentg,j,0);
					Elements[eleId]->addToElementalElasticSystemForces(indexOnElement,j,valueToAdd );
					if (std::isnan(valueToAdd)){
						cout<<" element: "<<eleId<<" g dimention: "<<indexOnElement<<" "<<j<<" is NaN after addition in lumen forces: "<<valueToAdd<<endl;
					}
				}
				else{
					//for fiexd pos clear the corresponding g as well?
					gsl_matrix_set(currentg,j,0,0);
				}
			}
		}
	}
}


void 	Lumen::calculateLumengFromElementalResiduals(gsl_matrix* g, const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	gsl_matrix* tmpg;
	for (int eleIndex= 0; eleIndex<nTriangleSize; ++eleIndex){
		vector <int> apicalNodeIds;
		int eleId = encapsulatingElementIds[eleIndex];
		Elements[eleId]->getApicalNodeIds(apicalNodeIds);
		for (int nodeIterator = 0; nodeIterator<3; nodeIterator++){
			if 		(nodeIterator==0)	{tmpg =g1[eleIndex];}
			else if (nodeIterator==1)	{tmpg =g2[eleIndex];}
			else if (nodeIterator==2)	{tmpg =g3[eleIndex];}
			int currNodeId = apicalNodeIds[nodeIterator];
			std::vector<int>::iterator iteratorToCurrNode = find(nodeIdsList.begin(), nodeIdsList.end(),currNodeId);
			int indexOfCurrNodeOnnodeIdsList = std::distance(nodeIdsList.begin(),iteratorToCurrNode);
			for (int k=0; k<Dim; ++k){
				//cout<<" size of node list: "<<nodeIdsList.size()<<" size of g: "<<Dim*nNode<<" indexOfCurrNodeOnnodeIdsList: "<<indexOfCurrNodeOnnodeIdsList<<endl;
				double valueToAdd = gsl_matrix_get(tmpg,k,0);
				double CurrValue = gsl_matrix_get(g,indexOfCurrNodeOnnodeIdsList*Dim+k,0);
				gsl_matrix_set(g,indexOfCurrNodeOnnodeIdsList*Dim+k,0,CurrValue+valueToAdd);
			}
		}
	}
	//cout<<"finished calculateLumengFromElementalResiduals "<<endl;
}

void	Lumen::calculateJacobian(const std::vector <std::unique_ptr<ShapeBase>>& Elements){
	gsl_matrix_set_zero(Kv);
	double KbulkrVover6V0 =  bulkModulus*rV /6.0 / currentIdealVolume;
	double KbulkoneOver36V0Sq = bulkModulus/currentIdealVolume/currentIdealVolume/36;
	//assemble global g:
	//get number of Nodes, make node list
	int  nNode = nodeIdsList.size();
	gsl_matrix* g = gsl_matrix_calloc(Dim*nNode,1);
	calculateLumengFromElementalResiduals(g, Elements);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,KbulkoneOver36V0Sq, g, g,0.0, Kv);
	for (int eleIndex=0; eleIndex<nTriangleSize; ++eleIndex){
		gsl_matrix* scaledxij = gsl_matrix_calloc(3,3);
		vector <int> apicalNodeIndices;
		int eleId= encapsulatingElementIds[eleIndex];
		Elements[eleId]->getApicalNodeIndicesOnElement(apicalNodeIndices);
		for (int i=0; i<Dim; ++i){
			for (int j=0; j<Dim; ++j){
				if (i==j){
					gsl_matrix_set_zero(scaledxij);
					continue;
				}
				else if (i==0 && j==1){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap3[eleIndex]);
					gsl_matrix_scale(scaledxij,(-1.0)*KbulkrVover6V0);
				}
				else if (i==0 && j==2){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap2[eleIndex]);
					gsl_matrix_scale(scaledxij,KbulkrVover6V0);
				}
				else if (i==1 && j==0){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap3[eleIndex]);
					gsl_matrix_scale(scaledxij,KbulkrVover6V0);
				}
				else if (i==1 && j==2){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap1[eleIndex]);
					gsl_matrix_scale(scaledxij,(-1.0)*KbulkrVover6V0);
				}
				else if (i==2 && j==0){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap2[eleIndex]);
					gsl_matrix_scale(scaledxij,(-1.0)*KbulkrVover6V0);
				}
				else if (i==2 && j==1){
					Elements[eleId]->createMatrixCopy(scaledxij,  xCap1[eleIndex]);
					gsl_matrix_scale(scaledxij,KbulkrVover6V0);
				}
				//Now I have the tile ij that is to be written on the elemental K.
				//The corresponding tile on elemental K will depend on the element node index, depending on where the apical surface lies
				//I will get the index of nodes 1-3 of the triangle, they are either 0-2 or 3-5 for prisims
				//int indexOfiOnElement = apicalNodeIndices[i];
				//int indexOfjOnElement = apicalNodeIndices[j];

				vector <int> apicalNodeIds;
				Elements[eleId]->getApicalNodeIds(apicalNodeIds);
				int currNodeIdi = apicalNodeIds[i];
				std::vector<int>::iterator iteratorToCurrNodei = find(nodeIdsList.begin(), nodeIdsList.end(),currNodeIdi);
				int indexOfCurrNodeiOnnodeIdsList = std::distance(nodeIdsList.begin(),iteratorToCurrNodei);
				int currNodeIdj = apicalNodeIds[j];
				std::vector<int>::iterator iteratorToCurrNodej = find(nodeIdsList.begin(), nodeIdsList.end(),currNodeIdj);
				int indexOfCurrNodejOnnodeIdsList = std::distance(nodeIdsList.begin(),iteratorToCurrNodej);
				for (int k=0; k <Dim ; ++k){
					for (int m=0; m <Dim ; ++m){
						//cout<<"k & m "<<k<<" "<<m<<endl;
						double valueToAdd = gsl_matrix_get(scaledxij,k, m);
						//Get the node position on nodelist, then obtain the position of its elements on Kv acordingly.
						//Kv size is not based on number of elements
						//it is based on total number of nodes
						int indexiOnKv = indexOfCurrNodeiOnnodeIdsList*3 +m;
						int indexjOnKv = indexOfCurrNodejOnnodeIdsList*3 +k;
						double valueFromKv = gsl_matrix_get(Kv,indexiOnKv,indexjOnKv);
						double newValue = valueFromKv - valueToAdd;
						//if (indexiOnKv == 2 && indexjOnKv == 4){
						//	cout<<"["<<indexiOnKv<<", "<<indexjOnKv<<"]: value to add"<<valueToAdd<<" valueFromKv "<<valueFromKv<<endl;
						//}
						gsl_matrix_set(Kv,indexiOnKv,indexjOnKv,newValue);
					}
				}
				//cout<<" outside K loop for element index: "<<eleIndex<<endl;
			}
		}
		gsl_matrix_free(scaledxij);
	}
	//cout<<" outside loop"<<endl;
	gsl_matrix_free(g);
}

void Lumen::writeLumenJacobianToSystemJacobian(gsl_matrix* K, const std::vector <std::unique_ptr<Node>>& Nodes){
	int nNode=nodeIdsList.size();
	for (int nodeiIndex = 0; nodeiIndex<nNode; ++nodeiIndex){
		for (int nodejIndex = 0; nodejIndex<nNode; ++nodejIndex){
			int nodeIdi = nodeIdsList[nodeiIndex];
			int nodeIdj = nodeIdsList[nodejIndex];
			for (int dimi = 0; dimi <3 ; dimi++){
				if (Nodes[nodeIdi]->FixedPos[dimi] ){
					continue;
				}
				for (int dimj = 0; dimj <3 ; dimj++){
					if (Nodes[nodeIdj]->FixedPos[dimj] ){
						continue;
					}
					int indexiOnSystemK = nodeIdi*Dim + dimi;
					int indexjOnSystemK = nodeIdj*Dim + dimj;
					int indexiOnLumenK = nodeiIndex*Dim + dimi;
					int indexjOnLumenK = nodejIndex*Dim + dimj;
					double valueOnSystemK = gsl_matrix_get(K,indexiOnSystemK,indexjOnSystemK);
					double valueOnLumenK = gsl_matrix_get(Kv,indexiOnLumenK,indexjOnLumenK);
					double newValue = valueOnSystemK + valueOnLumenK;
					gsl_matrix_set(K,indexiOnSystemK,indexjOnSystemK,newValue);
				}
			}
		}
	}
}

void Lumen::growLumen(double currentTimeInSec){
	currentIdealVolume = initialIdealVolume * exp(growthRate*currentTimeInSec);
	cout<<" lumen ideal volume: "<<currentIdealVolume<<endl;
}
