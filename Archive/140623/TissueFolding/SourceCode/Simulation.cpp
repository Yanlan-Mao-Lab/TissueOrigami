
#include "Simulation.h"
#include "Prism.h"
#include <string.h>

using namespace std;

Simulation::Simulation(){
	currElementId = 0;
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	time =0.0;
};

Simulation::~Simulation(){
};

void Simulation::initiateNodes(){
	cout<<"initiated"<<endl;
	//initiateSinglePrismNodes();
	initiateNodesByRowAndColumn(5,3,2,2);
}

void Simulation::initiateElements(){
	//initiateSinglePrismElement();
	initiateElementsByRowAndColumn(5,3);
}

void Simulation::initiateSystemForces(){
	const int n = Nodes.size();
	SystemForces = new double*[n];
	for (int i=0;i<n;++i){
		SystemForces[i] = new double[3];
		SystemForces[i][0]=0.0;
		SystemForces[i][1]=0.0;
		SystemForces[i][2]=0.0;
		//cout<<"systemforces[i][j]: "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<endl;
	}
}

void Simulation::initiateSinglePrismNodes(){
	for (int i =0;i<6;++i){
		double* tmp_dp = new double[3];
		Nodes.push_back(tmp_dp);
	}
	Nodes[0][0]=0;Nodes[0][1]=1;Nodes[0][2]=0;
	Nodes[1][0]=1;Nodes[1][1]=0;Nodes[1][2]=0;
	Nodes[2][0]=0;Nodes[2][1]=0;Nodes[2][2]=0;
	Nodes[3][0]=0;Nodes[3][1]=1;Nodes[3][2]=1;
	Nodes[4][0]=1;Nodes[4][1]=0;Nodes[4][2]=1;
	Nodes[5][0]=0;Nodes[5][1]=0;Nodes[5][2]=1;
}

void Simulation::initiateSinglePrismElement(){
	int* NodeIds;
	NodeIds = new int[6];
	for (int i = 0; i < 6 ; i++){
		NodeIds[i]=i;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
	Elements.push_back(PrismPnt01);
	currElementId++;
}

void Simulation::initiateNodesByRowAndColumn(int Row, int Column, float SideLength, float zHeight){
	//The height of the equilateral triangle with side length: SideLength
	double sqrt3 = 1.7321;
	float h = sqrt3/2*SideLength;
	vector <double> xPos, yPos;
	for (int ColCount = 0; ColCount < Column+1; ++ColCount){
		double CurrY = ColCount*h;
		int CurrRowNum = Row - ColCount;
		double RowOffset = 0.5*SideLength*ColCount;
		for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
			double CurrX = RowOffset + RowCount * SideLength;
			xPos.push_back(CurrX);
			yPos.push_back(CurrY);
		}
		if (ColCount>0){
			CurrY = (-1.0)*CurrY;
			for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
				double CurrX = RowOffset + RowCount * SideLength;
				xPos.push_back(CurrX);
				yPos.push_back(CurrY);
			}
		}
	}
	int n =  xPos.size();
	//Adding the lower level of nodes:
	for (int i =0; i< n; ++i){
		double* tmp_dp = new double[3];
		tmp_dp[0] = xPos[i];
		tmp_dp[1] = yPos[i];
		tmp_dp[2] = 0.0;
		Nodes.push_back(tmp_dp);
	}
	//Adding the upper level:
	for (int i =0; i< n; ++i){
			double* tmp_dp = new double[3];
			tmp_dp[0] = xPos[i];
			tmp_dp[1] = yPos[i];
			tmp_dp[2] = zHeight;
			Nodes.push_back(tmp_dp);
	}
}

void Simulation::initiateElementsByRowAndColumn(int Row, int Column){
	int xinit1 = 0;
	int xinit2 = xinit1+Row;
	int xinit3 = 0;
	int xinit4 = xinit1+2*Row-1;
	int n = Nodes.size()/2.0;

	for (int ColCount = 0; ColCount < Column; ++ColCount){
		int CurrRowNum = Row - ColCount;
		for (int RowCount = 0; RowCount<CurrRowNum-1; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit1+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			currElementId++;
		}
		for (int RowCount = 0; RowCount<CurrRowNum-2; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit2+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId);
			Elements.push_back(PrismPnt01);
			currElementId++;
		}
		xinit1 = xinit2;
		xinit2 = xinit4 + CurrRowNum-1;
		xinit3 = xinit4;
		xinit4 = xinit2 + CurrRowNum-2;
	}
	//for (int i = 0; i< Elements.size(); ++i){
	//	cout<<"Element : "<<i<<" positions: "<<endl;
	//	Elements[i]->displayPositions();
	//	Elements[i]->displayIdentifierColour();
	//}
}

void Simulation::calculateSystemCentre(){
	int n = Nodes.size();
	for (int i = 0; i< n; ++i){
		SystemCentre[0] += Nodes[i][0];
		SystemCentre[1] += Nodes[i][1];
		SystemCentre[2] += Nodes[i][2];
	}
	SystemCentre[0]= SystemCentre[0]/n;
	SystemCentre[1]= SystemCentre[1]/n;
	SystemCentre[2]= SystemCentre[2]/n;
}

void Simulation::RunOneStep(){
	double viscosity = 100;
	double dt = 0.01;
	//cout<<"time: "<<time<<endl;
	if (time == 0.0){
		growSystem(2.0);
	}

	const int n = Nodes.size();
	const int dim = 3;
	//memset(SystemForces,0.0,sizeof(SystemForces[0][0])*n*dim);
	for (int i=0;i<n;++i){
		for (int j=0;j<dim;++j){
			SystemForces[i][j]=0.0;
		}
	}
	//aligning the centres:
	for (int i=0;i<Elements.size();++i){
		Elements[i]->alignReference();
		Elements[i]->calculateForces(1,SystemForces);
	}

	double velocities[n][dim];
	for (int i=0;i<n;++i){
		for (int j=0; j<dim; ++j){
			velocities[i][j] = SystemForces[i][j]/viscosity;
			Nodes[i][j] += (velocities[i][j])*dt;
		}
	}
	for (int i=0;i<Elements.size();++i){
		Elements[i]->updatePositions(Nodes);
	}
	/*
	cout<<"velocities:"<<endl;
	for (int i=0;i<n;++i){
		for (int j=0;j<dim;++j){
			cout<<velocities[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	cout<<"SystemForces:"<<endl;
	for (int i=1;i<n;++i){
		for (int j=0;j<dim;++j){
			cout<<SystemForces[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	*/
	time +=dt;
}

void Simulation::growSystem(float scale){
	cout<<"growing the system with scale factor: "<<scale<<endl;
	for (int i=0;i<Elements.size();++i){
		Elements[i]->growShape(scale);
	}
}
