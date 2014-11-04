using namespace std;

#include "Simulation.h"
#include <vector>

int main(int argc, char * argv[])
{
	Simulation* Sim01 = new Simulation();
	Sim01->displayIsOn = false;
	if (argc<2){
		Sim01->DisplaySave = false;
		cerr<<"Using default settings"<<endl;
	}
	else{
		bool Success = Sim01->readExecutableInputs(argc, argv);
		if (!Success){
			cerr<<"Error in input to executable"<<endl;
			return 1;
		}
	}
	if (Sim01->DisplaySave){
		cerr<<"This is the executable for running the model without display"<<endl;
		return true;
	}
	else{
		Sim01->initiateSystem();
		cout<<"Initiating simulation in the background"<<endl;
		while (Sim01->timestep < Sim01->SimLength){
			//cout<<"running step: "<<Sim01.timestep<<", this is time: "<<Sim01.timestep*Sim01.dt<<" sec"<<endl;
			Sim01->runOneStep();
		}
	}
	delete Sim01;
	cout<<"Finished Simulation"<<endl;
}

