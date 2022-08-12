using namespace std;

#include "Simulation.h"
#include "ArgumentParserSimulation.h"
#include <vector>

int main(int argc, char * argv[])
{
	// read command-line arguments
	ArgumentSpace simArgs = ArgumentReader::readInput(argc, argv);

	// command-line arguments are valid in the sense that they all exist,
	// now we need to pass this information into the simulation
	Simulation* Sim01 = new Simulation();
	Sim01->displayIsOn = false;

	if (simArgs.getSimulationMode()==Default) {
		// default settings only require DisplaySave set to false
		Sim01->DisplaySave = false;
	}
	else 
	{
		// we will need to read in the inputs
		Sim01->readExecutableInputs(argc, argv);
	}

	if (Sim01->DisplaySave){
		cerr<<"This is the executable for running the model without display"<<endl;
		return true;
	}
	else{
		Sim01->initiateSystem();
		int n = Sim01->Elements.size();
		for (int i=0; i<n; ++i){
			Sim01->Elements[i]->updatePositions(Sim01->Nodes);
		}
		cout<<"Initiating simulation in the background"<<endl;
		while (Sim01->currSimTimeSec <= Sim01->SimLength){
			//cout<<"running step: "<<Sim01.timestep<<", this is time: "<<Sim01.timestep*Sim01.dt<<" sec"<<endl;
			bool Success = Sim01->runOneStep();
			if (Success == false ){
				break;
			}
		}
		Sim01->wrapUpAtTheEndOfSimulation();
		Sim01->writeRelaxedMeshFromCurrentState();
	}

	delete Sim01;
	cout<<"Finished Simulation"<<endl;
}

