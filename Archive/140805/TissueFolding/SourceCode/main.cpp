//#include <iostream>
using namespace std;

#include "Simulation.h"
#include <vector>

int main(int argc, char * argv[])
{
	Simulation Sim01;
	Sim01.initiateSystem();
	Sim01.Elements[0]->displayName();
	Sim01.Elements[0]->displayPositions();
}

