#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>

using namespace std;


/* 
 g++ -o Writer_JobSubmission Writer_JobSubmission.cpp
 g++ -o Writer_ContinueJobSubmission Writer_ContinueJobSubmission.cpp
 g++ -o Writer_DirectoryMakerCopier Writer_DirectoryMakerCopier.cpp
 ./Writer_JobSubmission
 ./Writer_ContinueJobSubmission
 ./Writer_DirectoryMakerCopier
*/
 
int main()
{
        //-----------set the inputs:----------------------------------------------------------------------//
        int initialindex=7005;
   	int sizeofset=5;
        int numberofrepetitions=1;
	
	int numThreads = 4;
        bool usingscratch=false;
        string mainpath = "/home/ucbpnkh/Scratch/bulkruns/";
      
        string executablenames[2000], modelinputnames[2000], ECMnames[2000];
	

	//executablenames[0]="/home/ucgamto/Scratch/Projects/TissueFolding/bin/TissueFoldingPeripodialRelease";	
	executablenames[0]="/home/ucbpnkh/Scratch/bin/TissueFolding";	

	
        for(int i=1;i<sizeofset;i++){executablenames[i]=executablenames[0];}
	for(int i=0;i<sizeofset;i++){
		stringstream ss;                   
                ss.fill('0');
                ss.width(5);
                ss<<i+initialindex;
		string s =(ss.str());
		modelinputnames[i]="/home/ucbpnkh/Scratch/bulkruns/ModelInputs/modelinput"+s;
		cout<<modelinputnames[i]<<endl;
	}
        ofstream tempjobsubmit;
        string head,counter;
        const char *name;

        int indexcounter=initialindex;
        string indexcounter_str;
        for (int i=0;i<numberofrepetitions;i++){
                for (int j=0;j<sizeofset;j++){
                        stringstream inter;
                        inter.fill('0');
                        inter.width(5);
                        inter<<indexcounter;
                        indexcounter_str=(inter.str());
                        string head="Job"+indexcounter_str;
                        name=head.c_str();
                        cout<<"head: "<<head<<" name: "<<name<<endl;
                        tempjobsubmit.open(name,ofstream::out);
			tempjobsubmit<<"#!/bin/bash -l"<<endl<<endl;
			tempjobsubmit<<"#the \"-l\" in bash setting is necessary to indicate login shell (it also must be on top of the script)s"<<endl;
			tempjobsubmit<<"#ulimit command sets the output limit to zerro, so that none of the spitted out error text will be recorded. This is to avoid filling up scratch under an unexpected error load "<<endl;
			tempjobsubmit<<"#the first 2 modules change the compiler to gcc"<<endl;
			tempjobsubmit<<"#the 3rd module loads gsl compatible with selected gsl version"<<endl;
			tempjobsubmit<<"#the 4th module loads openblas, so that lapack will work "<<endl;
			tempjobsubmit<<"#the 5th&6th modules are necessary to use ublas. 5th is blas itself, and it depends on pyhton hence 7th)"<<endl<<endl;
			tempjobsubmit<<"ulimit -c 0"<<endl;
			tempjobsubmit<<"module unload compilers mpi mkl"<<endl;
			tempjobsubmit<<"module load compilers/gnu/4.9.2"<<endl;
			//tempjobsubmit<<"module load mpi/openmpi/1.10.1/gnu-4.9.2"<<endl;
			tempjobsubmit<<"module load gsl/1.16/gnu-4.9.2 "<<endl;
			tempjobsubmit<<"module load openblas/0.2.14/gnu-4.9.2"<<endl;
			tempjobsubmit<<"module load python/2.7.9"<<endl;
			tempjobsubmit<<"module load boost/1_54_0/gnu-4.9.2"<<endl<<endl;
			tempjobsubmit<<"#the path settings are necessary for pardiso"<<endl;
			tempjobsubmit<<"#the OMP_NUM_THREADS is necessary for setting the maximum number of parallel threads"<<endl;
			tempjobsubmit<<"#the OPENBLAS_NUM_THREADS is necessary for setting the thread # of openBLAS separately."<<endl;								tempjobsubmit<<"#   the linked library is not compatible with openMP, if this env var is not set, it creates warning messages "<<endl;
			tempjobsubmit<<"#   for each call to an openBLAS function. "<<endl;

			tempjobsubmit<<"export PATH=/home/ucbpnkh/LocalLibs/:$PATH"<<endl;
			tempjobsubmit<<"export LD_LIBRARY_PATH=/home/ucbpnkh/LocalLibs/:$LD_LIBRARY_PATH"<<endl;
			tempjobsubmit<<"export LIBRARY_PATH=/home/ucbpnkh/LocalLibs/:$LIBRARY_PATH"<<endl;
			tempjobsubmit<<"export OMP_NUM_THREADS="<<numThreads<<endl;
			tempjobsubmit<<"export OPENBLAS_NUM_THREADS=1"<<endl;
			
			tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/Out"<<endl;
                        tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/Save*"<<endl;
			tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/tmp "<<endl;
			tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/cont/Out"<<endl;
                        tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/cont/Save*"<<endl;
                        tempjobsubmit<<endl;
			
                        tempjobsubmit<<executablenames[j]<<" -mode SimulationOnTheGo -i "<<modelinputnames[j]<<" -od "<<mainpath<<"Run"<<indexcounter_str<<"/ > "<<mainpath<<"Run"<<indexcounter_str<<"/tmp"<<endl;
                        //tempjobsubmit<<"rm "<<mainpath<<"Run"<<indexcounter_str<<"/tmp "<<endl;
                        tempjobsubmit.close();
                        indexcounter++;
                }
        }
}
   
