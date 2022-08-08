#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>


using namespace std;

/* 
 g++ -o Writer_ContinueJobSubmission Writer_ContinueJobSubmission.cpp 
 g++ -o Writer_JobSubmission Writer_JobSubmission.cpp
 g++ -o Writer_DirectoryMakerCopier Writer_DirectoryMakerCopier.cpp
 ./Writer_JobSubmission
 ./Writer_ContinueJobSubmission
 ./Writer_DirectoryMakerCopier
*/

int main()
{
        ofstream copierscript;
        int initialindex=7005;
   	int sizeofset=5;
        int numberofrepetitions=1;
	
	
	
    	string mainpath = "/home/ucbpnkh/Scratch/bulkruns/";
	bool makedirectories=true;

        ofstream outputfile;
        outputfile.open("copierscript",ofstream::out);
        int indexcounter=initialindex;
        string indexcounter_str;
        for (int i=0;i<numberofrepetitions;i++){
                for (int j=0;j<sizeofset;j++){
                        stringstream inter;
                        inter.fill('0');
                        inter.width(5);
                        inter<<indexcounter;
                        indexcounter_str=(inter.str());
                        if(makedirectories){
                                outputfile<<"mkdir "<<mainpath<<"Run"<<indexcounter_str<<endl;
                        }
                        outputfile<<"mv Job"<<indexcounter_str<<" "<<mainpath<<"Run"<<indexcounter_str<<"/"<<endl;
			outputfile<<"mv JobContinue"<<indexcounter_str<<" "<<mainpath<<"Run"<<indexcounter_str<<"/"<<endl;
                        indexcounter++;
                }
        }
}
