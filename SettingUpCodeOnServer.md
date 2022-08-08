Setting up directories on Myriad:

The main directory includes the following files/folders:
1.    LocalLibs folder: This currently includes an "emptyFile" but should contain the pardiso package (.so files)
2.    Pardiso.lic file: This file should include the license key for the pardiso package
3.    Scratch folder: This includes the following sub-folders:
  -   bin folder: This currently contains an "emptyFile" but the executable will be copied here automatically.
  -   bulkruns: I have copied the contents of this folder to github. This is where are the inputs and outputs of the code are stored.
  -   TissueFolding: This folder contains the actual source code. The version I currently have on the Server has all the sparse matrix solver references enabled and needs pardiso to run. I have copied this folder on github under the Myriad branch, but for the purpose of testing without pardiso, you can copy the version on ARCCollaboration branch where all references to pardiso have been commented out.

Creating executables:

1. Go to ./Scratch/TissueFolding and open TissueFolding.pro.
2. Update the file paths to direct to your local directory (currently they point to mine).
3. Details of the steps to follow to make executables when running the code for the first time should be in TissueFolding.pro. In brief, these commands should be run in order: 
- module load qt/4.8.6/gnu-4.9.2
- qmake TissueFolding.pro
- module unload compilers mpi mkl
- module load compilers/gnu/4.9.2
- module load gsl/1.16/gnu-4.9.2
- module load openblas/0.2.14/gnu-4.9.2
- module load python/2.7.9
- module load boost/1_54_0/gnu-4.9.2
- export PATH=/home/UCLUSERID/LocalLibs/:$PATH
- export LD_LIBRARY_PATH=/home/UCLUSERID/LocalLibs/:$LD_LIBRARY_PATH
- export LIBRARY_PATH=/home/UCLUSERID/LocalLibs/:$LIBRARY_PATH
- export OMP_NUM_THREADS=4
- export OPENBLAS_NUM_THREADS=1
- make
- cd ../bin
- rm -r TissueFolding
- scp -r ../TissueFolding/Debug/TissueFolding ./

Automatic job submssion:

1.    Go to ./Scratch/bulkruns/ScriptWriters.
2.    In the three .cpp files below modify initialindex and sizeofset to include the modelinput ids that you want to run. For example, to run simulations using modelinput07003 to modelinput07006, initialindex should be 7003 and sizeofset should be 3.
-    Writer_ContinueJobSubmission.cpp
-    Writer_DirectoryMakerCopier.cpp
-    Writer_JobSubmission.cpp
3.    Run these commands:
- g++ -o Writer_JobSubmission Writer_JobSubmission.cpp
- g++ -o Writer_ContinueJobSubmission Writer_ContinueJobSubmission.cpp
- g++ -o Writer_DirectoryMakerCopier Writer_DirectoryMakerCopier.cpp
- ./Writer_JobSubmission
- ./Writer_ContinueJobSubmission
- ./Writer_DirectoryMakerCopier
- bash copierscript
4.    cd to bulkruns and modify submitter.sh to include the files that you want to run. Then:
- bash submitter.sh

Note: Please note that paths to the code's inputs should also be updated in modelinput files in Scratch/bulkruns/ModelInputs/. These files currently point to my Myriad directory.


