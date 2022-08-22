For small rectangle and sphere, the inputs for the cpp files which were generated in Matlab are in the folder "/TissueOrigami/ToolBox/MeshGeneration/2DEllipsecppInputs/cppInputs". 
For small wing disc, the original input is the "48hrDiscSymmetricOutline" which is stored in /TissueOrigami/ToolBox/MeshGeneration/2DEllipse/inputOutlines. From this "Points.1.ele" and "Points.1.node" are created (please see details below). I have stored these in "/TissueOrigami/ToolBox/MeshGeneration/2DEllipse/cppInputs".
All .mesh outputs are in folder "/TissueOrigami/ToolBox/MeshGeneration/2DEllipse/cppOutputs".

Small Rectangle:
- the "smallRectangle.ele" and "smallRectangle.node" ccpInput files have to be moved to the folder /TissueOrigami/ToolBox/MeshGeneration/2DEllipse and renamed to "Points.1.ele" and "Points.1.node"

- In EllipseFromOutline.cpp, manually change the following parameters:
    selectTissueType=1
    then in condition: else if(selectTissueType == 1)
    actinHeight = 2.0;
    ECMHeight = 0.2; 
    
- cd to 2DEllipse folder and run the following commands:
    g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
    cp ./src/EllipseFromOutline ./
    ./EllipseFromOutline -1 5.2 2 3 0

This will create the outout "MeshFile.out". I have renamed this to "smallRectangle.mesh" and put it in the cppOutput folder.


SphericalTriangulation:
- copy the "SphericalTriangulation" cppInput to the folder  /TissueOrigami/ToolBox/MeshGeneration/2DEllipse 

-In EllipseFromOutline.cpp, manually change the following parameters:
selectTissueType=5
then in condition: else if(selectTissueType == 1)
actinHeight = 1.0;
ECMHeight = 0.2; 

- cd to 2DEllipse folder and run the following commands:
g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
cp ./src/EllipseFromOutline ./
./EllipseFromOutline -2 4.2 2 3 0

I have renamed the output of this to "smallSphere.mesh"

Small wing disc:
- In EllipseFromOutline.cpp, manually change the following parameters:
    -  selectTissueType=1
    - then in condition: else if(selectTissueType == 1)
        actinHeight = 2.0;
        ECMHeight = 0.2; 
    - Also, all folder paths need to be manually changed. They currently point to my local folder.


- cd to 2DEllipse folder and run the following commands:
g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
cp ./src/EllipseFromOutline ./
./EllipseFromOutline 1 35.56 50 27.3 27.3 12.5 9 2 0  ./inputOutlines/48hrDiscSymmetricOutline
this will triangulate the wing disc outline and create Points.1.node and Points.1.ele. I put these in "cppInputs" folder.
Then run:
./EllipseFromOutline -1 5.2 1.0 3 1 

I have renamed the output of this to "smallWingDisc.mesh"
