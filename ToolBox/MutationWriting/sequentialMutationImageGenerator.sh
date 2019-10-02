#!/bin/bash

#changes in paintGL:
#     DisplayPysProp = true;
#     PysPropToDisplay = 4;
#     PerspectiveView = false;
#     obj_pos[0] =  15.0f;
#     orthoViewLimits[0] = -250;
#     orthoViewLimits[1] =  250;
#     orthoViewLimits[2] = -130;//-130;//-60;
#     orthoViewLimits[3] =  130;//130;// 60;
#     orthoViewLimits[4] = -500;
#     orthoViewLimits[5] =  500;
#     orthoViewLimits[2] += 5.0*17;
#     orthoViewLimits[3] -= 5.0*17;
#     orthoViewLimits[0] = orthoViewLimits[2]*aspectratio;
#     orthoViewLimits[1] = orthoViewLimits[3]*aspectratio;

iRun=25
n=60

while [ $iRun -le $n ]
do
	printf -v CurrNo "%05d" $iRun
	~/Documents/TissueFolding/UserInterface/Debug/TissueFoldingUI -mode SimulationOnTheGo -i /home/melda/Documents/TissueFolding/ToolBox/MutationWriting/ModelInputs/modelinput$CurrNo -od /home/melda/Documents/TissueFolding/ToolBox/MutationWriting/runFolder/
	cp /home/melda/Documents/TissueFolding/ToolBox/MutationWriting/runFolder/ScreenShots/frame000002-2sec.png  ~/Pictures/$iRun.png
	iRun=$(( $iRun + 1 ))
done
