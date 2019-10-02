#!/bin/bash

baseModeliInputFile="./modelinputMutationBase"
radius=2.5
useAbsoluteGrowth=0
#growthFoldORRatePerHour=0.1493
growthFoldORRatePerHour=3
iRun=25

declare -A relativePosRange

relativePosRange[0,0]=0.05  #x0
relativePosRange[0,1]=0.50  #y0
relativePosRange[1,0]=0.15  #x1
relativePosRange[1,1]=0.50  #y1
relativePosRange[2,0]=0.25  #.
relativePosRange[2,1]=0.50  #.
relativePosRange[3,0]=0.35  #.
relativePosRange[3,1]=0.50  
relativePosRange[4,0]=0.45  
relativePosRange[4,1]=0.50  
relativePosRange[5,0]=0.55  
relativePosRange[5,1]=0.50  
relativePosRange[6,0]=0.65  
relativePosRange[6,1]=0.50  
relativePosRange[7,0]=0.75  
relativePosRange[7,1]=0.50  
relativePosRange[8,0]=0.85  
relativePosRange[8,1]=0.50  
relativePosRange[9,0]=0.95  
relativePosRange[9,1]=0.50 
 
relativePosRange[10,0]=0.15  
relativePosRange[10,1]=0.60  
relativePosRange[11,0]=0.25  
relativePosRange[11,1]=0.60  
relativePosRange[12,0]=0.35  
relativePosRange[12,1]=0.60 
relativePosRange[13,0]=0.45  
relativePosRange[13,1]=0.60  
relativePosRange[14,0]=0.55  
relativePosRange[14,1]=0.60  
relativePosRange[15,0]=0.65  
relativePosRange[15,1]=0.60 
relativePosRange[16,0]=0.75  
relativePosRange[16,1]=0.60  
relativePosRange[17,0]=0.85  
relativePosRange[17,1]=0.60  
relativePosRange[18,0]=0.95  
relativePosRange[18,1]=0.60 

relativePosRange[19,0]=0.25  
relativePosRange[19,1]=0.70  
relativePosRange[20,0]=0.35  
relativePosRange[20,1]=0.70 
relativePosRange[21,0]=0.45  
relativePosRange[21,1]=0.70  
relativePosRange[22,0]=0.55  
relativePosRange[22,1]=0.70  
relativePosRange[23,0]=0.65  
relativePosRange[23,1]=0.70 
relativePosRange[24,0]=0.75  
relativePosRange[24,1]=0.70  
relativePosRange[25,0]=0.85  
relativePosRange[25,1]=0.70  

relativePosRange[26,0]=0.35  
relativePosRange[26,1]=0.80 
relativePosRange[27,0]=0.45  
relativePosRange[27,1]=0.80  
relativePosRange[28,0]=0.55  
relativePosRange[28,1]=0.80  
relativePosRange[29,0]=0.65  
relativePosRange[29,1]=0.80 
relativePosRange[30,0]=0.75  
relativePosRange[30,1]=0.80  
relativePosRange[31,0]=0.85  
relativePosRange[31,1]=0.80

relativePosRange[32,0]=0.45  
relativePosRange[32,1]=0.90  
relativePosRange[33,0]=0.55  
relativePosRange[33,1]=0.90  
relativePosRange[34,0]=0.65  
relativePosRange[34,1]=0.90 
relativePosRange[35,0]=0.75  
relativePosRange[35,1]=0.90  

num_columns=35

iMutation=0
while [ $iMutation -le $num_columns ]
do
	printf -v CurrNo "%05d" $iRun
	cp  $baseModeliInputFile ./ModelInputs/modelinput$CurrNo
	echo "" >>./ModelInputs/modelinput$CurrNo
	echo "MutationOptions:" >>./ModelInputs/modelinput$CurrNo
	echo "  numberOfClones(int): 1" >>./ModelInputs/modelinput$CurrNo
	echo "  cloneInformation(double-relativeX,relativeY,micronRadius,usingAbsoluteGrowth(bool),growthRatePerHour_OR_growthFoldIncrease):" >>./ModelInputs/modelinput$CurrNo
	echo "  "${relativePosRange[$iMutation,0]}  ${relativePosRange[$iMutation,1]} $radius $useAbsoluteGrowth $growthFoldORRatePerHour >>./ModelInputs/modelinput$CurrNo
	echo $iRun ${relativePosRange[$iMutation,0]}-${relativePosRange[$iMutation,1]}, $growthFoldORRatePerHour fold
	iRun=$(( $iRun + 1 ))
	iMutation=$(( $iMutation + 1 ))
done
