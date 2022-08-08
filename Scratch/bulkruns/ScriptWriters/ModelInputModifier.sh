function addPlasticity {
	CurrNo=$1
	sed -i 's/ThereIsPlasticDeformation(bool): 0/ThereIsPlasticDeformation(bool): 1/' ../ModelInputs/modelinput$CurrNo
}

function addSoftPeripodial {
	CurrNo=$1
	sed -i 's/PeripodialMembraneYoungsModulus: 1000.0/PeripodialMembraneYoungsModulus: 250/' ../ModelInputs/modelinput$CurrNo
}

function addViscosity {
	CurrNo=$1
	sed -i 's/PeripodialMembraneApicalViscosity: 0.0/PeripodialMembraneApicalViscosity: 10.0/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/PeripodialMembraneBasalViscosity: 0.0/PeripodialMembraneBasalViscosity: 100.0/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ ApicalViscosity: 0.0/ ApicalViscosity: 10.0/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ BasalViscosity: 0.0/ BasalViscosity: 100.0/' ../ModelInputs/modelinput$CurrNo
}

function addSoftEdges {
	CurrNo=$1
	sed -i 's/AddSoftPeriphery(bool): 0/AddSoftPeriphery(bool): 1/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/SoftPeripheryRange(double-microns): 30.0/SoftPeripheryRange(double-microns): 10.0/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/SoftnessFraction(double-fraction): 0.1/SoftnessFraction(double-fraction): 0.5/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ApplyToApicalSurface(bool): 1/ApplyToApicalSurface(bool): 1/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ApplyToBasalSurface(bool): 0/ApplyToBasalSurface(bool): 0/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ApplyToColumnarLayer(bool): 1/ApplyToColumnarLayer(bool): 1/' ../ModelInputs/modelinput$CurrNo
	sed -i 's/ApplyToPeripodialMembrane(bool): 1/ApplyToPeripodialMembrane(bool): 1/' ../ModelInputs/modelinput$CurrNo
}

function createNewSet {
	ModelInputBase=$1
	setSize=$2
	next=$3
	iterator=0
	while [ $iterator -lt $setSize ]
	do
      	 	printf -v CurrNoBase "%05d" $ModelInputBase
      	 	printf -v CurrNoNew "%05d" $next
		#echo writing new modelinput modelinput$CurrNoNew from modelinput$CurrNoBase
 		cp ../ModelInputs/modelinput$CurrNoBase ../ModelInputs/modelinput$CurrNoNew
		iterator=$(( $iterator + 1 ))
		ModelInputBase=$(( $ModelInputBase + 1 ))
		next=$(( $next + 1 ))
	done
}

function displaySet {
	setSize=$1
	next=$2
	iterator=0
	while [ $iterator -lt $setSize ]
	do
      	 	printf -v CurrNoNew "%05d" $next
		echo modelinput$CurrNoNew $3 $4 $5 $6
 		next=$(( $next + 1 ))
		iterator=$(( $iterator + 1 ))
	done	
}

function addChangeToSet {
	ModelInputNew=$1
	setSize=$2
	#Adding plasticity and  soft peripodial membrane
	iterator=0
	while [ $iterator -lt $setSize ]
	do
      	 	printf -v CurrNoNew "%05d" $ModelInputNew
		#echo adding change: $3 to modelinput$CurrNoNew
		$3 $CurrNoNew
		iterator=$(( $iterator + 1 ))
		ModelInputNew=$(( $ModelInputNew + 1 ))
	done
}

ModelInputNo=53
setSize=12
nextNo=66
n=15   

    plasticityIndex=(1 0 0 0 1 1 1 1 1 1 0 0 0 0 1)
softPeripodialIndex=(0 1 0 0 1 0 0 1 1 0 1 1 1 0 1)
     viscosityIndex=(0 0 1 0 0 1 0 1 0 1 1 1 0 1 1)
     softEdgesIndex=(0 0 0 1 0 0 1 0 1 1 1 0 1 1 1)  
      
     
echo ${plasticityIndex[*]}
echo ${softPeripodialIndex[*]}
echo ${viscosityIndex[*]}
echo ${softEdgesIndex[*]}


i=0
while [ $i -lt $n ]
do
	#creating set to add change:
	createNewSet $ModelInputNo $setSize $nextNo
	if [ ${plasticityIndex[$i]} -eq 1 ]
	then
		addChangeToSet $nextNo $setSize addPlasticity
	fi
	
	if [  ${softPeripodialIndex[$i]}  -eq 1 ]
	then
		addChangeToSet $nextNo $setSize addSoftPeripodial
	fi
	
	if [  ${viscosityIndex[$i]}  -eq 1 ]
	then
		addChangeToSet $nextNo $setSize addViscosity
	fi
	
	if [  ${softEdgesIndex[$i]}  -eq 1 ]
	then
		addChangeToSet $nextNo $setSize addSoftEdges
	fi
	displaySet $setSize $nextNo ${plasticityIndex[$i]} ${softPeripodialIndex[$i]} ${viscosityIndex[$i]} ${softEdgesIndex[$i]}
	i=$(( $i + 1 ))
	nextNo=$(( $nextNo + $setSize ))
done


    plasticityIndex=(1 0 0 0 1 1 1 1 1 1 0 0 0 0 1)
softPeripodialIndex=(0 1 0 0 1 0 0 1 1 0 1 1 1 0 1)
     viscosityIndex=(0 0 1 0 0 1 0 1 0 1 1 1 0 1 1)
     softEdgesIndex=(0 0 0 1 0 0 1 0 1 1 1 0 1 1 1)  
     

