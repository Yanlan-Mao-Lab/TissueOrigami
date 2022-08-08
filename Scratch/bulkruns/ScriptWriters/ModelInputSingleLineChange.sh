#sleep  14400
i=28
n=43

while [ $i -le $n ]
do
	printf -v CurrNo "%05d" $i
	
	sed -i 's/LumenBulkModulus(Pa): 0/LumenBulkModulus(Pa): 0 \n  LumenGrowthRate(foldPer24hr): 0/' ../ModelInputs/modelinput$CurrNo
 	sed -i 's/LumenBulkModulus(Pa): 9600/LumenBulkModulus(Pa): 9600 \n  LumenGrowthRate(foldPer24hr): 0/' ../ModelInputs/modelinput$CurrNo
        sed -i 's/LumenBulkModulus(Pa): 19200/LumenBulkModulus(Pa): 19200 \n  LumenGrowthRate(foldPer24hr): 0/' ../ModelInputs/modelinput$CurrNo
        sed -i 's/LumenBulkModulus(Pa): 32000/LumenBulkModulus(Pa): 32000 \n  LumenGrowthRate(foldPer24hr): 0/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/timeOfStiffnessChange(hr): 30 36/timeOfStiffnessChange(hr): 30 32/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/melda\/Documents\/TissueFolding\/ToolBox/ucgamto\/Scratch\/Projects\/TissueFolding\/bulkruns/' ../ModelInputs/modelinput$CurrNo 
	#sed -i 's/MeshGeneration\/2DEllipse/MeshFiles/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/3-0/3-3/' ../ModelInputs/modelinput$CurrNo 	
	#sed -i 's/3-0/3-3/' ../ModelInputs/modelinput$CurrNo 	
	#sed -i 's/Grid20-11_Col_48-72_1-1/Grid20-11_Col_48-72_3-1/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/Grid20-11_Col_56-80_1-1/Grid20-11_Col_56-80_3-1/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/Grid20-11_Col_64-88_1-0/Grid20-11_Col_64-88_3-0/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/Grid20-11_Col_72-96_2-0/Grid20-11_Col_72-96new_2-0/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/Grid20-11_Col_56-88_1-0/Grid20-11_Col_56-88_3-0/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/Per_48-72_6-0/Per_48-72_7-0/' ../ModelInputs/modelinput$CurrNo
	
	#sed -i 's/Exp48hr-noPeri02ECM_10.mesh/Exp48hr-noPeri-ThinECM_11.mesh/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/Exp48hr-noPeri_10.mesh/Exp48hr-noPeri_11.mesh/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/TimeStep(sec): 600/TimeStep(sec): 1800/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/TimeStep(sec): 300/TimeStep(sec): 1800/' ../ModelInputs/modelinput$CurrNo
        #sed -i 's/TimeStep(sec): 900/TimeStep(sec): 1800/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/HingeMutation25/HingeMutation10/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/HeightAt96hr_/HeightAt96from80hr_/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/DataSaveInterval(sec):  3600/DataSaveInterval(sec):  1800/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ECMRemodellingHalfLife(hour): 12.0/ECMRemodellingHalfLife(hour): 8.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ECMColumnarYoungsModulus:  6400/ECMColumnarYoungsModulus:  1600/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/HeightAt72hr_Growth_Z/HeightAt80hr_Growth_Z/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/HeightAt96from72hr_Growth_Z-only/HeightAt96from80hr_Growth_Z-only/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/HeightAt80hr_Growth_Z/HeightAt72hr_Growth_Z/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/HeightAt96hr_Growth_Z-only/HeightAt96from72hr_Growth_Z-only/' ../ModelInputs/modelinput$CurrNo

	
	#sed -i 's/SimulationLength(sec): 21600/SimulationLength(sec): 180000/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SimulationLength(sec): 180000/SimulationLength(sec): 21600/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/ThereIsAdhesion(bool): 0/ThereIsAdhesion(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ThereIsNodeCollapse(bool): 0/ThereIsNodeCollapse(bool): 1/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/TimeStep(sec): 600/TimeStep(sec): 300/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SimulationLength(sec): 86400/SimulationLength(sec): 180000/' ../ModelInputs/modelinput$CurrNo
		
	#sed -i 's/TimeStep(sec): 900/TimeStep(sec): 600/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SimulationLength(sec): 21600/SimulationLength(sec): 57600/' ../ModelInputs/modelinput$CurrNo
	
	
	#sed -i 's/TimeStep(sec): 900/TimeStep(sec): 600/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SimulationLength(sec): 28800/SimulationLength(sec): 86400/' ../ModelInputs/modelinput$CurrNo
	
	#sed -i 's/TimeStep(sec): 600/TimeStep(sec): 150/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SimulationLength(sec): 86400/SimulationLength(sec): 180000/' ../ModelInputs/modelinput$CurrNo
	
	sed -i 's/AddLateralECM(bool): 1/AddLateralECM(bool): 0/' ../ModelInputs/modelinput$CurrNo
	
	#sed -i 's/DataSaveInterval(sec):  600/DataSaveInterval(sec):  2000/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/Per_/Rate70_Per_/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ThereIsPlasticDeformation(bool): 0/ThereIsPlasticDeformation(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ThereIsECMRemodellinbg(bool): 0/ThereIsECMRemodellinbg(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/GridGrowthsPinnedOnInitialMesh(bool):	0/GridGrowthsPinnedOnInitialMesh(bool):	1/' ../ModelInputs/modelinput$CurrNo
        #sed -i 's/AddCurvature(bool): 1/AddCurvature(bool): 0/' ../ModelInputs/modelinput$CurrNo   
	#sed -i 's/CurvatureDepthAtCentre(double-microns): -2.0/CurvatureDepthAtCentre(double-microns): 2.0/' ../ModelInputs/modelinput$CurrNo 
	#sed -i 's/PeripodialMembraneYoungsModulus: 250.0/PeripodialMembraneYoungsModulus: 750.0/' ../ModelInputs/modelinput$CurrNo  
	#sed -i 's/LinkerZoneApicalYoungsModulus: 1000.0/LinkerZoneApicalYoungsModulus: 2000.0/' ../ModelInputs/modelinput$CurrNo  
	#sed -i 's/LinkerZoneBasalYoungsModulus: 450/LinkerZoneBasalYoungsModulus: 3600/' ../ModelInputs/modelinput$CurrNo  
	
	#sed -i 's/YoungsModulusApical: 100.0/YoungsModulusApical: 25.0/' ../ModelInputs/modelinput$CurrNo
        #sed -i 's/YoungsModulusApical: 50.0/YoungsModulusApical: 25.0/' ../ModelInputs/modelinput$CurrNo
        #sed -i 's/YoungsModulusApical: 26.0/YoungsModulusApical: 25.0/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/YoungsModulusBasal: 25.0/YoungsModulusBasal: 100.0/' ../ModelInputs/modelinput$CurrNo  
	#sed -i 's/YoungsModulusMid: 25.0/YoungsModulusMid: 26.0/' ../ModelInputs/modelinput$CurrNo        
	#sed -i 's/PoissonsRatio: 0.3 /PoissonsRatio: 0.45/' ../ModelInputs/modelinput$CurrNo   
	#sed -i 's/LinkerBasalCircumferenceFix(bool-x,y,z): 0 0 1/LinkerBasalCircumferenceFix(bool-x,y,z): 0 0 0/' ../ModelInputs/modelinput$CurrNo  
	#sed -i 's/LinkerZoneBasalViscosity: 1000.0/LinkerZoneBasalViscosity: 100.0/' ../ModelInputs/modelinput$CurrNo 
	#sed -i 's/LinkerZoneBasalViscosity: 10000.0/LinkerZoneBasalViscosity: 100.0/' ../ModelInputs/modelinput$CurrNo 
	#sed -i 's/LinkerZoneApicalViscosity: 1000.0/LinkerZoneApicalViscosity: 100.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ThereIsMyosinFeedback(bool): 0/ThereIsMyosinFeedback(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/MyosinFeedbackCap: 200/MyosinFeedbackCap: 400/' ../ModelInputs/modelinput$CurrNo
	
	#sed -i 's/DiscProperApicalExternalViscosity: 10.0/DiscProperApicalExternalViscosity: 16000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/DiscProperBasalExternalViscosity: 0.0/DiscProperBasalExternalViscosity: 8000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/PeripodialMembraneApicalExternalViscosity: 5000.0/PeripodialMembraneApicalExternalViscosity: 10000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/PeripodialMembraneBasalExternalViscosity: 5000.0/PeripodialMembraneBasalExternalViscosity: 10000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/LinkerZoneApicalExternalViscosity: 5000.0/LinkerZoneApicalExternalViscosity: 10000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/LinkerZoneBasalExternalViscosity: 5000.0/LinkerZoneBasalExternalViscosity: 10000.0/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/CollapseNodesOnAdhesion(bool): 0/CollapseNodesOnAdhesion(bool): 1/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/ApicalViscosity:  0.0/ApicalViscosity:  100000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/BasalViscosity:   0.0/BasalViscosity:   100000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/MidLineViscosity: 0.0/MidLineViscosity: 100000.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/Pouch50percentZgrowth/PouchAndNotum50percentZgrowth/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/FixWithHighExternalViscosity(bool): 0/FixWithHighExternalViscosity(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/FixVis(x,y,z): 100000   100000  100000/FixVis(x,y,z): 0   0  250000/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/  CircumferenceFix(bool-x,y,z): 0 0 1/  CircumferenceFix(bool-x,y,z): 0 0 0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/bindCircumferenceXYToBasal(bool): 0/bindCircumferenceXYToBasal(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/  BasalCircumferenceFix(bool-x,y,z): 0 0 1/  BasalCircumferenceFix(bool-x,y,z): 0 0 0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ExtendToWholeTissue: 1/ExtendToWholeTissue: 0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/AddSoftPeriphery(bool): 0/AddSoftPeriphery(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ApplyToBasalSurface(bool): 0/ApplyToBasalSurface(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SoftPeripheryRange(double-microns): 10.0/SoftPeripheryRange(double-microns): 6.0/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/SoftnessFraction(double-fraction): 0.2/SoftnessFraction(double-fraction): 3.0/' ../ModelInputs/modelinput$CurrNo 
  	#sed -i 's/ThereIsECMSoftening(bool): 0/ThereIsECMSoftening(bool): 1/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/timeOfSoftening(hr): 30  36.0/timeOfSoftening(hr): 1  48.0/' ../ModelInputs/modelinput$CurrNo
  	#sed -i 's/2 0.55 0.65 0.45 0.50/1 0.0 1.05/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/ThereIsStiffnessPerturbation(bool): 0/ThereIsStiffnessPerturbation(bool): 1/' ../ModelInputs/modelinput$CurrNo	
	#sed -i 's/ThereIsECMSoftening(bool): 0/ThereIsECMSoftening(bool): 1/' ../ModelInputs/modelinput$CurrNo

	#sed -i 's/ApplyToPeripodialMembrane(bool): 1/ApplyToPeripodialMembrane(bool): 0/' ../ModelInputs/modelinput$CurrNo	
	 
	
	#sed -i 's/stiffnessChangeAppliedToEllipses(number,\[ellipseId\]\[ellipseId\]): 1 0/stiffnessChangeAppliedToEllipses(number,\[ellipseId\]\[ellipseId\]): 2 1 2/' ../ModelInputs/modelinput$CurrNo
	#sed -i 's/stiffnessPerturbationAppliedToEllipses(number,\[ellipseId\]\[ellipseId\]): 1 0/stiffnessPerturbationAppliedToEllipses(number,\[ellipseId\]\[ellipseId\]): 2 1 2/' ../ModelInputs/modelinput$CurrNo
	echo ../ModelInputs/modelinput$CurrNo

	i=$(( $i + 1 ))
done
