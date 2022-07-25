1. Open UserInterface/TissueFoldingUI.pro and modify according to your OS.
2. Build
3.
3.1. To run code with SimulationOnTheGo mode:
3.1.1. Open terminal
3.1.2. cd to the git folder
3.1.3. NAMEOFBUILDFOLDER/Debug/TissueFoldingUI.app/Contents/MacOS/TissueFoldingUI -mode SimulationOnTheGo -i ToolBox/SampleModelInputs/modelinputTest -od ./Samples/SimulationOnTheGo > ./Samples/SimulationOnTheGo/tm

Note: The outputs will be saves in the /Samples/SimulationOnTheGo. For Linux, update path to the build folder accordingly.

3.2. To run code with DisplaySave mode:
3.2.1. Open terminal
3.2.2. cd to the git folder
3.2.3. NAMEOFBUILDFOLDER/Debug/TissueFoldingUI.app/Contents/MacOS/TissueFoldingUI -mode DisplaySave -dInput ./Samples/DisplaySave/ > ./Samples/DisplaySave/tmp

Note: For Linux, update path to the build folder accordingly.
