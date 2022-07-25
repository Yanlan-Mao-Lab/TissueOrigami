- Open UserInterface/TissueFoldingUI.pro and modify according to your OS.
- Build

To run code with SimulationOnTheGo mode:
- Open terminal
- cd to the git folder
- NAMEOFBUILDFOLDER/Debug/TissueFoldingUI.app/Contents/MacOS/TissueFoldingUI -mode SimulationOnTheGo -i ToolBox/SampleModelInputs/modelinputTest -od ./Samples/SimulationOnTheGo > ./Samples/SimulationOnTheGo/tm

Note: The outputs will be saves in the /Samples/SimulationOnTheGo. For Linux, update path to the build folder accordingly.

To run code with DisplaySave mode:
- Open terminal
- cd to the git folder
- NAMEOFBUILDFOLDER/Debug/TissueFoldingUI.app/Contents/MacOS/TissueFoldingUI -mode DisplaySave -dInput ./Samples/DisplaySave/ > ./Samples/DisplaySave/tmp

Note: For Linux, update path to the build folder accordingly.
