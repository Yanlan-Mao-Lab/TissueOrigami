## Run 07007

#### Hard-paths changelog

Input file changes have been conducted, so as to remove hard-coded paths.

```bash
InputMeshParameters:
    MeshFile(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/MeshGeneration/2DEllipse/smallRectangle.mesh

YoungsModulusTimeseriesGrids:
    Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/StiffnessTimeSeries/Stiffness96hrRectangleWing_Reduction_0
    WhenToApplyInput(sec): 0
    Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/StiffnessTimeSeries/Stiffness96hrRectangleWing_Reduction_4
    WhenToApplyInput(sec): 600

GrowthFunctionType(int-seeDocumentation):
    InitialTime(sec): 120
    Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/GrowthRates/ShapeChangeRate96hrRectangleWingZ_Reduction_3

GrowthFunctionType(int-seeDocumentation): 3
    InitialTime(sec): 60
    Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/GrowthRates/ShapeChangeRate96hrRectangleWingXY_Reduction_3
```