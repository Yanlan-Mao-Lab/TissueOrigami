## Conducting the tests

The inputs for each test case `x` can be found in the respective `run0700x` directory.
In order to avoid hard-coded paths, the input files have had any absolute paths in them replaced by relative paths, _assuming that these files are now present in the same directory as the executable_ `TissueFolding`, when it is executed.

The `TissueFolding` executable itself is assumed to be in the `TissueFolding/` subdirectory, which has relative path `../../TissueFolding/` from this readme.

### Run 07007

#### File dependencies:

This run requires access to the following files:
- `smallRectangle.mesh`
- `Stiffness96hrRectangleWing_Reduction_0`
- `Stiffness96hrRectangleWing_Reduction_4`
- `ShapeChangeRate96hrRectangleWingZ_Reduction_3`
- `ShapeChangeRate96hrRectangleWingXY_Reduction_3`

The `.mesh` file will appear as a reference output for the mesh generation regression tests - we can later save storage space by pointing to this file, rather than having two copies of said file in the repository if we so wish.

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