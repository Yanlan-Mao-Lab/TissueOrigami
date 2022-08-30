# Regression Tests for Non-Pardiso Build

When building without Pardiso, the system of evolutionary equations is not solved.
Whilst this means that the computations are of little interest, it is an important benchmark in establishing that, up to linking to Pardiso, the `TissueFolding` executable has been successfully built.

Test runs are given a unique number, the inputs to each run being organised into subfoldered named with this unique identifier. Any generated outputs when running the regression tests are placed into the relevant folder.
That is, for `run0700x`, the following steps are conducted:
- Input files from the `run0700x` folder are copied to the directory containing the `TissueFolding` executable.
- The `TissueFolding` executable (without Pardiso) is called to produce outputs
- The outputs are moved to the `run07007/genertd_out` subfolder
- The contents of the `run07007/genertd_out` directory are compared to their homonymous counterparts in `run0700x/expected_out`.

The comparison of the output files uses the Python function `filecmp.cmp`, which can act on binary files and non-binary files.
Non-binary output files must satisfy `filecmp.cmp(reference, generated, shallow=False)==True` to pass.
The simulation outputs that are produced as binary files are:
 - Force
 - Growth
 - GrowthRate
 - GrowthRedistribution
 - Packing
 - PhysicalProp
 - SpecificElementAndNodeTypes
 - TensionCompression
When comparing these generated outputs to their reference counterparts, we allow for subtle differences in the generated binary files.
In the event that the binary files _do not_ match (according to `filecmp.cmp`), the following is done:
1. Reading the values stored in the binary back into the simulation (for example, if we were resuming a previously stopped simulation)
1. Writing them to a `.txt` file
1. Applying the same proceedure to the reference binary file
Obtaining `True` from acting `filecmp.cmp` on the resulting text files will allow the test to pass.
This event indicates that there is a difference between a binary file and its counterpart, but _not_ between the data that is read into the simulation. 
`pytest` will throw a warning detailing the file which differed.


## Conducting the tests

The inputs for each test case `x` can be found in the respective `run0700x` directory.
In order to avoid hard-coded paths, the input files have had any absolute paths in them replaced by relative paths, _assuming that these files are now present in the same directory as the executable_ `TissueFolding`, when it is executed.
The tests are conducted through the `test_SimulationNoPardiso.py` file.

The `TissueFolding` executable itself is assumed to be in the `TissueFolding/` subdirectory, which has relative path `../../TissueFolding/` from this readme.

The testing criteria may also require us to convert the outputs of the simulation from binary files to text files, in order to assert whether the simulation will read in identical values, and the binary-file difference is superflous.
With regards to this, the `simOutputsToTxt.o` executable must be located in the `tests/` directory, and have been built prior to running `pytest`.

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

### Run 07008

#### File dependencies:

This run requires access to the following files:
- `smallWingDisc.mesh`

The `.mesh` file will appear as a reference output for the mesh generation regression tests - we can later save storage space by pointing to this file, rather than having two copies of said file in the repository if we so wish.

#### Hard-paths changelog

Input file changes have been conducted, so as to remove hard-coded paths.

```bash
InputMeshParameters:
  MeshFile(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/MeshGeneration/2DEllipse/smallWingDisc.mesh
```

### Run 07009

#### File dependencies:

This run requires access to the following files:
- `smallSphere.mesh`
- `StiffnessMatrix_200_200_homogeneous.txt`
- `StiffnessMatrix_200_200_20_new.txt` - Filename changed from `StiffnessMatrix_200_200_20%_new.txt`

The `.mesh` file will appear as a reference output for the mesh generation regression tests - we can later save storage space by pointing to this file, rather than having two copies of said file in the repository if we so wish.

#### Hard-paths changelog

Input file changes have been conducted, so as to remove hard-coded paths.

```bash
InputMeshParameters:
  MeshFile(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/MeshGeneration/2DEllipse/smallSphere.mesh

YoungsModulusTimeseriesGrids:
  Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/StiffnessTimeSeries/StiffnessMatrix_200_200_homogeneous.txt
  WhenToApplyInput(sec): 0
  Filename(full-path): /Users/nkhalilgharibi/TissueFoldingProject/YanlanMaoLabRepo/TissueOrigami/ToolBox/StiffnessTimeSeries/StiffnessMatrix_200_200_20%_new.txt
  WhenToApplyInput(sec): 6000
```
