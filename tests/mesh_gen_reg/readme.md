# Regression Tests for mesh generation executable

The `mesh_gen_reg` folder contains the scripts and files for running regression tests on the mesh-generation executable, currently named `EllipseFromOutline`.
Below is an outline of the _original_ mesh generation process in each of the test cases that are provided.
As the codebase is developed, the tests themselves will adapt to take new input formats and adjust to code refactorings that take place, but ultimately must reproduce the following outputs from the given inputs at all times.

Testing is conducted using `pytest`, with python's `filecmp` handling the comparision between the reference outputs and those produced by the latest version of the code.

## `ref_outputs`

This folder contains the reference, or expected, outputs from the mesh generation executable.
These outputs were produced prior to modifications to the codebase (made after 2022-08-08), and future modifications to the codebase should be able to replicate the reference outputs, given the reference inputs.

## `mesh_inputs`

Contains the original inputs that were used to generate the reference outputs.

## The Regression Tests

There are three (regression) tests, for meshes of differing shapes:
- (small) Rectangle
- (small) Spherical triangulation
- (small) Wing-disc

Instructions for _the original_ production of each mesh file (given it's inputs) are as follows:

### Small Rectangle

- The `smallRectangle.ele` and `smallRectangle.node` input files have to be moved to the folder `TissueOrigami/ToolBox/MeshGeneration/2DEllipse` (although moving the executable to the same directory as the files will likely also work) and renamed to `Points.1.ele` and `Points.1.node` respectively.
- In `EllipseFromOutline.cpp`, manually change `selectTissueType=1` and then in the condition `else if(selectTissueType == 1)`, set the following variables:
```cpp
    actinHeight = 2.0;
    ECMHeight = 0.2;
```
- `cd` to the `2DEllipse` folder and run the following commands:
```bash
    g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
    cp ./src/EllipseFromOutline ./
    ./EllipseFromOutline -1 5.2 2 3 0
```

This will create the outout `MeshFile.out`. This is what has been renamed to `smallRectangle.mesh` and placed in the `ref_outputs` folder.

## Spherical Triangulation

- `SphericalTriangulation.txt` must be copied to the folder `TissueOrigami/ToolBox/MeshGeneration/2DEllipse` and renamed `SphericalTriangulation` (no file extension).
- In `EllipseFromOutline.cpp` manually change the parameter `selectTissueType=5`, and then in the condition `else if(selectTissueType == 5)` set
```cpp
    actinHeight = 1.0;
    ECMHeight = 0.2;
```
- `cd` to the `2DEllipse` folder and run the following commands:
```bash
g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
cp ./src/EllipseFromOutline ./
./EllipseFromOutline -2 4.2 2 3 0
```

This produces `MeshFile.out`, which has been renamed to `smallSphere.mesh` and placed in the `ref_outputs` folder.

## Wing disc

- In `EllipseFromOutline.cpp`, manually change `selectTissueType=1`, and in the condition `else if(selectTissueType == 1)` set the parameters
```cpp
    actinHeight = 2.0;
    ECMHeight = 0.2; 
```
- Ensure that folder paths have been manually changed to point to the correct directories.
- `cd` to `2DEllipse` and run the following commands:
```bash
g++ -std=c++11 -o ./src/EllipseFromOutline ./src/EllipseFromOutline.cpp
cp ./src/EllipseFromOutline ./
./EllipseFromOutline 1 35.56 50 27.3 27.3 12.5 9 2 0  ./inputOutlines/48hrDiscSymmetricOutline
```
This will triangulate the wing disc outline and create `Points.1.node` and `Points.1.ele`. Since these are technically outputs of the meshing process, a copy of these files under the names `WD_Points.1.node` and `WD_Points.1.ele` (respectively) is included in `ref_outputs`.
- Then run:
```bash
./EllipseFromOutline -1 5.2 1.0 3 1 
```

This should produce an output identical to `smallWingDisc.mesh`.