# Regression Tests for mesh generation executable

The `mesh_gen_reg` folder contains the scripts and files for running regression tests on the mesh-generation executable.
Below is an outline of the _original_ mesh generation process in each of the test cases that are provided.
As the codebase is developed, the tests themselves will adapt to take new input formats and adjust to code refactorings that take place, but ultimately must reproduce the following outputs from the given inputs at all times.

Testing is conducted using `pytest`, with python's `filecmp` handling the comparision between the reference outputs and those produced by the latest version of the code.

## The Regression Tests

Test cases are separated into folders containing the resources and reference files for one or more tests:
- `smallRectangle`, containing 1 test for meshing a rectangular tissue in 2d from a pre-built mesh
- `smallSphere`, containing 1 test for meshing a spherical tissue in 3d from a pre-built mesh
- `smallWingDisc`, containing 2 tests; the first for generating a tesselation from an outline file and via `triangle`, the second to generate the corresponding 2d mesh.

For each test, the following steps are run:
1. The necessary input files are copied to the mesh generation executable's directory
2. The executable is run
3. The desired outputs are moved back to the testing folder
4. These outputs are compared to the reference outputs
The test fails if an output file differs from its corresponding reference file.

### Adding additional reference tests

The file `test_MeshGeneration.py` is structured in such a way as to seemlessly allow for the introduction of additional regression tests.
To add a new test:
- Create a folder `new_test` in the `mesh_gen_reg` directory
- Populate `new_test` with the input files that the mesh generation executable requires, and the reference output files
- Create a new entry (`my_test`) in the `test_cases` attribute of the `Test_MeshGeneration` class in `test_MeshGeneration.py`
- Create a new instance of a `MeshGenRun` object (`my_test_obj = MeshGenRun(...)`) corresponding to the test (see the "Creating a test via `MeshGenRun`, below)
- Assign `test_info[my_test] = my_test_obj`.
The next time pytest runs, your test will also be included.
For concrete examples of this setup process, one can examine the `Test_MeshGeneration` class, and the relevant docstrings in `test_MeshGeneration.py`.

#### Creating a test via `MeshGenRun`

The `MeshGenRun` class in `test_MeshGeneration.py` provides a framework for processing all the information needed to run the mesh generation executable with a given set of inputs, extract (possibly a subset of) the outputs, and compare them to some given reference files.
Let us assume you have created a new folder for your test, named `my_test`, and have populated it with the following files:
- `cmd_line_input` : The input file that the mesh generation executable takes from the command line
- `extra_input_1, extra_input_2, ..., extra_input_N` : Additional input files the mesh generation executable needs for this run, for example outline files, pre-built meshes or tesselations, etc.
- `reference_1, reference_2, ..., reference_M` : The reference outputs - whatever the mesh generation executable produces given the input files above should match these files.
And let us assume that you are interested in outputs `out_1, out_2, ..., out_M` from the executable - note that you can have fewer reference files than output files, if for example two outputs should be equal.
However, we will assume a 1-to-1 mapping between reference files and outputs, so we want to make sure `out_j` is the same as `reference_j` for `j=1,...M`.

Next, open the `test_MeshGeneration.py` file and navigate to the `Test_MeshGeneration` class.
Choose a unique identifier for the test, `my_test_ID`; current convention is to match this ID with the name of the folder containing the relevant files, IE to set `my_test_ID = my_test`.
Append this ID to the `test_cases` variable;
```python
test_cases = [ "existing", "test", "cases", "my_test_ID" ]
```

Still within the `Test_MeshGeneration` class, create a new instance of a `MeshGenRun` object via
```python
myMeshGenObj = MeshGenRun("cmd_line_input", "my_test")
```
This informs pytest that for your test; the mesh generation executable takes `cmd_line_input` from the command line, and the folder containing the reference and input files is `my_test`.
Next, inform pytest about the additional input files required for your test, using
```python
myMeshGenObj.addExeInput("extra_input_i", "optional_new_name_i")
```
for each `i=1,...,N`.

The `optional_new_name_i` argument is optional; if provided, pytest will copy _and rename_ `extra_input_i` to `optional_new_name_i` before running mesh generation.
This is useful for giving your input files understandable and clear names, but for when the mesh generation executable is looking for a specific file name during runtime.

Finally, inform pytest which outputs should be recovered, and to which reference files they should be compared, using
```python
myMeshGenObj.addComparisonFile("out_j", "reference_j", rename_to="intermediate_name_j", trim=True/False)
```
for each `j=1,...,M`.

The `trim` parameter defaults to `False`, however should be set to `True` in the event that `out_j` and `reference_j` are produced from `triangle`. 
In such cases, triangle appends the command passed to `triangle.o` the end of the output file, which involves absolute paths and does not contain any useful meshing information, whilst simultaneously causing file-comparison to fail.

The `rename_to` argument is also optional; if provided, `out_j` will be copied back to the testing area with the name `intermediate_name_j` rather than `out_j`. This is particularly useful if `reference_j` is saved under the same name as `out_j`, as it prevents file overwriting and an error occuring.

## [Legacy - Running the Executable] Executable Runs

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

**This is also where the 1st part of the WingDisc regression tests ends, and part 2 begins.**

Then run:
```bash
./EllipseFromOutline -1 5.2 1.0 3 1 
```
This should produce an output identical to `smallWingDisc.mesh`.