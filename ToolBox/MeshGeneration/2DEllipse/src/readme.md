# Mesh Generation via the `EllipseFromOutline` executable

(flavour text here)

## Build Instructions for `EllipseFromOutline`

### [`triangle`]((https://www.cs.cmu.edu/~quake/triangle.html)) Dependency

[`triangle`](https://www.cs.cmu.edu/~quake/triangle.html) is a dependency of the meshing process.
The mesh generation executable, `EllipseFromOutline`, needs to be able to call the `triangle` executable during execution in order to process and write the meshfiles.
To do so, the path to the `triangle` executable must be included in the system `PATH` variable.
If `triangle` is not included in `PATH`, one can edit the `find_program` command in `CMakeLists.txt` to include a path to the `triangle` executable:
```
find_program(TRIANGLE_PATH
            triangle triangle.o
            "manual_path_to_triangle")
```

Bear in mind however, that `EllipseFromOutline` will continue to look in the provided path for `triangle` even if the location of the executable is changed.
Recompiling/rebuilding `EllipseFromOutline` will be neccessary - along with providing the new path to `triangle` - in such circumstances.

### Build via CMake

After installing `triangle`, change directory into the `Toolbox/MeshGeneration/2DEllipse/src` directory, and run the following commands:
```
mkdir build; cd build
cmake ../src/
cmake --build .
```
This will compile the `EllipseFromOutline` binary and place it into the `build` directory.
You can then move it to the desired location on your machine.

## Usage

`EllipseFromOutline` is designed to generate meshes for a number of different tissue samples, of various shapes, and in both two dimensions (for the planar layers of the tissue) and three dimensions (for linking the two dimensional meshes for each layer together). To accomodate the variance in the required inputs in each of these cases, the executable reads from user specifications contained in an input file, and writes the resulting mesh to an output file. For details on the syntax of the input file, see below.

Command-line usage is as follows:
```bash
EllipseFromOutline [options] inputFile [outputFile]
```
Required inputs:
- `inputFile`: Path to the input file
Optional inputs:
- `outputFile`: Path to write the output mesh to. Defaults to `MeshFile.out` if not provided.
Options:
- `-h, --help`: Prints the command-line help message.

### Input File Syntax

Input files follow a `.yml` style syntax of variable names, followed by colons, then by the value assigned to the variable, separated by newlines.
Certain variable values may require further variables to be provided, which are indicated by indentation by _two character spaces_.
The general syntax for the input files is thus as follows:
```yml
input 1: value
  input 1 dependency 1: value
  input 1 dependency 2: value
input 2: value
  input 2 dependency 1: value
input 3: value
```

The first variable name should be `meshing_mode`, and should take one of the following values:
- `rec`: The tissue to be meshed is of rectangular shape
- `wgd`: The tissue to be meshed is of wing-disc shape - this includes elliptical (and circular) shapes
- `2d` : The tissue is to be meshed in 2D using a pre-built tesselation
- `3d` : The tissue is to be meshed in 3D using a pre-built tesselation
- `3d_cyl` : The tissue is to be meshed in 3D using a pre-built cylindrical tesselation
In the event that the value `wdg` is supplied, four dependencies must be provided (need diagram and details to be filled in by someone with contextual knowledge!)
- `length1`:
- `length2`:
- `width1` :
- `width2` :
- `outline`: Path to an outline file that traces out the boundary of the tissue. For more information see the "shape outline files" section.

Following the `meshing_mode`, the following arguments are required (in order):
- `ABHeight`: The combined height (length in the $z$-direction) of the actin and basal tissue layers. When ECM is present this is technically the combined height of the actin, ECM, and basal tissue layers.
- `PrismSideLength`: The target side length for the (prism shaped) finite volume elements.
- `nzLayers`: The number of tissue layers in the $z$-direction.
- `symY`: Toggles whether the mesh should be made symmetric in y.

Finally, one should specify a `tissueType`:
- `tissueType`: Determines the type of tissue, and thus the material properties etc, that is to be meshed.
Options for the value of this variable are:
- 0: Wing disc after 48 hours (no ECM)
- 1: Wing disc after 48 hours (ECM)
- 2: Wing disc after 72 hours
- 3: Optic cup
- 4: Half-disc / quarter circle
- 5: Spherical organoid
- 6: Tubular organoid
- 7: Rectangle (ECM)
- 8: Rectangle (no ECM)
**NOTE**: There is cause to revise the inputs that are provided, especially since the tissue types provided are only then ultimately used to set a number of parameters to hard-coded values.

### Shape Outline Files

Not sure what these are (or rather, where they come from - MATLAB?). Need contextual input.