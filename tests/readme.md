# `tests`

This directory contains all regression tests for the TissueOrigami pipeline.
These regression tests cover:
- Building and running the mesh generation executable `EllipseFromOutline`, tests are found in the folder `mesh_gen_reg`
- Building and running the non-GUI `TissueFolding` executable, tests are found in the folder `sim_no_pardiso_reg`

For more information, see the `readme`s in the relevant folders:
- `mesh_gen_reg`: Regression tests concerning mesh generation
- `sim_no_pardiso_reg`: Regression tests concerning `TissueFolding`, built without Pardiso.

## `simOutputsToTxt`

It is sometimes necessary for us to translate the (binary file) outputs of `TissueFolding` into text files, usually for testing whether the data _read into the simulation_ from these files matches the reference files.

This executable can be compiled via
```bash
g++ simOutputsToTxt.cpp -o ./simOutputsToTxt.o
```
**NOTE:** The regression tests assume that `simOutputsToTxt.o` is located in the `tests/` directory, and automated testing compiles this code and places the executable at this location.
Each test folder will explicity mention whether it depends on this executable being present.