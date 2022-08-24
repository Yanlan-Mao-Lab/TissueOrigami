# `TissueFolding`

This executable conducts the simulation (add biological context - expertise required)

## Requirements

`TissueFolding` has the following dependencies:

- [OpenMP](https://www.openmp.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [OpenBLAS and LAPACK](https://www.openblas.net/)
- [boost](https://www.boost.org/)

These libraries must be included in your `C++` compiler's include path.

Additionally, `TissueFolding` relies on [Pardiso](https://www.pardiso-project.org/) to solve the evolutionary equations. Please visit the Pardiso website to obtain the Pardiso library and a license key, and the appropriate location on your machine to save these. Then follow the "Building with Pardiso" instructions below. `TissueFolding` can be built without Pardiso - in such a case, functionality is limited to visualising the results of a previously-completed simulation, and validating input files.

## Building without Pardiso

To build `TissueFolding` without Pardiso, navigate to the `TissueFolding/SourceCode` directory and run the following commands:
```bash
mkdir build; cd build;
cmake ..
cmake --build .
```
Then run
```bash
mv ./TissueFolding ../../
```
to move the `TissueFolding` executable out of the build directory and into the `TissueFolding` directory.