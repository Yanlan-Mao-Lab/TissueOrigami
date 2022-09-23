# `TissueFolding`

This executable conducts the simulation (add biological context - expertise required)

## Requirements

`TissueFolding` has the following dependencies:

- [OpenMP](https://www.openmp.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [OpenBLAS and LAPACK](https://www.openblas.net/)
- [boost](https://www.boost.org/)

These libraries must be included in your `C++` compiler's include path.

Additionally, `TissueFolding` relies on [PARDISO](https://www.pardiso-project.org/) to solve the evolutionary equations. Please visit the Pardiso website to obtain the Pardiso library and a license key, and the appropriate location on your machine to save these. Then follow the "Building with Pardiso" instructions below. `TissueFolding` can be built without Pardiso - in such a case, functionality is limited to visualising the results of a previously-completed simulation, and validating input files.

## Building the Executable

The `TissueFolding` executable is built via `cmake`; by default, it will assume you want to build the executable _without_ the PARDISO library.
If you wish to build the executable with PARDISO, please follow the setup instructions in the relevent section.

#### Building on Mac

MacOS' default compiler, `clang`, does not support compiling with/against `OpenMP`, which is a dependency of `TissueFolding`. To install on MacOS, we recommend that you use `HomeBrew` to install `g++`:
```bash
brew install gcc
```
Then pass the `CMAKE_CXX_COMPILER` flag to `cmake` when configuring the build - that is, add the flag
```bash
-DCMAKE_CXX_COMPILER=/path/to/g++
```
to the call to `cmake ..` in the instructions below.
You can execute the bash command `which g++` in a terminal to print the `/path/to/g++` that `HomeBrew` has installed `g++` to, or simply pass in
```bash
-DCMAKE_CXX_COMPILER=$(which g++)
```
if you are certain that the version of `g++` you wish to use to compile is on your `${PATH}`.

### Building without PARDISO

Building `TissueFolding` without PARDISO is the same procedure across all platforms.
Navigate to the `TissueFolding/SourceCode` directory and run the following commands:
```bash
mkdir build; cd build
cmake ..
cmake --build .
```

This will create the `TissueFolding` executable and place it into the `build/` directory.
For testing purposes, we recommend you move the executable out of the build directory. 
The automated tests assume that the `TissueFolding` executable is located in the `/TissueOrigami/TissueFolding/` directory, so we recommend you move it here.

Once the executable has been moved, you can safely delete the `build` directory.

### Building with PARDISO

Building `TissueFolding` with PARDISO requires some additional setup to accomodate the PARDISO library.
Additional dependencies are:
- `pthread`,
- `gfortran`,
- Ensure that you have the math library `<math.h>` on your compiler's include path.

For the specific PARDISO instructions, please see [the manual available on the PARDISO website](https://pardiso-project.org/manual/manual.pdf).
The essential steps to carry out are:
- Ensure that you have downloaded the shared (`.so`) PARDISO library file. Make a note of the path to the directory it is saved in, and a note of the _exact_ filename (including the `.so` extension).
- Ensure that your license key is either in your home directory, or in the same directory as the working directory you wish to run `TissueFolding` from.

The build instructions then only differ from the build without PARDISO in the options that must be passed to the call to `cmake`:
```bash
cd TissueOrigami/TissueFolding/SourceCode
mkdir build; cd build
cmake .. -DPARDISO=ON -DPLOC=pardiso_path -DPNAME=pardiso_so_name
cmake --build .
```
Here, `-DPARDISO` is the option that informs `cmake` that you wish to build using PARDISO. _This must be set to `ON`._
The remaining flags `-DPLOC` and `DPNAME` are optional, although in practice the user will have to specfiy them in order to aid `cmake` in locating the PARDISO library.
`pardiso_path` should be set to the path at which the PARDISO `.so` file can be found (either a complete path to this file, or to the directory containing this file).
`pardiso_so_name` should be set to the full name (including the lib prefix if present, and the `.so` extension) of the PARDISO library you wish to build against.
If you only supply one (or neither) of these parameters, `cmake` will attempt to locate PARDISO on your system by itself - but this attempt is likely to be futile due.
`cmake` will let you know whether it was able to find PARDISO, and if so _where_ it found PARDISO so you can check that the correct library has been found.

Once the build completes, you should move the `TissueFolding` executable out of the `build/` directory, as detailed in the "build without PARDISO" section.

## Usage

`TissueFolding` takes in a maximum of 4 command-line arguments, indicated by flags:
```bash
./TissueFolding [-i input_file] [REQUIRED -mode sim_mode] [-od output_dir] [-dInput prev_save_dir]
```
The arguments can be passed in any order, however the value of the argument must be passed immediately after the corresponding flag, and separated by a single whitespace.
- `-i input_file` : Path to the input file for this simulation. See "input files to the simulation", below.
- `-mode sim_mode`: One of `SimulationOnTheGo`, `ContinueFromSave`, or `DisplaySave`
    - `SimulationOnTheGo`   : Starts a new simulation from scratch.
    - `ContinueFromSave`    : Reads in the outputs of a finished simulation from the directory specified by `-dInput`, and resumes the simulation from this endpoint.
    - `DisplaySave`         : Reads in the outputs of a finished simulation from the directory specified by `-dInput`, for the Qt-GUI visualisation. **NOTE:** Surely this should only be used when compiling with the GUI?!?!?!!?!
- `od output_dir`           : The directory to write the outputs to. Will overwrite any previous outputs in this directory.
- `dInput prev_save_dir`    : Specifies a directory containing the output of a finished simulation, which will be read in by the program. Note required if `SimulationOnTheGo` is set with `-mode`.

#### NOTE:

It should not be necessary to `export LD_LIBRARY_PATH=path/to/pardiso.so` when running the compiled executable - more research incoming.
It is likely still necessary to export `OMP_NUM_THREADS` however before running.


### Input files to the simulation

(Domain knowledge explaination of the contents of the input files, and how much is optional/not optional, whether an order is expected, etc.)

### Outputs

(Domain knowledge explaination of the names and contents of the output files)