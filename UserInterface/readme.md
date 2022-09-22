# `VisualiseTissueFolding`

This executable visualises the simulation results (add biological context - expertise required)

## Requirements

`VisualiseTissueFolding` has the following dependencies:

- [OpenMP](https://www.openmp.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [boost](https://www.boost.org/)
- [Qt with OpenGL]()

These libraries must be included in your `C++` compiler's include path.

## Building the Executable

The `VisualiseTissueFolding` executable is built via `cmake`;

Navigate to the `UserInterface/SourceCode` directory and run the following commands:
```bash
mkdir build; cd build
cmake ..
cmake --build .
```

This will create the `VisualiseTissueFolding` executable and place it into the `build/` directory.

## Usage

```bash
./VisualiseTissueFolding [-dInput input_folder] [-mode sim_mode]
```
For example, you could
- Open the terminal
- cd to the git folder
- ./UserInterface/Sourcecode/build/VisualiseTissueFolding -mode DisplaySave -dInput ./Samples/DisplaySave/ > ./Samples/DisplaySave/tmp
to display the example provided in this repository.