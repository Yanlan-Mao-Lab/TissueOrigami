# `VisualiseTissueFolding`

This executable visualises the simulation results (add biological context - expertise required)

## Requirements

`VisualiseTissueFolding` has the following dependencies:

- [OpenMP](https://www.openmp.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [boost](https://www.boost.org/)
- [Qt](https://www.qt.io/), in particular its [OpenGL and Widget components](https://doc.qt.io/qt-6/qtopengl-index.html#qt-opengl-and-qt-widgets). If installation difficulties are encountered, there is an alternative [command line installation tool](https://github.com/miurahr/aqtinstall) (which we use in the CI).
- [OpenGL](https://github.com/miurahr/aqtinstall)

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
./VisualiseTissueFolding [-dInput input_folder]
```
`input_folder` is a required argument, and must point to a directory in which the output of a completed simulation is saved.

For example, you could
- Open the terminal
- cd to the git folder
- ./UserInterface/Sourcecode/build/VisualiseTissueFolding -mode DisplaySave -dInput ./Samples/DisplaySave/ > ./Samples/DisplaySave/tmp
to display the example provided in this repository.