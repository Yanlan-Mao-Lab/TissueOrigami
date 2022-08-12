# TissueFolding

Introduction to the purpose of the executable

## Build Instructions

Instructions for building the executable on a user's machine

## Usage

`TissueFolding` is designed to be called from the command-line.

#### Command-Line Functionality

The `TissueFolding` executable, supports command-line help by passing the `-h` flag.

Recognised patterns are:
`TissueFolding -d`
Runs the TissueFolding using the default parameters. Passing this option will cause the executable to ignore all other inputs passed, including paths to the desired input and output directories.
`TissueFolding onthego input_file output_dir`
Starts a fresh TissueFolding simulation, setting up the initial state from the `input_file` and placing the results in the `output_dir`
`TissueFolding continue input_file output_dir input_dir`
Resumes a previously finished simulation by reading its initial input file from `input_file`, and final state from the `input_dir`. Writes the results to the `output_dir`.
`TissueFolding displaysave input_file output_dir input_dir`
Displays the results of a finished simulation by reading its initial input file from `input_file`, and final state from the `input_dir`. Soley intended for use with the GUI. Results can optionally be written again to `output_dir`.

#### Preparing the `input_file`

Information on how to setup the input should go here

#### Output

Information on what each of the output files are should go here