# include "SimulationArgumentReader.h"
# include <stdexcept>
# include <fstream>
using namespace std;

SimulationArgumentReader::SimulationArgumentReader(int argc, char **argv) {
    // there should be an odd number of arguments passed in (each flag must be followed by its value, plus the executable name)
    if (argc % 2 == 0) {
        printHelp(true);
    }
    // cycle through the input arguments and update the corresponding attributes
    int i=1;
    while (i<argc) {
        string argument = argv[i];
        if (argument=="-mode") {
            // always safe to check next argument: argc is odd and i is odd and LESS than argc, so i+1 is a valid index
            string mode_type = argv[i+1];
            // set the simulation run type - for the non-GUI we DON'T allow DisplaySave
            if (mode_type=="SimulationOnTheGo" || mode_type=="ContinueFromSave") { mode=mode_type; }
            else if (mode_type=="DisplaySave") {
                // we don't permit DisplaySave runs when just using the executable,
                // but provide a hint that the UI can be used instead
                fprintf(stdout, "Error: DisplaySave is not accepted as a mode for the non-UI executable.\n"
                "Did you mean to run the GUI executable instead?\n");
                printHelp(true);
            }
            else { throw runtime_error("Improper usage of executable."); }
        }
        else if (argument=="-i") {
            input_file = argv[i+1];     // set the input file path
        }
        else if (argument=="-od") {
            output_dir = argv[i+1];     // set the output directory
        }
        else if (argument=="-dInput") {
            read_in_dir = argv[i+1];    // set the directory to read in from
        }
        else {
            // unrecognised argument
            throw runtime_error("Invalid use of executable");
        }
        i += 2; // incriment 2, to the next "flag"
    }
}

void SimulationArgumentReader::printHelp(bool bad_input) {
    if (bad_input) {
        fprintf(stdout, "Error - invalid usage of TissueFolding executable. See help below.\n\n");
    }
    fprintf(stdout, "Runs the TissueFolding simulation (without GUI) from the command line.\n"
                    "Usage:\n"
                    "TissueFolding [-i inputFile] [-mode SimMode] [-od outputDirectory=.] [-dInput readInDirectory] \n"
                    "Required:\n"
                    "\t inputFile:\t Input file to read from. See the readme for format specifications.\n"
                    "\t mode:\t Either SimulationOnTheGo or ContinueFromSave. See the readme for run mode specifications.\n"
                    "\t outputDirectory:\t [Default .] Directory to write output files to.\n"
                    "Conditional Requirements:\n"
                    "\t readInDirectory:\t If ContinueFromSave is selected for the mode, this must be the path to the directory containing the output of the simulation from which to continue.\n");
    if (bad_input) {
        exit(-1);
    }
    else {
        exit(0);
    }
}