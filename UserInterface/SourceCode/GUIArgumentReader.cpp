# include "GUIArgumentReader.h"
# include <stdexcept>
# include <fstream>
using namespace std;

GUIArgumentReader::GUIArgumentReader(int argc, char **argv) {
    // there should be an odd number of arguments passed in (each flag must be followed by its value, plus the executable name)
    if (argc % 2 == 0) {
        printHelp(true);
    }
    // cycle through the input arguments and update the corresponding attributes
    int i=1;
    while (i<argc) {
        string argument = argv[i];
        if (argument=="-mode") {
            fprintf(stdout, "Error: VisualiseTissueFolding is a visualisation tool only.\n"
            "Did you mean to run the simulation executable (TissueFolding) instead?\n");
            printHelp(true);
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
            printHelp(true);
        }
        i += 2; // increment 2, to the next "flag"
    }
}

void GUIArgumentReader::printHelp(bool bad_input) {
    if (bad_input) {
        fprintf(stdout, "Error - invalid usage of VisualiseTissueFolding executable. See help below.\n\n");
    }
    fprintf(stdout, "Visual display executable for the TissueFolding simulation.\n"
                    "Usage:\n"
                    "VisualiseTissueFolding [-i inputFile] [-dInput readInDirectory] [-od outputDirectory=.] \n"
                    "Required:\n"
                    "\t inputFile:\t Input file to read from. See the readme for format specifications.\n"
                    "\t readInDirectory:\t The directory containing the output of the simulation to visualise.\n"
                    "Optional:\n"
                    "\t outputDirectory:\t [Default .] Directory to write any output files to.\n");
    if (bad_input) {
        exit(-1);
    }
    else {
        exit(0);
    }
}