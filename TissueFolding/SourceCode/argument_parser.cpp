# include <stdexcept>
# include "argument_parser.h"
# include <iostream>
# include <sys/stat.h>

using namespace std;

/**
 * @brief Checks whether path points to a directory that is accessible
 * 
 * @param path Path to a directory on the filesystem
 * @return true The path points to a directory
 * @return false The path does not exist, nor point to a directory
 */
bool DirectoryExists(string path) {

    struct stat info;

    if (stat(path.c_str(), &info) != 0) {
        // directory does not exist
        return false;
    }
    else if (info.st_mode && S_IFDIR)
    {
        // path exists on the filesystem, but it's not pointing to a directory
        return false;
    }
    else {
        return true;
    }
}

/**
 * @brief Read arguments from the command line and place them into an ArgumentSpace instance.
 * 
 * The ArgumentSpace instance interprets the command-line inputs, and determines the mode the simulation will run in, the location of the input file (and directory, if appropriate), and the output directory.
 * These are then passed to the simulation instance that the main program creates.
 * 
 * @param argc 
 * @param argv 
 * @return ArgumentSpace Instance storing parsed input values and determining type of mesh generation scheme.
 */
ArgumentSpace ArgumentReader::readInput(int argc, char **argv) {

    auto args = ArgumentSpace(argc, argv);

    // if the user has requested help we should provide it here
    if (args.hasFlag("-h"))
    {
        printHelp(args.getSimulationMode());
        exit(0);
    }
    // if no simulation mode was specified, but the user did not ask for help, exit with error
    else if (args.getSimulationMode()==NotSet) {
        fprintf(stdout, "Error - no simulation mode specified. See the help below for usage information.\n");
        printHelp(NotSet);
        exit(-1);
    }
    // if we get to here then a run type was specified, and the user didn't ask for help
    // we thus check that the correct input arguments were given, and if so, process them into the members of args
    else if ((args.getSimulationMode()==OnTheGo) && args.nNonFlags==3 ) {
        /* We require an input file path and output file directory. Input directory is NOT required.
        We have also been blessed with the correct number of inputs */
        args.pathToInputFile = args.getNonFlag(1);
        args.pathToOutputDir = args.getNonFlag(2);
    }
    else if (!(args.getSimulationMode() == OnTheGo) && args.nNonFlags == 4)
    {
        /* We are running DisplaySave or Continue, but need the same number of inputs regardless, of which we have enough.
        We require an input file path, output file directory, and the input file directory to continue from */
        args.pathToInputFile = args.getNonFlag(1);
        args.pathToOutputDir = args.getNonFlag(2);
        args.pathToInputDir  = args.getNonFlag(3);

        if (!DirectoryExists(args.pathToInputDir))
        {
            cout << "Error; (input) directory not found at path: " + args.pathToInputDir + "\n";
            exit(-1);
        }
    }
    // if we got to here, the input was invalid (too few arguments and didn't ask for help)
    else {
        fprintf(stdout, "Error - incorrect number of inputs. "
                        "See function help below.\n");
        printHelp(args.getSimulationMode());
        exit(-1);
    }

    // check input file and output directory are where the user claims they are
    // these checks are common to all simulation modes
    ifstream inputFile(args.pathToInputFile.c_str());
    if (!inputFile.good())
    {
        cout << "Error; input file not found at path provided: " + args.pathToInputFile + "\n";
        exit(-1);
    }
    else if (!DirectoryExists(args.pathToOutputDir))
    {
        cout << "Error; (output) directory not found at path: " + args.pathToOutputDir + "\n";
        exit(-1);
    }
    // finally we can return args, which contains the information we need for the run to execute
    return args;
};

/**
 * @brief Prints the help message for the executable.
 *
 * Update this function as additional functionality is added to the ArgumentReader class
 */
void ArgumentReader::printHelp(SimMode mode)
{
    switch (mode)
    {
    case NotSet:
        fprintf(stdout, "Usage:\n"
                        "TissueFolding [options] mode \n"
                        "Required:\n"
                        "mode:\t Determines the simulation mode.\n"
                        "\t Recognised values are:\n"
                        "\t \t onthego: Start a new simulation from an input file\n"
                        "\t \t continue: Resume a simulation that has previously finished\n"
                        "\t \t display: Read data from the output of a finished simulation\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case OnTheGo:
        fprintf(stdout, "Usage:\n"
                        "TissueFolding onthego [options] input_file output_dir \n"
                        "Starts and runs a new simulation using the input_file, placing the results in output_dir\n"
                        "Required:\n"
                        "input_file: \t Path to the input file the simulation should read from\n"
                        "output_dir: \t Path to the directory to write simulation results to\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case Continue:
        fprintf(stdout, "Usage:\n"
                        "TissueFolding continue [options] input_file output_dir input_dir \n"
                        "Resumes a simulation by reading from the directory input_dir, and input file input_file \n"
                        "Required:\n"
                        "input_file: \t Path to the input file the simulation should read from\n"
                        "output_dir: \t Path to the directory to write simulation results to\n"
                        "input_dir : \t Path to the directory to read the starting state from\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case DisplaySave:
        fprintf(stdout, "Usage:\n"
                        "TissueFolding display [options] input_file output_dir input_dir \n"
                        "Reads simulation output from input_dir and displays the information \n"
                        "Required:\n"
                        "input_file: \t Path to the input file the simulation should read from\n"
                        "output_dir: \t Path to the directory to write simulation results to\n"
                        "input_dir : \t Path to the directory to read the starting state from\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    }
};

/**
 * @brief Construct a new Argument Space:: Argument Space object.
 * 
 * Reads in the flags provided to the executable, and sorts them into flags and non-flags.
 * 
 * @param argc Number of command-line arguments passed.
 * @param argv Strings of command-line arguments passed.
 */
ArgumentSpace::ArgumentSpace(int argc, char **argv) {

    allArgs = vector<string>(argv+1, argv+argc);
    nNonFlags = 0;
    nFlags = 0;
    for(const auto &a : allArgs) {
        if(!argIsFlag(a)) {
            nonFlags.push_back(a);
            nNonFlags++;
        }
        else {
            flags.push_back(a);
            nFlags++;
        }
    }
    nArgs = nNonFlags + nFlags;
    // attempt to identify the simulation mode. If it cannot be identified, leave it as 0 to force an exit later
    simulationMode = NotSet;
    if (nNonFlags >= 1)
    {
        if (nonFlags[0]=="onthego") {simulationMode = OnTheGo;}
        else if (nonFlags[0]=="continue") {simulationMode = Continue;}
        else if (nonFlags[0]=="displaysave") {simulationMode = DisplaySave;}
    }
};

/**
 * @brief Prints the specifications provided by the user to stdout, for this meshing process
 * 
 */
void ArgumentSpace::printModeSpecs() {
    cout << "Run type: " + nonFlags[0] + "\n";
    cout << "Variable values as follows:\n"
            "pathToInputFile = " + pathToInputFile + "\n"
            "pathToOutputDir = " + pathToOutputDir + "\n"
            "pathToInputDir  = " + pathToInputDir  + "\n";
}

/**
 * @brief Get the mode that the simulation will run
 * 
 * @return SimMode Mode of simulation
 */
SimMode ArgumentSpace::getSimulationMode() {
    return simulationMode;
}

/**
 * @brief Returns the i-th non-flag argument passed from the command line
 * 
 * @param i Index of argument passed from command line
 * @return string Value of required argument i passed from command line
 */
string ArgumentSpace::getNonFlag(int i) {
    return nonFlags[i];
}

/**
 * @brief Check whether a given command line input is a flag.
 *
 * Flags are identified by the first character being a single dash, '-'.
 *
 * @param[in] argument The command line argument passed
 * @return true The argument is a flag (by the above definition)
 * @return false The argument is not a flag
 */
bool ArgumentSpace::argIsFlag(string argument)
{
    return argument[0] == '-';
};

/**
 * @brief Checks whether a given flag was passed to the executable.
 * 
 * If any string in args matches the pattern flag, then that flag has been passed to the executable.
 * 
 * @param flag Flag pattern to search for
 * @return true Flag has been passed to executable
 * @return false Flag has not been passed to executable
 */
bool ArgumentSpace::hasFlag(string const &flag) {
    for(const auto &a : allArgs) {
        if (a == flag) {
            return true;
        }
    }
    return false;
};