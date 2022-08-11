# include <stdexcept>
# include "argument_parser.hpp"
# include <iostream>

using namespace std;

/**
 * @brief Read arguments from the command line and place them into an ArgumentSpace instance.
 * 
 * The ArgumentSpace instance stores the input values, and determines the nature of, the mesh generation process.
 * 
 * @param argc 
 * @param argv 
 * @return ArgumentSpace Instance storing parsed input values and determining type of mesh generation scheme.
 */
ArgumentSpace ArgumentReader::read_input(int argc, char **argv) {

    auto args = ArgumentSpace(argc, argv);

    // if the user has requested help we should provide it here
    if (args.has_flag("-h"))
    {
        print_help(args.get_run_type());
        exit(0);
    }
    // if no run type was specified, but the user did not ask for help, exit with error
    else if (args.get_run_type()==NO_RUN_SET) {
        fprintf(stdout, "Error - no run type specified. See the help below for usage information.\n");
        print_help(NO_RUN_SET);
        exit(-1);
    }
    // if we get to here then a run type was specified, and the user didn't ask for help
    // we thus check that the correct input arguments were given, and if so, process them into the members of args
    switch (args.get_run_type())
    {
    case WDC_RUN:
        /* We require l1, l2, w1, w2, abh, slp, nlz, path_to_input to be present in the input */
        if (args.num_non_flags == 9) {
            // there are the correct number of input arguments, place them into the appropriate variables
            args.length1 = stof(args.get_non_flag(1));
            args.length2 = stof(args.get_non_flag(2));
            args.width1 = stof(args.get_non_flag(3));
            args.width2 = stof(args.get_non_flag(4));
            args.ABHeight = stof(args.get_non_flag(5));
            args.prism_side_len = stof(args.get_non_flag(6));
            args.n_layers_z = stof(args.get_non_flag(7));
            args.path_to_input = args.get_non_flag(8); // path to file needs to stay of type string
            // check input file is where the user claims it is
            ifstream inputOutline(args.path_to_input.c_str());
            if (!inputOutline.good()) {
                cout << "Error; input file not found at path provided: " + args.path_to_input + "\n";
                exit(-1);
            }
        }
        else {
            fprintf(stdout, "Error - incorrect number of inputs. "
                            "See function help below.\n");
            print_help(WDC_RUN);
            exit(-1);
        }
        break;
    case REC_RUN:
        /* We require abh, slp, nlz */
        if (args.num_non_flags == 4) {
            args.ABHeight = stof(args.get_non_flag(1));
            args.prism_side_len = stof(args.get_non_flag(2));
            args.n_layers_z = stof(args.get_non_flag(3));
        }
        else {
            fprintf(stdout, "Error - incorrect number of inputs. "
                            "See function help below.\n");
            print_help(REC_RUN);
            exit(-1);
        }
        break;
    case TESS_2D:
        /* We require abh, slp, nlz */
        if (args.num_non_flags == 4) {
            args.ABHeight = stof(args.get_non_flag(1));
            args.prism_side_len = stof(args.get_non_flag(2));
            args.n_layers_z = stof(args.get_non_flag(3));
        }
        else {
            fprintf(stdout, "Error - incorrect number of inputs. "
                            "See function help below.\n");
            print_help(TESS_2D);
            exit(-1);
        }
        break;
    case TESS_3D:
        /* We require abh, slp, nlz */
        if (args.num_non_flags == 4) {
            args.ABHeight = stof(args.get_non_flag(1));
            args.prism_side_len = stof(args.get_non_flag(2));
            args.n_layers_z = stof(args.get_non_flag(3));
        }
        else {
            fprintf(stdout, "Error - incorrect number of inputs. "
                            "See function help below.\n");
            print_help(TESS_3D);
            exit(-1);
        }
        break;
    }

    // also, if -sym is flagged, toggle this to true
    if (args.has_flag("-sym")) {args.sym_y = 1.;}

    // finally we can return args, which contains the information we need for the run to execute
    return args;
};

/**
 * @brief Prints the help message for the executable.
 *
 * Update this function as additional functionality is added to the ArgumentReader class
 */
void ArgumentReader::print_help(int mode)
{
    switch (mode)
    {
    case NO_RUN_SET:
        fprintf(stdout, "Usage:\n"
                        "EllipseFromOutline [options] type \n"
                        "Required:\n"
                        "type:\t Determines the mode of execuation and mesh generation.\n"
                        "\t Recognised values are:\n"
                        "\t \t wdc: Outline of tissue to be meshed is shaped like a wing disc or cicle\n"
                        "\t \t rec: Outline of tissue to be meshed is shaped like a rectangle\n"
                        "\t \t 2d : Providing a pre-built 2D tesselation file\n"
                        "\t \t 3d : Providing a pre-built 3D tesselation file\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case WDC_RUN:
        fprintf(stdout, "Usage:\n"
                        "EllipseFromOutline wdc [options] [-sym] l1 l2 w1 w2 abh slp nlz path_to_input \n"
                        "Generates a mesh from an input file containing a wingdisc or circular outline.\n"
                        "Required:\n"
                        "l1, l2, w1, w2: \t Respectively length 1, length 2, width 1 and width 2 of the input outline\n"
                        "abh: \t ABHeight\n"
                        "slp: \t Target side length of prisms\n"
                        "nlz: \t Number of layers in z\n"
                        "path_to_input: \t Path to input file containing tissue outline\n"
                        "Flags:\n"
                        "-sym: \t Flags symmetry in y to be true\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case REC_RUN:
        fprintf(stdout, "Usage:\n"
                        "EllipseFromOutline rec [options] [-sym] abh slp nlz\n"
                        "Generates a rectangular mesh.\n"
                        "Required:\n"
                        "abh: \t ABHeight\n"
                        "slp: \t Target side length of prisms\n"
                        "nlz: \t Number of layers in z\n"
                        "Flags:\n"
                        "-sym: \t Flags symmetry in y to be true\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case TESS_2D:
        fprintf(stdout, "Usage:\n"
                        "EllipseFromOutline rec [options] [-sym] abh slp nlz\n"
                        "Generates a 2D mesh from a pre-existing tesselation.\n"
                        "Tesselation is ASSUMED to be in the fileS Points.1.node and Points.1.ele\n"
                        "Required:\n"
                        "abh: \t ABHeight\n"
                        "slp: \t Target side length of prisms\n"
                        "nlz: \t Number of layers in z\n"
                        "Flags:\n"
                        "-sym: \t Flags symmetry in y to be true\n"
                        "Options:\n"
                        "-h:\t Display this help message\n");
        break;
    case TESS_3D:
        fprintf(stdout, "Usage:\n"
                        "EllipseFromOutline rec [options] [-sym] abh slp nlz\n"
                        "Generates a 3D mesh from a pre-existing tesselation.\n"
                        "Tesselation is ASSUMED to be in the file SphericalMeshTesselation.\n"
                        "Required:\n"
                        "abh: \t ABHeight\n"
                        "slp: \t Target side length of prisms\n"
                        "nlz: \t Number of layers in z\n"
                        "Flags:\n"
                        "-sym: \t Flags symmetry in y to be true\n"
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

    all_args = vector<string>(argv+1, argv+argc);
    num_non_flags = 0;
    num_flags = 0;
    for(const auto &a : all_args) {
        if(!arg_is_flag(a)) {
            non_flags.push_back(a);
            num_non_flags++;
        }
        else {
            flags.push_back(a);
            num_flags++;
        }
    }
    num_args = num_non_flags + num_flags;
    // attempt to identify the run type. If it cannot be identified, leave it as 0 to force an exit later
    run_type = NO_RUN_SET;
    if (num_non_flags >= 1)
    {
        if (non_flags[0]=="wdc") {run_type = WDC_RUN;}
        else if (non_flags[0]=="rec") {run_type = REC_RUN;}
        else if (non_flags[0]=="2d") {run_type = TESS_2D;}
        else if (non_flags[0]=="3d") {run_type = TESS_3D;}
    }
};

/**
 * @brief Prints the specifications provided by the user to stdout, for this meshing process
 * 
 */
void ArgumentSpace::print_mode_specs() {
    cout << "Run type: " + non_flags[0] + "\n";
    cout << "Variable values as follows:\n"
            "length1 = " + to_string(length1) + "\n"
            "length2 = " + to_string(length2) + "\n"
            "width1 = " + to_string(width1) + "\n"
            "width2 = " + to_string(width2) + "\n"
            "ABHeight = " + to_string(ABHeight) + "\n"
            "prism_side_len = " + to_string(prism_side_len) + "\n"
            "n_layers_z = " + to_string(n_layers_z) + "\n"
            "sym_y = " + to_string(sym_y) + "\n"
            "path_to_input = " + path_to_input + "\n";
}

/**
 * @brief Get the run_type for mesh generation
 * 
 * @return int Mode of mesh generation
 */
int ArgumentSpace::get_run_type() {
    return run_type;
}

/**
 * @brief Returns the i-th non-flag argument passed from the command line
 * 
 * @param i Index of argument passed from command line
 * @return string Value of required argument i passed from command line
 */
string ArgumentSpace::get_non_flag(int i) {
    return non_flags[i];
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
bool ArgumentSpace::arg_is_flag(string argument)
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
bool ArgumentSpace::has_flag(string const &flag) {
    for(const auto &a : all_args) {
        if (a == flag) {
            return true;
        }
    }
    return false;
};