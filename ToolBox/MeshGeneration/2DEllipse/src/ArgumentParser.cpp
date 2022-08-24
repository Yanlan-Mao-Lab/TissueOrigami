# include <stdexcept>
# include "ArgumentParser.hpp"
# include <iostream>
# include <fstream>

using namespace std;

explicit ArgumentReader::ArgumentReader(int argc, char **argv) {
    // print command-line help if it was requested
    all_args = vector<string>(argv + 1, argv + argc);
    num_non_flags = 0;
    num_flags = 0;
    for (const auto &a : all_args)
    {
        if (!arg_is_flag(a))
        {
            non_flags.push_back(a);
            num_non_flags++;
        }
        else
        {
            flags.push_back(a);
            num_flags++;
        }
    }
    num_args = num_non_flags + num_flags;
    validate_cmdline_inputs();
    return;
}

void ArgumentReader::validate_cmdline_inputs() {
    // if the -h or --help flags were passed, return the inline help
    if (has_flag("-h") || has_flag("--help")) {
        print_help(DEFAULT);
        exit(0);
    }
    // check that the user didn't provide too many input arguments
    else if (num_non_flags>2) {
        print_help(TOO_MANY_INPUTS);
        exit(1);
    }
    // validate inputs to executable
    else {
        // is the input file visible?
        ifstream inputFile;
        // non_flags[0] should contain the input file location
        inputFile.open(non_flags[0], ifstream::in);
        if (!inputFile.good()) {
            throw runtime_error("Error: could not locate input file " + non_flags[0]);
        }
        path_to_input = non_flags[0];
        // if provided, non_flags[1] should contain the output file location
        if (num_non_flags>1) {
            ofstream outputFile;
            outputFile.open(non_flags[1]);
            if (!outputFile.good()) {
                throw runtime_error("Error: could not write to output " + non_flags[1]);
            }
            path_to_output = non_flags[1];
        }
        // if it was not provided, use the default output name
        else {
            path_to_output = "MeshFile.out";
        }
    }
}

void ArgumentReader::print_help(helpOpts mode)
{
    switch (mode) {
        case DEFAULT:
            fprintf(stdout, "Usage:\n"
                            "EllipseFromOutline [options] inputFile [outputFile] \n"
                            "Required:\n"
                            "inputFile:\t Input file to read from. See the readme for format specifications.\n"
                            "Optional:\n"
                            "outputFile: [Default MeshFile.out] Output file to write mesh to"
                            "Options:\n"
                            "-h, --help:\t Display command-line help message.\n");
            break;
        case TOO_MANY_INPUTS:
            fprintf(stdout, "Error: too many inputs provided, see usage below.\n");
            print_help(DEFAULT);
            break;
    }
};

ArgumentSpace::ArgumentSpace(string input_file) {
    // we have already checked that the input exists, so we can open it
    ifstream input;
    input.open(input_file, ifstream::in);

    // now read the information from the input file, and place it into the attributes of this class
    // ASSIGNING SOME ATTRIBUTES HERE AND READING THROUGH THE FILE IN A MEANINGFUL WAY
};


void ArgumentSpace::print_mode_specs() {

}

bool ArgumentReader::arg_is_flag(string argument)
{
    return argument[0] == '-';
};
bool ArgumentReader::has_flag(string const &flag) {
    for(const auto &a : all_args) {
        if (a == flag) {
            return true;
        }
    }
    return false;
};