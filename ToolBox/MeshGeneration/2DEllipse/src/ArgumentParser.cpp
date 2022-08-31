# include <stdexcept>
# include "ArgumentParser.hpp"
# include <iostream>
# include <fstream>
# include <sstream>

using namespace std;

mesh_mode interpretMeshMode(string mode) {
    if (mode=="rec") {
        return REC;
    }
    else if (mode=="wgd") {
        return WGD;
    }
    else if (mode=="2d") {
        return T2D;
    }
    else if (mode=="3d") {
        return T3D;
    }
    else if (mode=="3d_cyl") {
        return T3D_CYL;
    }
    else {
        throw runtime_error("Invalid meshing_mode read from input: " + mode);
    }
}

tissue_type interpretTissueType(string tt) {
    int t = stoi(tt);
    if(t<0 || t>9) {
        throw runtime_error("Error: invalid tissue type given (" + tt + ")");
    }
    else {
        return tissue_type(t);
    }
}

explicit ArgumentReader::ArgumentReader(int argc, char **argv) {
    // parse command-line arguments
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
    // print help if it was requested
    // otherwise, assign the path_to_input and path_to_output variables
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
    ifstream input(input_file, ifstream::in);

    // now read the information from the input file, and place it into the attributes of this class
    string line;
    int expecting_variable = 0;
    while(getline(input, line)) {
        auto colon_pos = line.find(":");
        // does this line contain a colon?
        if (colon_pos == string::npos)
        {
            // no colon found, move on to next line
            // is accounts for additional whitespace in the inputs
            continue;
        }
        else if (colon_pos >= line.length()-1) {
            // colon is the last character in the line, invalid input!
            throw runtime_error("Error - no variable value provided on input line: " + line);
        }
        else {
            // there is indeed a colon - read up to it to identify the variable
            string variable = line.substr(0,colon_pos);
            // and then continue reading to identify the value that was provided
            string value = line.substr(colon_pos+1);
            // now assign the value provided to the variable
            assign_input(variable, value);
            cout << "I want to assign " << value << "to the variable " << variable;
        }
    }
    // having read all that we can from the input file, confirm we have been given a valid set of inputs
    validate_input_file_contents();
};

void ArgumentSpace::assign_input(string variable, string value) {
    if (variable == "meshing_mode") {
        meshing_mode = interpretMeshMode(value);
        vars_set.meshing_mode = true;
    }
    else if (variable == "tissueType") {
        selectTissueType = interpretTissueType(value);
        vars_set.selectTissueType = true;
    }
    else if (variable == "ABHeight") {
        ABHeight = stod(value);
        vars_set.ABHeight = true;
    }
    else if (variable == "PrismSideLength") {
        prismSideLen = stod(value);
        vars_set.PrismSideLen = true;
    }
    else if (variable == "nzLayers") {
        nzLayers = stoi(value);
        vars_set.nzLayers = true;
    }
    // potentially non-present arguments
    else if (variable == "length1") {
        length[0] = stod(value);
        vars_set.length[0] = true;
    }
    else if (variable == "length2") {
        length[1] = stod(value);
        vars_set.length[1] = true;
    }
    else if (variable == "width1") {
        width[0] = stod(value);
        vars_set.width[0] = true;
    }
    else if (variable == "width2") {
        width[1] = stod(value);
        vars_set.width[1] = true;
    }
    else if (variable == "outline") {
        outline_file_path = value;
        vars_set.outline = true;
    }
    // purely optional input arguments
    else if (variable == "symY") {
        symY = !(value == "0" || value == "false" || value == "False");
    }
    // unrecognised field, throw an error
    else {
        throw runtime_error("Error - did not recognise the variable name: " + variable);
    }
}

void ArgumentSpace::validate_input_file_contents() {
    // we always require the following to be set:
    // meshing_mode, selectTissueType, ABHeight, PrismSideLen, nzLayers
    if (!vars_set.meshing_mode) { throw runtime_error("Error - no meshing_mode variable defined in input\n");}
    else if (!vars_set.selectTissueType) { throw runtime_error("Error - no tissueType variable defined in input\n");}
    else if (!vars_set.ABHeight) { throw runtime_error("Error - no ABHeight variable defined in input\n");}
    else if (!vars_set.PrismSideLen) { throw runtime_error("Error - no PrismSideLen variable defined in input\n");}
    else if (!vars_set.nzLayers) { throw runtime_error("Error - no nzLayers variable defined in input\n");}

    // if we are meshing a wing disc, check that optional arguments were passed
    if (meshing_mode == WGD) {
        if (!vars_set.length[0]) { throw runtime_error("Error - meshing wing-disc, but no length1 variable provided\n");}
        else if (!vars_set.length[1]) { throw runtime_error("Error - meshing wing-disc, but no length2 variable provided\n");}
        else if (!vars_set.width[0]) { throw runtime_error("Error - meshing wing-disc, but no width1 variable provided\n");}
        else if (!vars_set.width[1]) { throw runtime_error("Error - meshing wing-disc, but no width2 variable provided\n");}
        else if (!vars_set.outline) { throw runtime_error("Error - meshing wing-disc, but no outline variable provided\n");}
    }
}

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