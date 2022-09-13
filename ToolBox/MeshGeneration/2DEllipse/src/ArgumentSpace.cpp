# include "ArgumentSpace.hpp"
# include <iostream>
# include <fstream>
# include <algorithm>

using namespace std;

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
            // is accounts for additional newlines in the inputs
            continue;
        }
        else if (colon_pos >= line.length()-1) {
            // colon is the last character in the line
            // continue reading in case this argument was optional (errors will be caught when flags for required variables are checked)
            cout << "No value provided on line: " << line << endl;
        }
        else {
            // there is indeed a colon - read up to it to identify the variable
            string variable = line.substr(0,colon_pos);
            // strip character spaces
            variable.erase( remove_if( variable.begin(), variable.end(), ::isspace), variable.end() );
            // continue reading to identify the value that was provided
            string value = line.substr(colon_pos+1);
            // strip character spaces
            value.erase(remove_if(value.begin(), value.end(), ::isspace), value.end());
            // check that, prior to attempting to assign, we have been left with non-empty strings value and variable
            if (!(variable.empty() || value.empty())) {
                // assign the value provided to the variable
                assign_input(variable, value);
            }
        }
    }
    // having read all that we can from the input file, confirm we have been given a valid set of inputs
    validate_input_file_contents();
};

void ArgumentSpace::assign_input(string variable, string value) {
    if (variable == "meshingMode") {
        interpretMeshMode(value);
        vars_set.meshingMode = true;
    }
    else if (variable == "tissueType") {
        interpretTissueType(value);
        vars_set.selectTissueType = true;
    }
    else if (variable == "ABHeight") {
        ABHeight = stod(value);
        vars_set.ABHeight = true;
    }
    else if (variable == "prismSideLength") {
        prismSideLen = stod(value);
        vars_set.prismSideLen = true;
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

void ArgumentSpace::interpretMeshMode(string mode)
{
    if (mode == "rec")
    {
        meshing_mode = MeshMode::RECTANGLE;
    }
    else if (mode == "wgd")
    {
        meshing_mode = MeshMode::WINGDISC;
    }
    else if (mode == "2d")
    {
        meshing_mode = MeshMode::TESSELATION2D;
    }
    else if (mode == "3d")
    {
        meshing_mode = MeshMode::TESSELATION3D;
    }
    else if (mode == "3d_cyl")
    {
        meshing_mode = MeshMode::TESSELATION3D_CYLINDER;
    }
    else
    {
        throw runtime_error("Invalid meshing_mode read from input: " + mode);
    }
}

void ArgumentSpace::interpretTissueType(string tt)
{
    int t = stoi(tt);
    if (t < 0 || t > 9)
    {
        throw runtime_error("Error: invalid tissue type given (" + tt + ")");
    }
    else
    {
        selectTissueType = TissueType(t);
    }
}

void ArgumentSpace::validate_input_file_contents() {
    // we always require the following to be set:
    // meshing_mode, selectTissueType, ABHeight, PrismSideLen, nzLayers
    if (!vars_set.meshingMode) { throw runtime_error("Error - no meshing_mode variable defined in input\n");}
    else if (!vars_set.selectTissueType) { throw runtime_error("Error - no tissueType variable defined in input\n");}
    else if (!vars_set.ABHeight) { throw runtime_error("Error - no ABHeight variable defined in input\n");}
    else if (!vars_set.prismSideLen) { throw runtime_error("Error - no PrismSideLen variable defined in input\n");}
    else if (!vars_set.nzLayers) { throw runtime_error("Error - no nzLayers variable defined in input\n");}

    // if we are meshing a wing disc, check that optional arguments were passed
    if (meshing_mode == MeshMode::WINGDISC) {
        if (!vars_set.length[0]) { throw runtime_error("Error - meshing wing-disc, but no length1 variable provided\n");}
        else if (!vars_set.length[1]) { throw runtime_error("Error - meshing wing-disc, but no length2 variable provided\n");}
        else if (!vars_set.width[0]) { throw runtime_error("Error - meshing wing-disc, but no width1 variable provided\n");}
        else if (!vars_set.width[1]) { throw runtime_error("Error - meshing wing-disc, but no width2 variable provided\n");}
        else if (!vars_set.outline) { throw runtime_error("Error - meshing wing-disc, but no outline variable provided\n");}
    }
}

void ArgumentSpace::print_mode_specs() {
    cout << "Printing meshing specifications:" << endl;
    cout << "Meshing mode: ";
    switch (meshing_mode)
    {
        case MeshMode::WINGDISC:
            cout << "Wing-disc / elliptical (WGD)" << endl;
            cout << "Length = {" << length[0] << " , " << length[1] << "}" << endl;
            cout << "Width = {" << width[0] << " , " << width[1] << "}" << endl;
            cout << "outline file: " << outline_file_path << endl;
            break;
        case MeshMode::RECTANGLE:
            cout << "Rectangular (REC)" << endl;
            break;
        case MeshMode::TESSELATION2D:
            cout << "Pre-built 2D tesselation (2D)" << endl;
            break;
        case MeshMode::TESSELATION3D:
            cout << "Pre-built 3D tesselation (3D)" << endl;
            break;
        case MeshMode::TESSELATION3D_CYLINDER:
            cout << "Pre-built cylindrical 3D tesselation (3D_CYL)" << endl;
            break;
    }
    // common variables
    cout << "ABHeight: " << ABHeight << endl;
    cout << "Prism side length: " << prismSideLen << endl;
    cout << "No. z layers: " << nzLayers << endl;
}
