// compiling this file:
// g++ read_sim_output.cpp -o ./reader.o

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// converts the output file Save_Force to txt from binary - functionally identical to read_GrowthRate
void read_Force(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 3; ++j)
            {
                double sysForceijk;
                input_file.read((char *)&sysForceijk, sizeof sysForceijk);
                output_file << sysForceijk << ",";
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_Growth to txt from binary
void read_Growth(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    double Fgjk;
                    input_file.read((char *)&Fgjk, sizeof Fgjk);
                    output_file << Fgjk << ",";
                }
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_GrowthRate to txt from binary
void read_GrowthRate(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 3; ++j)
            {
                double rxyz;
                input_file.read((char *)&rxyz, sizeof rxyz);
                output_file << rxyz << ",";
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_GrowthRedistribution to txt from binary
void read_GrowthRedistribution(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 3; ++j)
            {
                bool readbool;
                input_file.read((char *)&readbool, sizeof readbool);
                output_file << readbool << ",";
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_Packing to txt from binary
void read_Packing(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        int n;
        input_file.read((char *) &n, sizeof n);
        output_file << n << "\n";

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            // reads the element ids (genuinely done twice each time)
            for (int j = 0; j < n; ++j)
            {
                int elementId;
                input_file.read((char *)&elementId, sizeof elementId);
                output_file << elementId << ",";
                input_file.read((char *)&elementId, sizeof elementId);
                output_file << elementId << ",";
            }
            // reads the forces (genuinely done twice each time)
            for (int j = 0; j < n; ++j)
            {
                double Fx, Fy, Fz;
                input_file.read((char *)&Fx, sizeof Fx);
                input_file.read((char *)&Fy, sizeof Fy);
                input_file.read((char *)&Fz, sizeof Fz);
                output_file << Fx << "," << Fy << "," << Fz;
                input_file.read((char *)&Fx, sizeof Fx);
                input_file.read((char *)&Fy, sizeof Fy);
                input_file.read((char *)&Fz, sizeof Fz);
                output_file << Fx << "," << Fy << "," << Fz;
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_PhysicalProp to txt from binary - functionally the same as read_GrowthRate
void read_PhysicalProp(int argc, char **argv) {
    /** see line 2624 of Simulation.cpp
     * For each element, the format is:
     * [Young's modulus] [internal viscosity] [zRemodelling]
     * Then for each node:
     * [external viscosity_x] [external viscosity_y] [external viscosity_z]
     */
    // we can thus "cheat" this read-in and just read 3 values at a time until we reach the end-of-file, as all inputs are read as doubles
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 3; ++j)
            {
                double variable;
                input_file.read((char *)&variable, sizeof variable);
                output_file << variable << ",";
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_SpecificElementAndNodeTypes to txt from binary
void read_SpecificElementAndNodeTypes(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        // actin layer marks
        int elementId, actinCounter;
        input_file.read((char *) &actinCounter, sizeof actinCounter);
        output_file << actinCounter << "\n";
        for (int j=0; j<actinCounter; j++) {
            input_file.read((char *) &elementId, sizeof elementId);
            output_file << elementId << ",";
        }
        output_file << ",";

        // ECM layer marks
        int ecmCounter;
        input_file.read((char *)&ecmCounter, sizeof ecmCounter);
        output_file << ecmCounter << "\n";
        for (int j = 0; j < actinCounter; j++)
        {
            input_file.read((char *)&elementId, sizeof elementId);
            output_file << elementId << ",";
        }
        output_file << ",";

        // marker ellipses for elements
        int ellipseCounter, ellipseBandID;
        input_file.read((char *) &ellipseCounter, sizeof ellipseCounter);
        output_file << ellipseCounter << "\n";
        for (int j = 0; j<ellipseCounter; j++) {
            input_file.read((char *) &elementId, sizeof elementId);
            input_file.read((char *) &ellipseBandID, sizeof ellipseBandID);
            output_file << elementId << "," << ellipseBandID << ",";
        }
        output_file << "\n";

        // marker ellipses to display for nodes
        int nodesCounter, nodeID;
        input_file.read((char *) &nodesCounter, sizeof nodesCounter);
        output_file << nodesCounter << "\n";
        for (int j = 0; j < nodesCounter; j++) {
            input_file.read((char *) &nodeID, sizeof nodeID);
            input_file.read((char *) &ellipseBandID, sizeof ellipseBandID);
            output_file << nodeID << "," << ellipseBandID << ",";
        }
        output_file << "\n";
        
        input_file.close();
        output_file.close();
    }
}

// converts the output file Save_TensionCompression to txt from binary
void read_TensionCompression(int argc, char **argv) {
    for (int fnumber = 2; fnumber < argc; fnumber++)
    {
        string read_from = argv[fnumber];
        ifstream input_file(read_from, ifstream::in);
        ofstream output_file(read_from + "-read.txt");
        cout << "Now reading " << read_from << " (into) -> " << read_from << "-read.txt" << endl;

        while (input_file.good())
        {
            if (input_file.eof())
            {
                break;
            }
            for (int j = 0; j < 6; ++j)
            {
                double S_element;
                input_file.read((char *)&S_element, sizeof S_element);
                output_file << S_element << ",";
            }
            output_file << "\n";
        }
        input_file.close();
        output_file.close();
    }
}

/**
 * @brief Add a readme file when this is complete!
 *
 * Reads the binary outputs of the simulation in, and then writes them out again as txt files.
 * This enables visual scanning of the file contents, and direct comparison of the values being read into the GUI from the outputs.
 *
 * Command-line use is:
 * EXE_NAME outputfile_type [pass files in]
 *
 * where outputfile_type is one of:
 * - Force
 * - Growth
 * - GrowthRate
 * - GrowthRedistribution
 * - Packing
 * - PhysicalProp
 * - SpecificElementAndNodeTypes
 * - TensionCompression
 *
 * Outputs are written to files of the same name as the inputs, appended with "-read.txt"
 *
 * The Frames and NodeBinding files are already saved as text files.
 */
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        throw runtime_error("Not enough inputs provided!");
        exit(1);
    }

    string mode = string(argv[1]);
    if (mode == "Force")
    {
        read_Force(argc, argv);
    }
    else if (mode == "Growth")
    {
        read_Growth(argc, argv);
    }
    else if (mode == "GrowthRate")
    {
        read_GrowthRate(argc, argv);
    }
    else if (mode == "GrowthRedistribution")
    {
        read_GrowthRedistribution(argc, argv);
    }
    else if (mode == "Packing")
    {
        read_Packing(argc, argv);
    }
    else if (mode == "PhysicalProp")
    {
        read_PhysicalProp(argc, argv);
    }
    else if (mode == "SpecificElementAndNodeTypes")
    {
        read_SpecificElementAndNodeTypes(argc, argv);
    }
    else if (mode == "TensionCompression")
    {
        read_TensionCompression(argc, argv);
    }
    else
    {
        //throw runtime_error("Error: please provide a valid output file type to read\n");
        exit(1);
    }
}