#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {

    int n = 250;
    string run_mode, fname;

    if (argc<3) {
        throw runtime_error("Error: not enough input arguments:\n Provide [-r read, -w write] filename\n");
    }
    else {
        run_mode = argv[1];
        fname = argv[2];
    }
    
    if (run_mode=="-r") {
        // want to read in the binary then dump the output
        double read_in[n];
        ifstream infile(fname, ifstream::in);

        string f_out = fname + "-out.txt";
        ofstream outfile(f_out, ofstream::out);

        for (int i = 0; i < n; i++)
        {
            infile.read((char *)&read_in[i], sizeof read_in[i]);
            outfile << read_in[i] << ",";
        }

        infile.close();
        outfile.close();
    }
    else if (run_mode=="-w") {
        // want to write the binary file
        ofstream outfile(fname, ofstream::binary);
        
        for (int i=0; i<n; i++)
        {
            double to_write = 1. / ((double)i + 1.);
            outfile.write((char*)&to_write, sizeof to_write);
        }
        outfile.close();
    }
    else {
        throw runtime_error("Error: unrecognised run mode (provide -r [read] or -w [write])");
    }
}
