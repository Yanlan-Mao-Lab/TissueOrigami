# include <string>
# include <vector>
# include <fstream>
# include <sys/stat.h>

using namespace std;

/**
 * @brief Define aliases for the different modes the simulation can run in.
 * 
 */
enum SimMode {
    NotSet = 0,
    OnTheGo = 1,
    Continue = 2,
    DisplaySave = 3
};

bool DirectoryExists(string path);

class ArgumentSpace{
    private:
        vector<string> allArgs;
        vector<string> flags;
        vector<string> nonFlags;

        SimMode simulationMode; // replace parameters[0]

    public:
        int nArgs;
        int nFlags;
        int nNonFlags;

        // will contain the path to the input file
        string pathToInputFile;
        // will contain the path to the output directory
        string pathToOutputDir;
        // will contain the path to the input directory (Continue or DisplaySave)
        string pathToInputDir;

        explicit ArgumentSpace(int argc, char **argv);
        void printModeSpecs();

        SimMode getSimulationMode();
        string getNonFlag(int i);
        
        bool argIsFlag(string argument);
        bool hasFlag(string const &flag);
};

class ArgumentReader{
    private:
        static void printHelp(SimMode mode);

    public:
        explicit ArgumentReader();
        static ArgumentSpace readInput(int argc, char **argv);
};