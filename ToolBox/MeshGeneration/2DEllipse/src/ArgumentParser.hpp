# include <string>
# include <vector>

using namespace std;

enum MeshMode {
    RECTANGLE = 0,
    WINGDISC = 1,
    TESSELATION2D  = 2,
    TESSELATION3D  = 3,
    TESSELATION3D_CYLINDER = 4
};
/**
 * @brief Translates values for meshing_mode from input files into enums
 * 
 * @param mode meshing_mode value read from the input file
 * @return mesh_mode The type of meshing to be performed
 */
MeshMode interpretMeshMode(string mode);

enum TissueType {
    WINGDISC_48HR = 0,
    WINGDISC_48HR_ECM = 1,
    WINGDISC_72HR = 2,
    OPTIC_CUP = 3,
    HALF_DISC = 4,
    SPHERICAL_ORGANOID = 5,
    TUBULAR_ORGANOID = 6,
    RECTANGLE_WITH_ECM = 7,
    RECTANGLE_WITHOUT_ECM = 8
};
TissueType interpretTissueType(string tt);

// this structure stores the flags that we can mark when certain inputs are read
struct argument_flags
{
    bool meshing_mode = false;
    bool selectTissueType = false;
    bool ABHeight = false;
    bool PrismSideLen = false;
    bool nzLayers = false;
    bool length[2] = {false};
    bool width[2] = {false};
    bool outline = false;
};

// the maximum number of fields that can appear in an input file
const int max_number_inputs = 11;

class ArgumentSpace{
    private:
        // private flags for input validation
        argument_flags vars_set;

    public:
        // signals whether we are meshing a wingdisc, rectangle, or pre-built 2d or 3d tesselation
        MeshMode meshing_mode;
        // determines the type of tissue that we are meshing
        TissueType selectTissueType;

        // wing disc specific required parameters

        double length[2], width[2];
        string outline_file_path;

        // parameters required for all meshing methods

        double ABHeight, prismSideLen;
        int nzLayers;

        // optional parameters - set default values here

        bool symY = false;

        /**
         * @brief Construct a new Argument Space object from the input file, simultaneously validating its entries
         * 
         * @param input_file Path to the input file to read from
         */
        explicit ArgumentSpace(string input_file);

        /**
         * @brief Given the name of an input variable, assign value to the relevent class atrribute
         *
         * @param variable Name of the variable to assign to
         * @param value Value to assign
         */
        void assign_input(string variable, string value);
        /**
         * @brief Check that the combination of parsed inputs are sufficient to execute mesh generation
         * 
         */
        void validate_input_file_contents();

        /**
         * @brief Print to stdout the information that has been read into the executable
         * 
         */
        void print_mode_specs();
};

/**
 * @brief Case-toggle for the different "help" messages that can be displayed to the screen
 * 
 */
enum helpOpts {
    DEFAULT = 0,
    TOO_MANY_INPUTS = 1,
    NO_INPUTS = 2
};

class ArgumentReader{
    private:
        // contains all input arguments passed from the command line
        vector<string> all_args;
        // contains all flag arguments (those beginning - or --)
        vector<string> flags;
        // contains all non-flag arguments
        vector<string> non_flags;

    public:
        // the number of arguments recieved from the command line
        int num_args;
        // the number of flags recieved from the command line
        int num_flags;
        // the number of non-flag arguments recieved from the command line
        int num_non_flags;

        // path to the input file
        string path_to_input;
        // path to the output file
        string path_to_output;

        /**
         * @brief Read the inputs passed to the executable
         *
         * @param argc Number of command-line arguments passed
         * @param argv Command line arguments
         * @return ArgumentSpace
         */
        explicit ArgumentReader(int argc, char **argv);

        /**
         * @brief Command-line help interface, printed when -h flag is passed
         *
         */
        static void print_help(helpOpts mode);

        /**
         * @brief Validate the inputs provided to the executable
         * 
         */
        void validate_cmdline_inputs();

        /**
         * @brief Determines whether a given flag was passed to the executable
         *
         * @param flag Flag pattern to check for
         * @return true The flag was passed
         * @return false The flag was not passed
         */
        bool has_flag(string const &flag);
        /**
         * @brief Check whether a given command line input is a flag.
         *
         * Flags are identified by the first character being a single dash, '-'.
         *
         * @param[in] argument The command line argument passed
         * @return true The argument is a flag (by the above definition)
         * @return false The argument is not a flag
         */
        bool arg_is_flag(string argument);
};