# include <string>
# include <vector>

using namespace std;

class ArgumentSpace{
    private:

    public:

        /**
         * @brief Construct a new Argument Space object from the input file, simultaneously validating its entries
         * 
         * @param input_file Path to the input file to read from
         */
        explicit ArgumentSpace(string input_file);

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
    TOO_MANY_INPUTS = 1
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