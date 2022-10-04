# ifndef GUI_ARGUMENT_PARSER_H
# define GUI_ARGUMENT_PARSER_H

# include <string>

class GUIArgumentReader
{
public:
    GUIArgumentReader(int argc, char **argv);

    std::string input_file;         // path to the modelinput file
    std::string output_dir = ".";   // path to the output directory
    std::string read_in_dir;        // path to directory containing the output of a previous simulation, to be read in


private:
    /**
     * @brief Print help message for the executable, then exit.
     * 
     * @param bad_input If true an invalid input was passed to the executable, and the printed message and exit status are changed accordingly.
     */
    void printHelp(bool bad_input);
};

# endif