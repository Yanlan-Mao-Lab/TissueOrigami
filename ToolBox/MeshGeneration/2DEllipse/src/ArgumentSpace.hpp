# ifndef ARGUMENT_SPACE_H
# define ARGUMENT_SPACE_H

# include <string>
# include <vector>

enum class MeshMode {
    RECTANGLE = 0,
    WINGDISC = 1,
    TESSELATION2D  = 2,
    TESSELATION3D  = 3,
    TESSELATION3D_CYLINDER = 4
};

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

// this structure stores the flags that we can mark when certain inputs are read
struct ArgumentFlags
{
    bool meshingMode = false;
    bool selectTissueType = false;
    bool ABHeight = false;
    bool prismSideLen = false;
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
        ArgumentFlags var_flags;

    public:
        // signals whether we are meshing a wingdisc, rectangle, or pre-built 2d or 3d tesselation
        MeshMode meshing_mode;
        // determines the type of tissue that we are meshing
        TissueType selectTissueType;

        // wing disc specific required parameters

        double length[2], width[2];
        std::string outline_file_path;

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
        explicit ArgumentSpace(std::string input_file);

        /**
         * @brief Given the name of an input variable, assign value to the relevent class atrribute
         *
         * @param variable Name of the variable to assign to
         * @param value Value to assign
         */
        void assignInput(std::string variable, std::string value);
        /**
         * @brief Translates values for meshing_mode from input files into enums
         *
         * @param mode meshing_mode value read from the input file
         * @return mesh_mode The type of meshing to be performed
         */
        void interpretMeshMode(std::string mode);
        /**
         * @brief Translates values for tissue_type from input files into enums
         *
         * @param mode tissue_type value read from the input file
         * @return TissueType The type of tissue being meshed
         */
        void interpretTissueType(std::string tt);
        /**
         * @brief Check that the combination of parsed inputs are sufficient to execute mesh generation
         * 
         */
        void validateInputFileContents();

        /**
         * @brief Print to stdout the information that has been read into the executable
         * 
         */
        void printModeSpecs();
};

# endif // ARGUMENT_SPACE_H
