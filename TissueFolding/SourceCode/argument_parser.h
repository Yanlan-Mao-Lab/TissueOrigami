# include <string>
# include <vector>
# include <fstream>

using namespace std;

// Define masks for the type of mesh generation
#define NO_RUN_SET 0
#define WDC_RUN 1
#define REC_RUN 2
#define TESS_2D -1
#define TESS_3D -2

class ArgumentSpace{
    private:
        vector<string> all_args;
        vector<string> flags;
        vector<string> non_flags;

        int run_type; // replace parameters[0]

    public:
        int num_args;
        int num_flags;
        int num_non_flags;

        double length1 = 0., length2 = 0., width1 = 0., width2 = 0.;     // replace parameters[1,2,3,4]
        double ABHeight = 0., prism_side_len = 0., n_layers_z = 0.; // replace parameters[5,6,7]. types could be improved
        double sym_y = 0.;                           // replace parameters[8], could be bool?
        string path_to_input = ".";                  // path to inputOutline file

        explicit ArgumentSpace(int argc, char **argv);
        void print_mode_specs();

        int get_run_type();
        string get_non_flag(int i);
        
        bool arg_is_flag(string argument);
        bool has_flag(string const &flag);
};

class ArgumentReader{
    private:
        static void print_help(int mode);

    public:
        explicit ArgumentReader();
        static ArgumentSpace read_input(int argc, char **argv);
};