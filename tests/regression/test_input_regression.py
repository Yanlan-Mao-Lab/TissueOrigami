## test_input_regression
## Created by: William Graham, 2022-08-04

## Regression test for input files to TissueFolding simulation
## Tests whether input mesh files to the Finite Volume simulation remain identical to those generated prior to beginning edits to the codebase.
## The three test cases are:
## - Rectangular tissue
## - Spherical tissue
## - Wingdisc tissue
## Additionally, these simulations will encompass tissue evolution with the following mechanics:
## - Growth map specified
## - Stiffness specified
## - shape map specified

## The reference (ground-truth) input files for three test cases are saved to the (relative) directory ./regression_inputs
## Filenames should be structured as regression_input_XXXX, where XXXX is replaced by rect, sphere, or wingd.
## These files are commited to the git repository and should not be edited under any circumstances.
## Relative paths will need to be updated if these files are relocated.

## The test will call the mesh generation functions for each of the tissue test cases listed above, and should save the generated input files to ./generated_inputs
## Filenames should be structured as generated_input_XXXX, under the same conventions for XXXX as above.
## The test will then check that the contents of the generated input file matches the contents of the reference input files for the corresponding simulation.

from cgi import test
import glob, os
import filecmp

def test_input_regression():

    # test case labels
    test_cases = ['rect', 'sphere', 'wingd']
    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # PLACEHOLDER FOR CALL TO CODE FOR GENERATION OF INPUT MESHES
    
    # compare matching filenames
    for t_case in test_cases:       
        # obtain the corresponding files from reference_inputs and generated_inputs
        r_input = glob.glob(dir_path + "/regression_inputs/regression_input_" + t_case + ".txt")
        g_input = glob.glob(dir_path + "/generated_inputs/generated_input_" + t_case + ".txt")

        # r_input and g_input should be lists of length 1, containing strings. 
        # Fail test and throw error (test not working!) if this is not the case.
        if len(g_input)!=1:
            raise(RuntimeError("Non-unique generated input file for test_input_regression in test case: " + t_case + "\n Got: %d, Expected 1. Did you clear previously generated input files?" % len(g_input)))
        elif len(r_input)!=1:
            raise(RuntimeError("Non unique reference input file for test_input_regression in test case: " + t_case + "\n Got: %d, Expected 1. Has the location of the reference input files been modified?" % len(r_input)))
        else:
            # proceed with test
            r_input_file = r_input[0]
            g_input_file = g_input[0]
            # assert file contents are identical, print error in test case if they are not
            assert filecmp.cmp(r_input_file, g_input_file, shallow=False), "Input file mismatch in " + t_case + " case"