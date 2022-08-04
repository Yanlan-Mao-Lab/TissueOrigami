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

class Test_RegressionClass():

    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # test case labels
    test_cases = ["rect", "sphere", "wingd"]

    # reference file location, relative to this file location
    ref_location = "ref_files"
    # generated file location, relative to this file location
    gen_location = "gen_files"

    # fetch the names of the reference files now and place them into dicts
    ref_mesh_input_to_sim = dict(); gen_mesh_input_to_sim = dict()
    ref_sim_output_pre_qt = dict(); gen_sim_output_pre_qt = dict()

    for t_case in test_cases:
        ref_mesh_input_to_sim[t_case] = dir_path + "/" + ref_location + "/ref_mesh_input_to_sim_" + t_case + ".txt"
        ref_sim_output_pre_qt[t_case] = dir_path + "/" + ref_location + "/ref_sim_output_pre_qt_" + t_case + ".txt"
        gen_mesh_input_to_sim[t_case] = dir_path + "/" + gen_location + "/gen_mesh_input_to_sim_" + t_case + ".txt"
        gen_sim_output_pre_qt[t_case] = dir_path + "/" + gen_location + "/gen_sim_output_pre_qt_" + t_case + ".txt"

    def test_mesh_input_to_simulation_regression(self):

        for t_case in self.test_cases:
            # PLACEHOLDER FOR CALL TO CODE FOR GENERATION OF INPUT MESHES
            # SAVE THESE INPUTS TO THE NAMES PROVIDED IN mesh_input_to_sim_fdict
            # This should look something like:
            # call executable that creates the mesh for the test case t_case, and save it to the file gen_mesh_input_to_sim[t_case]

            # now assert that the input file generated matches the reference input file, for this t_case
            # if assert fails, print out test case that caused failure
            assert filecmp.cmp(self.ref_mesh_input_to_sim[t_case], self.gen_mesh_input_to_sim[t_case], shallow=False), "Input file mismatch in " + t_case + " case"              

    def test_simulation_output_regression(self):

        for t_case in self.test_cases:
            # PLACEHOLDER FOR CALL TO SIMULATION RUN
            # This should look something like:
            # run simulation using input file in ref_mesh_input_to_sim[t_case], saving it to file gen_sim_output_pre_qt[t_case]
        
            # assert file contents are identical, print error in test case if they are not
            assert filecmp.cmp(self.ref_sim_output_pre_qt[t_case], self.gen_sim_output_pre_qt[t_case], shallow=False), "Input file mismatch in " + t_case + " case"