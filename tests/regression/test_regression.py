## test_regression.py
## Created by: William Graham, 2022-08-04

## Regression test for TissueFolding simulation
## 
## The pipeline for the TissueFolding simulation can be broken down into "modules" as follows:
## mesh generation -> simulation -> qt visualisation
## These tests exist to ensure that changes to the codebase of any module does not alter the functionality of said module.
## This is done by storing reference files in the /ref_files directory - these are files generated by the unmodified codebase, storing inputs and outputs of a given module for known data.
## These are then checked against what the current codebase generates for any discrepencies.
## 
## There are 7 sample runs that have been produced to use as test cases:
## Run Number   Mesh type           StiffnessChange     Growth/ShapeChange  LumenGrowth 
## 07000        Rect                Y                   Y (xy, z)           N
## 07001        Sphere (small)      Y                   N                   Y
## 07002        Sphere (fine)       Y                   N                   Y
## 07003        WingDisc (small)    N                   Y (uniform)         N
## 07004        WingDisc (fine)     N                   Y (uniform)         N
## 07005        Sphere (small)      Y                   N                   Y
## 07006        Sphere (fine)       Y                   N                   Y
##
## The mesh files generated for these runs are as follows:
## 07000    96hrRectangleWing_SmallMesh_posCoord.mesh
## 07001    HalfSphere_50_mesh_withECM.mesh
## 07002    HalfSphere_20_mesh_withECM.mesh
## 07003    SmallExp48hNoPeriWithECM.mesh
## 07004    FineMeshExp48hrNoPeriWithECM.mesh
## 07005    HalfSphere_50_mesh_withECM.mesh
## 07006    HalfSphere_20_mesh_withECM.mesh
## The input to the simulation (note: one of the fields is the path to the mesh file!) is saved in the style of ref_input_to_sim_XXXX where XXXX is replaced with the corresponding run number.

## For the output of a simulation given an input mesh, there are further cases

import os, filecmp

class Test_RegressionClass():
    '''
    Class container grouping together regression tests.

    Attributes
    ----------
    dir_path : str
        string containing path to the test directory (where this source file is stored)
    test_cases : ["rect", "sphere", "wingd"]
        list containing string identifiers for the regression test cases
    ref_location : str
        relative path to directory where reference files are stored
    gen_location : str
        relative path to directory where generated files (from testing) are stored
    ref_XXXX : dict
        keys are entries of test_cases, values are paths to reference files for the test XXXX
    gen_XXXX : dict
        as ref_XXXX, but for files generated by the codebase

    Tests
    -----
    test_mesh_input_to_sim(self)
        Tests whether meshes generated for each of the test cases match the reference inputs
    test_sim_output_pre_qt(self)
        Tests whether the simulation outputs for each of the test cases match the reference outputs
    '''

    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # test case labels
    test_cases = ["07000", "07001", "07002", "07003", "07004", "07005", "07006"]

    # reference file location, relative to this file location
    ref_location = "ref_files"
    # generated file location, relative to this file location
    gen_location = "gen_files"

    # fetch the names of the reference files now and place them into dicts
    ref_mesh_input_to_sim = dict(); gen_mesh_input_to_sim = dict()
    ref_input_to_sim = dict()
    ref_sim_output_pre_qt = dict(); gen_sim_output_pre_qt = dict()
    for t_case in test_cases:
        ref_mesh_input_to_sim[t_case] = dir_path + "/" + ref_location + "/ref_mesh_input_to_sim_" + t_case + ".txt"
        ref_sim_output_pre_qt[t_case] = dir_path + "/" + ref_location + "/ref_sim_output_pre_qt_" + t_case + ".txt"
        ref_input_to_sim[t_case] = dir_path + "/" + ref_location + "/ref_input_to_sim_" + t_case + ".txt"
        gen_mesh_input_to_sim[t_case] = dir_path + "/" + gen_location + "/gen_mesh_input_to_sim_" + t_case + ".txt"
        gen_sim_output_pre_qt[t_case] = dir_path + "/" + gen_location + "/gen_sim_output_pre_qt_" + t_case + ".txt"

    def test_ref_files_exist(self):
        '''Check that all reference files can be found before beginning tests!
        
        Parameters
        ----------
        '''

        ## only check for existence of reference files - if these are missing or have been renamed we are in trouble!
        for t_case in self.test_cases:
            assert os.path.exists(self.ref_mesh_input_to_sim[t_case]), "Could not find reference file: " + self.ref_mesh_input_to_sim[t_case]
            assert os.path.exists(self.ref_sim_output_pre_qt[t_case]), "Could not find reference file: " + self.ref_sim_output_pre_qt[t_case]
            assert os.path.exists(self.ref_input_to_sim[t_case]), "Could not find reference file: " + self.ref_input_to_sim[t_case]

    def test_mesh_input_to_sim(self):
        '''For each test case, checks whether a generated mesh matches the reference mesh

        Parameters
        ----------
        '''
        for t_case in self.test_cases:
            # PLACEHOLDER FOR CALL TO CODE FOR GENERATION OF INPUT MESHES
            # SAVE THESE INPUTS TO THE NAMES PROVIDED IN mesh_input_to_sim_fdict
            # This should look something like:
            # call executable that creates the mesh for the test case t_case, and save it to the file gen_mesh_input_to_sim[t_case]

            # now assert that the input file generated matches the reference input file, for this t_case
            # if assert fails, print out test case that caused failure
            assert filecmp.cmp(self.ref_mesh_input_to_sim[t_case], self.gen_mesh_input_to_sim[t_case], shallow=False), "Input file mismatch in " + t_case + " case"              

    def test_sim_output_pre_qt(self):
        '''For each test case, checks whether a simulation output (prior to Qt processing) matches the reference output
        
        Parameters
        ----------
        '''
        for t_case in self.test_cases:
            # PLACEHOLDER FOR CALL TO SIMULATION RUN
            # This should look something like:
            # run simulation using input file in ref_mesh_input_to_sim[t_case], saving it to file gen_sim_output_pre_qt[t_case]
        
            # assert file contents are identical, print error in test case if they are not
            assert filecmp.cmp(self.ref_sim_output_pre_qt[t_case], self.gen_sim_output_pre_qt[t_case], shallow=False), "Input file mismatch in " + t_case + " case"