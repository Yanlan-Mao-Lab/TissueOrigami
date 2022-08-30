import os, filecmp, shutil, subprocess
import pytest
from pytest_check import check
import warnings

def cleanup(fnames):
    '''
    Removes files whose names match those in fnames
    Parameters
    ----------
    fnames: \t List of strings, files to be removed
    '''

    for f in fnames:
        os.remove(f)
    return

class Test_SimulationNoPardiso():
    '''
    Tests the output of the simulation executable, built without Pardiso, to reference outputs.
    '''

    # path to preappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))
    # subfolder of run directories that reference results are saved to
    ref_res_subdir = "expected_out"
    # subfolder of run directories that generated results are saved to
    gen_res_subdir = "generatd_out"

    # location of the executable (relative to this file)
    exe_loc = dir_path + "/../../TissueFolding"
    # name of the executable
    exe_name = "TissueFolding"

    # location of the binary-to-text cpp source
    converter_exe = dir_path + "/../simOutputsToTxt.o"

    # variables that will store the information each test needs to run on

    # the run numbers
    run_numbers = [7, 8, 9]

    # files that the run needs to be able to see
    files_needed = dict()
    # run 07007 needs access to the following (non-modelinput) files
    files_needed[7] = [ "smallRectangle.mesh", \
                        "Stiffness96hrRectangleWing_Reduction_0", \
                        "Stiffness96hrRectangleWing_Reduction_4", \
                        "ShapeChangeRate96hrRectangleWingZ_Reduction_3", \
                        "ShapeChangeRate96hrRectangleWingXY_Reduction_3"
    ]
    # run 07008 needs access to the following (non-modelinput) files
    files_needed[8] = [ "smallWingDisc.mesh" ]
    # run 07009 needs access to the following (non-modelinput) files
    files_needed[9] = [ "smallSphere.mesh", \
                        "StiffnessMatrix_200_200_homogeneous.txt", \
                        "StiffnessMatrix_200_200_20_new.txt"
    ]

    # output file names that the simulation produces
    # we will have to convert these before comparing
    bin_outputs = [ "Force", \
                    "Growth", \
                    "GrowthRate", \
                    "GrowthRedistribution", \
                    "Packing", \
                    "PhysicalProp", \
                    "SpecificElementAndNodeTypes", \
                    "TensionCompression"
                    ]
    # these outputs can be compared directly, if we wish
    oth_outputs = [ "MeshFromEndPoint.mesh", \
                    "Save_Frame", \
                    "Save_NodeBinding", \
                    "Out"
                    ]
    # these outputs should be ignored - they are logs, dependent on the machine, etc
    ign_outputs = [ "tmp" ]

    @pytest.mark.parametrize('run_number', run_numbers)
    def test_run0700x(self, run_number):
        '''
        Runs the simulation run0700{run_number} (without Pardiso) and compares the results to the reference outputs.

        Output non-binary files are directly compared using filecmp.
        Output binary files are first compared using filecmp:
            In the event that these differ from the reference outputs, the two files are then _read in_ as they are done via the simulation, using the simOutputsToTxt.o executable.
            This produces two text files whose (comma-separated) entries are the values read in by the simulation from each file, with newlines placed at points where the behaviour in reading the file buffer changes.
            If the converted text files match, the test passes but a warning is thrown.
            Otherwise, the test fails.

        Parameters
        ----------
        run_number: \t Element of self.run_numbers, specifies the run number to execute.
        '''

        # this is the folder which we will be copying from/to, and which contains the reference outputs
        run_folder = self.dir_path + "/run0700" + str(run_number)
        model_input_file = "modelinput0700" + str(run_number)
        ref_model_input_file = run_folder + "/" + model_input_file

        # copy required files across to the location of the executable
        # model input file
        shutil.copyfile(ref_model_input_file, self.exe_loc + "/" + model_input_file)
        # auxillary files: growthmaps, stiffness maps, etc
        for file in self.files_needed[run_number]:
            shutil.copyfile(run_folder + "/" + file, self.exe_loc + "/" + file)

        # now that files have been copied across, run the executable with the appropriate command-line options
        command = "./" + self.exe_name + " -i ./" + model_input_file + " -mode SimulationOnTheGo -od ."
        # the tmp output is the result of piping stdout to the file tmp when running the simulation executable
        # to do this with subprocess, we can just map stdout to a file of the same name
        subprocess.run(command.split(), cwd=self.exe_loc, stdout = open(self.exe_loc + "/tmp", "w") )

        # once the run is over, move the output files back to the testing area, in the generated results location
        gen_res_loc = run_folder + "/" + self.gen_res_subdir
        # first check that the folder to move them to exists (and create it if not)
        if (not os.path.exists(gen_res_loc)):
            os.mkdir(gen_res_loc)
        # move output files back to testing area
        for file in self.oth_outputs:
            shutil.move(self.exe_loc + "/" + file, gen_res_loc + "/" + file)
        for f_type in self.bin_outputs:
            shutil.move(self.exe_loc + "/Save_" + f_type, gen_res_loc + "/Save_" + f_type)
        # move our ignored outputs to save trying to clean them up later
        for file in self.ign_outputs:
            shutil.move(self.exe_loc + "/" + file, gen_res_loc + "/" + file)

        # cleanup the input files that we copied across, leaving the exe_loc clean
        os.remove(self.exe_loc + "/" + model_input_file)
        for file in self.files_needed[run_number]:
            os.remove(self.exe_loc + "/" + file)

        # run comparisions between the reference outputs and the generated outputs
        # where the reference outputs are stored
        ref_res_loc = run_folder + "/" + self.ref_res_subdir

        # directly compare the outputs that do not require conversion
        for file in self.oth_outputs:
            ref_op = ref_res_loc + "/" + file
            gen_op = gen_res_loc + "/" + file
            # use nonfatal assert so that we always compare every file
            with check:
                # make the errors (if printed) more readable by not passing in filecmp.cmp into assert (
                # (avoids contextual expansion, as this info is printed in the error anyway)
                assert filecmp.cmp(ref_op, gen_op, shallow=False), "Output mismatch between: " + ref_op + " and " + gen_op

        # outputs that require conversion need to be converted, compared, and cleaned up
        # check that the converter exists!
        if not os.path.exists(self.converter_exe):
            raise RuntimeError("Error: output binary to txt converter not found at location: %s" % self.converter_exe)
        for f_type in self.bin_outputs:
            # we need to append the "Save_" to these because of how we stored them, so we could pass to the converter
            ref_op = ref_res_loc + "/Save_" + f_type
            gen_op = gen_res_loc + "/Save_" + f_type

            # attempt direct comparison of binaries - not recommended
            bin_compare = filecmp.cmp(ref_op, gen_op, shallow=False)
            txt_compare = False
            # if bin_compare is false, the binaries do not match,
            # but if we read them into the simulation, do they match now?
            if not bin_compare:
                # raise a warning that we are having to binary convert!
                warnings.warn("Warning: having to convert on output: %s" % f_type)
                # run the converter on these output files, the mode to pass to the converter is given by f_type
                command = self.converter_exe + " " + f_type + " " + ref_op + " " + gen_op
                subprocess.run(command.split(), cwd=self.dir_path)
                # compare the converted outputs
                txt_compare = filecmp.cmp(ref_op+"-read.txt", gen_op+"-read.txt", shallow=False)
                # cleanup the converted files that we made
                os.remove(ref_op+"-read.txt")
                os.remove(gen_op+"-read.txt")
            # test passes IFF bin_compare = True or (bin_compare = False and txt_compare = True)
            # since we do not assign to txt_compare unless bin_compare is false, we can simply check
            # bin_compare || txt_compare
            with check:
                assert (bin_compare or txt_compare), "Output mismatch between: " + ref_op + " and " + gen_op
        return