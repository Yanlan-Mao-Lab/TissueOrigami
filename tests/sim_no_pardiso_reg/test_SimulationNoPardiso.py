import os, filecmp, shutil, subprocess, glob
from weakref import ref

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
    outputs = [ "MeshFromEndPoint.mesh", \
                "Out", \
                "Save_Force", \
                "Save_Frame", \
                "Save_Growth", \
                "Save_GrowthRate", \
                "Save_GrowthRedistribution", \
                "Save_NodeBinding", \
                "Save_Packing", \
                "Save_PhysicalProp", \
                "Save_SpecificElementAndNodeTypes", \
                "Save_TensionCompression", \
                "tmp"
    ]

    def launch_run0700x(self, run_number):
        '''
        
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
        for file in self.outputs:
            shutil.move(self.exe_loc + "/" + file, gen_res_loc + "/" + file)

        # cleanup the input files that we copied across, leaving the exe_loc clean
        os.remove(self.exe_loc + "/" + model_input_file)
        for file in self.files_needed[run_number]:
            os.remove(self.exe_loc + "/" + file)

        # run comparisions between the reference outputs and the generated outputs
        # where the reference outputs are stored
        ref_res_loc = run_folder + "/" + self.ref_res_subdir
        for file in self.outputs:
            # tmp contains absolute paths - there is no way we can get these to match the reference files
            # as such, skip over these files
            if file == "tmp":
                continue
            else:
                ref_op = ref_res_loc + "/" + file
                gen_op = gen_res_loc + "/" + file
                assert filecmp.cmp(ref_op, gen_op, shallow=False), "Output mismatch between " + ref_op + " and " + gen_op
        
        return

    def test_allruns(self):
        for rn in self.run_numbers:
            self.launch_run0700x(rn)
        return