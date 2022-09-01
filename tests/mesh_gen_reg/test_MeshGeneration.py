## test_regression.py
## Created by: William Graham, 2022-08-04

import pytest
from pytest_check import check
import os, filecmp, shutil, subprocess

def cleanup(fnames):
    '''
    Removes files whose names match those in fnames, skipping over entries that do not match a path

    Parameters
    ----------
    fnames: \t List of strings, files to be removed
    '''

    for f in fnames:
        if os.path.exists(f):
            os.remove(f)
        else:
            print("Could not find file %s, skipping removal" % f)
    return

def trim_final_line(fnames):
    '''
    Removes the final line from a file

     Parameters
    ----------
    fnames: \t List of strings, files to have final line removed
    '''

    for fname in fnames:
        with open(fname, "r+", encoding = "utf-8") as file:
            # Move pointer to end of file
            file.seek(0, os.SEEK_END)
            # In the case that the last line is null, delete it
            pos = file.tell() - 1

            # Read characters backwards until hitting a newline character
            while pos > 0 and file.read(1) != "\n":
                pos -= 1
                file.seek(pos, os.SEEK_SET)
            
            # So long was we're not at the start of the file, delete everything FORWARD of the current position
            if pos > 0:
                file.seek(pos, os.SEEK_SET)
                file.truncate()
    return

class OutputFileRelocation():
    '''
    
    '''
    def __init__(self, out_name, gen_name, ref_name, trim=False):
        '''
        
        '''
        self.output_name = out_name
        self.rename_to = gen_name
        self.compare_to = ref_name
        self.trim = trim
        return

class InputFileRelocation():
    '''
    '''
    def __init__(self, ref_name, reloc_name=""):
        '''
        '''
        self.ref_name = ref_name
        if reloc_name=="":
            self.reloc_name = ref_name
        else:
            self.reloc_name = reloc_name
        return

class MeshGenRun():
    '''
    
    '''
    def __init__(self, exe_input_file, run_folder, exe_mesh_out="MeshFile.out"):
        '''
        
        '''
        self.exe_in = exe_input_file
        self.folder = run_folder

        self.oth_inputs = []

        self.exe_mesh_out = exe_mesh_out

        self.outputs_to_compare = []
        return

    def addComparisonFile(self, output_name, gen_name, ref_name, trim=False):
        '''
        
        '''
        self.outputs_to_compare.append(OutputFileRelocation(output_name, gen_name, ref_name, trim))
        return
    
    def addExeInput(self, ref_name, reloc_name=""):
        '''
        
        '''
        self.oth_inputs.append(InputFileRelocation(ref_name, reloc_name))
        return
    
    def copy_inputs_to_loc(self, path):
        '''
        '''
        # move the required input file first - do not allow renaming!
        shutil.copyfile(self.folder + "/" + self.exe_in, path + "/" + self.exe_in)
        # move any other required inputs we have, renaming them if specified
        for f in self.oth_inputs:
            shutil.copyfile(self.folder + "/" + f.ref_name, path + "/" + f.reloc_name)
        return

    def run_executable(self, exe_loc, exe_name):
        '''
        '''
        # create the command
        command = "./" + exe_name + " " + self.exe_in + " " + self.exe_mesh_out
        # executate mesh generation executable,
        # supplying the input file self.exe_in and output name self.exe_mesh_out as command-line args
        subprocess.Popen(command.split(), cwd=exe_loc)
        return
    
    def move_gen_outputs(self, exe_loc):
        '''
        '''
        # for each instance in outputs_to_compare, move the file to self.folder
        for f in self.outputs_to_compare:
            src = exe_loc + "/" + f.out_name
            dst = self.folder + "/" + f.rename_to
            shutil.move(src, dst)
        return
    
    def cleanup_moved_inputs(self, path):
        '''
        '''
        # do not allow destruction of the reference files!
        if os.path.abspath(path)==os.path.abspath(self.folder):
            raise RuntimeError("Error: will not delete reference folder files!")
        # clean up the copied input file, error if not found since we could be looking in the wrong directory
        if os.path.exists(path + "/" + self.exe_in):
            os.remove(path + "/" + self.exe_in)
        else:
            raise RuntimeError("Error: %s not found - are you looking in the correct directory?" % path + "/" + self.exe_in)
        # clean up any other input files that we copied across
        for f in self.oth_inputs:
            if os.path.exists(path + "/" + f.reloc_name):
                os.remove(path + "/" + f.reloc_name)
            else:
                raise RuntimeError("Error: %s not found - are you looking in the correct directory?" % path + "/" + f.reloc_name)
        return

    def compare_output_files(self):
        '''
        '''
        # perform _all_ comparisons, even if one pair do not match
        with check:
            # for each instance in self.outputs_to_compare, run filecmp.cmp on f.rename_to and f.compare_to
            for f in self.outputs_to_compare:
                ref_file = self.folder + "/" + f.compare_to
                gen_file = self.folder + "/" + f.rename_to
                # trim final line if this is flagged, before comparing
                if f.trim:
                    trim_final_line([ref_file, gen_file])
                assert filecmp(ref_file, gen_file), "Error: difference between reference %s and generated %s" % (ref_file, gen_file)
        return

class Test_MeshGeneration():
    '''
    Test that the functionality of the mesh generation executable has not been signficantly altered by changes to the codebase (regression test).

    These tests assume that the mesh generation executable is located in Toolbox/MeshGeneraion/2DEllipse/ and has been built using the CMake instructions.

    Attributes
    ----------

    Tests
    -----

    '''

    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # identifiers for the regression tests that we have
    # these should be the folder names containing the inputs to the executable and the reference outputs
    test_cases = [ "smallRectangle", "smallSphere" ]

    test_info = dict()
    smallRecRun = MeshGenRun("inputFile_smallRectangle", dir_path + "/smallRectangle")
    smallRecRun.addExeInput("smallRectangle.ele", "Points.1.ele")
    smallRecRun.addExeInput("smallRectangle.node", "Points.1.node")
    smallRecRun.addComparisonFile("MeshFile.out", "smallRectangle-gen.mesh", "smallRectangle.mesh")
    test_info["smallRectangle"] = smallRecRun

    smallSphRun = MeshGenRun("inputFile_smallSphere", dir_path + "/smallSphere")
    smallSphRun.addExeInput("SphericalTriangulation.txt", "SphericalTriangulation")
    smallSphRun.addComparisonFile("MeshFile.out", "smallSphere-gen.mesh", "smallSphere.mesh")
    test_info["smallSphere"] = smallSphRun

    # location of the folder containing the mesh generation executable
    executable_loc = dir_path + "/../../ToolBox/MeshGeneration/2DEllipse"
    # name of the mesh generation executable
    executable_name = "EllipseFromOutline"

    # output files that are produced by the executable that need to be cleaned up every time
    always_cleanup = [ executable_loc + "/NodesPostTesselation.out", executable_loc + "/VectorsPostTesselation.out", executable_loc + "/Points.node" ]

    @pytest.mark.parametrize("tc", test_cases)
    def test_mesh_generation(self, tc):
        '''
        
        '''
        # easily aliased information for this test
        mesh_gen_test = self.test_info[tc]

        # copy the required input files to the executable directory
        mesh_gen_test.copy_inputs_to_loc(self.executable_loc)
        
        # having copied the input files across, run the executable
        mesh_gen_test.run_executable(self.executable_loc, self.executable_name)

        # move the outputs (that we care about) to the testing area
        mesh_gen_test.move_gen_outputs(self.executable_loc)

        # cleanup the input file copies that we made
        mesh_gen_test.cleanup_moved_inputs(self.executable_loc)
        # cleanup any other files that always need to be removed from the directory after running the executable
        cleanup(self.always_cleanup)

        # now perform the file comparisons
        # ASSERT will be performed in here
        mesh_gen_test.compare_output_files()
        return