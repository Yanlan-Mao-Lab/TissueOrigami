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
    Removes the final line from each member of a list of files

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
    The purpose of this class is to aid in the extraction of output files from the mesh generation executable from the executable's directory back to the testing area, and then to track which reference file this particular generated output should be compared against.

    An output from the executable will be identified by output_name, it will be saved to the testing area under rename_to, and the corresponding reference file to compare it to is compare_to.
    If trim is true, then this is an output produced by triangle, and so the final line needs to be removed (due to containing hard-coded paths) before comparison happens.

    Attributes
    ----------
    output_name : Name of an output file produced by the mesh generation exe
    compare_to : Name of the reference file to compare this output to
    rename_to : [Default output_name] Name to save this output to in the testing area
    trim : [Default false] If true, remove the finl line from this output file before comparing

    Construction
    ------------
    OutputFileRelocation(output_name, rename_to, compare_to, trim=False)

    Methods
    -------
    '''
    def __init__(self, output_name, compare_to, rename_to="", trim=False):
        '''
        Creates a new instance of OutputFileRelocation with the assigned attributes.
        '''
        self.output_name = output_name
        if rename_to=="":
            self.rename_to = output_name
        else:
            self.rename_to = rename_to
        self.compare_to = compare_to
        self.trim = trim
        return

class InputFileRelocation():
    '''
    The purpose of this class is to aid in the placement of auxillary input files for each regression test. That is, to move additional input files (such as outline files, pre-built meshes, etc) to the mesh generation executable's directory.

    A required input to the executable is named ref_name, and will be copied to the executable's location under the name reloc_name. If reloc_name is not provided, it will default to ref_name.

    Attributes
    ----------
    ref_name : Name of an input file to be passed to the mesh generation executable
    reloc_name : [Default ref_name] Name under which to save a copy of ref_name to the executable's directory

    Construction
    ------------
    InputFileRelocation(ref_name, reloc_name="")

    Methods
    -------
    '''
    def __init__(self, ref_name, reloc_name=""):
        '''
        Creates a new instance of InputFileRelocation with the assigned attributes.
        '''
        self.ref_name = ref_name
        if reloc_name=="":
            self.reloc_name = ref_name
        else:
            self.reloc_name = reloc_name
        return

class MeshGenRun():
    '''
    Class that handles all information associated with a given regression test case for the mesh generation process.
    
    Attributes
    ----------
    exe_in : Name of the input file that is to be passed (on the command line) to the mesh generation executable
    folder : Name of the folder (within mesh_gen_reg) that this test's files are stored in

    oth_inputs : List of InputFileRelocation, additional files that need to be visible to the mesh generation executable (and the appropriate name for them) when running this test

    exe_mesh_out : [Default MeshFile.out] If a meshfile is being produced, pass this on the command line to rename the meshfile via the executable.

    outputs_to_compare : List of OutputFileRelocation, output files to pull from the executable, and the corresponding reference files to compare them to

    Construction
    ------------
    MeshGenRun(exe_in, folder, exe_mesh_out="MeshFile.out")
    addComparisonFile(output_name, rename_to, compare_to, trim=False)
    addExeInput(ref_name, reloc_name="")

    Methods
    -------
    copy_inputs_to_loc(path)
    run_executable(exe_loc, exe_name)
    move_gen_outputs(exe_loc)
    cleanup_moved_inputs(path)
    compare_output_files()
    cleanup_output_files()
    '''

    # Construction methods

    def __init__(self, exe_in, folder, exe_mesh_out="MeshFile.out"):
        '''
        Creates a new MeshGenRun instance with the exe_in, folder, and exe_mesh_out fields set.
        '''
        self.exe_in = exe_in
        self.folder = folder

        self.oth_inputs = []

        self.exe_mesh_out = exe_mesh_out

        self.outputs_to_compare = []
        return

    def addComparisonFile(self, output_name, compare_to, rename_to="", trim=False):
        '''
        Appends a new instance of OutputFileRelocation to the list of outputs to extract for this test.
        See OutputFileRelocation for parameter details.
        '''
        self.outputs_to_compare.append(OutputFileRelocation(output_name, compare_to, rename_to, trim))
        return
    
    def addExeInput(self, ref_name, reloc_name=""):
        '''
        Appends a new instance of InputFileRelocation to the list of inputs to copy for this test.
        See InputFileRelocation for parameter details.
        '''
        self.oth_inputs.append(InputFileRelocation(ref_name, reloc_name))
        return
    
    # Usage methods

    def copy_inputs_to_loc(self, path):
        '''
        Copies the input file self.exe_in, as well as any additional input files in self.oth_inputs, to the path

        Parameters
        ----------
        path : Location to copy input files to (usually the executable location)
        '''
        # move the required input file first - do not allow renaming!
        shutil.copyfile(self.folder + "/" + self.exe_in, path + "/" + self.exe_in)
        # move any other required inputs we have, renaming them if specified
        for f in self.oth_inputs:
            shutil.copyfile(self.folder + "/" + f.ref_name, path + "/" + f.reloc_name)
        return

    def run_executable(self, exe_loc, exe_name):
        '''
        Runs the executable exe_name in the working directory exe_loc, providing the input file self.exe_in and (optional) output name for the meshfile self.exe_mesh_out on the command line.

        Parameters
        ----------
        exe_loc : Path to the directory in which the (mesh generation) executable is stored
        exe_name : Name of the mesh generation executable
        '''
        # create the command
        command = "./" + exe_name + " " + self.exe_in + " " + self.exe_mesh_out
        # executate mesh generation executable,
        # supplying the input file self.exe_in and output name self.exe_mesh_out as command-line args
        running_exe = subprocess.Popen(command.split(), cwd=exe_loc)
        # we MUST wait for this subprocess to complete before we allow the program to continue
        running_exe.wait()
        return
    
    def move_gen_outputs(self, exe_loc):
        '''
        Moves the outputs self.outputs_to_compare from the exe_loc to the testing area.

        Parameters
        ----------
        exe_loc : Directory in which the output files are stored (usually the directory containing the executable)
        '''
        # for each instance in outputs_to_compare, move the file to self.folder
        for f in self.outputs_to_compare:
            src = exe_loc + "/" + f.output_name
            dst = self.folder + "/" + f.rename_to
            shutil.move(src, dst)
        return
    
    def cleanup_moved_inputs(self, path):
        '''
        Removes the copies of the input files self.exe_in and self.oth_inputs from the directory path.
        Throws an error if path matches the testing directory, so that reference files are not deleted.

        Parameters
        ----------
        path : Directory in which the input file copies have been saved to
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
        Compare all output files in self.outputs_to_compare to their reference-file counterparts, ASSERTing if the files match.
        '''
        # we should track if there is a failed check, so that we don't delete the outputs after we are done
        all_checks_passed = True
        # perform _all_ comparisons, even if one pair do not match
        with check:
            # for each instance in self.outputs_to_compare, run filecmp.cmp on f.rename_to and f.compare_to
            for f in self.outputs_to_compare:
                ref_file = self.folder + "/" + f.compare_to
                gen_file = self.folder + "/" + f.rename_to
                # trim final line if this is flagged, before comparing
                # however, we DON'T want to alter the reference files, so we need to make temporary copies, then clean them up
                if f.trim:
                    # create temporary filenames
                    tmp_ref = self.folder + "/tmp-trimmedline-" + f.compare_to
                    tmp_gen = self.folder + "/tmp-trimmedline-" + f.rename_to
                    # copy outputs to temporary files
                    shutil.copyfile(ref_file, tmp_ref)
                    shutil.copyfile(gen_file, tmp_gen)
                    # trim final lines off the temporary files
                    trim_final_line([tmp_ref, tmp_gen])
                    # compare the now-trimmed files
                    files_match = filecmp.cmp(tmp_ref, tmp_gen)
                    if not files_match:
                        all_checks_passed = False
                    assert files_match, "Error: difference between reference %s and generated %s" % (ref_file, gen_file)
                    # cleanup temporary files
                    os.remove(tmp_ref)
                    os.remove(tmp_gen)
                else:
                    # no need for trimming - begin checking if the files match
                    files_match = filecmp.cmp(ref_file, gen_file)
                    if not files_match:
                        all_checks_passed = False
                    assert files_match, "Error: difference between reference %s and generated %s" % (ref_file, gen_file)
        return all_checks_passed

    def cleanup_output_files(self):
        '''
        Deletes any output files that were copied into the testing area.
        Only to be used after successful tests.
        '''
        # for each instance in self.outputs_to_compare, delete f.rename_to, which is the name under which this generated output is saved
        for f in self.outputs_to_compare:
            os.remove(self.folder + "/" + f.rename_to)
        return

class Test_MeshGeneration():
    '''
    Test that the functionality of the mesh generation executable has not been signficantly altered by changes to the codebase (regression test).

    These tests assume that the mesh generation executable is located in the Toolbox/MeshGeneraion/2DEllipse/ directory, and has been built using the CMake instructions.

    Each test is given a unique ID and the necessary reference and input files are saved to a subdirectory of mesh_gen_reg/ with that ID as the folder name.

    Attributes
    ----------
    dir_path : Absolute path to the directory containing this file

    test_cases : Identifiers for the regression test cases that we have, also coincide with the folder naming convention
    test_info : dictionary whose keys are the entries of test_cases. Each entry is an instance of MeshGenRun providing the information needed to execute the test whose ID coincides with the key

    executable_loc : Location of the mesh generation executable to be called
    executable_name : Name of the mesh generation executable to be called

    always_cleanup : List of filenames, these are auxillary files that are produced by the mesh generation executable that need to be cleaned up after each call.

    Tests
    -----
    test_mesh_generation(TEST_CASE):

    For each test ID;
    - Copy the input files to the executable directory
    - Run the executable
    - Move the output files back to the testing directory
    - Compare the outputs to the reference files

    '''

    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # identifiers for the regression tests that we have
    # these should be the folder names containing the inputs to the executable and the reference outputs
    test_cases = [ "smallRectangle", "smallSphere", "smallWingDisc-pt1", "smallWingDisc-pt2" ]

    test_info = dict()
    # create the smallRectangle test information
    smallRecRun = MeshGenRun("inputFile_smallRectangle", dir_path + "/smallRectangle")
    smallRecRun.addExeInput("smallRectangle.ele", "Points.1.ele")
    smallRecRun.addExeInput("smallRectangle.node", "Points.1.node")
    smallRecRun.addComparisonFile("MeshFile.out", "smallRectangle.mesh", rename_to="smallRectangle-gen.mesh", )
    test_info["smallRectangle"] = smallRecRun

    # create the smallSphere test information
    smallSphRun = MeshGenRun("inputFile_smallSphere", dir_path + "/smallSphere")
    smallSphRun.addExeInput("SphericalTriangulation.txt", "SphericalTriangulation")
    smallSphRun.addComparisonFile("MeshFile.out", "smallSphere.mesh", rename_to="smallSphere-gen.mesh", )
    test_info["smallSphere"] = smallSphRun

    # create the smallWingDisc (pt1, 2d meshing from outline file) information
    smallWgdRun1 = MeshGenRun("inputFile_smallWingDisc-1", dir_path + "/smallWingDisc")
    smallWgdRun1.addExeInput("48hrDiscSymmetricOutline")
    smallWgdRun1.addComparisonFile("Points.1.ele", "WD_Points1.ele", rename_to="WD_Points1-gen.ele", trim=True)
    smallWgdRun1.addComparisonFile("Points.1.node", "WD_Points1.node", rename_to="WD_Points1-gen.node", trim=True)
    test_info["smallWingDisc-pt1"] = smallWgdRun1

    # create the smallWingDisc (pt2, 3d meshing) information
    smallWgdRun2 = MeshGenRun("inputFile_smallWingDisc-2", dir_path + "/smallWingDisc")
    smallWgdRun2.addExeInput("WD_Points1.ele", "Points.1.ele")
    smallWgdRun2.addExeInput("WD_Points1.node", "Points.1.node")
    smallWgdRun2.addComparisonFile("MeshFile.out", "smallWingDisc.mesh", rename_to="smallWingDisc-gen.mesh")
    test_info["smallWingDisc-pt2"] = smallWgdRun2

    # location of the folder containing the mesh generation executable
    executable_loc = dir_path + "/../../ToolBox/MeshGeneration/2DEllipse"
    executable_loc = os.path.abspath(executable_loc) # makes debugging slightly nicer
    # name of the mesh generation executable
    executable_name = "EllipseFromOutline"

    # output files that are produced by the executable that need to be cleaned up every time
    always_cleanup = [ executable_loc + "/NodesPostTesselation.out", executable_loc + "/VectorsPostTesselation.out", executable_loc + "/Points.node", executable_loc + "/NodesPreTesselation.out", executable_loc + "/MeshFile.out" ]

    @pytest.mark.parametrize("tc", test_cases)
    def test_mesh_generation(self, tc):
        '''
        For the given test ID, stored in tc:
        - Copy the input files to the executable directory
        - Run the executable
        - Move the output files back to the testing directory
        - Compare the outputs to the reference files

        Parameters
        ----------
        tc : test ID, matching a value in test_cases
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
        all_passed = mesh_gen_test.compare_output_files()

        # finally, cleanup the generated output files that we created, provided NO checks came back negative
        if all_passed:
            mesh_gen_test.cleanup_output_files()
        else:
            print("Some checks failed, preserving outputs")
        return