import os, filecmp, shutil, subprocess
import numpy as np
import pytest
from pytest_check import check
import warnings

def CompareTxtElements(e1, e2, line_num, ele_num, f_name1, f_name2, tol=1e-8):
    '''
    Compares two strings to determine whether they are identical, under the following assumptions:
    - If the strings are equal via string comparison, return True
    - Otherwise, attempt to interpret the strings are doubles. Return True if these doubles are within tol of each other.
    If the above conditions are not satisfied, the elements are deemed to be different and False is returned.
    '''

    if e1==e2:
        # either these are the same text fields, or the same number has been written out
        return True
    else:
        # more drastic measures - can be convert to doubles?
        try:
            double1 = np.double(e1)
            double2 = np.double(e2)
        except ValueError:
            # could not convert strings to doubles, different files
            print('String difference at line %d, element %d:' % (line_num, ele_num))
            print('\t (%s) %s' % (f_name1, e1))
            print('\t (%s) %s' % (f_name2, e2))
            return False
        except:
            # unexpected error thrown, print error information
            raise RuntimeError("Error: unexpected error whilst handling string -> double conversion")
        # otherwise, these were actually numbers, is the difference less than tol?
        double_diff = np.abs( double1 - double2 )
        # are the elements different by a significant amount?
        if double_diff >= tol:
            print('Line %d, element %d, values [ %.8e , %.8e ] have difference greater than %.1e' % (line_num, ele_num, double1, double2), tol)
            return False
        else:
            # upon conversion to doubles, elements match
            return True
    # we always return in the above if statement, no need for concluding return here

def CompareLinesInRawTxt(line1, line2, line_num, f_name1, f_name2, tol=1e-8):
    '''
    Compare the contents of two lines from text files.
    The lines are split into elements, which are separated by whitespace.
    Elements are deemed equal if:
        - They are both strings, and these strings are identical in length and content
        - They are both doubles, and these doubles are within tol of each other
    If all corresponding pairs of elements are equal, the lines are equal.

    The lines passed in should be matching lines from a TEXT output of the simulation, and the corresponding reference file.
    '''
    values1 = list(filter(None, [x.strip() for x in line1.split()])); n_v1 = len(values1)
    values2 = list(filter(None, [x.strip() for x in line2.split()])); n_v2 = len(values1)
    if n_v1 != n_v2:
        # different number of fields on this line, not the same file
        return False
    else:
        # compare field-by-field
        # compare as strings first (as these might be text characters)
        # otherwise, try to convert to double and compare - catch errors if these really are strings that can't be converted to doubles, and are genuinely different
        for ele_num in range(n_v1):
            e1 = values1[ele_num]; e2 = values2[ele_num]
            # compare the two elements; either as strings, or failing this as doubles saved as strings
            # CompareTxtElements returns False if the elements are different
            if not CompareTxtElements(e1, e2, line_num, ele_num, f_name1, f_name2, tol):
                return False
    # if we don't escape early, the lines are the same as all elements are equal
    return True

def CompareLinesInConvertedTxt(line1, line2, line_num, tol=1e-8):
    '''
    Compare the contents of two lines of (comma-separated) values.
    Values are compared element-wise, and are deemed to be equal when they coincide to a difference of tol.

    The lines passed in should be matching lines from a BINARY output of the simulation, converted to a text file via the converter executable, and the corresponding line from (the converted) reference file.
    '''

    values1 = list(filter(None, [x.strip() for x in line1.split(sep=',')])); n_v1 = len(values1)
    values2 = list(filter(None, [x.strip() for x in line2.split(sep=',')])); n_v2 = len(values2)

    if n_v1 != n_v2:
        # different number of fields on this line, not the same file
        return False
    else:
        # compare field-by-field, these are doubles so we should be able to cast
        for ele_num in range(n_v1):
            double1 = np.double(values1[ele_num])
            double2 = np.double(values2[ele_num])
            double_diff = np.abs( double1 - double2 )
            # are the elements different by a significant amount?
            if double_diff >= tol:
                print('Line %d, element %d, values [ %.8e , %.8e ] have difference greater than %.1e' % (line_num, ele_num, np.double(values1[ele_num]), np.double(values2[ele_num]), tol))
                return False
    # if we didn't return early, all elements must have the same value, so the lines are equal
    return True

def CompareTxt(f_name1, f_name2, tol=1e-8, type='raw_txt'):
    '''
    Compares the contents of the files f_name1 and f_name2, returning True if the contents match according to the following criteria:
    - all non-numerical text must match exactly between files
    - numerical values must agree to within tol

    For converted binaries, we can just skip straight to comparing the numerical values, and this can be toggled on by passing type as anything EXCEPT its default value.
    '''

    # diagnostics
    print('CompareTxt running on %s vs %s' % (f_name1, f_name2))
    # open files, get the data, then close them
    f1 = open(f_name1, 'r'); 
    f2 = open(f_name2, 'r'); 
    lines1 = f1.readlines(); f1.close()
    lines2 = f2.readlines(); f2.close()

    n_lines_f1 = len(lines1)
    n_lines_f2 = len(lines2)
    # immediately flag a difference if number of lines is different
    if n_lines_f1!=n_lines_f2:
        return False
    elif type=='raw_txt':
        # comparing raw text files spat out by the program
        # these might have odd formating patterns, so we need to be careful as we compare
        # compare line-by-line
        for line_num in range(n_lines_f1):
            # CompareLinesInRawTxt will return True if the files are the same
            # Thus, we return False if this function also returned False
            line1 = lines1[line_num].strip('\n'); line2 = lines2[line_num].strip('\n')
            if not CompareLinesInRawTxt(line1, line2, line_num, f_name1, f_name2, tol):
                return False
    else:
        # comparing converted binaries, which are just .csv files
        # compare field-by-field
        for line_num in range(n_lines_f1):
            line1 = lines1[line_num]; line2 = lines2[line_num]
            # CompareLinesInConvertedTxt will return True if the files are the same
            # Thus, we return False if this function also returned False
            if not CompareLinesInConvertedTxt(lines1[line_num], lines2[line_num], line_num, tol):
                return False
    # if we get to here, we did not exit in the above comparisons, so the files must be the same
    return True

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

    # the value x in the test runs; run0700x. Tests iterate over these values
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

    # tolerance to accept differences in doubles to
    tolerance = 1e-8

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
            dif_compare = False
            # if bin_compare is false, the binaries do not match,
            # but if we read them into the simulation, do they match now?
            if not bin_compare:
                # raise a warning that we are having to binary convert!
                warnings.warn("Warning: having to examine via CompareTxt on: %s" % f_type)
                # run the converter on these output files, the mode to pass to the converter is given by f_type
                command = self.converter_exe + " " + f_type + " " + ref_op + " " + gen_op
                subprocess.run(command.split(), cwd=self.dir_path)
                # compare the converted outputs
                txt_compare = filecmp.cmp(ref_op+"-read.txt", gen_op+"-read.txt", shallow=False)
                if not txt_compare:
                    # need to actually read in the values now
                    dif_compare = CompareTxt(ref_op+"-read.txt", gen_op+"-read.txt", type="bin_compare", tol=self.tolerance)
                # cleanup the converted files that we made
                os.remove(ref_op+"-read.txt")
                os.remove(gen_op+"-read.txt")
            # test passes IFF bin_compare = True or (bin_compare = False and txt_compare = True)
            # since we do not assign to txt_compare unless bin_compare is false, we can simply check
            # bin_compare || txt_compare
            with check:
                assert (bin_compare or txt_compare or dif_compare), "Output mismatch between: " + ref_op + " and " + gen_op
        return