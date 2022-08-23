import os, filecmp, shutil, subprocess

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

    # location of the executable (relative to this file)
    exe_loc = "../../TissueFolding/"
    # name of the executable
    exe_name = "TissueFolding"

    def test_run07007():
        '''
        
        '''

        # create a list of the files that this run is dependant on
        files_needed = [ "smallRectangle.mesh", \
                        "Stiffness96hrRectangleWing_Reduction_0", \
                        "Stiffness96hrRectangleWing_Reduction_4", \
                        "ShapeChangeRate96hrRectangleWingZ_Reduction_3", \
                        "ShapeChangeRate96hrRectangleWingXY_Reduction_3"
        ]
