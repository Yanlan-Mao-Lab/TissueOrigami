## test_regression.py
## Created by: William Graham, 2022-08-04

from cgi import test
import os, filecmp, shutil, subprocess

class Test_MeshGeneration():
    '''

    Attributes
    ----------

    Tests
    -----

    '''

    # path to reappend in order to find files to compare
    dir_path = os.path.dirname(os.path.abspath(__file__))

    # location of the folder containing the reference outputs
    ref_location = dir_path + "/ref_outputs"
    # location of the folder containing the inputs to mesh generation executable
    input_loc = dir_path + "/mesh_inputs"
    # location of the folder to place the generated outputs into
    gen_location = dir_path + "/gen_putputs"

    # location of the folder containing the mesh generation executable
    executable_loc = dir_path + "/../../ToolBox/MeshGeneration/2DEllipse"
    # name of the mesh generation executable
    executable_name = "EllipseFromOutline"

    # at present, the executable just dumps MeshFile.out to the directory
    raw_output_src = executable_loc + "/MeshFile.out"           

    def test_smallRectangle(self):
        '''
        
        '''

        # reference output to compare to
        ref_output = self.ref_location + "/smallRectangle.mesh"
        # check this file can be found, fail if not
        assert os.path.exists(ref_output), "Could not find reference file: " + ref_output
        # where to place the generated output
        gen_output = self.gen_location + "/smallRectangle.mesh"

        # copy input files to executable directory
        src_nodes = self.input_loc + "/smallRectangle.node"
        src_ele = self.input_loc + "/smallRectangle.ele"
        dst_nodes = self.executable_loc + "/Points.1.nodes"
        dst_ele = self.executable_loc + "/Points.1.ele"
        shutil.copyfile(src_nodes, dst_nodes)
        shutil.copyfile(src_ele, dst_ele)

        # ./EllipseFromOutline -1 5.2 2 3 0,
        command = "./" + self.executable_name + " -1 5.2 2 3 0"
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # the output should then be moved (and renamed) to gen_output, in case we wish to inspect it later
        shutil.move(self.raw_output_src, gen_output)

        # now compare the contents of the generated output and the reference output
        assert filecmp.cmp(ref_output, gen_output, shallow=False), "Output meshfile mismatch between " + ref_output + " and " + gen_output

        # clean up copied files
        os.remove(dst_nodes)
        os.remove(dst_ele)