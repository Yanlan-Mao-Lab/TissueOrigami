## test_regression.py
## Created by: William Graham, 2022-08-04

from cgi import test
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
    gen_location = dir_path + "/gen_outputs"

    # location of the folder containing the mesh generation executable
    executable_loc = dir_path + "/../../ToolBox/MeshGeneration/2DEllipse"
    # name of the mesh generation executable
    executable_name = "EllipseFromOutline"

    # at present, the executable just dumps MeshFile.out to the directory
    raw_output_src = executable_loc + "/MeshFile.out"           

    def test_smallRectangle(self):
        '''
        Regression test for the smallRectangle mesh.
        Reference output: smallRectangle.mesh
        '''

        # create the output directory if it doesn't already exist
        if (not os.path.exists(self.gen_location)):
            os.mkdir(self.gen_location)

        # reference output to compare to
        ref_output = self.ref_location + "/smallRectangle.mesh"
        # check this file can be found, fail if not
        assert os.path.exists(ref_output), "Could not find reference file: " + ref_output
        # where to place the generated output
        gen_output = self.gen_location + "/smallRectangle.mesh"

        # copy input files to executable directory
        src_nodes = self.input_loc + "/smallRectangle.node"
        src_ele = self.input_loc + "/smallRectangle.ele"
        dst_nodes = self.executable_loc + "/Points.1.node"
        dst_ele = self.executable_loc + "/Points.1.ele"
        shutil.copyfile(src_nodes, dst_nodes)
        shutil.copyfile(src_ele, dst_ele)

        # ./EllipseFromOutline -1 5.2 2 3 0, additional 1 for TissueType input allowance
        command = "./" + self.executable_name + " -1 5.2 2 3 0 1"
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # the output should then be moved (and renamed) to gen_output, in case we wish to inspect it later
        shutil.move(self.raw_output_src, gen_output)

        # now compare the contents of the generated output and the reference output
        assert filecmp.cmp(ref_output, gen_output, shallow=False), "Output meshfile mismatch between " + ref_output + " and " + gen_output

        # clean up auxillary and copied files
        cleanup([dst_nodes, dst_ele, self.executable_loc + "/NodesPostTesselation.out", self.executable_loc + "/VectorsPostTesselation.out"])

        return

    def test_smallSphere(self):
        '''
        Regression test for the smallSphere mesh.
        Reference output: smallSphere.mesh
        '''

        # create the output directory if it doesn't already exist
        if (not os.path.exists(self.gen_location)):
            os.mkdir(self.gen_location)

        # reference output to compare to
        ref_output = self.ref_location + "/smallSphere.mesh"
        # check this file can be found, fail if not
        assert os.path.exists(ref_output), "Could not find reference file: " + ref_output
        # where to place the generated output
        gen_output = self.gen_location + "/smallSphere.mesh"

        # copy input files to executable directory
        src_tri = self.input_loc + "/SphericalTriangulation.txt"
        dst_tri = self.executable_loc + "/SphericalTriangulation"
        shutil.copyfile(src_tri, dst_tri)

        # ./EllipseFromOutline -2 4.2 2 3 0, additional 5 for TissueType input allowance
        command = "./" + self.executable_name + " -2 4.2 2 3 0 5"
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # the output should then be moved (and renamed) to gen_output, in case we wish to inspect it later
        shutil.move(self.raw_output_src, gen_output)

        # now compare the contents of the generated output and the reference output
        assert filecmp.cmp(ref_output, gen_output, shallow=False), "Output meshfile mismatch between " + ref_output + " and " + gen_output

        # clean up auxillary and copied files
        cleanup([dst_tri, self.executable_loc + "/NodesPostTesselation.out", self.executable_loc + "/VectorsPostTesselation.out"])

        return

    def test_smallWingDisc(self):
        '''
        Regression test for the smallWingDisc mesh.
        Reference output: smallWingDisc.mesh
        '''

        # create the output directory if it doesn't already exist
        if (not os.path.exists(self.gen_location)):
            os.mkdir(self.gen_location)

        # reference output to compare to
        ref_output = self.ref_location + "/smallWingDisc.mesh"
        # check this file can be found, fail if not
        assert os.path.exists(ref_output), "Could not find reference file: " + ref_output
        # where to place the generated output
        gen_output = self.gen_location + "/smallWingDisc.mesh"

        # copy input files to executable directory
        src_tri = self.input_loc + "/SphericalTriangulation.txt"
        dst_tri = self.executable_loc + "/SphericalTriangulation"
        shutil.copyfile(src_tri, dst_tri)

        # ./EllipseFromOutline -2 4.2 2 3 0, additional 5 for TissueType input allowance
        command = "./" + self.executable_name + " -2 4.2 2 3 0 5"
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # the output should then be moved (and renamed) to gen_output, in case we wish to inspect it later
        shutil.move(self.raw_output_src, gen_output)

        # now compare the contents of the generated output and the reference output
        assert filecmp.cmp(ref_output, gen_output, shallow=False), "Output meshfile mismatch between " + ref_output + " and " + gen_output

        # clean up auxillary and copied files
        cleanup([dst_tri, self.executable_loc + "/NodesPostTesselation.out", self.executable_loc + "/VectorsPostTesselation.out"])

        return