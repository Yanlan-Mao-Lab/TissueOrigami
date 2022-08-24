## test_regression.py
## Created by: William Graham, 2022-08-04

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
        Intermediaries: WD_Points.1.node, WD_Points.1.ele
        '''

        # create the output directory if it doesn't already exist
        if (not os.path.exists(self.gen_location)):
            os.mkdir(self.gen_location)

        # reference output to compare to
        ref_output = self.ref_location + "/smallWingDisc.mesh"
        # check this file can be found, fail if not
        assert os.path.exists(ref_output), "Could not find reference file: " + ref_output

        # reference intermediary output files
        ref_im_nodes = self.ref_location + "/WD_Points.1.node"
        ref_im_ele = self.ref_location + "/WD_Points.1.ele"
        # check the files exist
        assert os.path.exists(ref_im_nodes), "Could not find reference file: " + ref_im_nodes
        assert os.path.exists(ref_im_ele), "Could not find reference file: " + ref_im_ele

        # where to place the generated output and intermediary files
        gen_output = self.gen_location + "/smallWingDisc.mesh"
        gen_im_nodes = self.gen_location + "/WD_Points.1.node"
        gen_im_ele = self.gen_location + "/WD_Points.1.ele"

        # expected names of the intermediary files
        raw_im_nodes = self.executable_loc + "/Points.1.node"
        raw_im_ele = self.executable_loc + "/Points.1.ele"

        # PHASE 1: Copy outline file and create 2D mesh

        # copy files across to target directories
        src_outline = self.input_loc + "/48hrDiscSymmetricOutline"
        # create the target directory if it doesn't already exist
        dst_outline = self.executable_loc + "/48hrDiscSymmetricOutline"
        shutil.copyfile(src_outline, dst_outline)

        # ./EllipseFromOutline 1 35.56 50 27.3 27.3 12.5 9 2 0  ./48hrDiscSymmetricOutline, insert 1 before the path to pass the TissueType
        command = "./" + self.executable_name + " 1 35.56 50 27.3 27.3 12.5 9 2 0 1 " + dst_outline
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # we should have produced the points.1.{ele, node} files
        # copy these to the output directory (to preserve them for comparison later)
        shutil.move(raw_im_nodes, gen_im_nodes)
        shutil.move(raw_im_ele, gen_im_ele)

        # now cleanup the files we produced and copied across
        cleanup([dst_outline, self.executable_loc + "/Points.node", self.executable_loc + "/NodesPreTesselation.out"])

        # PHASE 2: Extrude to 3 dimensions

        # copy across the reference input files to the locations they are expected to be at
        shutil.copyfile(ref_im_nodes, raw_im_nodes)
        shutil.copyfile(ref_im_ele, raw_im_ele)

        # ./EllipseFromOutline -1 5.2 1.0 3 1, extra 1 to pass the TissueType
        command = "./" + self.executable_name + " -1 5.2 1.0 3 1 1"
        # run the executable...
        subprocess.run(command.split(), cwd=self.executable_loc)
        # save the output
        shutil.move(self.raw_output_src, gen_output)
        
        # cleanup auxillary files
        cleanup([raw_im_ele, raw_im_nodes, self.executable_loc + "/NodesPostTesselation.out", self.executable_loc + "/VectorsPostTesselation.out"])

        # # PHASE 3: Compare intermediaries and output files

        assert filecmp.cmp(ref_output, gen_output, shallow=False), "Output meshfile mismatch between " + ref_output + " and " + gen_output

        # triangle-generated files append the ABSOLUTE path to the end of the file
        # this means that the final line of the reference intermediary- and the generated intermediary files will always be different, as the call to triangle uses a different path
        # the solution is to trim this final line from both files and compare the result
        # make temporary directory for this
        tmp_dir = self.dir_path + "/temp_im_comp"
        if (not os.path.exists(tmp_dir)):
            os.mkdir(tmp_dir)
        # for ease of making new files in the directory
        tmp_dir += "/"
        # copy files across
        shutil.copyfile(ref_im_nodes, tmp_dir + "ref.node")
        shutil.copyfile(ref_im_ele, tmp_dir + "ref.ele")
        shutil.copyfile(gen_im_nodes, tmp_dir + "gen.node")
        shutil.copyfile(gen_im_ele, tmp_dir + "gen.ele")
        # trim final lines
        trim_final_line([tmp_dir + "ref.node", tmp_dir + "ref.ele", tmp_dir + "gen.node", tmp_dir + "gen.ele"])
        # compare remaining parts of files, which should match
        assert filecmp.cmp(tmp_dir + "ref.node", tmp_dir + "gen.node", shallow=False), "Intermediary file mismatch between " + ref_im_nodes + " and " + gen_im_nodes
        assert filecmp.cmp(tmp_dir + "ref.ele", tmp_dir + "gen.ele", shallow=False), "Intermediary file mismatch between " + ref_im_ele + " and " + gen_im_ele
        # cleanup the temporary directory
        shutil.rmtree(tmp_dir)

        return