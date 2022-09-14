import os, glob, subprocess

#os.path.dirname(os.path.abspath(__file__))

non_pardiso_code_dir = "/home/ccaegra/Documents/TissueOrigami/TissueFolding/SourceCode"
pardiso_code_dir = "/home/ccaegra/Documents/TF-SecondRepoForCodeCompare/TissueOrigami/Scratch/TissueFolding/SourceCode"
log_file_dir = "/home/ccaegra/Documents/TissueOrigami"

# get .cpp files to compare
np_cpp_files = glob.glob(non_pardiso_code_dir + "/*.cpp")
for i,f in enumerate(np_cpp_files):
    np_cpp_files[i] = os.path.basename(f)
p_cpp_files = glob.glob(pardiso_code_dir + "/*.cpp")
for i,f in enumerate(p_cpp_files):
    p_cpp_files[i] = os.path.basename(f)
# validate: should have same number of files in each here!
if len(np_cpp_files)!=len(p_cpp_files):
    print("non-pardiso cpp files (%d):\n" % len(np_cpp_files))
    print(np_cpp_files)
    print("pardiso cpp files (%d):\n" % len(p_cpp_files))
    print(p_cpp_files)   
    raise ValueError("Number of .cpp source files is different!")

# get header files to compare
np_h_files = glob.glob(non_pardiso_code_dir + "/*.h")
for i,f in enumerate(np_h_files):
    np_h_files[i] = os.path.basename(f)
p_h_files = glob.glob(pardiso_code_dir + "/*.h")
for i,f in enumerate(p_h_files):
    p_h_files[i] = os.path.basename(f)
# validate: should have same number of file in each here!
if len(np_h_files)!=len(p_h_files):
    print("non-pardiso h files (%d):\n" % len(np_h_files))
    print(np_h_files)
    print("pardiso h files (%d):\n" % len(p_h_files))
    print(p_h_files) 
    raise ValueError("Number of /h files is different!")

print("I found %d .h files and %d .cpp files" % (len(np_cpp_files), len(np_h_files)))

# run diff on each, and pipe the result to an output file
output_file = open(log_file_dir + "/difference-log.txt", 'a')
for cpp_file in np_cpp_files:
    # files to run diff on
    npf = non_pardiso_code_dir + "/" + cpp_file
    pf = pardiso_code_dir + "/" + cpp_file
    # note in the log that these files are being X-examined
    output_file.write('\n=======! Comparison of: ' + cpp_file + " !=======\n")
    output_file.flush()
    # run diff on them
    command = "diff " + npf + " " + pf
    run_diff = subprocess.Popen(command.split(), stdout=output_file, stderr=output_file)
    run_diff.wait()

for h_file in np_h_files:
    # files to run diff on
    npf = non_pardiso_code_dir + "/" + h_file
    pf = pardiso_code_dir + "/" + h_file
    # note in the log that these files are being X-examined
    output_file.write('\n=======! Comparison of: ' + h_file + " !=======\n")
    output_file.flush()
    # run diff on them
    command = "diff " + npf + " " + pf
    run_diff = subprocess.Popen(command.split(), stdout=output_file, stderr=output_file)
    run_diff.wait()

# close log file we were using
output_file.close()