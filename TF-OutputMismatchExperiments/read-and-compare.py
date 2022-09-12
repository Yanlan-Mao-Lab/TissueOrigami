import subprocess, os, filecmp

# where am I?
this_dir_path = os.path.dirname(os.path.abspath(__file__))
print("Working in: %s" % this_dir_path)
# prefixes to use in comparisons
prefixes = ["RF-", "GUI-"]
# suffix when read by reader.o
suffix = "-read.txt"
# name of the files to compare
fnames = ["Force", "Growth", "GrowthRate", "GrowthRedistribution", "Packing", "PhysicalProp", "SpecificElementAndNodeTypes", "TensionCompression"]
# name of the reader executables
executables = dict()
executables[fnames[0]] = "./growth_reader.o"
executables[fnames[1]] = "./growthrate_reader.o"

for i,p1 in enumerate(prefixes):
    for p2 in prefixes[i+1:]:
        # compare files pairwise
        for f in fnames:
            f1 = p1+f; f2 = p2+f
            f1read = "./" + f1 + suffix; f2read = "./" + f2 + suffix

            # read from binary to text so I can comprehend these things
            command = executables[f] + " ./" + f1 + " ./" + f2
            reading = subprocess.Popen(command.split(), cwd=this_dir_path)
            reading.wait()

            # compare f1read and f2read to find the differences...
            with open(f1read, 'r') as F1, open(f2read, 'r') as F2:
                lines1 = F1.readlines(); lines2 = F2.readlines()

                # if there are a different number of lines, we have different outputs!
                if len(lines1)!=len(lines2):
                    print("Uneven number of lines! %s (%d) vs %s (%d)" % (f1, len(lines1), f2, len(lines2)))
                else:
                    # compare line by line
                    mismatches = []
                    nLines = len(lines1)
                    for nl, l1 in enumerate(lines1):
                        l2 = lines2[nl]
                        if l1 != l2:
                            # lines do not match!
                            mismatches.append(nl)
                    # how many lines did not match?
                    if mismatches:
                        # if mismatches is non-empty, the if statement evaluates to true
                        print("Counted %d different lines, now printing:" % len(mismatches))
                        for n in mismatches:
                            print("[Line %d] =========" % n)
                            print("%s : %s" % (f1, lines1[n]))
                            print("%s : %s" % (f2, lines2[n]))
                            print("-------------------")
                    else:
                        # they're the same - does filecmp think they're the same?
                        print("Contents match. Running filecmp.cmp on these files returns: %d" % filecmp.cmp(f1read, f2read, shallow=False))

            # cleanup afterwards
            os.remove(f1read); os.remove(f2read)