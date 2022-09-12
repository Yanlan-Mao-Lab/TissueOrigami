import filecmp

# ReFerence, MYriad, GEnerated (locally), 
# NargessK 2nd local run, DocKerfile w/ qmake build, DockerFile w/ CMake build
# Dockerfile GUI build
prefixes = ["RF-", "NK2-", "GE-", "MY-", "DK-", "DFCM-", "GUI-"]
N = len(prefixes)
files = ["Save_Growth", "Save_GrowthRate"]

matches = 0; differences = 0
exp_comparisons = N*(N-1) # (N(N-1)/2 * 2 files per comparison)

for i,p1 in enumerate(prefixes):
    for p2 in prefixes[i+1:]:
        # for good measure, let's also check if files are the same as themselves... what could go wrong?
        for f in files:
            if (not filecmp.cmp(p1+f, p2+f, shallow=False)):
                # these files were different from each other!
                print("Difference: \t" + p1 + f + " , " + p2 + f)
                differences += 1
            else:
                # these files were the same
                print("Match: \t" + p1 + f + " , " + p2 + f)
                matches += 1

print("COMPARE BINARIES: Found %d differences and %d matches (expected %d total comparisons)" % (differences, matches, exp_comparisons))