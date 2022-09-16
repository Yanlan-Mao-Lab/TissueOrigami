# Regression Tests for Pardiso Build

**IMPORTANT NOTE:** We maintain a testing framework that runs on Linux.
We are aware of potential small (order machine-precision) discrepencies between the outputs of the `TissueFolding` when built on other operating systems, such as MacOS.

Building with PARDISO requires a license key that is tied to the username of the account running the tests, which makes running these tests through CI impossible.
Instead the tests need to be run locally, with the user:
- Following the build instructions on a machine that has a PARDISO license,
- Running the simulations with the reference inputs,
- Copying the output files to the respective `generatd_out` directories,
- Running `test_SimulationPardiso.py` via `pytest`.

Besides the outputs, the inputs to each simulation `run0700x` are the same as those to the non-PARDISO build.
The pass conditions for the tests are largely identical to the `sim_no_pardiso_reg` case, however exceptions are made due to the possibility of numerical imprecision when time-evolving the system.
For a given output `out` and its corresponding reference file `ref`, the following checks are performed:

If `python`'s `filecmp.cmp(ref, out)` method returns `True`, the files are identical and the test passes.
If `out` and `ref` are binary files, they are converted to `txt` files using the `simOutputsToTxt` executable, and compared as text files with `filecmp.cmp` again - the test passes if these files are identical.
Otherwise, the (now possibly converted to text) files `out` and `ref` are compared using the `CompareTxt` function:
- All strings within the files are directly compared against each other. If these differ, then the test _fails_.
- All saved numerical values in the files are checked to be within a set tolerance of each other. If any pair of values differs, then the test _fails_.