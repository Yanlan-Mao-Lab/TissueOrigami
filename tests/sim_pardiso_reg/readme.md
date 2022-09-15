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
The pass conditions for the tests, and the contextual setup for each, is the same - please see the readme file in the `sim_no_pardiso_reg` for details.