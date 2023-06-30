# Reduction Utilities

Included are a few algorithms for doing reaction network reduction.

## Running Example Script

Running the example script requires supplying a library to load the network from. To generate the
default library, one can navigate to `../../examples/rp-process/` and run the following command:
```bash
./rp_process.py te108 --write_lib rp-process-lib
```
The script takes a range of command line arguments -- run `./reduction.py --help` to see all of
them. To run with a smaller network (*ni56* endpoint), a small dataset (64 points), and no MPI, do:
```bash
./reduction.py -e ni56
```
To run with 2 MPI processes (assuming *mpi4py* is installed), do:
```bash
mpiexec -n 2 ./reduction.py -e ni56 --use_mpi
```
To run with the full network (up to *te108*), omit the `-e ni56`. The `--use_numpy` argument tells
it to use the NumPy version of the algorithm.
