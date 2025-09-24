# Example Reduction

This is an example of performing a reduction on a large network.

## Generating a starting library

To start the example, we want to create a reduced library that has a
set endpoint, and eliminates duplicate rates (in particular, there are
some heavy nuclei for which 2 different rates are provided).

To generate the library, run the following command:

```bash
./rp_process.py te108 --write_lib rp-process-lib
```

This will put the reduced library in the default pynucastro
library directory (so you won't see it via `ls` in the directory
where you ran this command.

## Reducing

An example of reduction is provided by the `reduction_driver.py`
script.

The script takes a range of command line arguments -- run:

```bash
./reduction_driver.py --help
```

to see all of them. To run with a smaller network (*ni56* endpoint), a
small dataset (64 points), do

```bash
./reduction_driver.py -e ni56
```

To run with 2 MPI processes (assuming *mpi4py* is installed), do:

```bash
mpiexec -n 2 ./reduction_driver.py -e ni56 --use_mpi
```

To run with the full network (up to *te108*), omit the `-e ni56`.

The `--use_numpy` argument can be used to have it using a NumPy
version of the algorithm.
