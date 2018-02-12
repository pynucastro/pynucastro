# Network Examples

Included are a few examples using pynucastro to generate reaction
networks in Python and Fortran.

For Fortran code generation, pynucastro offers two options.

- Standalone: Generate right hand side routines and an integration
  driver program using the included VODE integrator in
  `util/VODE/`. These networks implement only species and energy
  integration at constant temperature, as pynucastro does not include
  an equation of state.

- Microphysics: Generate network files to copy directly into the
  StarKiller Microphysics repository for use with astrophysical
  simulation codes. These networks implement species, temperature, and
  energy integration.

## All networks

The first few steps to setting up any of the available network options
are as follows. See the following section for details about the
particular kind of network you are interested in creating.

* Create a working directory under `examples/`

* Obtain Reaclib rate files (in Reaclib 1 or 2 formats) for your problem and
  either place them in your working directory or `pynucastro/library`.

* Write a short python script to generate your network,
  e.g. `mynet.py`, based on the following example scripts:
    - For Python networks, `CNO/cno.py`
    - For standalone Fortran networks, `urca-23_f90/urca.py`
    - For StarKiller Microphysics networks, `urca-23_starkiller/urca.py`

* Run your python script (pynucastro must be in your PYTHONPATH).

```
$ python mynet.py
```

## Python network

To integrate a Python network using the SciPy integration routines,
customize `CNO/burn.py` to initialize and run your network using the
right hand side module you generated above.

## Standalone Fortran network

The `urca-23_f90` example builds a Fortran 90 network together with a
GNU Makefile and integration driver program using the included VODE
package for ODE integration.

* In the steps for all networks above, pynucastro will create several
  Fortran 90 files as well as a GNU Makefile set up to compile the
  integrator using gfortran. So next do:

```
$ make
```

* Edit the generated `inputs_integration` file to specify your initial
  conditions as well as integration tolerances, stop time, and output
  sampling interval.

* Run the integrator as:

```
$ main.Linux.gfortran.exe inputs_integration
```

* The integration results will be stored in a text file named by
  default, `integration_history.dat`

Notes on the build system:

* If you would like to compile with debugging symbols, do:

```
$ make NDEBUG=
```

* To remove the executable and object files produced during compilation, do:

```
$ make realclean
```

## StarKiller Microphysics network

The `urca-23_starkiller` example builds the right hand side, jacobian,
and helper Fortran modules to copy into the `networks/` subdirectory
of the StarKiller Microphysics repository.

No additional customization is required after running the steps for
all networks above.

## Tabular Rates

Tabular rates for reactions of the form `A -> B` are supported by the
standalone Fortran and StarKiller Microphysics network outputs. If you
would like to include tabular rates, for now they must be in the form
of, eg. `23Na-23Ne_electroncapture.dat` in `urca-23_f90/`, indexed by
density*ye and temperature, with the same number and order of
variables. Prepare a rate file to describe how to read the table in
the form of, eg. `urca-23_f90/na23--ne23-toki`:

```
t
       [parent nuclide]  [daughter nuclide]
[rate table file name]
[number of header lines to eat before the first line of data]
[number of density*ye values]
[number of temperature values]
```
