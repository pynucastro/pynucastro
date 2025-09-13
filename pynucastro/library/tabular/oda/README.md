# oda

This contains the rates from Oda, T., Hino, M., Muto, K., Takahara,
M., & Sato, K. (1994)

These cover nuclei masses 17 <= A <= 39, and a thermodynamic range:

* 0.01 < T_9 < 100
* 1 < log(rho Y_e) ?< 11

The actual data table was downloaded from the NSCL Charge-Exchange Group
https://groups.frib.msu.edu/charge_exchange/weakrate_tables/

## preparing the data

There are 2 python files which manage the data:

* `extract.py` : this splits the original `odarates.dat` file into
  separate files---one for each parent daughter pair.

* `weak.py` : this is run after `extract.py`loops over the
  tables of parent-daughter rates and outputs 2 files in the format
  that pynucastro usesm one for its effective beta decay (positron capture
  and beta-decay), and electron-capture (beta+ decay and
  electron capture).

  Note that by construction, this will combine the electron-capture
  and beta+-decay into a single rate and likewise, the e+-capture and
  beta-decay into a single rate.

  Everything is also put into CGS units.
