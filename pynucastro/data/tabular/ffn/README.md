# ffn

This contains the rates from Fuller, Fowler, and Newman (1982)
https://ui.adsabs.harvard.edu/abs/1982ApJ...252..715F/abstract

These cover nuclei masses 21 <= A <= 60, and a thermodynamic range:

* 0.01 < T_9 < 100
* 1 < log(rho Y_e) ?< 11

The actual data table was downloaded from the NSCL Charge-Exchange Group
https://groups.frib.msu.edu/charge_exchange/weakrates.html

## preparing the data

There are 2 python files which manage the data:

* `extract.py` : this splits the original `ffnrates.dat` file into
  separate files---one for each parent daughter pair.

* `beta_rates.py` : this is run after `extract.py`loops over the
  tables of parent-daughter rates and outputs 4 files in the format
  that pynucastro uses (a "toki" files and the data file each for
  electron-capture and beta-decay.

  Note that by construction, this will combine the electron-capture
  and beta+-decay into a single rate and likewise, the e+-capture and
  beta-decay into a single rate.

  Everything is also put into CGS units.
  
