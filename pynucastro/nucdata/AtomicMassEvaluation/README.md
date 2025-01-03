# Atomic Mass Evaluation / Nubase

The Nubase / 2020 Atomic Mass Evaluations data were downloaded from:
https://www-nds.iaea.org/amdc/

The Atomic Mass Evaluation 2020 is published as Chinese Phys. C 45 (2021) 030002, 030003

The Nubase evaluation is published as Chinese Physics C 45 (2021) 030001

Note: these tables use the CODATA 2018 evaluation, so we need to use
those constants when converting between mass excess and binding energy
for consistency.

## Scripts

* `extract_spin.py` : this extracts the spins from the Nubase version
  of the data.  This outputs `spins2020.txt`.

* `extract_mass_excess.py` : this extracts the mass excess from the
  Nubase version of the data.  This outputs `mass_excess2020.txt`.

* `extract_halflife.py` : this extracts the half life from the Nubase
  version of the data.  This outputs `halflife2020.txt`.

For each of these, you provide the name of the table as an argument, e.g.:

```
python extract_mass_excess.py nubase_4.mas20.txt
```
