# Atomic Mass Evaluation

The script `extract_table.py` along with the classes in `ame_table.py`
and `ame_nuclide.py` can be used to read the Atomic Mass Evaluation
tables.

For example, to generate the `binding_2016.txt` table, download the
AME 2016 table file `mass16.txt` and run:

```
./extract_table.py mass16.txt -o binding_2016.txt
```

The 2020 Atomic Mass Evaluations were retrieved from
https://www-nds.iaea.org/amdc/

The Atomic Mass Evaluation 2020 is published as  Chinese Phys. C 45 (2021) 030002, 030003

The Nubase evaluation is published as Chinese Physics C 45 (2021) 030001


## Scripts

* `extract_spin_1.py` : this extracts the spins from the Nubase version of the data,
  currently using `nubase_4.mas20`.

* `extract_mass_excess.py` : this extracts the mass excess from the Nubase version of
  the data, currently using `nubase_4.mas20`.

