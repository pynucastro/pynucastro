# Atomic Mass Evaluation

The script `extract_table.py` along with the classes in `ame_table.py`
and `ame_nuclide.py` can be used to read the Atomic Mass Evaluation
tables.

For example, to generate the `binding_2016.txt` table, download the
AME 2016 table file `mass16.txt` and run:

```
./extract_table.py mass16.txt -o binding_2016.txt
```

The 2012 and 2016 Atomic Mass Evaluations were retrieved from
https://www-nds.iaea.org/amdc.

The Atomic Mass Evaluation 2016 is published as Chinese Physics C 41 (2017) 030002, 030003.

The Atomic Mass Evaluation 2012 is published as Chinese Physics C 36 (2012) 1287-1602.
