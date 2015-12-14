# pyreaclib

Python interfaces to the JINA reaclib nuclear reactions database.  The
main goal is to provide a parser for the reaclib format to produce
either a callable python function for a rate or output the python code
for a function that can be incorporated into a rate module.

For performance reasons, the common temperature factors needed by a
rate are stored in a `Tfactors` class---a `Tfactors` object is passed
into each rate function to get the rate.


# example

```
r = Rate("c13-pg-n14-nacr")

# output python code for this rate
print r.function_string()

# evaluate this rate at T = 1.e9 K
print r.eval(1.0e9)

```
