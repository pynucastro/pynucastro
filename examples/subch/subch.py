# C-burning with A=23 URCA rate module generator

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180228default2"
mylibrary = pyna.rates.Library(library_file)

subCh = pyna.rates.Library()

all_reactants = [("He4", "He4", "He4"),
                 ("C12", "He4"),
                 ("N14", "He4"),
                 ("F18", "He4"),
                 ("C12", "p"),
                 ("N13", "He4"),
                 ("O16", "Ne20")]

print(mylibrary)

for r in all_reactants:
    react = [pyna.rates.Nucleus(q) for q in r]
    print(react)
    rfilter = pyna.rates.RateFilter(reactants=react)
    print(rfilter)
    _library = mylibrary.filter(rfilter)
    print(_library)
    print("Adding")
    subCh += _library

print(subCh)




