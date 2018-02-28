# C-burning with A=23 URCA rate module generator

import pynucastro as pyna

library_file = "20180201ReaclibV2.22"
mylibrary = pyna.rates.Library(library_file)

subCh = pyna.rates.Library()

reactants = [['he4', 'he4', 'he4'],
             ['c12', 'he4']]

for r in reactants:
    rf = pyna.rates.RateFilter(reactants=r)
    rsc = mylibrary.filter(rf)
    print(rsc)
    subCh += rsc

print('constructed the library:')
print(subCh)
