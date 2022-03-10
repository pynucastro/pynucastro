import pynucastro as pyna

library_file = "20180319default2"

reaclib_library = pyna.rates.Library(library_file)

#all_reactants = ["n", "p", "o16", "he4", "c12", "ne20", "si28", "s31", "s32", "p31"]

all_reactants = ["n", "p",
                 "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                 "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                 "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                 "c14", "n13", "n14", "o18", "f18", "ne21", "mg23", "na23", "si27", "s31"]

reduced_library = reaclib_library.linking_nuclei(all_reactants)

reduced_library.validate(reaclib_library)

print(reduced_library)




