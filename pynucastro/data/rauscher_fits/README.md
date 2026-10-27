# Rauscher SMARAGD Rate Fits

These files contain the Rauscher 2025 (version 2) fits to stellar reactivities in the REACLIB seven parameter formats. Calculations can be found here https://doi.org/10.1140/epja/s10050-026-01809-4.

Data can be found at https://doi.org/10.5281/zenodo.18075461.

`fits_p.dat`: Fits to astrophysical reactivities for proton-induced reactions

`fits_a.dat`: Fits to astrophysical reactivities for alpha-induced reactions

`fits_n.dat`: Fits to astrophysical reactivities for neutron-induced reactions

# File Format

Each reaction entry is made of three lines.

Line 1: Identification line with the corresponding reaction.

Line 2: 7 parameters for forward rate.

Line 3: 7 parameters for reverse rate.

**Note: reverse reactivity obtained from these parameters has to be multiplied by the ratio of the normalised partition functions of target and final nuclide.**

