# CNO reactions

This example reads in the JINA reaclib rates for the main CNO cycle
and a portion of the hot CNO cycle (proton capture onto N13) and
outputs a python module that has each of the rates in its own
function as well as a RHS routine for the system of ODEs (in terms
of molar masses, Y).

To use:

 * run `cno.py` to generate the `cno_rhs.py` routine -- this requires
   that `reaclib.py` be in your path.

 * run `burn.py` -- this will use the `scipy` ODE routines to
   integrate the system.
   

