The URCA_f90 example integrates a reaction network which implements
simple C12-burning along with the A=23 URCA reactions until either
23Ne or 23Na are depleted.

To run:

1) Generate the reaction network using urca.py (pyreaclib directory
must be in your PATH):

$ python urca.py

2) Edit net_rates.f90 and fill in the binding energies (MeV) in the
init_net_info subroutine: (relative to proton BE)
    
    self%ebind_per_nucleon(self%in)   = 0.782347d0
    self%ebind_per_nucleon(self%ip)   = 0.0d0
    self%ebind_per_nucleon(self%ihe4)   = 7.073915d0
    self%ebind_per_nucleon(self%ic12)   = 7.680144d0
    self%ebind_per_nucleon(self%ine20)   = 8.032240d0
    self%ebind_per_nucleon(self%ine23)   = 7.955255d0
    self%ebind_per_nucleon(self%ina23)   = 8.111493d0
    self%ebind_per_nucleon(self%img23)   = 7.901104d0

3) Add your SUNDIALS libraries to the Makefile, eg.

SUNINCDIR = /home/eugene/local/sundials/instdir/include
SUNLIBDIR = /home/eugene/local/sundials/instdir/lib

4) Make and run

$ make
$ ./intnet.exe


