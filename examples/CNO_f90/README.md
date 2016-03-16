The CNO_f90 example uses the CVODE implementation of SUNDIALS to
integrate the CNO reaction network corresponding to the python-based
CNO example. The integration is carried out until H is depleted.

To run:

1) Generate the reaction network using cno.py (pyreaclib directory
must be in your PATH):

$ python cno.py

2) Edit net_rates.f90 and fill in the binding energies (MeV) in the
init_net_info subroutine:
    
    self%ebind_per_nucleon(self%ip) = 0.0d0
    self%ebind_per_nucleon(self%ihe4) = 7.073915d0
    self%ebind_per_nucleon(self%ic12) = 7.680144d0
    self%ebind_per_nucleon(self%ic13) = 7.469849d0
    self%ebind_per_nucleon(self%in13) = 7.238863d0
    self%ebind_per_nucleon(self%in14) = 7.475614d0
    self%ebind_per_nucleon(self%in15) = 7.699460d0
    self%ebind_per_nucleon(self%io14) = 7.052301d0
    self%ebind_per_nucleon(self%io15) = 7.46369d0

3) Add your SUNDIALS libraries to the Makefile, eg.

SUNINCDIR = /home/eugene/local/sundials/instdir/include
SUNLIBDIR = /home/eugene/local/sundials/instdir/lib

4) Make and run

$ make
$ ./intnet.exe
