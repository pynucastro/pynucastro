# Common Imports
from __future__ import print_function

import os
import re
import glob

from pynucastro.networks import BaseFortranNetwork

class SundialsNetwork(BaseFortranNetwork):
    def __init__(self, *args, **kwargs):
        # Initialize BaseFortranNetwork parent class
        super(SundialsNetwork, self).__init__(*args, **kwargs)
        
        # Set up some directories
        self.sundials_dir = os.path.join(self.pynucastro_dir,
                                    'templates',
                                    'sundials-cvode')
        self.template_file_select = os.path.join(self.sundials_dir,
                                            '*.template')
        self.template_files = glob.glob(self.template_file_select)

        # Tags specific to this network
        self.ftags['<net_ynuc>'] = self.ynuc
        self.ftags['<cvodeneq>'] = self.cvodeneq
        self.ftags['<net_ymass_init>'] = self.net_ymass_init

        # Initialize values specific to this network
        self.name_rate_data = 'screened_rates'
        self.name_reactvec = 'reactvec'
        self.name_y = 'Y'
        self.name_ydot = 'YDOT'
        self.name_ydot_nuc = 'YDOT'
        self.name_jacobian = 'DJAC'
        self.name_jacobian_nuc = 'DJAC'

    def headerline(self, n_indent, of):
        of.write('{}write(2, fmt=hfmt) '.format(self.indent*n_indent))
        for nuc in self.unique_nuclei:
            of.write("'Y_{}', ".format(nuc))
        of.write("'E_nuc', 'Time'\n")

    def ynuc(self, n_indent, of):
        for nuc in self.unique_nuclei:
            of.write('{}double precision :: y{}\n'.format(
                self.indent*n_indent, nuc))

    def yinit_nuc(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write("{}cv_data%Y0(j{})   = net_initial_abundances%y{}\n".format(
                self.indent*n_indent, n, n))

    def final_net_print(self, n_indent, of):
        of.write('{}write(*,*) "MASS FRACTIONS:"\n'.format(self.indent*n_indent))
        for n in self.unique_nuclei:
            of.write("{}write(*,'(A,ES25.14)') '{}: ', cv_data%Y(j{})*aion(j{})\n".format(self.indent*n_indent, n, n, n))

    def cvodeneq(self, n_indent, of):
        of.write('{} '.format(self.indent*n_indent))
        of.write('integer*8 :: NEQ = {} ! Size of ODE system\n'.format(
            len(self.unique_nuclei)+1))

    def net_ymass_init(self, n_indent, of):
        for n in self.unique_nuclei:
            of.write('{}net_initial_abundances%y{} = 0.0d0\n'.format(
                self.indent*n_indent, n))
