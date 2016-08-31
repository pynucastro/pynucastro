# Common Imports
from __future__ import print_function

import glob
import os
import sys
import shutil
import re
import sympy

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from periodictable import elements

# Import Network_f90
from pyreaclib.networks.net_fortran import Network_f90

class Network_boxlib(Network_f90):
    def __init__(self, parent_instance_object=None):
        # Inherit all the instance attributes of
        # the parent_instance_object if it's passed.
        rem = re.compile('__(.*)__')
        if parent_instance_object:
            for d in dir(parent_instance_object):
                if not rem.match(d):
                    setattr(self, d, getattr(parent_instance_object, d))

        # Initialize Network_f90 stuff
        Network_f90.__init__(self)
                    
        # Set up some directories
        self.boxlib_dir = os.path.join(self.pyreaclib_dir,
                                       'templates',
                                       'boxlib-microphysics')
        self.template_file_select = os.path.join(self.boxlib_dir,
                                                 '*.template')
        self.template_files = glob.glob(self.template_file_select)

        # Initialize values specific to this network
        self.name_rate_data = 'screened_rates'
        self.name_y         = 'Y'
        self.name_ydot      = 'state%ydot'
        self.name_ydot_nuc  = 'ydot_nuc'
        self.name_jacobian  = 'state%jac'
        self.name_jacobian_nuc  = 'dfdy_nuc'

    def compute_tabular_rates_rhs(self, n_indent, of):
        if len(self.tabular_rates) > 0:
            of.write('{}! Included only if there are tabular rates\n'.format(self.indent*n_indent))
            of.write('{}do i = 1, nrat_tabular\n'.format(self.indent*n_indent))
            of.write('{}call tabular_evaluate(table_meta(i), rhoy, temp, reactvec)\n'.format(self.indent*(n_indent+1)))
            of.write('{}j = i + nrat_reaclib\n'.format(self.indent*(n_indent+1)))
            of.write('{}unscreened_rates(:,j) = reactvec(1:4)\n'.format(self.indent*(n_indent+1)))
            of.write('{}dqweak(i) = reactvec(5)\n'.format(self.indent*(n_indent+1)))
            of.write('{}epart(i)  = reactvec(6)\n'.format(self.indent*(n_indent+1)))
            of.write('{}end do\n'.format(self.indent*n_indent))

    def compute_tabular_rates_jac(self, n_indent, of):
        if len(self.tabular_rates) > 0:
            of.write('{}! Included only if there are tabular rates\n'.format(self.indent*n_indent))
            of.write('{}do i = 1, nrat_tabular\n'.format(self.indent*n_indent))
            of.write('{}call tabular_evaluate(table_meta(i), rhoy, temp, reactvec)\n'.format(self.indent*(n_indent+1)))
            of.write('{}j = i + nrat_reaclib\n'.format(self.indent*(n_indent+1)))
            of.write('{}unscreened_rates(:,j) = reactvec(1:4)\n'.format(self.indent*(n_indent+1)))
            of.write('{}end do\n'.format(self.indent*n_indent))

    def enuc_dqweak(self, n_indent, of):
        # Add tabular dQ corrections to the energy generation rate
        for nr, r in enumerate(self.rates):
            if nr in self.tabular_rates:
                if len(r.reactants) != 1:
                    print('ERROR: Unknown tabular dQ corrections for a reaction where the number of reactants is not 1.')
                    exit()
                else:
                    reactant = r.reactants[0]
                    of.write('{}enuc = enuc + N_AVO * {}(j{}) * dqweak(j_{})\n'.format(self.indent*n_indent, self.name_ydot, reactant, r.fname))
        
    def enuc_epart(self, n_indent, of):
        # Add particle energy generation rates (gamma heating and neutrino loss from decays)
        # to the energy generation rate (doesn't include plasma neutrino losses)
        for nr, r in enumerate(self.rates):
            if nr in self.tabular_rates:
                if len(r.reactants) != 1:
                    print('ERROR: Unknown particle energy corrections for a reaction where the number of reactants is not 1.')
                    exit()
                else:
                    reactant = r.reactants[0]
                    of.write('{}enuc = enuc + N_AVO * {}(j{}) * epart(j_{})\n'.format(self.indent*n_indent, self.name_y, reactant, r.fname))
