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

class Network_sundials(Network_f90):
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
        self.sundials_dir = os.path.join(self.pyreaclib_dir,
                                    'templates',
                                    'sundials-cvode')
        self.template_file_select = os.path.join(self.sundials_dir,
                                            '*.template')
        self.template_files = glob.glob(self.template_file_select)

        # Initialize values specific to this network
        self.name_rate_data = 'screened_rates'
        self.name_reactvec = 'reactvec'
        self.name_y         = 'Y'
        self.name_ydot      = 'YDOT'
        self.name_ydot_nuc      = 'YDOT'
        self.name_jacobian  = 'DJAC'
        self.name_jacobian_nuc  = 'DJAC'

    def enuc_dqweak(self, n_indent, of):
        # Add tabular dQ corrections to the energy generation rate
        for nr, r in enumerate(self.rates):
            if nr in self.tabular_rates:
                if len(r.reactants) != 1:
                    print('ERROR: Unknown tabular dQ corrections for a reaction where the number of reactants is not 1.')
                    exit()
                else:
                    reactant = r.reactants[0]
                    of.write('{}{}(net_ienuc) = {}(net_ienuc) + N_AVO * {}(j{}) * {}(i_dqweak, k_{})\n'.format(self.indent*n_indent, self.name_ydot, self.name_ydot, self.name_ydot, reactant, self.name_reactvec, r.fname))
        
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
                    of.write('{}{}(net_ienuc) = {}(net_ienuc) + N_AVO * {}(j{}) * {}(i_epart, k_{})\n'.format(self.indent*n_indent, self.name_ydot, self.name_ydot, self.name_y, reactant, self.name_reactvec, r.fname))
