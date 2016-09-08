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

# Import pyreaclib modules
import pyreaclib.util
import pyreaclib.amemass
import pyreaclib.networks
import pyreaclib.rates


def make_network(file_list, net_type, use_cse=False):
    rc = pyreaclib.networks.RateCollection(file_list, use_cse)
    rc.make_network(net_type)

