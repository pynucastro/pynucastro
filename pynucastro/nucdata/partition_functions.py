"""
Classes for reading a partition function table and storing
partition functions in PartitionFunction objects.

Supported tables at present are the following:

(RT2000) Thomas Rauscher and Friedrich-Karl Thielemann, Atomic Data and Nuclear
Data Tables, 75:1–351, 2000

(R2003) Thomas Rauscher, Astrophysical Journal Supplement Series, 147:403-408,
2003

For R2003, both Table 2 (FRDM) and Table 3 (ETFSI-Q) are supported. By
default, pynucastro uses the FRDM values from R2003 for consistency
with RT2000.
"""

import os
import numpy as np
import pynucastro
from scipy.interpolate import InterpolatedUnivariateSpline

class PartitionFunction(object):
    """
    Holds information for a single partition function table
    corresponding to a single nucleus.
    
    Units of temperature are Kelvins.
    """
    def __init__(self, name=None, temperature=None, partition_function=None):
        self.name = name
        self.temperature = temperature
        self.partition_function = partition_function
        self.interpolant = None
        self.interpolant_order = None
        if (type(self.temperature) == np.ndarray and
            type(self.partition_function) == np.ndarray and
            len(self.temperature) == len(self.partition_function)):
            self.construct_spline_interpolant()

    def lower_partition(self):
        return self.partition_function[0]

    def upper_partition(self):
        return self.partition_function[-1]

    def lower_temperature(self):
        return self.temperature[0]

    def upper_temperature(self):
        return self.temperature[-1]

    def __add__(self, other):
        """
        Adding two PartitionFunction objects is implemented by simply
        appending the temperature and partition function values of the
        higher-temperature partition function to those of the lower-temperature
        partition function. If the temperature ranges overlap,
        however, an exception is generated.

        If either the PartitionFunction objects added have already had
        a spline interpolant constructed, then construct a spline
        interpolant for the returned PartitionFunction of order equal
        to the maximum order of the added PartitionFunction objects.
        """
        assert(self.upper_temperature() < other.lower_temperature() or
               self.lower_temperature() > other.upper_temperature())
        if self.upper_temperature() < other.lower_temperature():
            lower = self
            upper = other
        else:
            lower = other
            upper = self
        temperature = np.array(list(lower.temperature) +
                               list(upper.temperature))
        partition_function = np.array(list(lower.partition_function) +
                                      list(upper.partition_function))
        name = '{} + {}'.format(lower.name, upper.name)
        newpf = PartitionFunction(name=name, temperature=temperature,
                                  partition_function=partition_function)

        if self.interpolant_order and other.interpolant_order:
            order = max(self.interpolant_order, other.interpolant_order)
        elif self.interpolant_order:
            order = self.interpolant_order
        elif other.interpolant_order:
            order = other.interpolant_order
        else:
            order = None
        if order:
            newpf.construct_spline_interpolant(order=order)

        return newpf

    def construct_spline_interpolant(self, order=3):
        """
        Construct an interpolating univariate spline of order >= 1 and
        order <= 5 using the scipy InterpolatedUnivariateSpline
        implementation. 
        """
        self.interpolant = InterpolatedUnivariateSpline(self.temperature,
                                                        self.partition_function,
                                                        k=order)
        self.interpolant_order = order

class PartitionFunctionTable(object):
    """ 
    Class for reading a partition function table file. A
    PartitionFunction object is constructed for each nucleus and
    stored in a dictionary keyed by the lowercase nucleus name in the
    form, e.g. "ni56". 
    """
    def __init__(self, file_name):
        self._partition_functions = {}
        self.name = None
        self.read_table(file_name)

    def _add_nuclide_pfun(self, nuc, pfun):
        """
        Given a Nucleus object nuc and partition function pfun,
        check to make sure nuc does not already have a partition function
        and if not then add it.
        """
        assert(not nuc in self._partition_functions)
        self._partition_functions[str(nuc)] = pfun

    def get_partition_function(self, nuc):
        """
        Given a nucleus object or string, return its partition function or None
        if it is not in the table.
        """
        if str(nuc) in self._partition_functions.keys():
            return self._partition_functions[str(nuc)]
        else:
            return None

    def read_table(self, file_name):
        fin = open(file_name, 'r')

        # Get name from header
        fhead = fin.readline()
        hsplit = fhead.strip().split('name: ')
        self.name = hsplit[-1]

        # Throw away 6 lines
        for _ in range(6):
            fin.readline()

        # Get the temperatures
        temp_strings = fin.readline().strip().split()
        temperatures = np.array([float(t) for t in temp_strings])

        # Read the rest of the lines
        lines = []
        for line in fin:
            ls = line.strip()
            if ls:
                lines.append(ls)
        fin.close()

        # Get partition functions for each nucleus
        while lines:
            nuc = pynucastro.rates.Nucleus(lines.pop(0))
            pfun_strings = lines.pop(0).split()
            partitionfun = np.array([float(pf) for pf in pfun_strings])
            pfun = PartitionFunction(name=self.name,
                                     temperature=temperatures,
                                     partition_function=partitionfun)
            self._add_nuclide_pfun(nuc, pfun)

class PartitionFunctionCollection(object):
    """ The PartitionFunctionCollection holds a set of partition function tables. """
    def __init__(self):
        self.partition_function_tables = {}
        self._read_collection()

    def _add_table(self, table):
        """ Checks to ensure table isn't already in the collection and adds it. """
        assert(not table.name in self.partition_function_tables)
        self.partition_function_tables[table.name] = table

    def _read_collection(self):
        """ Read the partition function tables included in pynucastro. """
        nucdata_dir = os.path.dirname(os.path.realpath(__file__))
        partition_function_dir = os.path.join(nucdata_dir, 'PartitionFunctions')

        pft = PartitionFunctionTable(os.path.join(partition_function_dir,
                                                  'partition_functions_rathpf.txt'))
        self._add_table(pft)

        pft = PartitionFunctionTable(os.path.join(partition_function_dir,
                                                  'partition_functions_rauscher2003_frdm.txt'))
        self._add_table(pft)

        pft = PartitionFunctionTable(os.path.join(partition_function_dir,
                                                  'partition_functions_rauscher2003_etfsiq.txt'))
        self._add_table(pft)

    def get_nuc_partition_function(self, nuc, high_temperature_partition_functions="rauscher2003_FRDM"):
        """
        Given a Nucleus object nuc or string representation, return its partition function. 
        If no partition function is located for this nuclide, return None.

        The argument high_temperature_partition_functions may be
        supplied with one of two options: "rauscher2003_ETFSIQ" or
        "rauscher2003_FRDM". The latter is the default.

        These correspond to partition functions from:

        Thomas Rauscher and Friedrich-Karl Thielemann, Atomic Data and Nuclear
        Data Tables, 75:1–351, 2000 (RT 2000)

        together with the high temperature partition functions from:

        Thomas Rauscher, Astrophysical Journal Supplement Series, 147:403-408,
        2003. (R 2003)

        The option high_temperature_partition_functions="rauscher2003_FRDM" 
        will combine the partition functions from RT 2000 with the partition 
        functions from Table 2 of R 2003.

        The option high_temperature_partition_functions="rauscher2003_ETFSIQ" 
        will combine the partition functions from RT 2000 with the partition 
        functions from Table 3 of R 2003.
        """
        assert(type(nuc) == pynucastro.rates.Nucleus or type(nuc) == str)

        pf_lo_temp = self.partition_function_tables['rathpf']
        pf_hi_temp = self.partition_function_tables[high_temperature_partition_functions]

        pflo = pf_lo_temp.get_partition_function(nuc)
        pfhi = pf_hi_temp.get_partition_function(nuc)

        pf = None
        if pflo and pfhi:
            pf = pflo + pfhi
        elif pflo:
            pf = pflo
        elif pfhi:
            pf = pfhi
        return pf
