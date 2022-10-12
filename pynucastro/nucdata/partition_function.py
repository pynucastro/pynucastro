import os

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline


class PartitionFunction:
    """
    The necessary class public variables of PartitionFunction(nucleus, name, temperature, partition_function)
    are characterized as follows:

    nucleus            : a string variable composed by an element and the atomic number, e.g ni56.
    name               : the name of the table on which the nucleus is read.
    temperature        : a list with all the temperatures involved in the table named in the previous variable.
    partition_function : a list with all the partition values given in the same order of the previos list.
    interpolant        : stores the interpolant function.
    interpolant_order  : stores the interpolation spline order.

    The public methods of this class are:

    lower_partition()   : returns the lowest temperature value of the partition_function list.
    upper_partition()   : returns the highest temperature value of the partition_function list.
    lower_temperature() : returns the lowest value temperature value of the temperature list.
    upper_temperature() : returns the lowest value temperature value of the temperature list.
    construct_spline_interpolant(order) : interpolates temperature vs log(partition_function), using the
    spline interpolation of order=3 by default, returning the function of T.

    The dunder methods of this class are

    __add__ : if two partition functions do not overlap their temperatures, we define the addition at the incorporation
          of all the temperatures and their partition function values, respectively.
    __call__: This object allow us to treat the class object as function of T, returning the appropiate value of
          the partition function.

    The purpose of this oject is to encompass all the nucleus partition function values into a single object, on which +
    is defined.
    """

    def __init__(self, nucleus, name, temperature, partition_function):

        assert isinstance(nucleus, str)

        self.nucleus = nucleus
        self.name = name
        self.temperature = temperature
        self.partition_function = partition_function
        self.interpolant = None
        self.interpolant_order = None

        if (isinstance(temperature, np.ndarray) and
            isinstance(partition_function, np.ndarray) and
            len(temperature) == len(partition_function)):
            self.construct_spline_interpolant()
        else:
            self.interpolant_order = 0
            self.interpolant = lambda x: 0.0

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

        assert self.nucleus == other.nucleus

        if self.upper_temperature() < other.lower_temperature():
            lower = self
            upper = other
        else:
            lower = other
            upper = self
        if lower.upper_temperature() >= upper.lower_temperature():
            raise ValueError("temperature ranges cannot overlap")

        temperature = np.array(list(lower.temperature) +
                               list(upper.temperature))

        partition_function = np.array(list(lower.partition_function) +
                             list(upper.partition_function))

        name = f'{lower.name}+{upper.name}'

        newpf = PartitionFunction(nucleus=self.nucleus, name=name,
                                  temperature=temperature, partition_function=partition_function)

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

    def __eq__(self, other):

        return (np.all(self.partition_function == other.partition_function) and
                np.all(self.temperature == other.temperature))

    def construct_spline_interpolant(self, order=3):
        """
        Construct an interpolating univariate spline of order >= 1 and
        order <= 5 using the scipy InterpolatedUnivariateSpline
        implementation.

        Interpolate in log space for the partition function and in GK
        for temperature.
        """

        self.interpolant = InterpolatedUnivariateSpline(self.temperature/1.0e9,
                                                        np.log10(self.partition_function),
                                                        k=order)

        self.interpolant_order = order

    def __call__(self, T):

        assert self.interpolant
        try:
            T = float(T)/1.0e9
        except ValueError:
            print("invalid temperature")
            raise
        else:
            if self.interpolant_order == 0:
                return 10**self.interpolant(T)
            return 10**self.interpolant(T, ext='const')  # extrapolates keeping the boundaries fixed.


class PartitionFunctionTable:
    """
    Class for reading a partition function table file. A
    PartitionFunction object is constructed for each nucleus and
    stored in a dictionary keyed by the lowercase nucleus name in the
    form, e.g. "ni56". The table files are stored in the PartitionFunction
    sub directory.

    The class PartitionFunctionTable(file_name) is characterized by the public variable self.name,
    which stores the name of the table. The private variable self._partition_function collects all
    the tables we have previously converted bu using their scripts.

    """

    def __init__(self, file_name):
        self._partition_function = {}
        self.name = None
        self.temperatures = None
        self._read_table(file_name)

    def _add_nuclide_pfun(self, nuc, pfun):
        assert isinstance(nuc, str)
        assert nuc not in self._partition_function
        self._partition_function[nuc] = pfun

    def get_nuclei(self):
        return list(self._partition_function)

    def get_partition_function(self, nuc):
        assert isinstance(nuc, str)
        if nuc in self._partition_function:
            return self._partition_function[nuc]
        return None

    def _read_table(self, file_name):
        with open(file_name, 'r') as fin:

            # get headers name
            fhead = fin.readline()
            hsplit = fhead.split('name: ')
            self.name = hsplit[-1].strip('\n')

            # throw away the six subsequent lines
            for _ in range(6):
                fin.readline()

            # Now, we want to read the lines of the file where
            # the temperatures are located
            temp_strings = fin.readline().strip().split()
            self.temperatures = np.array([float(t) for t in temp_strings])

            # Now, we append on the array lines = [] all the remaining file, the structure
            # 1. The nucleus
            # 2. The partition value relative to the nucleus defined in 1.

            lines = []
            for line in fin:
                ls = line.strip()
                if ls:
                    lines.append(ls)

        # Using .pop(0) twice we construct each nucleus partition function.
        while lines:
            nuc = lines.pop(0)
            pfun_strings = lines.pop(0).split()
            partitionfun = np.array([float(pf) for pf in pfun_strings])
            pfun = PartitionFunction(nuc, self.name, self.temperatures, partitionfun)
            self._add_nuclide_pfun(nuc, pfun)


class PartitionFunctionCollection:

    """ The PartitionFunctionCollection holds a collection of PartitionFunctionTable objects in a dictionary keyed
    by the name of the tables

    In our discussion we have two different set of tables"""

    def __init__(self, use_high_temperatures=True, use_set='frdm'):
        self._partition_function_tables = {}
        self.use_high_temperatures = use_high_temperatures
        self.use_set = use_set
        self._read_collection()

    def _add_table(self, table):
        """
        This private function appends a PartitionFunctionTable object to each key characterized by a file_name.
        """
        assert table.name not in self._partition_function_tables
        self._partition_function_tables[table.name] = table

    def _read_collection(self):

        """
        This private function construct the whole collection of tables
        """

        nucdata_dir = os.path.dirname(os.path.realpath(__file__))
        partition_function_dir = os.path.join(nucdata_dir, 'PartitionFunction')

        pft = PartitionFunctionTable(file_name=os.path.join(partition_function_dir, 'etfsiq_low.txt'))
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=os.path.join(partition_function_dir, 'frdm_low.txt'))
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=os.path.join(partition_function_dir, 'etfsiq_high.txt'))
        self._add_table(pft)

        pft = PartitionFunctionTable(file_name=os.path.join(partition_function_dir, 'frdm_high.txt'))
        self._add_table(pft)

    def get_nuclei(self):

        nuclei = set()
        for table in self._partition_function_tables.values():
            nuclei.update(table.get_nuclei())
        return nuclei

    def __iter__(self):
        for nuc in self.get_nuclei():
            yield self.get_partition_function(nuc)

    def get_partition_function(self, nuc):

        """This function access to the partition function for a given nucleus"""
        assert isinstance(nuc, str)

        if self.use_set == 'frdm':
            pf_lo_table = self._partition_function_tables['frdm_low']
            pf_lo = pf_lo_table.get_partition_function(nuc)

            pf_hi_table = self._partition_function_tables['frdm_high']
            pf_hi = pf_hi_table.get_partition_function(nuc)

        elif self.use_set == 'etfsiq':
            pf_lo_table = self._partition_function_tables['etfsiq_low']
            pf_lo = pf_lo_table.get_partition_function(nuc)

            pf_hi_table = self._partition_function_tables['etfsiq_high']
            pf_hi = pf_hi_table.get_partition_function(nuc)
        else:
            raise ValueError("invalid partition function type")

        if self.use_high_temperatures:
            if pf_lo and pf_hi:
                pf = pf_lo + pf_hi
            elif pf_lo:
                pf = pf_lo
            elif pf_hi:
                pf = pf_hi
            else:
                #name = 'default'
                #pf_default = PartitionFunction(nuc, name, pf_lo_table.temperatures, np.ones_like(pf_lo_table.temperatures))
                #pf = pf_default
                raise ValueError
                #if str(nuc) != 'h1' and str(nuc) != 'n' and str(nuc) != 'he4':
                #    print(f'WARNING: {nuc} partition function is not supported: set pf = 1.0 by default')

        else:
            if pf_lo:
                pf = pf_lo
            else:
                raise ValueError
                #name = 'default'
                #pf_default = PartitionFunction(nuc, name, pf_lo_table.temperatures, np.ones_like(pf_lo_table.temperatures))
                #pf = pf_default
                #if str(nuc) != 'p' and str(nuc) != 'n' and str(nuc) != 'he4':
                #    print(f'WARNING: {nuc} partition function is not supported: set pf = 1.0 by default')

        return pf
