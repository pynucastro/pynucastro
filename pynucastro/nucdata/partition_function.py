"""
Class holding a partition function for a given nuclide.

The partition function is stored as a table of partition function
values evaluated at a sequence of temperatures.

The minimum and maximum valid temperatures for this partition function
are also stored.
"""

from scipy.interpolate import InterpolatedUnivariateSpline

class PartitionFunction(object):
    """
    Holds information for a single partition function table
    corresponding to a single nucleus.
    
    Units of temperature are Kelvins.
    """
    def __init__(self, name=None, temperature=None, partition_function=None):
        self.name = name
        self.temperature = temp
        self.partition_function = pfun
        self.interpolant = None
        self.interpolant_order = None

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
        tmin = lower.minimum_temperature
        tmax = upper.maximum_temperature
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
