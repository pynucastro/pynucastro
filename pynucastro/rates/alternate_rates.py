"""A collection of alternate rates that are not contained in the other
sources.

"""

import numpy as np

from pynucastro.nucdata import Nucleus
from pynucastro.rates.reaclib_rate import ReacLibRate, SingleSet
from pynucastro.rates.temperature_tabular_rate import TemperatureTabularRate


class DeBoerC12agO16(ReacLibRate):
    """A drop-in replacement for "C12(a,g)O16" from ReacLib that uses
    the formulation from DeBoer et al. 2017

    """

    def __init__(self):

        labelprops = "deboer"

        reactants = [Nucleus("c12"), Nucleus("he4")]
        products = [Nucleus("o16")]

        Q = 0
        for n in reactants:
            Q += -n.A * n.nucbind
        for n in products:
            Q += n.A * n.nucbind

        # we'll use the same chapter notation as ReacLib
        chapter = 4

        super().__init__(reactants=reactants, products=products,
                         chapter=chapter, Q=Q,
                         labelprops=labelprops)

        # from DeBoer et al. 2017
        # https://arxiv.org/pdf/1709.03144
        # table XXVI

        # non-resonant
        a = [24.1, 0.0, -32.0, -5.9, 1.8, -0.17, -2./3.]
        self.sets.append(SingleSet(a, labelprops=labelprops))

        # resonance
        a = [7.4, -30.0, 0.0, 0.0, 0.0, 0.0, -1.5]
        self.sets.append(SingleSet(a, labelprops=labelprops))


class IliadisO16pgF17(TemperatureTabularRate):
    """The O16(p,g)F17 rate from Iliadis et al. 2022.  This uses
    the median value for the rate from the table at the end of the
    paper.
    """

    def __init__(self):

        T9 = np.array([3.00e-03, 4.00e-03, 5.00e-03, 6.00e-03, 7.00e-03, 8.00e-03,
                       9.00e-03, 1.00e-02, 1.10e-02, 1.20e-02, 1.30e-02, 1.40e-02,
                       1.50e-02, 1.60e-02, 1.80e-02, 2.00e-02, 2.50e-02, 3.00e-02,
                       4.00e-02, 5.00e-02, 6.00e-02, 7.00e-02, 8.00e-02, 9.00e-02,
                       1.00e-01, 1.10e-01, 1.20e-01, 1.30e-01, 1.40e-01, 1.50e-01,
                       1.60e-01, 1.80e-01, 2.00e-01, 2.50e-01, 3.00e-01, 3.50e-01,
                       4.00e-01, 4.50e-01, 5.00e-01, 6.00e-01, 7.00e-01, 8.00e-01,
                       9.00e-01, 1.00e+00, 1.25e+00, 1.50e+00, 1.75e+00, 2.00e+00,
                       2.50e+00, 3.00e+00, 3.50e+00])

        rate = np.array([4.338e-41, 1.399e-36, 2.244e-33, 6.230e-31, 5.555e-29, 2.253e-27,
                         5.142e-26, 7.595e-25, 7.988e-24, 6.405e-23, 4.117e-22, 2.203e-21,
                         1.010e-20, 4.067e-20, 4.770e-19, 3.968e-18, 2.749e-16, 6.919e-15,
                         7.513e-13, 2.082e-11, 2.601e-10, 1.943e-09, 1.016e-08, 4.098e-08,
                         1.359e-07, 3.865e-07, 9.731e-07, 2.219e-06, 4.659e-06, 9.131e-06,
                         1.688e-05, 4.991e-05, 1.266e-04, 8.106e-04, 3.309e-03, 1.011e-02,
                         2.528e-02, 5.464e-02, 1.058e-01, 3.126e-01, 7.381e-01, 1.493e+00,
                         2.699e+00, 4.482e+00, 1.229e+01, 2.632e+01, 4.806e+01, 7.863e+01,
                         1.686e+02, 2.971e+02, 4.591e+02])

        reactants = [Nucleus("o16"), Nucleus("p")]
        products = [Nucleus("f17")]

        Q = 0
        for n in reactants:
            Q += -n.A * n.nucbind
        for n in products:
            Q += n.A * n.nucbind

        super().__init__(np.log10(T9), np.log10(rate),
                         reactants=reactants, products=products, Q=Q)
