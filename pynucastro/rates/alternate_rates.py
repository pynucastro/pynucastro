from pynucastro.nucdata import Nucleus
from pynucastro.rates import ReacLibRate, SingleSet


class DeBoerC12agO16(ReacLibRate):

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
