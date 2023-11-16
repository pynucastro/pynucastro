from pynucastro.networks.rate_collection import RateCollection



class EosState:

    def __init__(self, rho, T, composition):
        self.rho = rho
        self.T = T
        self.ye = composition.eval_ye()
        self.abar = composition.eval_abar()

        self.p = 0.0
        self.dpdr = 0.0
        self.dpdT = 0.0

        self.s = 0.0
        self.dsdr = 0.0
        self.dsdT = 0.0

        self.e = 0.0
        self.dedr = 0.0
        self.dedT = 0.0

        self.eta = 0.0
        self.xn = 0.0



# def initialize(self):

#     state = EosState

#     self.p = 0.0
#     self.dpdr = 0.0
#     self.dpdt = 0.0

#     self.s = 0.0
#     self.dsdr = 0.0
#     self.dsdT = 0.0

#     self.e = 0.0
#     self.dedr = 0.0
#     self.dedT = 0.0

#     self.eta = 0.0
#     self.xn = 0.0

#     self.abar = 0.0

#     return