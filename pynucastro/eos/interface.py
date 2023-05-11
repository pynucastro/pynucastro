class EosState:

    def __init__(self, rho, T, ye):
        self.rho = rho
        self.T = T
        self.ye = ye

        self.p = 0.0
        self.s = 0.0
        self.eint = 0.0

        self.dpdr = 0.0
        self.eta = 0.0
        self.xn = 0.0