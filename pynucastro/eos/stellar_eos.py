"""Classes a methods for a full stellar equation of state (ions +
electrons + radiation)

"""


class EOSState:
    """A container to hold a thermodynamic state that depends on
    density, temperature, and composition

    """

    def __init__(self,
                 eta=0.0,
                 n=0.0, p=0.0, e=0.0,
                 dn_drho=0.0, dn_dT=0.0,
                 dp_drho=0.0, dp_dT=0.0,
                 de_drho=0.0, de_dT=0.0):

        self.eta = eta

        self.n = n
        self.p = p
        self.e = e

        self.dn_drho = dn_drho
        self.dn_dT = dn_dT

        self.dp_drho = dp_drho
        self.dp_dT = dp_dT

        self.de_drho = de_drho
        self.de_dT = de_dT
