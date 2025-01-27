import numpy
from pynucastro.eos import brentq_eta
from pynucastro.constants import constants

class ElectronPositronTable:

    def __init__(self, rho_bounds, temp_bounds, ye, rho_samples_per_decade=20, temp_samples_per_decade=20):

        rho_l, rho_h = rho_bounds
        temp_l, temp_h = temp_bounds

        self.rho_points = self._log_per_decade(rho_l, rho_h, rho_samples_per_decade)
        self.temp_points = self._log_per_decade(temp_l, temp_h, temp_samples_per_decade)
        self.y_e = y_e
        self.rho_e = self.rho * self.ye

        low_eta = -34
        max_eta = 25

        self.eta_points = []
        for rho in rho_points:
            for temp in t_points:
                beta = k*T/(m_e*c**2)
                eta = brentq_eta(self.rho, self.beta, low_eta, max_eta, max_iter=150)
                self.eta_points.append(eta)

    def _log_per_decade(x_init, x_end, samples_per_decade):
        ndecades = np.log10(x_end) - np.log10(x_init)
        npoints = int(ndecades) * samples_per_decade
        x_points = np.logspace(np.log10(x_init), np.log10(x_end), num=npoints, endpoint=True, base=10)
        return x_points

    def pressure(self):
        pass

    def entropy(self):
        pass

    def free_energy(self):
        pass
