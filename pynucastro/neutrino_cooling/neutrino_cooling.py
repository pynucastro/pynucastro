"""Classes and methods for exploring thermal neutrinos."""

import matplotlib.pyplot as plt
import numpy as np

from pynucastro.neutrino_cooling.sneut5_mod import sneut5


class NeutrinoCooling:
    """A class to provide an interface to explore neutrino cooling.
    This calls a specific implementation of the cooling that includes
    contributions from pairs, plasma, recombination, bremsstrahlung,
    and photoneutrinos.

    """

    def __init__(self, neutrino_function=sneut5):

        self.func = neutrino_function

    def plot(self, *, Tmin=1.e7, Tmax=1.e10, rhomin=1.e3, rhomax=1.e10,
             abar=20, zbar=10, npts_temp=50, npts_rho=50):
        """Plot the cooling term in as a function of density and
        temperature given a fixed abar/zbar.

        """

        Ts = np.logspace(np.log10(Tmin), np.log10(Tmax), npts_temp)
        rhos = np.logspace(np.log10(rhomin), np.log10(rhomax), npts_rho)

        data = np.zeros((npts_temp, npts_rho))

        for i, T in enumerate(Ts):
            for j, rho in enumerate(rhos):
                data[i, j] = np.log10(self.func(rho, T, abar=abar, zbar=zbar))

        fig, ax = plt.subplots()

        # by default, imshow puts the first index on the y axis
        # and the second on the x axis, so label accordingly

        extent = (np.log10(rhomin), np.log10(rhomax),
                  np.log10(Tmin), np.log10(Tmax))

        im = ax.imshow(data, cmap="viridis",
                       extent=extent, aspect="auto", origin="lower")
        fig.colorbar(im, ax=ax)

        ax.set_ylabel(r"$\log(T)$ [K]")
        ax.set_xlabel(r"$\log(\rho)$ [g/cm$^3$]")

        return fig
