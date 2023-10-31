import copy
import warnings

import numpy as np
from scipy import constants
from scipy.optimize import fsolve

from pynucastro._version import version
from pynucastro.networks.rate_collection import Composition, RateCollection
from pynucastro.nucdata import Nucleus
from pynucastro.rates import TabularRate
from pynucastro.screening import NseState, potekhin_1998


class NSETableEntry:
    def __init__(self, rho, T, Ye, *,
                 comp=None, ydots=None, enu=None,
                 comp_reduction_func=None):
        """a simple container to hold a single entry in the NSE table.

        Here, comp_reduction_func(comp) is a function that converts
        the NSE composition into a smaller set of nuclei.  It takes a
        Composition object and returns a dictionary with the nucleus
        name (like "Ni56") as the key and the corresponding mass fraction
        as the value.  It should be ordered in the way you want the nuclei
        output into the NSE table file.

        """

        self.rho = rho
        self.T = T
        self.Ye = Ye

        self.comp = comp
        self.ydots = ydots
        self.enu = enu

        # compute the bits we need for the table

        if comp:
            # mean molecular weight of the full NSE state
            self.abar = comp.eval_abar()

            # average binding energy / nucleon for the full NSE state
            self.bea = sum(q.nucbind * self.comp.X[q] for q in self.comp.X)

            # evolution of the electron fraction from weak rates alone
            self.dYedt = sum(q.Z * self.ydots[q] for q in self.comp.X)

            # evolution of Abar from weak rates alone
            self.dabardt = -self.abar**2 * sum(self.ydots[q] for q in self.comp.X)

            self.X = None
            if comp_reduction_func:
                self.X = comp_reduction_func(self.comp)

    def __str__(self):
        return f"({self.rho:12.6g}, {self.T:12.6g}, {self.Ye:6.4f}): {self.abar:6.4f}  {self.bea:6.4f}  {self.dYedt:12.6g}  {self.enu:12.6g}"

    def value(self):
        """a simple integer used for sorting.  This has the format
        (logrho)(logT)(1-Ye)"""

        return int(f"{np.log10(self.rho):08.5f}{np.log10(self.T):08.5f}{1-self.Ye:08.5f}".replace(".", ""))

    def __lt__(self, other):
        return self.value() < other.value()


class NSENetwork(RateCollection):
    """a network for solving for the NSE composition and outputting
    tabulated NSE quantities"""

    def _evaluate_n_e(self, state, Xs):

        m_u = constants.value("unified atomic mass unit") * 1.0e3  # atomic unit mass in g

        n_e = 0.0
        for nuc in self.unique_nuclei:
            n_e += nuc.Z * state.dens * Xs[nuc] / (nuc.A * m_u)

        return n_e

    def _evaluate_mu_c(self, n_e, state, use_coulomb_corr=True):
        """ A helper equation that finds the mass fraction of each nuclide in NSE state,
        u[0] is chemical potential of proton  while u[1] is chemical potential of neutron"""

        # If there is proton included in network, upper limit of ye is 1
        # And if neutron is included in network, lower limit of ye is 0.
        # However, there are networks where either of them are included
        # So here I add a general check to find the upper and lower limit of ye
        # so the input doesn't go outside of the scope and the solver won't be able to converge if it did
        ye_low = min(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        ye_max = max(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        assert state.ye >= ye_low and state.ye <= ye_max, "input electron fraction goes outside of scope for current network"

        # Setting up the constants needed to compute mu_c
        k = constants.value("Boltzmann constant") * 1.0e7          # boltzmann in erg/K
        Erg2MeV = 624151.0

        # u_c is the coulomb correction term for NSE
        # Calculate the composition at NSE, equations found in appendix of Calder paper
        u_c = {}

        for nuc in set(list(self.unique_nuclei) + [Nucleus("p")]):
            if use_coulomb_corr:
                # These are three constants for calculating coulomb corrections of chemical energy, see Calders paper: iopscience 510709, appendix
                A_1 = -0.9052
                A_2 = 0.6322
                A_3 = -0.5 * np.sqrt(3.0) - A_1 / np.sqrt(A_2)

                Gamma = state.gamma_e_fac * n_e ** (1.0/3.0) * nuc.Z ** (5.0 / 3.0) / state.temp
                u_c[nuc] = Erg2MeV * k * state.temp * (A_1 * (np.sqrt(Gamma * (A_2 + Gamma)) - A_2 * np.log(np.sqrt(Gamma / A_2) +
                                      np.sqrt(1.0 + Gamma / A_2))) + 2.0 * A_3 * (np.sqrt(Gamma) - np.arctan(np.sqrt(Gamma))))
            else:
                u_c[nuc] = 0.0

        return u_c

    def _nucleon_fraction_nse(self, u, u_c, state):

        # Define constants: amu, boltzmann, planck, and electron charge
        m_u = constants.value("unified atomic mass unit") * constants.kilo  # atomic unit mass in g
        k = constants.value("Boltzmann constant") / constants.erg           # boltzmann in erg/K
        h = constants.value("Planck constant") / constants.erg              # in cgs
        Erg2MeV = constants.erg / (constants.eV * constants.mega)

        Xs = {}
        up_c = u_c[Nucleus("p")]

        for nuc in self.unique_nuclei:
            if nuc.partition_function:
                pf = nuc.partition_function.eval(state.temp)
            else:
                pf = 1.0

            if not nuc.spin_states:
                raise ValueError(f"The spin of {nuc} is not implemented for now.")

            nse_exponent = (nuc.Z * u[0] + nuc.N * u[1] - u_c[nuc] + nuc.Z * up_c + nuc.nucbind * nuc.A) / (k * state.temp * Erg2MeV)
            nse_exponent = min(500.0, nse_exponent)

            Xs[nuc] = m_u * nuc.A_nuc * pf * nuc.spin_states / state.dens * (2.0 * np.pi * m_u * nuc.A_nuc * k * state.temp / h**2) ** 1.5 \
                    * np.exp(nse_exponent)

        return Xs

    def _constraint_eq(self, u, u_c, state):
        """ Constraint Equations used to evaluate chemical potential for proton and neutron,
        which is used when evaluating composition at NSE"""

        # Don't use eval_ye() since it does automatic mass fraction normalization.
        # However, we should force normalization through constraint eq1.

        Xs = self._nucleon_fraction_nse(u, u_c, state)

        eq1 = sum(Xs[nuc] for nuc in self.unique_nuclei) - 1.0
        eq2 = sum(Xs[nuc] * nuc.Z / nuc.A for nuc in self.unique_nuclei) - state.ye

        return [eq1, eq2]

    def get_comp_nse(self, rho, T, ye, init_guess=(-3.5, -15),
                     tol=1.0e-11, use_coulomb_corr=False, return_sol=False):
        """
        Returns the NSE composition given density, temperature and prescribed electron fraction
        using scipy.fsolve.

        Parameters:
        -------------------------------------
        rho: NSE state density

        T: NSE state Temperature

        ye: prescribed electron fraction

        init_guess: optional, initial guess of chemical potential of proton and neutron, [mu_p, mu_n]

        tol: optional, sets the tolerance of scipy.fsolve

        use_coulomb_corr: Whether to include coulomb correction terms

        return_sol: Whether to return the solution of the proton and neutron chemical potential.
        """

        # From here we convert the init_guess list into a np.array object:

        init_guess = np.array(init_guess)
        state = NseState(T, rho, ye)

        u_c = {}
        for nuc in set(list(self.unique_nuclei) + [Nucleus("p")]):
            u_c[nuc] = 0.0

        Xs = {}

        j = 0
        init_guess = np.array(init_guess)
        is_pos_old = False
        found_sol = False

        # Filter out runtimewarnings from fsolve, here we check convergence by np.isclose
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        # This nested loops should fine-tune the initial guess if fsolve is unable to find a solution
        while j < 20:
            i = 0
            guess = copy.deepcopy(init_guess)
            init_dx = 0.5

            while i < 20:
                u = fsolve(self._constraint_eq, guess, args=(u_c, state), xtol=tol, maxfev=800)
                Xs = self._nucleon_fraction_nse(u, u_c, state)
                n_e = self._evaluate_n_e(state, Xs)
                u_c = self._evaluate_mu_c(n_e, state, use_coulomb_corr)

                res = self._constraint_eq(u, u_c, state)
                is_pos_new = all(k > 0 for k in res)
                found_sol = np.all(np.isclose(res, [0.0, 0.0], rtol=1.0e-7, atol=1.0e-7))

                if found_sol:
                    Xs = self._nucleon_fraction_nse(u, u_c, state)
                    comp = Composition(self.unique_nuclei)
                    comp.X = Xs

                    if return_sol:
                        return comp, u

                    return comp

                if is_pos_old != is_pos_new:
                    init_dx *= 0.8

                if is_pos_new:
                    guess -= init_dx
                else:
                    guess += init_dx

                is_pos_old = is_pos_new
                i += 1

            j += 1
            init_guess[0] -= 0.5

        raise ValueError("Unable to find a solution, try to adjust initial guess manually")

    def generate_table(self, rho_values=None, T_values=None, Ye_values=None,
                       comp_reduction_func=None,
                       verbose=False, outfile="nse.tbl"):

        # initial guess
        mu_p0 = -3.5
        mu_n0 = -15.0

        # arrays to cache the chemical potentials as mu_p(rho, Ye)
        mu_p = np.ones((len(rho_values), len(Ye_values)), dtype=np.float64) * mu_p0
        mu_n = np.ones((len(rho_values), len(Ye_values)), dtype=np.float64) * mu_n0

        nse_states = []
        for T in reversed(T_values):
            for irho, rho in enumerate(reversed(rho_values)):
                for iye, ye in enumerate(reversed(Ye_values)):
                    initial_guess = (mu_p[irho, iye], mu_n[irho, iye])
                    try:
                        comp, sol = self.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                                      init_guess=initial_guess,
                                                      return_sol=True)
                    except ValueError:
                        initial_guess = (-3.5, -15)
                        comp, sol = self.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                                      init_guess=initial_guess,
                                                      return_sol=True)

                    mu_p[irho, iye] = sol[0]
                    mu_n[irho, iye] = sol[1]

                    # get the dY/dt for just the weak rates
                    ydots = self.evaluate_ydots(rho, T, comp,
                                                screen_func=potekhin_1998,
                                                rate_filter=lambda r: isinstance(r, TabularRate))

                    _, enu = self.evaluate_energy_generation(rho, T, comp,
                                                             screen_func=potekhin_1998,
                                                             return_enu=True)

                    nse_states.append(NSETableEntry(rho, T, ye,
                                                    comp=comp, ydots=ydots, enu=enu,
                                                    comp_reduction_func=comp_reduction_func))
                    if verbose:
                        print(nse_states[-1])

        with open(outfile, "w") as of:

            # write the header
            of.write(f"# NSE table generated by pynucastro {version}\n")
            of.write(f"# original NSENetwork had {len(self.unique_nuclei)} nuclei\n")
            of.write("#\n")
            of.write(f"# {'log10(rho)':^15} {'log10(T)':^15} {'Ye':^15} ")
            of.write(f"{'Abar':^15} {'<B/A>':^15} {'dYe/dt':^15} {'dAbar/dt':^15} {'e_nu':^15} ")

            if nse_states[0].X:
                for nuc, _ in nse_states[0].X:
                    _tmp = f"X({nuc})"
                    of.write(f"{_tmp:^15} ")

            of.write("\n")

            for entry in sorted(nse_states):
                of.write(f"{np.log10(entry.rho):15.10f} {np.log10(entry.T):15.10f} {entry.Ye:15.10f} ")
                of.write(f"{entry.abar:15.10f} {entry.bea:15.10f} {entry.dYedt:15.8g} {entry.dabardt:15.8g} {entry.enu:15.8g} ")

                if entry.X:
                    for _, val in entry.X:
                        of.write(f"{val:15.10g} ")
                of.write("\n")