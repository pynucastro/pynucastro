"""A collection of classes and methods for working with nuclear
statistical equilibrium.

"""

import warnings

import numpy as np
from scipy.optimize import fsolve

from pynucastro._version import version
from pynucastro.constants import constants
from pynucastro.networks.rate_collection import Composition, RateCollection
from pynucastro.rates import TabularRate
from pynucastro.screening import NseState, potekhin_1998


class NSETableEntry:
    """A simple container to hold a single entry in the NSE table.

    Parameters
    ----------
    rho : float
        the density of the NSE state
    T : float
        the temperature of the NSE state
    Ye : float
        the electron fraction
    comp : Composition
        the NSE composition
    ydots : dict
        a dictionary of dY/dt keyed by Nucleus.  This is meant to be
        the weak nuclear rates only, since those affect the NSE state.
    enu : float
        the weak rate neutrino energy loss
    comp_reduction_function : Callable
        a function that converts the NSE composition into a smaller set
        of nuclei.  It takes a Composition object and returns a dictionary
        with the nucleus name (like "Ni56") as the key and the corresponding
        mass fraction as the value.  It should be ordered in the way you
        want the nuclei output into the NSE table file.
    """

    def __init__(self, rho, T, Ye, *,
                 comp=None, ydots=None, enu=None,
                 comp_reduction_func=None):

        self.rho = rho
        self.T = T
        self.Ye = Ye

        self.comp = comp
        self.ydots = ydots
        self.enu = enu

        # compute the bits we need for the table

        if comp:
            # mean molecular weight of the full NSE state
            self.abar = comp.abar

            # average binding energy / nucleon for the full NSE state
            self.bea = sum(q.nucbind * self.comp.X[q] for q in self.comp.X)

            # evolution of the electron fraction from weak rates alone
            self.dYedt = sum(q.Z * self.ydots[q] for q in self.comp.X)

            # evolution of Abar from weak rates alone
            self.dabardt = -self.abar**2 * sum(self.ydots[q] for q in self.comp.X)

            # evolution of B/A from weak rates alone
            self.dbeadt = sum(self.ydots[q] * q.nucbind for q in self.comp.X)

            self.X = None
            if comp_reduction_func:
                self.X = comp_reduction_func(self.comp)

    def __str__(self):
        return f"({self.rho:12.6g}, {self.T:12.6g}, {self.Ye:6.4f}): {self.abar:6.4f}  {self.bea:6.4f}  {self.dYedt:12.6g}  {self.enu:12.6g}"

    def value(self):
        """Return an integer used for sorting.  This has the format
        (logrho)(logT)(1-Ye)

        """

        return int(f"{np.log10(self.rho):08.5f}{np.log10(self.T):08.5f}{1-self.Ye:08.5f}".replace(".", ""))

    def __lt__(self, other):
        return self.value() < other.value()


class NSENetwork(RateCollection):
    """A network for solving the NSE composition and outputting
    tabulated NSE quantities

    """

    def __init__(self, *args, use_unreliable_spins=True, **kwargs):
        """Parameters
        ----------
        use_unreliable_spins : bool
            whether to allow nuclei with weakly supported spin data
        """

        self.use_unreliable_spins = use_unreliable_spins
        super().__init__(*args, **kwargs)

    def _evaluate_mu_c(self, state, use_coulomb_corr=True):
        """Find the Coulomb potential of each nuclide in NSE state.
        The Coulomb potential is evaluated using the Helmholtz free
        energy following Eq. 28 from Chabrier & Potekhin 1998 DOI:
        10.1103/PhysRevE.58.4941 Also see appendix in Calder 2007 DOI:
        10.1086/510709

        A1 = -0.9052
        A2 = 0.6322
        A3 = -sqrt(3) / 2 - A1 / sqrt(A2)
        Γ = Z^{5/3} q_e^2 (4 π n_e / 3)^{1/3} / (k_B T)
        u^c / (k_B T)=   A1 [ sqrt(Γ (A2 + Γ)) - A2 ln(sqrt(Γ / A2) + sqrt(1 + Γ / A2)) ]
                       + 2 A3 [ sqrt(Γ) - atan(sqrt(Γ)) ]

        Parameters
        ----------
        state : NseState
            NSE state
        use_coulomb_corr: bool
            whether to include coulomb correction terms

        Returns
        -------
        u_c : dict
            A dictionary of coulomb potential keyed by Nucleus.

        """

        if not use_coulomb_corr:
            u_c = {n: 0 for n in self.unique_nuclei}
            return u_c

        u_c = {}
        n_e = state.ye * state.dens / constants.m_u_C18
        A_1 = -0.9052
        A_2 = 0.6322
        A_3 = -0.5 * np.sqrt(3.0) - A_1 / np.sqrt(A_2)
        Gamma_e = state.gamma_e_fac * n_e**(1.0/3.0) / state.temp

        for nuc in self.unique_nuclei:
            Gamma = Gamma_e * nuc.Z**(5.0 / 3.0)
            u_c[nuc] = constants.erg2MeV * constants.k * state.temp * \
                (A_1 * (np.sqrt(Gamma * (A_2 + Gamma)) - A_2 * np.log(np.sqrt(Gamma / A_2) +
                np.sqrt(1.0 + Gamma / A_2))) + 2.0 * A_3 * (np.sqrt(Gamma) - np.arctan(np.sqrt(Gamma))))

        return u_c

    def _nucleon_fraction_nse(self, u, u_c, state):
        """Compute the NSE mass fraction for a given NSE state.

        NSE condition says that the chemical potential of the i-th nuclei
        is equal to the chemical potential of the free nucleon.

        i.e. μ_i = Z_i μ_p + N_i μ_n

        where μ_p and μ_n are the chemical potential of free proton and neutron,
        which can be decomposed into the kinetic part μ^{id}, coulomb potential
        part μ^c and its rest energy.

        μ_p = μ^{id}_p + μ^{c}_p + m_p c^2
        μ_n = μ^{id}_n           + m_n c^2

        Assuming Boltzmann statistics and the NSE condition, we get:

        X_i = m_i / ρ (2 J_i + 1) G_i (m_i k_B T / (2 π ℏ^2))^{3/2}
            * exp{ Z_i (μ^{id}_p + μ^{c}_p) + N_i μ^{id}_n + B_i - μ^{c}_i}

        See full derivation in appendix Smith Clark 2023 DOI: 10.3847/1538-4357/acbaff

        Parameters
        ----------
        u : (ndarray, tuple, list)
            chemical potentials of proton and neutron, where we choose
            u[0] = μ^{id}_p + μ^{c}_p and u[1] = μ^{id}_n

            One can make the choice of setting u[0] = μ^{id}_p only and
            add μ^{c}_p as a separate constant. However it won't affect the end result.
        u_c : dict
            A dictionary of coulomb potential keyed by Nucleus.
        state : NseState
            NSE state.

        Returns
        -------
        Xs : dict
            A dictionary of NSE mass fractions keyed by Nucleus.

        """

        Xs = {}
        for nuc in self.unique_nuclei:
            nse_exponent = 0.0
            if nuc.partition_function is not None:
                nse_exponent += nuc.partition_function.eval(state.temp)

            if not nuc.spin_states:
                raise ValueError(f"The spin of {nuc} is not implemented for now.")
            if not self.use_unreliable_spins and not nuc.spin_reliable:
                raise ValueError(f"The spin of {nuc} is determined by a weak experimental or theoretical argument. "
                                 "Pass in use_unreliable_spins=True as a parameter to NSENetwork() to override.")

            nse_exponent += (nuc.Z * u[0] + nuc.N * u[1] - u_c[nuc] + nuc.nucbind * nuc.A) / (constants.k_MeV * state.temp)
            nse_exponent = min(500.0, nse_exponent)

            Xs[nuc] = (nuc.A_nuc * constants.m_u_C18)**2.5 * pf * nuc.spin_states / state.dens * \
                (constants.k * state.temp / (2.0 * np.pi * constants.hbar**2))**1.5 * np.exp(nse_exponent)

        return Xs

    def _constraint_eq(self, u, u_c, state):
        """Implement the constraints for our system to evaluate
        chemical potential for proton and neutron, which is used when
        evaluating composition at NSE.

        1) Conservation of Mass:   Σ_k X_k - 1 = 0
        2) Conservation of Charge: Σ_k Z_k X_k / A_k - Y_e = 0

        Parameters
        ----------
        u : (ndarray, tuple, list)
            chemical potentials of proton and neutron.
            u[0] = μ^{id}_p + μ^{c}_p and u[1] = μ^{id}_n
        u_c : dict
            A dictionary of coulomb potential keyed by Nucleus.
        state : NseState
            NSE state

        Returns
        -------
        constraint_eqs : List
            Constraint equations for NSE.

        """

        Xs = self._nucleon_fraction_nse(u, u_c, state)
        constraint_eqs = [sum(Xs[nuc] for nuc in self.unique_nuclei) - 1.0,
                          sum(Xs[nuc] * nuc.Z / nuc.A for nuc in self.unique_nuclei) - state.ye]

        return constraint_eqs

    def get_comp_nse(self, rho, T, ye, init_guess=(-3.5, -15),
                     tol=1.0e-11, use_coulomb_corr=False,
                     return_sol=False):
        """Return the NSE composition given density, temperature and
        prescribed electron fraction using scipy.fsolve.

        Parameters
        ----------
        rho : float
            NSE state density
        T : float
            NSE state Temperature
        ye : float
            prescribed electron fraction
        init_guess : (tuple, list)
            initial guess of chemical potential of proton and neutron, [μ^{id}_p + μ^{c}_p, μ^{id}_n]
        tol : float
            tolerance of scipy.fsolve
        use_coulomb_corr : bool
            whether to include coulomb correction terms
        return_sol : bool
            whether to return the solution of the proton and neutron chemical potential.

        Returns
        -------
        comp : Composition
            the NSE composition
        u : ndarray
            the chemical potentials found by solving the nonlinear system.

        """

        # Determine the upper and lower bound of possible ye for the current network.
        # Make sure the prescribed ye are within this range.
        ye_low = min(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        ye_max = max(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        assert ye_low <= ye <= ye_max, "input electron fraction goes outside of scope for current network"

        init_guess = np.array(init_guess)
        state = NseState(T, rho, ye)

        # Evaluate coulomb potential
        u_c = self._evaluate_mu_c(state, use_coulomb_corr)

        j = 0
        is_pos_old = False
        found_sol = False

        # This nested loops should fine-tune the initial guess if fsolve is unable to find a solution
        while j < 20:
            i = 0
            guess = init_guess.copy()
            init_dx = 0.5

            while i < 20:
                # Filter out runtimewarnings from fsolve, here we check convergence by np.isclose
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=RuntimeWarning)
                    u = fsolve(self._constraint_eq, guess, args=(u_c, state), xtol=tol, maxfev=800)
                Xs = self._nucleon_fraction_nse(u, u_c, state)
                res = self._constraint_eq(u, u_c, state)
                is_pos_new = all(k > 0 for k in res)
                found_sol = np.all(np.isclose(res, [0.0, 0.0], rtol=1.0e-11, atol=1.0e-11))

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
        """Generate a table of NSE properties.  For every combination
        of density, temperature, and Ye, we solve for the NSE state
        and output composition properties to a file.

        Parameters
        ----------
        rho_values : numpy.ndarray
            values of density to use in the tabulation
        T_values : numpy.ndarray
            values of temperature to use in the tabulation
        Ye_values : numpy.ndarray
            values of electron fraction to use in the tabulation
        comp_reduction_func : Callable
            a function that takes the NSE composition and return a reduced
            composition
        verbose : bool
            output progress on creating the table as we go along
        outfile : str
            filename for the NSE table output

        """

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
            of.write(f"{'Abar':^15} {'<B/A>':^15} {'dYe/dt':^15} {'dAbar/dt':^15} {'d<B/A>/dt':^15} {'e_nu':^15} ")

            if nse_states[0].X:
                for nuc, _ in nse_states[0].X:
                    _tmp = f"X({nuc})"
                    of.write(f"{_tmp:^15} ")

            of.write("\n")

            for entry in sorted(nse_states):
                of.write(f"{np.log10(entry.rho):15.10f} {np.log10(entry.T):15.10f} {entry.Ye:15.10f} ")
                of.write(f"{entry.abar:15.10f} {entry.bea:15.10f} {entry.dYedt:15.8g} {entry.dabardt:15.8g} {entry.dbeadt:15.8g} {entry.enu:15.8g} ")

                if entry.X:
                    for _, val in entry.X:
                        of.write(f"{val:15.10g} ")
                of.write("\n")
