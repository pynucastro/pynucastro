"""
Python implementations of screening routines.
"""
import numpy as np
from scipy import constants

# use the jitclass placeholder from rate.py
from pynucastro.rates.rate import jitclass, numba

if numba is not None:
    njit = numba.njit
else:
    def njit(func):
        return func

__all__ = ["PlasmaState", "ScreenFactors", "chugunov_2007", "chugunov_2009",
           "make_plasma_state", "make_screen_factors", "potekhin_1998"]


amu = constants.value("atomic mass constant") / constants.gram  # kg to g
q_e = constants.value("elementary charge") * (constants.c * 100) / 10  # C to statC (esu)
hbar = constants.value("reduced Planck constant") / constants.erg  # J*s to erg*s
k_B = constants.value("Boltzmann constant") / constants.erg  # J/K to erg/K


@jitclass()
class PlasmaState:
    """
    Stores precomputed values that are reused for all screening correction
    factor calculations.

    :var temp:        temperature in K
    :var dens:        density in g/cm^3
    :var abar:        average atomic mass
    :var zbar:        average ion charge
    :var z2bar:       average (ion charge)^2
    :var n_e:         electron number density
    :var gamma_e_fac: temperature-independent part of Gamma_e
    """
    temp: float
    dens: float
    abar: float
    zbar: float
    z2bar: float
    n_e: float
    gamma_e_fac: float

    def __init__(self, temp, dens, Ys, Zs):
        """
        :param temp: temperature in K
        :param dens: density in g/cm^3
        :param Ys:   molar fractions of each ion
        :type Ys: numpy ndarray
        :param Zs:   charge of each ion, in the same order as Ys
        :type Zs: numpy ndarray
        """
        self.temp = temp
        self.dens = dens
        ytot = np.sum(Ys)
        self.abar = 1 / ytot
        self.zbar = np.sum(Zs * Ys) / ytot
        self.z2bar = np.sum(Zs ** 2 * Ys) / ytot

        # Average mass and total number density
        mbar = self.abar * amu
        ntot = self.dens / mbar
        # Electron number density
        # zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
        self.n_e = self.zbar * ntot

        # temperature-independent part of Gamma_e, from Chugunov 2009 eq. 6
        self.gamma_e_fac = q_e ** 2 / k_B * np.cbrt(4 * np.pi / 3) * np.cbrt(self.n_e)


@jitclass()
class NseState:
    """
    Stores precomputed values that are reused in the NSE state screening
    calculations

    :var temp:        temperature in K
    :var dens:        density in g/cm^3
    :var ye:          electron molar fraction
    :var n_e:         electron number density
    :var gamma_e_fac: temperature-independent part of Gamma_e
    """

    temp: float
    dens: float
    ye: float
    gamma_e_fac: float

    def __init__(self, temp, dens, ye):

        """
        :param temp:        temperature in K
        :param dens:        density in g/cm^3
        :param ye:          electron molar fraction
        :param Xs:          nucleon fraction of each ion
        :type Xs: numpy ndarray
        :param As:          atomic mass number of each ion
        :type As: numpy ndarray
        :param Zs:          atomic number of each ion
        :type Zs: numpy ndarray
        """

        self.temp = temp
        self.dens = dens
        self.ye = ye
        self.gamma_e_fac = q_e ** 2 / k_B * np.cbrt(4.0 * np.pi / 3.0)


def make_plasma_state(temp, dens, molar_fractions):
    """
    Construct a PlasmaState object from simulation data.

    :param temp:            temperature in K
    :param dens:            density in g/cm^3
    :param molar_fractions: dictionary of molar abundances for each ion,
                            as returned by :meth:`.Composition.get_molar`
    """
    nuclei = list(molar_fractions.keys())
    Ys = np.asarray([molar_fractions[n] for n in nuclei])
    Zs = np.asarray([n.Z for n in nuclei])
    return PlasmaState(temp, dens, Ys, Zs)


@jitclass()
class ScreenFactors:
    """
    Stores values that will be used to calculate the screening correction factor
    for a specific pair of nuclei.

    :var z1: atomic number of first nucleus
    :var z2: atomic number of second nucleus
    :var a1: atomic mass of first nucleus
    :var a2: atomic mass of second nucleus
    :var aznut: combination of a1, z1, a2, z2 raised to 1/3 power
    :var ztilde: effective ion radius factor for a MCP
    """
    z1: int
    z2: int
    a1: int
    a2: int
    aznut: float
    ztilde: float

    def __init__(self, z1, a1, z2, a2):
        self.z1 = z1
        self.z2 = z2
        self.a1 = a1
        self.a2 = a2
        self.aznut = np.cbrt(z1 ** 2 * z2 ** 2 * a1 * a2 / (a1 + a2))
        self.ztilde = 0.5 * (np.cbrt(z1) + np.cbrt(z2))


def make_screen_factors(n1, n2):
    """
    Construct a ScreenFactors object from a pair of nuclei.

    :param Nucleus n1: first nucleus
    :param Nucleus n2: second nucleus
    """
    return ScreenFactors(n1.Z, n1.A, n2.Z, n2.A)


@njit
def smooth_clip(x, limit, start):
    """Smoothly transition between y=limit and y=x with a half-cosine.

    Clips smaller values if limit < start and larger values if start < limit.

    :param x:     the value to clip
    :param limit: the constant value to clip x to
    :param start: the x-value at which to start the transition
    :returns: y
    """
    if limit < start:
        lower = limit
        upper = x
    else:
        lower = x
        upper = limit

    if x < min(limit, start):
        return lower
    if x > max(limit, start):
        return upper

    tmp = np.pi * (x - min(limit, start)) / (start - limit)
    f = (1 - np.cos(tmp)) / 2

    return (1 - f) * lower + f * upper


@njit
def chugunov_2007(state, scn_fac):
    """Calculates screening factors based on Chugunov et al. 2007.

    Follows the approach in Yakovlev 2006 to extend to a multi-component plasma.

    :param PlasmaState state:     the precomputed plasma state factors
    :param ScreenFactors scn_fac: the precomputed ion pair factors
    :returns: screening correction factor

    References:
        | Chugunov, DeWitt, and Yakovlev 2007, PhRvD, 76, 025028
        | Yakovlev, Gasques, Afanasjev, Beard, and Wiescher 2006, PhRvC, 74, 035803
        | Chugunov and DeWitt 2009, PhRvC, 80, 014611
    """
    # Plasma temperature T_p
    # This formula comes from working backwards from zeta_ij (Chugunov 2009 eq. 12)
    # through Chugunov 2007 eq. 3 to Chugunov 2007 eq. 2.
    # Ultimately, the changes from the expression in Chugunov 2007 are:
    #   Z^2 -> Z1 * Z2
    #   n_i -> n_e / ztilde^3, where ztilde = (Z1^(1/3) + Z2^(1/3)) / 2
    #   m_i -> 2 mu12 (reduced mass)
    # This prescription reduces to the expressions from Chugunov 2007 in the case
    # of an OCP, and to Chugunov 2009 in the case of a binary ionic mixture.
    # This also matches Yakovlev et al. 2006, eq. 10.
    #
    # For reference, MESA r21.12.1 does:
    #   Z^2 -> Z1 * Z2
    #   n_i -> n_e / zbar (=ntot)
    #   m_i -> m_u * abar
    #
    # Sam Jones' Fortran implementation does:
    #   Z^2 -> zbar^2
    #   n_i -> ntot
    #   m_i -> m_u * abar
    mu12 = scn_fac.a1 * scn_fac.a2 / (scn_fac.a1 + scn_fac.a2)
    z_factor = scn_fac.z1 * scn_fac.z2
    n_i = state.n_e / scn_fac.ztilde ** 3
    m_i = 2 * mu12 * amu

    T_p = hbar / k_B * q_e * np.sqrt(4 * np.pi * z_factor * n_i / m_i)

    # Normalized temperature
    T_norm = state.temp / T_p

    # The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    # should be nearly temperature-independent (in the pycnonuclear regime),
    # and we clip the temperature to 0.1 T_p at small T.
    # start the transition here
    T_norm_fade = 0.2
    # minimum value of T/T_p
    T_norm_min = 0.1

    T_norm = smooth_clip(T_norm, limit=T_norm_min, start=T_norm_fade)

    # Coulomb coupling parameter from Yakovlev 2006, eq. 10
    Gamma = state.gamma_e_fac * scn_fac.z1 * scn_fac.z2 / (scn_fac.ztilde * T_norm * T_p)

    # The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    Gamma_fade = 590
    Gamma_max = 600
    Gamma = smooth_clip(Gamma, limit=Gamma_max, start=Gamma_fade)

    # Chugunov 2007 eq. 3
    zeta = np.cbrt(4 / (3 * np.pi ** 2 * T_norm ** 2))

    # Gamma tilde from Chugunov 2007 eq. 21
    fit_alpha = 0.022
    fit_beta = 0.41 - 0.6 / Gamma
    fit_gamma = 0.06 + 2.2 / Gamma

    # Polynomial term in Gamma tilde
    poly = 1 + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta))

    gamtilde = Gamma / np.cbrt(poly)

    # fit parameters just after Chugunov 2007 eq. 19
    A1 = 2.7822
    A2 = 98.34
    A3 = np.sqrt(3) - A1 / np.sqrt(A2)
    B1 = -1.7476
    B2 = 66.07
    B3 = 1.12
    B4 = 65
    gamtilde2 = gamtilde ** 2

    # Chugunov 2007 eq. 19
    term1 = 1 / np.sqrt(A2 + gamtilde)
    term2 = 1 / (1 + gamtilde)
    term3 = gamtilde ** 2 / (B2 + gamtilde)
    term4 = gamtilde2 / (B4 + gamtilde2)

    h = gamtilde ** (3 / 2) * (A1 * term1 + A3 * term2) + B1 * term3 + B3 * term4

    # machine limit the output
    h_max = 300
    h = min(h, h_max)
    scor = np.exp(h)

    return scor


@njit
def f0(gamma):
    r"""Calculate the free energy per ion in a OCP from Chugunov & DeWitt 2009 eq. 24

    :param gamma: Coulomb coupling parameter
    :returns: free energy
    """
    A1 = -0.907
    A2 = 0.62954
    A3 = -np.sqrt(3) / 2 - A1 / np.sqrt(A2)
    B1 = 0.00456
    B2 = 211.6
    B3 = -1e-4
    B4 = 0.00462

    term1 = np.sqrt(gamma * (A2 + gamma))
    term2 = np.log(np.sqrt(gamma / A2) + np.sqrt(1 + gamma / A2))
    gamma_12 = np.sqrt(gamma)
    term3 = gamma_12 - np.arctan(gamma_12)
    term4 = np.log1p(gamma / B2)
    term5 = np.log1p(gamma ** 2 / B4)

    return (
        A1 * (term1 - A2 * term2) +
        2 * A3 * term3 +
        B1 * (gamma - B2 * term4) +
        B3 / 2 * term5
    )


@njit
def chugunov_2009(state, scn_fac):
    """Calculates screening factors based on Chugunov & DeWitt 2009.

    :param PlasmaState state:     the precomputed plasma state factors
    :param ScreenFactors scn_fac: the precomputed ion pair factors
    :returns: screening correction factor

    References:
        | Chugunov and DeWitt 2009, PhRvC, 80, 014611
    """
    z1z2 = scn_fac.z1 * scn_fac.z2
    zcomp = scn_fac.z1 + scn_fac.z2

    # Gamma_e from eq. 6
    Gamma_e = state.gamma_e_fac / state.temp

    # Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    Gamma_1 = Gamma_e * scn_fac.z1 ** (5 / 3)
    Gamma_2 = Gamma_e * scn_fac.z2 ** (5 / 3)
    Gamma_comp = Gamma_e * zcomp ** (5 / 3)

    Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde

    # Coulomb barrier penetrability, eq. 10
    tau_factor = np.cbrt(27 / 2 * (np.pi * q_e ** 2 / hbar) ** 2 * amu / k_B)
    tau_12 = tau_factor * scn_fac.aznut / np.cbrt(state.temp)

    # eq. 12
    zeta = 3 * Gamma_12 / tau_12

    # additional fit parameters, eq. 25
    y_12 = 4 * z1z2 / zcomp ** 2
    c1 = 0.013 * y_12 ** 2
    c2 = 0.406 * y_12 ** 0.14
    c3 = 0.062 * y_12 ** 0.19 + 1.8 / Gamma_12

    poly = 1 + zeta*(c1 + zeta*(c2 + c3*zeta))
    t_12 = np.cbrt(poly)

    # strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    # Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    # values similar to those from Chugunov 2007.
    term1 = f0(Gamma_1 / t_12)
    term2 = f0(Gamma_2 / t_12)
    term3 = f0(Gamma_comp / t_12)
    h_fit = term1 + term2 - term3

    # weak screening correction term, eq. A3
    corr_C = (
        3 * z1z2 * np.sqrt(state.z2bar / state.zbar) /
        (zcomp ** 2.5 - scn_fac.z1 ** 2.5 - scn_fac.z2 ** 2.5)
    )

    # corrected enhancement factor, eq. A4
    Gamma_12_2 = Gamma_12 ** 2
    numer = corr_C + Gamma_12_2
    denom = 1 + Gamma_12_2
    h12 = numer / denom * h_fit

    # machine limit the output
    h12_max = 300
    h12 = min(h12, h12_max)
    scor = np.exp(h12)

    return scor


@njit
def potekhin_1998(state, scn_fac):
    """Calculates screening factors based on Chabrier & Potekhin 1998.

    :param PlasmaState state:     the precomputed plasma state factors
    :param ScreenFactors scn_fac: the precomputed ion pair factors
    :returns: screening correction factor

    References:
        Chabrier and Potekhin 1998, PhRvE, 58, 4941
    """

    Gamma_e = state.gamma_e_fac / state.temp
    zcomp = scn_fac.z1 + scn_fac.z2

    Gamma_1 = Gamma_e * scn_fac.z1 ** (5 / 3)
    Gamma_2 = Gamma_e * scn_fac.z2 ** (5 / 3)
    Gamma_comp = Gamma_e * zcomp ** (5 / 3)

    A_1 = -0.9052
    A_2 = 0.6322
    A_3 = -0.5 * np.sqrt(3) - A_1/np.sqrt(A_2)

    f1 = A_1 * (np.sqrt(Gamma_1 * (A_2 + Gamma_1)) - A_2 * np.log(np.sqrt(Gamma_1 / A_2) +
                  np.sqrt(1.0 + Gamma_1/A_2))) + 2.0 * A_3 * (np.sqrt(Gamma_1) - np.arctan(np.sqrt(Gamma_1)))

    f2 = A_1 * (np.sqrt(Gamma_2 * (A_2 + Gamma_2)) - A_2 * np.log(np.sqrt(Gamma_2 / A_2) +
                  np.sqrt(1.0 + Gamma_2/A_2))) + 2.0 * A_3 * (np.sqrt(Gamma_2) - np.arctan(np.sqrt(Gamma_2)))

    f12 = A_1 * (np.sqrt(Gamma_comp * (A_2 + Gamma_comp)) - A_2 * np.log(np.sqrt(Gamma_comp / A_2) +
                  np.sqrt(1.0 + Gamma_comp/A_2))) + 2.0 * A_3 * (np.sqrt(Gamma_comp) - np.arctan(np.sqrt(Gamma_comp)))

    h12 = f1 + f2 - f12

    # machine limit the output
    h12_max = 300
    h12 = min(h12, h12_max)
    scor = np.exp(h12)

    return scor
