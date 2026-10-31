"""Python implementations of common electron screening routines."""

from collections import namedtuple

import numpy as np

import pynucastro.numba_util as numba
from pynucastro.autodiff import dual
from pynucastro.autodiff.derivative_pair import DerivativePair
from pynucastro.constants import constants
from pynucastro.nucdata import Composition, Nucleus
from pynucastro.numba_util import jitclass, njit

# jit-decorated functions should be listed here so they
# get documented by sphinx
__all__ = [
    "NseState",
    "PlasmaState",
    "ScreenFactors",
    "chugunov_2007",
    "chugunov_2009",
    "debye_huckel",
    "f0",
    "get_screening_func",
    "make_plasma_state",
    "make_plasma_state_with_derivatives",
    "make_screen_factors",
    "potekhin_1998",
    "screen5",
    "screening_check",
    "smooth_clip",
]


# A native tuple can carry either floats or dual-valued composition factors.
_PlasmaStateData = namedtuple("PlasmaStateData",
                             "temp dens qlam0z taufac aa abar zbar z2bar n_e gamma_e_fac")


@njit
def _plasma_state_from_moments(temp, dens, ytot, ye, s2):
    """Compute state factors at fixed total molar abundance."""
    abar = 1 / ytot
    zbar = ye / ytot
    z2bar = s2 / ytot

    # ntot
    rr = dens * ytot

    # Part version of Eq. 19 in Graboske:1973
    # pp = sqrt( \tilde{z}*(rho/u_I/T) )
    pp = dual.sqrt(rr/temp*(z2bar + zbar))
    qlam0z = 1.88e8 / temp * pp

    # Part of Eq.6 in Itoh:1979
    # 4.248719e3 = (27*pi^2*e^4*m_u/(2*k_B*hbar^2))^(1/3)
    # the extra (1/3) to make tau -> tau/3
    co2 = dual.cbrt(27*np.pi**2*constants.q_e**4*constants.m_u_C18 /
                    (2*constants.k*constants.hbar**2)) / 3
    taufac = co2 / dual.cbrt(temp)

    xni = dual.cbrt(rr * zbar)

    # Part of Eq.4 in Itoh:1979
    # 2.27493e5 = e^2 / ( (3*m_u/(4pi))^(1/3) *k_B )
    aa_factor = constants.q_e**2 / (dual.cbrt(3*constants.m_u_C18/(4*np.pi)) *
                                    constants.k)
    aa = aa_factor / temp * xni

    # Average mass and total number density
    mbar = abar * constants.m_u_C18
    ntot = dens / mbar
    # Electron number density
    # zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
    n_e = zbar * ntot

    # temperature-independent part of Gamma_e, from Chugunov 2009 eq. 6
    gamma_e_fac = constants.q_e**2 / constants.k * dual.cbrt(4 * np.pi / 3) * dual.cbrt(n_e)
    return _PlasmaStateData(temp, dens, qlam0z, taufac, aa, abar, zbar, z2bar, n_e, gamma_e_fac)


@jitclass([
    ("temp", numba.float64),
    ("dens", numba.float64),
    ("qlam0z", numba.float64),
    ("taufac", numba.float64),
    ("aa", numba.float64),
    ("abar", numba.float64),
    ("zbar", numba.float64),
    ("z2bar", numba.float64),
    ("n_e", numba.float64),
    ("gamma_e_fac", numba.float64)
])
class PlasmaState:
    """Store precomputed values that are reused for all screening
    correction factor calculations.

    Parameters
    ----------
    temp : float
        temperature in K
    dens : float
        density in g/cm^3
    Ys : Iterable
        molar fractions of the composition
    Zs : Iterable
         proton numbers of the composition

    Attributes
    ----------
    qlam0z : float
        a common factor from Graboske 1973
    taufac : float
        a common factor from Itoh 1979
    aa : float
        a common factor from Itoh 1979
    abar : float
        average atomic mass
    zbar : float
        average ion charge
    z2bar : float
       average (ion charge)^2
    n_e : float
         electron number density
    gamma_e_fac : float
         temperature-independent part of Gamma_e

    """

    def __init__(self, temp, dens, Ys, Zs):
        state = _plasma_state_from_moments(temp, dens, np.sum(Ys),
                                           np.sum(Zs * Ys), np.sum(Zs**2 * Ys))
        self.temp = state.temp
        self.dens = state.dens
        self.qlam0z = state.qlam0z
        self.taufac = state.taufac
        self.aa = state.aa
        self.abar = state.abar
        self.zbar = state.zbar
        self.z2bar = state.z2bar
        self.n_e = state.n_e
        self.gamma_e_fac = state.gamma_e_fac


@jitclass([
    ("temp", numba.float64),
    ("dens", numba.float64),
    ("ye", numba.float64),
    ("gamma_e_fac", numba.float64)
])
class NseState:
    """Store precomputed values that are reused in the NSE state
    screening calculations

    Parameters
    ----------
    temp : float
        Temperature in K
    dens : float
        Density in g/cm^3
    ye : float
        Electron molar fraction

    Attributes
    ----------
    gamma_e_fac : float
        Temperature-independent part of Gamma_e.

    """

    def __init__(self, temp, dens, ye):

        self.temp = temp
        self.dens = dens
        self.ye = ye
        self.gamma_e_fac = constants.q_e**2 / constants.k * dual.cbrt(4.0 * np.pi / 3.0)


def make_plasma_state(temp, dens, composition, *, compute_derivatives=False):
    """Construct a PlasmaState object from simulation data.

    Parameters
    ----------
    temp : float
        temperature in K
    dens : float
        density in g/cm^3
    composition : Composition or dict(Nucleus)
        a Composition object or a dictionary of molar abundances for each ion,
        as returned by :meth:`.Composition.get_molar`
    compute_derivatives : bool
        Seed derivatives with respect to Ye = sum(Z Y) and S2 = sum(Z**2 Y).
        The default is False.

    Returns
    -------
    PlasmaState or tuple
        With compute_derivatives=True, a named tuple whose
        composition-dependent fields carry two derivatives.

    """

    if isinstance(composition, Composition):
        Ys = composition.get_molar_array()
        Zs = np.asarray([nuc.Z for nuc in composition])
    else:
        nuclei = list(composition.keys())
        Ys = np.asarray([composition[n] for n in nuclei])
        Zs = np.asarray([n.Z for n in nuclei])
    if compute_derivatives:
        return make_plasma_state_with_derivatives(temp, dens, Ys, Zs)
    return PlasmaState(temp, dens, Ys, Zs)


@njit
def make_plasma_state_with_derivatives(temp, dens, Ys, Zs):
    """Construct a state carrying two composition derivatives.

    Parameters
    ----------
    temp : float
        Temperature in K, held fixed during differentiation.
    dens : float
        Density in g/cm^3, held fixed during differentiation.
    Ys : numpy.ndarray
        Molar abundances, Y_i = X_i / A_i.
    Zs : numpy.ndarray
        Nuclear charges in the same order as Ys.

    Returns
    -------
    tuple
        Named tuple with the same attributes as PlasmaState. Composition-dependent
        fields are Dual values with DerivativePair seeds: d1 corresponds to
        Ye = sum(Z_i Y_i), and d2 to S2 = sum(Z_i**2 Y_i).

    Notes
    -----
    The total molar abundance is held fixed when seeding these moments.
    It cancels out of the screening expressions. For a returned log screening
    factor h, dh/dY_i = Z_i * h.deriv.d1 + Z_i**2 * h.deriv.d2, without
    renormalizing the abundances. Reuse this state for different ion pairs.

    """
    ye = dual.Dual(np.sum(Zs * Ys), DerivativePair(1.0, 0.0))
    s2 = dual.Dual(np.sum(Zs**2 * Ys), DerivativePair(0.0, 1.0))
    return _plasma_state_from_moments(temp, dens, np.sum(Ys), ye, s2)


@jitclass([
    ("z1", numba.int32),
    ("z2", numba.int32),
    ("a1", numba.int32),
    ("a2", numba.int32),
    ("zs13", numba.float64),
    ("zhat", numba.float64),
    ("zhat2", numba.float64),
    ("lzav", numba.float64),
    ("aznut", numba.float64),
    ("ztilde", numba.float64)
])
class ScreenFactors:
    """Store values that will be used to calculate the screening
    correction factor for a specific pair of nuclei.

    Parameters
    ----------
    z1 : float
        atomic number of first nucleus
    z2 : float
        atomic number of second nucleus
    a1 : float
        atomic mass of first nucleus
    a2 : float
        atomic mass of second nucleus

    Attributes
    ----------
    zs13 : float
        (z1+z2)**(1/3)
    zhat : float
        combination of z1 and z2 raised to the 5/3 power
    zhat2 : float
        combination of z1 and z2 raised to the 5/12 power
    lzav : float
        log of effective charge
    aznut : float
        combination of a1, z1, a2, z2 raised to 1/3 power
    ztilde : float
        effective ion radius factor for a MCP

    """

    def __init__(self, z1, a1, z2, a2):
        self.z1 = z1
        self.z2 = z2
        self.a1 = a1
        self.a2 = a2
        self.zs13 = dual.cbrt(z1 + z2)
        self.zhat = (z1 + z2)**(5/3) - z1**(5/3) - z2**(5/3)
        self.zhat2 = (z1 + z2)**(5/12) - z1**(5/12) - z2**(5/12)
        self.lzav = (5/3) * dual.log(z1 * z2 / (z1 + z2))
        self.aznut = dual.cbrt(z1**2 * z2**2 * a1 * a2 / (a1 + a2))
        self.ztilde = 0.5 * (dual.cbrt(z1) + dual.cbrt(z2))


def make_screen_factors(n1, n2):
    """Construct a ScreenFactors object from a pair of nuclei.

    Parameters
    ----------
    n1 : Nucleus, str
        first nucleus
    n2 : Nucleus, str
        second nucleus

    Returns
    -------
    ScreenFactors

    """

    n1 = Nucleus.cast(n1)
    n2 = Nucleus.cast(n2)
    return ScreenFactors(n1.Z, n1.A, n2.Z, n2.A)


@njit
def debye_huckel(state, scn_fac):
    """Calculate the Debye-Huckel enhancement factor for weak Coloumb
    coupling, following the appendix of :cite:t:`chugunov:2009`.

    Parameters
    ----------
    state : PlasmaState
        the precomputed plasma state factors
    scn_fac : ScreenFactors
        the precomputed ion pair factors

    Returns
    -------
    float or Dual
        Logarithm of the screening enhancement, with derivatives when the
        state was constructed with make_plasma_state_with_derivatives.

    """
    z1z2 = scn_fac.z1 * scn_fac.z2

    # Gamma_e from eq. 6
    Gamma_e = state.gamma_e_fac / state.temp

    # eq. A1
    h_DH = z1z2 * dual.sqrt(3 * Gamma_e**3 * state.z2bar / state.zbar)

    return h_DH


@njit
def screen5(state: PlasmaState, scn_fac):
    """Calculate screening factors following the appendix of
    :cite:t:`Wallace:1982`.

    Based on :cite:t:`graboske:1973` for weak screening. Based on
    :cite:t:`alastuey:1978` with plasma parameters from :cite:t:`itoh:1979`,
    for strong screening.

    Parameters
    ----------
    state : PlasmaState
        the precomputed plasma state factors
    scn_fac : ScreenFactors
        the precomputed ion pair factors

    Returns
    -------
    float or Dual
        Logarithm of the screening enhancement, with derivatives when the
        state was constructed with make_plasma_state_with_derivatives.

    """
    fact = dual.cbrt(2)
    gamefx = 0.3e0  # lower gamma limit for intermediate screening
    gamefs = 0.8e0  # upper gamma limit for intermediate screening

    # Get the ion data based on the input index
    z1 = scn_fac.z1
    z2 = scn_fac.z2

    # calculate individual screening factors
    bb = z1 * z2
    gamp = state.aa

    # In Eq.4 in Itoh:1979, this term is 2*Z_1*Z_2/(Z_1^(1/3) + Z_2^(1/3))
    # However here we follow Wallace:1982 Eq. A13, which is Z_1*Z_2*(2/(Z_1+Z_2))^(1/3)

    qq = fact * bb / scn_fac.zs13

    # Full Equation of Wallace:1982 Eq. A13

    gamef = qq * gamp

    # Full version of Eq.6 in Itoh:1979 with extra 1/3 factor
    # the extra 1/3 factor is there for convenience.
    # tau12 = Eq.6 / 3

    tau12 = state.taufac * scn_fac.aznut

    # alph12 = 3*gamma_ij/tau_ij

    alph12 = gamef / tau12

    # limit alph12 to 1.6 to prevent unphysical behavior.
    # See Introduction in Alastuey:1978

    # this should really be replaced by a pycnonuclear reaction rate formula
    if alph12 > 1.6:
        alph12 = dual.constant_like(alph12, 1.6e0)

        # redetermine previous factors if 3*gamma_ij/tau_ij > 1.6

        gamef = dual.constant_like(gamef, 1.6e0 * tau12)

        gamp = gamef * scn_fac.zs13/(fact * bb)

    # weak screening regime
    # Full version of Eq. 19 in Graboske:1973 by considering weak regime
    # and Wallace:1982 Eq. A14. Here the degeneracy factor is assumed to be 1.

    h12w = bb * state.qlam0z

    h12 = h12w

    # intermediate and strong sceening regime

    if gamef > gamefx:

        # gamma_ij^(1/4)
        gamp14 = gamp**0.25

        # Here we follow Eq. A9 in Wallace:1982
        # See Eq. 25 Alastuey:1978, Eq. 16 and 17 in Jancovici:1977 for reference
        cc = (0.896434e0 * gamp * scn_fac.zhat +
              -3.44740e0 * gamp14 * scn_fac.zhat2 +
              -0.5551e0 * (dual.log(gamp) + scn_fac.lzav) +
              -2.996e0)

        # (3gamma_ij/tau_ij)^3
        a3 = alph12 * alph12 * alph12

        # Part of Eq. 28 in Alastuey:1978
        qq = 0.014e0 + 0.0128e0*alph12

        # Part of Eq. 28 in Alastuey:1978
        rr = (5.0/32.0) - alph12*qq

        # Part of Eq. 28 in Alastuey:1978
        ss = tau12*rr

        # Part of Eq. 31 in Alastuey:1978
        tt = -0.0098e0 + 0.0048e0*alph12

        # Part of Eq. 31 in Alastuey:1978
        uu = 0.0055e0 + alph12*tt

        # Part of Eq. 31 in Alastuey:1978
        vv = gamef * alph12 * uu

        # Exponent of Eq. 32 in Alastuey:1978, which uses Eq.28 and Eq.31
        # Strong screening factor
        h12 = cc - a3 * (ss + vv)

        # See conclusion and Eq. 34 in Alastuey:1978
        # This is an extra factor to account for quantum effects
        rr = 1.0 - 0.0562e0*a3

        # In extreme case, rr is 0.77, see conclusion in Alastuey:1978
        xlgfac = max(dual.constant_like(rr, 0.77), rr)

        # Include the extra factor that accounts for quantum effects
        h12 = h12 + dual.log(xlgfac)

        # If gamma_ij < upper limit of intermediate regime
        # then it is in the intermediate regime, else strong screening.
        if gamef <= gamefs:
            dgamma = 1.0e0/(gamefs - gamefx)

            rr = dgamma*(gamefs - gamef)

            ss = dgamma*(gamef - gamefx)

            # Then the screening factor is a combination
            # of the strong and weak screening factor.
            h12 = h12w*rr + h12*ss

        # end of intermediate and strong screening

    # limit to avoid the pycnonuclear regime
    log_scor = max(h12, dual.constant_like(h12, 0.0))

    return log_scor


@njit
def smooth_clip(x, limit, start):
    """Create a smooth transition between y=limit and y=x with a
    half-cosine.

    Clips smaller values if limit < start and larger values if start <
    limit.

    Parameters
    ----------
    x : float
        the value to clip
    limit : float
        the constant value to clip x to
    start : float
        the x-value at which to start the transition

    Returns
    -------
    float

    """
    bound = dual.constant_like(x, limit)
    if limit < start:
        lower = bound
        upper = x
    else:
        lower = x
        upper = bound

    if x < min(limit, start):
        return lower
    if x > max(limit, start):
        return upper

    tmp = np.pi * (x - min(limit, start)) / (start - limit)
    f = (1 - dual.cos(tmp)) / 2

    return (1 - f) * lower + f * upper


@njit
def chugunov_2007(state, scn_fac):
    """Calculate screening factors based on :cite:t:`chugunov:2007`.

    Follows the approach in :cite:t:`yakovlev:2006` to extend to a
    multi-component plasma.

    Parameters
    ----------
    state : PlasmaState
        the precomputed plasma state factors
    scn_fac : ScreenFactors
        the precomputed ion pair factors

    Returns
    -------
    float or Dual
        Logarithm of the screening enhancement, with derivatives when the
        state was constructed with make_plasma_state_with_derivatives.

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
    n_i = state.n_e / scn_fac.ztilde**3
    m_i = 2 * mu12 * constants.m_u_C18

    T_p = constants.hbar / constants.k * constants.q_e * dual.sqrt(4 * np.pi * z_factor * n_i / m_i)

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
    zeta = dual.cbrt(4 / (3 * np.pi**2 * T_norm**2))

    # Gamma tilde from Chugunov 2007 eq. 21
    fit_alpha = 0.022
    fit_beta = 0.41 - 0.6 / Gamma
    fit_gamma = 0.06 + 2.2 / Gamma

    # Polynomial term in Gamma tilde
    poly = 1 + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta))

    gamtilde = Gamma / dual.cbrt(poly)

    # fit parameters just after Chugunov 2007 eq. 19
    A1 = 2.7822
    A2 = 98.34
    A3 = dual.sqrt(3) - A1 / dual.sqrt(A2)
    B1 = -1.7476
    B2 = 66.07
    B3 = 1.12
    B4 = 65
    gamtilde2 = gamtilde**2

    # Chugunov 2007 eq. 19
    term1 = 1 / dual.sqrt(A2 + gamtilde)
    term2 = 1 / (1 + gamtilde)
    term3 = gamtilde**2 / (B2 + gamtilde)
    term4 = gamtilde2 / (B4 + gamtilde2)

    h = gamtilde**(3 / 2) * (A1 * term1 + A3 * term2) + B1 * term3 + B3 * term4

    # Assume h >= 0
    log_scor = max(h, dual.constant_like(h, 0.0))

    return log_scor


@njit
def f0(gamma):
    r"""Calculate the free energy per ion in a OCP from
    :cite:t:`chugunov:2009` eq. 24

    Parameters
    ----------
    gamma : float
        Coulomb coupling parameter

    Returns
    -------
    float

    """
    A1 = -0.907
    A2 = 0.62954
    A3 = -dual.sqrt(3) / 2 - A1 / dual.sqrt(A2)
    B1 = 0.00456
    B2 = 211.6
    B3 = -1e-4
    B4 = 0.00462

    term1 = dual.sqrt(gamma * (A2 + gamma))
    term2 = dual.log(dual.sqrt(gamma / A2) + dual.sqrt(1 + gamma / A2))
    gamma_12 = dual.sqrt(gamma)
    term3 = gamma_12 - dual.arctan(gamma_12)
    term4 = dual.log1p(gamma / B2)
    term5 = dual.log1p(gamma**2 / B4)

    return (
        A1 * (term1 - A2 * term2) +
        2 * A3 * term3 +
        B1 * (gamma - B2 * term4) +
        B3 / 2 * term5
    )


@njit
def chugunov_2009(state, scn_fac):
    """Calculate screening factors based on :cite:t:`chugunov:2009`.

    Parameters
    ----------
    state : PlasmaState
        the precomputed plasma state factors
    scn_fac : ScreenFactors
        the precomputed ion pair factors

    Returns
    -------
    float or Dual
        Logarithm of the screening enhancement, with derivatives when the
        state was constructed with make_plasma_state_with_derivatives.

    """
    z1z2 = scn_fac.z1 * scn_fac.z2
    zcomp = scn_fac.z1 + scn_fac.z2

    # Gamma_e from eq. 6
    Gamma_e = state.gamma_e_fac / state.temp

    # Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    Gamma_1 = Gamma_e * scn_fac.z1**(5 / 3)
    Gamma_2 = Gamma_e * scn_fac.z2**(5 / 3)
    Gamma_comp = Gamma_e * zcomp**(5 / 3)

    Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde

    # Coulomb barrier penetrability, eq. 10
    tau_factor = dual.cbrt(27 / 2 * (np.pi * constants.q_e**2 / constants.hbar)**2 * constants.m_u_C18 / constants.k)
    tau_12 = tau_factor * scn_fac.aznut / dual.cbrt(state.temp)

    # eq. 12
    zeta = 3 * Gamma_12 / tau_12

    # additional fit parameters, eq. 25
    y_12 = 4 * z1z2 / zcomp**2
    c1 = 0.013 * y_12**2
    c2 = 0.406 * y_12**0.14
    c3 = 0.062 * y_12**0.19 + 1.8 / Gamma_12

    poly = 1 + zeta*(c1 + zeta*(c2 + c3*zeta))
    t_12 = dual.cbrt(poly)

    # strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    # Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    # values similar to those from Chugunov 2007.
    term1 = f0(Gamma_1 / t_12)
    term2 = f0(Gamma_2 / t_12)
    term3 = f0(Gamma_comp / t_12)
    h_fit = term1 + term2 - term3

    # weak screening correction term, eq. A3
    corr_C = (
        3 * z1z2 * dual.sqrt(state.z2bar / state.zbar) /
        (zcomp**2.5 - scn_fac.z1**2.5 - scn_fac.z2**2.5)
    )

    # corrected enhancement factor, eq. A4
    Gamma_12_2 = Gamma_12**2
    numer = corr_C + Gamma_12_2
    denom = 1 + Gamma_12_2
    h12 = numer / denom * h_fit

    # Assume h >= 0
    log_scor = max(h12, dual.constant_like(h12, 0.0))

    return log_scor


@njit
def potekhin_1998(state, scn_fac):
    """Calculate screening factors based on
    :cite:t:`chabrier_potekhin:1998`.

    Parameters
    ----------
    state : PlasmaState
        the precomputed plasma state factors
    scn_fac : ScreenFactors
        the precomputed ion pair factors

    Returns
    -------
    float or Dual
        Logarithm of the screening enhancement, with derivatives when the
        state was constructed with make_plasma_state_with_derivatives.

    """

    Gamma_e = state.gamma_e_fac / state.temp
    zcomp = scn_fac.z1 + scn_fac.z2

    Gamma_1 = Gamma_e * scn_fac.z1**(5 / 3)
    Gamma_2 = Gamma_e * scn_fac.z2**(5 / 3)
    Gamma_comp = Gamma_e * zcomp**(5 / 3)

    A_1 = -0.9052
    A_2 = 0.6322
    A_3 = -0.5 * dual.sqrt(3) - A_1/dual.sqrt(A_2)

    f1 = A_1 * (dual.sqrt(Gamma_1 * (A_2 + Gamma_1)) - A_2 * dual.log(dual.sqrt(Gamma_1 / A_2) +
                  dual.sqrt(1.0 + Gamma_1/A_2))) + 2.0 * A_3 * (dual.sqrt(Gamma_1) - dual.arctan(dual.sqrt(Gamma_1)))

    f2 = A_1 * (dual.sqrt(Gamma_2 * (A_2 + Gamma_2)) - A_2 * dual.log(dual.sqrt(Gamma_2 / A_2) +
                  dual.sqrt(1.0 + Gamma_2/A_2))) + 2.0 * A_3 * (dual.sqrt(Gamma_2) - dual.arctan(dual.sqrt(Gamma_2)))

    f12 = A_1 * (dual.sqrt(Gamma_comp * (A_2 + Gamma_comp)) - A_2 * dual.log(dual.sqrt(Gamma_comp / A_2) +
                  dual.sqrt(1.0 + Gamma_comp/A_2))) + 2.0 * A_3 * (dual.sqrt(Gamma_comp) - dual.arctan(dual.sqrt(Gamma_comp)))

    h12 = f1 + f2 - f12

    # Assume h >= 0
    log_scor = max(h12, dual.constant_like(h12, 0.0))

    return log_scor


def screening_check(check_func=debye_huckel, threshold: float = 1e-2):
    """Create a decorator factory that wraps a screening function with
    a check that determines whether that function can be skipped for a
    given plasma state and screening pair.

    :param func: the function to check against the threshold
    :param threshold: the threshold to check against. If screen_check
                      is less than the threshold, skip screen_func
    :returns: a decorator for wrapping screening functions

    """

    def screening_decorator(screen_func):
        """Decorate a screening function with a computation
        that determines whether the check is skippable.

        :param screen_func: the screening function being wrapped
        :returns: a wrapped screening function that performs this check
        """

        @njit
        def screening_wrapper(state, scn_fac):
            F0 = check_func(state, scn_fac)
            if F0 <= threshold:
                return F0
            return screen_func(state, scn_fac)
        return screening_wrapper
    return screening_decorator


# A list of all different screening functions
SCREEN_METHODS = {
    "screen5": screen5,
    "chugunov_2007": chugunov_2007,
    "chugunov_2009": chugunov_2009,
    "potekhin_1998": potekhin_1998,
    "debye_huckel": debye_huckel
}


def get_screening_func(screen_method):
    """Return the screening function by its name

    Parameters
    ----------
    screen_method : str
        name of the screening function

    Returns
    -------
    Callable
    """

    if screen_method is None:
        return None

    try:
        return SCREEN_METHODS[screen_method]
    except KeyError as exc:
        raise ValueError(
            f"{screen_method} is not a valid screening method. "
            f"Choose from {list(SCREEN_METHODS)}"
        ) from exc
