"""
Python implementations of screening routines.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import constants

from pynucastro.nucleus import Nucleus

__all__ = ["PlasmaState", "ScreenFactors", "chugunov_2007", "chugunov_2009"]


amu = constants.value("atomic mass constant") / constants.gram  # kg to g
q_e = constants.value("elementary charge") * (constants.c * 100) / 10  # C to statC (esu)
hbar = constants.value("reduced Planck constant") / constants.erg  # J*s to erg*s
k_B = constants.value("Boltzmann constant") / constants.erg  # J/K to erg/K
n_A = constants.value("Avogadro constant")


@dataclass
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
    n_e: float = field(init=False)
    gamma_e_fac: float = field(init=False)

    def __post_init__(self):
        # Average mass and total number density
        mbar = self.abar * amu
        ntot = self.dens / mbar
        # Electron number density
        # zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
        self.n_e = self.zbar * ntot

        # temperature-independent part of Gamma_e, from Chugunov 2009 eq. 6
        self.gamma_e_fac = q_e ** 2 / k_B * np.cbrt(4 * np.pi / 3) * np.cbrt(self.n_e)

    @classmethod
    def fill(cls, temp, dens, molar_fractions):
        """Construct a PlasmaState object from simulation data.

        :param temp:            temperature in K
        :param dens:            density in g/cm^3
        :param molar_fractions: dictionary of molar abundances for each ion,
                                as returned by :meth:`.Composition.get_molar`
        """
        ytot = sum(y for y in molar_fractions.values())
        abar = 1 / ytot

        zbar = sum(n.Z * y for n, y in molar_fractions.items()) / ytot
        z2bar = sum(n.Z ** 2 * y for n, y in molar_fractions.items()) / ytot
        return cls(temp, dens, abar, zbar, z2bar)


@dataclass
class ScreenFactors:
    """
    Stores values that will be used to calculate the screening correction factor
    for a specific pair of nuclei.

    :var n1: first nucleus
    :var n2: second nucleus
    :var aznut: combination of a1, z1, a2, z2 raised to 1/3 power
    :var ztilde: effective ion radius factor for a MCP
    """
    n1: Nucleus
    n2: Nucleus
    aznut: float = field(init=False)
    ztilde: float = field(init=False)

    def __post_init__(self):
        self.aznut = np.cbrt(
            self.n1.Z ** 2 * self.n2.Z ** 2 * self.n1.A * self.n2.A / (self.n1.A + self.n2.A)
        )
        self.ztilde = 0.5 * (np.cbrt(self.n1.Z) + np.cbrt(self.n2.Z))


def smooth_clip(x, dx_dT, limit, start):
    """Smoothly transition between y=limit and y=x with a half-cosine.

    Clips smaller values if limit < start and larger values if start < limit.

    :param x:     the value to clip
    :param dx_dT: derivative of x w.r.t. temperature
    :param limit: the constant value to clip x to
    :param start: the x-value at which to start the transition
    :returns: (y, dy/dT)
    """
    if limit < start:
        lower = limit
        dlower_dT = 0.0
        upper = x
        dupper_dT = dx_dT
    else:
        lower = x
        dlower_dT = dx_dT
        upper = limit
        dupper_dT = 0.0

    if x < min(limit, start):
        return lower, dlower_dT
    if x > max(limit, start):
        return upper, dupper_dT

    delta = start - limit
    tmp = np.pi * (x - min(limit, start)) / delta
    dtmp_dT = np.pi / delta * dx_dT

    f = (1 - np.cos(tmp)) / 2
    df_dT = np.sin(tmp) * dtmp_dT / 2

    tmp = (1 - f) * lower + f * upper
    dtmp_dT = dlower_dT + df_dT * (upper - lower) + f * (dupper_dT - dlower_dT)
    return tmp, dtmp_dT


def chugunov_2007(state, scn_fac):
    """Calculates screening factors based on Chugunov et al. 2007.

    Follows the approach in Yakovlev 2006 to extend to a multi-component plasma.

    :param state:   the precomputed :class:`PlasmaState`
    :param scn_fac: a :class:`ScreenFactors` object
    :returns: (screening correction factor, derivative w.r.t. temperature)

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
    mu12 = scn_fac.n1.A * scn_fac.n2.A / (scn_fac.n1.A + scn_fac.n2.A)
    z_factor = scn_fac.n1.Z * scn_fac.n2.Z
    n_i = state.n_e / scn_fac.ztilde ** 3
    m_i = 2 * mu12 * amu

    T_p = hbar / k_B * q_e * np.sqrt(4 * np.pi * z_factor * n_i / m_i)

    # Normalized temperature
    T_norm = state.temp / T_p
    dT_norm_dT = 1 / T_p

    # The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    # should be nearly temperature-independent (in the pycnonuclear regime),
    # and we clip the temperature to 0.1 T_p at small T.
    # start the transition here
    T_norm_fade = 0.2
    # minimum value of T/T_p
    T_norm_min = 0.1

    T_norm, dT_norm_dT = smooth_clip(
        T_norm, dT_norm_dT, limit=T_norm_min, start=T_norm_fade
    )

    # Coulomb coupling parameter from Yakovlev 2006, eq. 10
    Gamma = state.gamma_e_fac * scn_fac.n1.Z * scn_fac.n2.Z / (scn_fac.ztilde * T_norm * T_p)
    dGamma_dT = -Gamma / T_norm * dT_norm_dT

    # The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    Gamma_fade = 590
    Gamma_max = 600
    Gamma, dGamma_dT = smooth_clip(Gamma, dGamma_dT, limit=Gamma_max, start=Gamma_fade)

    # Chugunov 2007 eq. 3
    zeta = np.cbrt(4 / (3 * np.pi ** 2 * T_norm ** 2))
    dzeta_dT = -2 / (3 * T_norm) * zeta * dT_norm_dT

    # Gamma tilde from Chugunov 2007 eq. 21
    fit_alpha = 0.022
    fit_beta = 0.41 - 0.6 / Gamma
    tmp = dGamma_dT / Gamma ** 2
    dfit_beta_dT = 0.6 * tmp
    fit_gamma = 0.06 + 2.2 / Gamma
    dfit_gamma_dT = -2.2 * tmp

    # Polynomial term in Gamma tilde
    poly = 1 + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta))
    dpoly_dT = (
        (fit_alpha + 2 * zeta * fit_beta + 3 * fit_gamma * zeta ** 2) * dzeta_dT +
        zeta ** 2 * (dfit_beta_dT + dfit_gamma_dT * zeta)
    )

    gamtilde = Gamma / np.cbrt(poly)
    # this is gamtilde * dlog(gamtilde)/dT
    dgamtilde_dT = gamtilde * (dGamma_dT / Gamma - dpoly_dT / poly / 3)

    # fit parameters just after Chugunov 2007 eq. 19
    A1 = 2.7822
    A2 = 98.34
    A3 = np.sqrt(3) - A1 / np.sqrt(A2)
    B1 = -1.7476
    B2 = 66.07
    B3 = 1.12
    B4 = 65
    gamtilde2 = gamtilde ** 2
    dgamtilde2_dT = 2 * gamtilde * dgamtilde_dT

    # Chugunov 2007 eq. 19
    term1 = 1 / np.sqrt(A2 + gamtilde)
    dterm1_dT = -term1 / 2 / (A2 + gamtilde) * dgamtilde_dT

    term2 = 1 / (1 + gamtilde)
    dterm2_dT = -(term2 ** 2) * dgamtilde_dT

    term3 = gamtilde ** 2 / (B2 + gamtilde)
    dterm3_dT = gamtilde * (2 * B2 + gamtilde) / (B2 + gamtilde) ** 2 * dgamtilde_dT

    term4 = gamtilde2 / (B4 + gamtilde2)
    dterm4_dT = B4 / (B4 + gamtilde2) ** 2 * dgamtilde2_dT

    inner = A1 * term1 + A3 * term2
    dinner_dT = A1 * dterm1_dT + A3 * dterm2_dT

    gamtilde32 = gamtilde ** (3 / 2)
    dgamtilde32_dT = 3 / 2 * np.sqrt(gamtilde) * dgamtilde_dT
    h = gamtilde32 * inner + B1 * term3 + B3 * term4
    dh_dT = (
        dgamtilde32_dT * inner +
        gamtilde32 * dinner_dT +
        B1 * dterm3_dT +
        B3 * dterm4_dT
    )

    # machine limit the output
    h_max = 300
    h = min(h, h_max)
    scor = np.exp(h)

    if h == h_max:
        dscor_dT = 0
    else:
        dscor_dT = scor * dh_dT

    return scor, dscor_dT


def f0(gamma, dlog_dT):
    r"""Calculate the free energy per ion in a OCP from Chugunov & DeWitt 2009 eq. 24

    :param gamma: Coulomb coupling parameter
    :param dlog_dT: log derivative of gamma w.r.t. temperature:
                    :math:`d(log(Gamma))/dT = dGamma/dT / Gamma`

    :returns: (free energy, derivative w.r.t. temperature)
    """
    A1 = -0.907
    A2 = 0.62954
    A3 = -np.sqrt(3) / 2 - A1 / np.sqrt(A2)
    B1 = 0.00456
    B2 = 211.6
    B3 = -1e-4
    B4 = 0.00462

    term1 = np.sqrt(gamma * (A2 + gamma))
    dterm1_dgamma = (A2 / 2 + gamma) / term1

    term2 = np.log(np.sqrt(gamma / A2) + np.sqrt(1 + gamma / A2))
    dterm2_dgamma = 1 / (2 * term1)

    gamma_12 = np.sqrt(gamma)
    term3 = gamma_12 - np.atan(gamma_12)
    dterm3_dgamma = gamma_12 / (1 + gamma) / 2

    term4 = np.log(1 + gamma / B2)
    dterm4_dgamma = 1 / (B2 + gamma)

    term5 = np.log(1 + gamma ** 2 / B4)
    dterm5_dgamma = 2 * gamma / (B4 + gamma ** 2)

    f = (
        A1 * (term1 - A2 * term2) +
        2 * A3 * term3 +
        B1 * (gamma - B2 * term4) +
        B3 / 2 * term5
    )
    df_dgamma = (
        A1 * (dterm1_dgamma - A2 * dterm2_dgamma) +
        2 * A3 * dterm3_dgamma +
        B1 * (1 - B2 * dterm4_dgamma) +
        B3 / 2 * dterm5_dgamma
    )
    df_dT = df_dgamma * dlog_dT * gamma
    return f, df_dT


def chugunov_2009(state, scn_fac):
    """Calculates screening factors based on Chugunov & DeWitt 2009.

    :param state:   the precomputed :class:`PlasmaState`
    :param scn_fac: a :class:`ScreenFactors` object
    :returns: (screening correction factor, derivative w.r.t. temperature)

    References:
        | Chugunov and DeWitt 2009, PhRvC, 80, 014611
    """
    z1z2 = scn_fac.n1.Z * scn_fac.n2.Z
    zcomp = scn_fac.n1.Z + scn_fac.n2.Z

    # Gamma_e from eq. 6
    Gamma_e = state.gamma_e_fac / state.temp
    # work in terms of log derivatives, since it's much simpler
    dlog_Gamma_dT = -1 / state.temp

    # Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    Gamma_1 = Gamma_e * scn_fac.n1.Z ** (5 / 3)
    Gamma_2 = Gamma_e * scn_fac.n2.Z ** (5 / 3)
    Gamma_comp = Gamma_e * zcomp ** (5 / 3)

    Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde

    # Coulomb barrier penetrability, eq. 10
    tau_factor = np.cbrt(27 / 2 * (np.pi * q_e ** 2 / hbar) ** 2 * amu / k_B)
    tau_12 = tau_factor * scn_fac.aznut / np.cbrt(state.temp)
    dlog_tau_12_dT = -1 / 3 / state.temp

    # eq. 12
    zeta = 3 * Gamma_12 / tau_12
    dzeta_dT = zeta * (dlog_Gamma_dT - dlog_tau_12_dT)

    # additional fit parameters, eq. 25
    y_12 = 4 * z1z2 / zcomp ** 2
    c1 = 0.013 * y_12 ** 2
    c2 = 0.406 * y_12 ** 0.14
    c3 = 0.062 * y_12 ** 0.19 + 1.8 / Gamma_12
    dc3_dT = -1.8 / Gamma_12 * dlog_Gamma_dT

    poly = 1 + zeta*(c1 + zeta*(c2 + c3*zeta))
    dpoly_dT = (c1 + zeta*(2*c2 + 3*c3*zeta)) * dzeta_dT + dc3_dT * zeta ** 3
    t_12 = np.cbrt(poly)
    dlog_t_12_dT = dpoly_dT / (3 * poly)

    # strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    # Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    # values similar to those from Chugunov 2007.
    dlog_dT = dlog_Gamma_dT - dlog_t_12_dT
    term1, dterm1_dT = f0(Gamma_1 / t_12, dlog_dT)
    term2, dterm2_dT = f0(Gamma_2 / t_12, dlog_dT)
    term3, dterm3_dT = f0(Gamma_comp / t_12, dlog_dT)
    h_fit = term1 + term2 - term3
    dh_fit_dT = dterm1_dT + dterm2_dT - dterm3_dT

    # weak screening correction term, eq. A3
    corr_C = (
        3 * z1z2 * np.sqrt(state.z2bar / state.zbar) /
        (zcomp ** 2.5 - scn_fac.n1.Z ** 2.5 - scn_fac.n2.Z ** 2.5)
    )

    # corrected enhancement factor, eq. A4
    Gamma_12_2 = Gamma_12 ** 2
    dGamma_12_2_dT = 2 * Gamma_12_2 * dlog_Gamma_dT
    numer = corr_C + Gamma_12_2
    denom = 1 + Gamma_12_2
    h12 = numer / denom * h_fit
    dh12_dT = h12 * (dGamma_12_2_dT/numer - dGamma_12_2_dT/denom + dh_fit_dT/h_fit)

    # machine limit the output
    h12_max = 300
    h12 = min(h12, h12_max)
    scor = np.exp(h12)
    if h12 == h12_max:
        dh12_dT = 0
    else:
        dscor_dT = scor * dh12_dT

    return scor, dscor_dT
