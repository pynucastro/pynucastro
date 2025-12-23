import numba
import numpy as np
from pynucastro.constants import constants
from numba.experimental import jitclass

from pynucastro.rates import (TableIndex, TableInterpolator, TabularRate,
                              TempTableInterpolator, TemperatureTabularRate,
                              Tfactors)
from pynucastro.screening import PlasmaState, ScreenFactors

jn = 0
jp = 1
jhe4 = 2
jc12 = 3
jo16 = 4
jne20 = 5
jne23 = 6
jna23 = 7
jmg23 = 8
nnuc = 9

A = np.zeros((nnuc), dtype=np.int32)

A[jn] = 1
A[jp] = 1
A[jhe4] = 4
A[jc12] = 12
A[jo16] = 16
A[jne20] = 20
A[jne23] = 23
A[jna23] = 23
A[jmg23] = 23

Z = np.zeros((nnuc), dtype=np.int32)

Z[jn] = 0
Z[jp] = 1
Z[jhe4] = 2
Z[jc12] = 6
Z[jo16] = 8
Z[jne20] = 10
Z[jne23] = 10
Z[jna23] = 11
Z[jmg23] = 12

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jn] = 0.001505349762871528
mass[jp] = 0.0015040963047307696
mass[jhe4] = 0.0059735574859708365
mass[jc12] = 0.017909017027273523
mass[jo16] = 0.023871099855618198
mass[jne20] = 0.029837079292893483
mass[jne23] = 0.03431735827046045
mass[jna23] = 0.03431034746033777
mass[jmg23] = 0.03431684618276469

names = []
names.append("n")
names.append("H1")
names.append("He4")
names.append("C12")
names.append("O16")
names.append("Ne20")
names.append("Ne23")
names.append("Na23")
names.append("Mg23")

def to_composition(Y):
    """Convert an array of molar fractions to a Composition object."""
    from pynucastro import Composition, Nucleus
    nuclei = [Nucleus.from_cache(name) for name in names]
    comp = Composition(nuclei)
    for i, nuc in enumerate(nuclei):
        comp.X[nuc] = Y[i] * A[i]
    return comp


def energy_release(dY):
    """return the energy release in erg/g (/s if dY is actually dY/dt)"""
    enuc = 0.0
    for i, y in enumerate(dY):
        enuc += y * mass[i]
    enuc *= -1*constants.N_A
    return enuc

@jitclass([
    ("C12_C12_to_He4_Ne20_reaclib", numba.float64),
    ("C12_C12_to_n_Mg23_reaclib", numba.float64),
    ("C12_C12_to_p_Na23_reaclib", numba.float64),
    ("He4_C12_to_O16_reaclib", numba.float64),
    ("n_to_p_reaclib", numba.float64),
    ("He4_He4_He4_to_C12_reaclib", numba.float64),
    ("Na23_to_Ne23_weaktab", numba.float64),
    ("Ne23_to_Na23_weaktab", numba.float64),
])
class RateEval:
    def __init__(self):
        self.C12_C12_to_He4_Ne20_reaclib = np.nan
        self.C12_C12_to_n_Mg23_reaclib = np.nan
        self.C12_C12_to_p_Na23_reaclib = np.nan
        self.He4_C12_to_O16_reaclib = np.nan
        self.n_to_p_reaclib = np.nan
        self.He4_He4_He4_to_C12_reaclib = np.nan
        self.Na23_to_Ne23_weaktab = np.nan
        self.Ne23_to_Na23_weaktab = np.nan

# note: we cannot make the TableInterpolator global, since numba doesn't like global jitclass
# load data for Na23 --> Ne23
Na23_to_Ne23_weaktab_rate = TabularRate(rfile='suzuki-23na-23ne_electroncapture.dat')
Na23_to_Ne23_weaktab_info = (Na23_to_Ne23_weaktab_rate.table_rhoy_lines,
                  Na23_to_Ne23_weaktab_rate.table_temp_lines,
                  Na23_to_Ne23_weaktab_rate.tabular_data_table)

# load data for Ne23 --> Na23
Ne23_to_Na23_weaktab_rate = TabularRate(rfile='suzuki-23ne-23na_betadecay.dat')
Ne23_to_Na23_weaktab_info = (Ne23_to_Na23_weaktab_rate.table_rhoy_lines,
                  Ne23_to_Na23_weaktab_rate.table_temp_lines,
                  Ne23_to_Na23_weaktab_rate.tabular_data_table)

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def C12_C12_to_He4_Ne20_reaclib(rate_eval, tf):
    # C12 + C12 --> He4 + Ne20
    rate = 0.0

    # cf88r
    rate += np.exp(  61.2863 + -84.165*tf.T913i + -1.56627*tf.T913
                  + -0.0736084*tf.T9 + -0.072797*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.C12_C12_to_He4_Ne20_reaclib = rate

@numba.njit()
def C12_C12_to_n_Mg23_reaclib(rate_eval, tf):
    # C12 + C12 --> n + Mg23
    rate = 0.0

    # cf88r
    rate += np.exp(  -12.8056 + -30.1498*tf.T9i + 11.4826*tf.T913
                  + 1.82849*tf.T9 + -0.34844*tf.T953)

    rate_eval.C12_C12_to_n_Mg23_reaclib = rate

@numba.njit()
def C12_C12_to_p_Na23_reaclib(rate_eval, tf):
    # C12 + C12 --> p + Na23
    rate = 0.0

    # cf88r
    rate += np.exp(  60.9649 + -84.165*tf.T913i + -1.4191*tf.T913
                  + -0.114619*tf.T9 + -0.070307*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.C12_C12_to_p_Na23_reaclib = rate

@numba.njit()
def He4_C12_to_O16_reaclib(rate_eval, tf):
    # C12 + He4 --> O16
    rate = 0.0

    # nac2 
    rate += np.exp(  254.634 + -1.84097*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 137.303*tf.lnT9)
    # nac2 
    rate += np.exp(  69.6526 + -1.39254*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 70.3554*tf.lnT9)

    rate_eval.He4_C12_to_O16_reaclib = rate

@numba.njit()
def n_to_p_reaclib(rate_eval, tf):
    # n --> p
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.78161)

    rate_eval.n_to_p_reaclib = rate

@numba.njit()
def He4_He4_He4_to_C12_reaclib(rate_eval, tf):
    # 3 He4 --> C12
    rate = 0.0

    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)
    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)
    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)

    rate_eval.He4_He4_He4_to_C12_reaclib = rate

@numba.njit()
def Na23_to_Ne23_weaktab(rate_eval, T, rho, Y):
    # Na23 --> Ne23
    rhoY = rho * ye(Y)
    Na23_to_Ne23_weaktab_interpolator = TableInterpolator(*Na23_to_Ne23_weaktab_info)
    r = Na23_to_Ne23_weaktab_interpolator.interpolate(np.log10(rhoY), np.log10(T), TableIndex.RATE.value)
    rate_eval.Na23_to_Ne23_weaktab = 10.0**r

@numba.njit()
def Ne23_to_Na23_weaktab(rate_eval, T, rho, Y):
    # Ne23 --> Na23
    rhoY = rho * ye(Y)
    Ne23_to_Na23_weaktab_interpolator = TableInterpolator(*Ne23_to_Na23_weaktab_info)
    r = Ne23_to_Na23_weaktab_interpolator.interpolate(np.log10(rhoY), np.log10(T), TableIndex.RATE.value)
    rate_eval.Ne23_to_Na23_weaktab = 10.0**r

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    C12_C12_to_He4_Ne20_reaclib(rate_eval, tf)
    C12_C12_to_n_Mg23_reaclib(rate_eval, tf)
    C12_C12_to_p_Na23_reaclib(rate_eval, tf)
    He4_C12_to_O16_reaclib(rate_eval, tf)
    n_to_p_reaclib(rate_eval, tf)
    He4_He4_He4_to_C12_reaclib(rate_eval, tf)

    # tabular rates
    Na23_to_Ne23_weaktab(rate_eval, T, rho=rho, Y=Y)
    Ne23_to_Na23_weaktab(rate_eval, T, rho=rho, Y=Y)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12_to_He4_Ne20_reaclib *= scor
        rate_eval.C12_C12_to_n_Mg23_reaclib *= scor
        rate_eval.C12_C12_to_p_Na23_reaclib *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12_to_O16_reaclib *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4_to_C12_reaclib *= scor * scor2

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jn] = (
          -Y[jn]*rate_eval.n_to_p_reaclib  +
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_n_Mg23_reaclib
       )

    dYdt[jp] = (
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_p_Na23_reaclib  +
          +Y[jn]*rate_eval.n_to_p_reaclib
       )

    dYdt[jhe4] = (
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_He4_Ne20_reaclib  +
          -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib  +
          + -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4_to_C12_reaclib
       )

    dYdt[jc12] = (
          + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_He4_Ne20_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_p_Na23_reaclib  +
          -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib  +
          +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4_to_C12_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_n_Mg23_reaclib
       )

    dYdt[jo16] = (
          +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib
       )

    dYdt[jne20] = (
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_He4_Ne20_reaclib
       )

    dYdt[jne23] = (
          ( -Y[jne23]*rate_eval.Ne23_to_Na23_weaktab +Y[jna23]*rate_eval.Na23_to_Ne23_weaktab )
       )

    dYdt[jna23] = (
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_p_Na23_reaclib  +
          ( +Y[jne23]*rate_eval.Ne23_to_Na23_weaktab -Y[jna23]*rate_eval.Na23_to_Ne23_weaktab )
       )

    dYdt[jmg23] = (
          +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12_to_n_Mg23_reaclib
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    C12_C12_to_He4_Ne20_reaclib(rate_eval, tf)
    C12_C12_to_n_Mg23_reaclib(rate_eval, tf)
    C12_C12_to_p_Na23_reaclib(rate_eval, tf)
    He4_C12_to_O16_reaclib(rate_eval, tf)
    n_to_p_reaclib(rate_eval, tf)
    He4_He4_He4_to_C12_reaclib(rate_eval, tf)

    # tabular rates
    Na23_to_Ne23_weaktab(rate_eval, T, rho=rho, Y=Y)
    Ne23_to_Na23_weaktab(rate_eval, T, rho=rho, Y=Y)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12_to_He4_Ne20_reaclib *= scor
        rate_eval.C12_C12_to_n_Mg23_reaclib *= scor
        rate_eval.C12_C12_to_p_Na23_reaclib *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12_to_O16_reaclib *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4_to_C12_reaclib *= scor * scor2

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jn, jn] = (
       -rate_eval.n_to_p_reaclib
       )

    jac[jn, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_n_Mg23_reaclib
       )

    jac[jp, jn] = (
       +rate_eval.n_to_p_reaclib
       )

    jac[jp, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_p_Na23_reaclib
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4_to_C12_reaclib
       )

    jac[jhe4, jc12] = (
       -rho*Y[jhe4]*rate_eval.He4_C12_to_O16_reaclib
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_He4_Ne20_reaclib
       )

    jac[jc12, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4_to_C12_reaclib
       )

    jac[jc12, jc12] = (
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_He4_Ne20_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_n_Mg23_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_p_Na23_reaclib
       -rho*Y[jhe4]*rate_eval.He4_C12_to_O16_reaclib
       )

    jac[jo16, jhe4] = (
       +rho*Y[jc12]*rate_eval.He4_C12_to_O16_reaclib
       )

    jac[jo16, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12_to_O16_reaclib
       )

    jac[jne20, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_He4_Ne20_reaclib
       )

    jac[jne23, jne23] = (
       -rate_eval.Ne23_to_Na23_weaktab
       )

    jac[jne23, jna23] = (
       +rate_eval.Na23_to_Ne23_weaktab
       )

    jac[jna23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_p_Na23_reaclib
       )

    jac[jna23, jne23] = (
       +rate_eval.Ne23_to_Na23_weaktab
       )

    jac[jna23, jna23] = (
       -rate_eval.Na23_to_Ne23_weaktab
       )

    jac[jmg23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12_to_n_Mg23_reaclib
       )

    return jac
