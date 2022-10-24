import numba
import numpy as np
from numba.experimental import jitclass

from pynucastro.rates import Tfactors, _find_rate_file
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

names = []
names.append("n")
names.append("h1")
names.append("he4")
names.append("c12")
names.append("o16")
names.append("ne20")
names.append("ne23")
names.append("na23")
names.append("mg23")

@jitclass([
    ("c12_c12__he4_ne20", numba.float64),
    ("c12_c12__n_mg23", numba.float64),
    ("c12_c12__p_na23", numba.float64),
    ("he4_c12__o16", numba.float64),
    ("n__p__weak__wc12", numba.float64),
    ("he4_he4_he4__c12", numba.float64),
    ("na23__ne23", numba.float64),
    ("ne23__na23", numba.float64),
])
class RateEval:
    def __init__(self):
        self.c12_c12__he4_ne20 = np.nan
        self.c12_c12__n_mg23 = np.nan
        self.c12_c12__p_na23 = np.nan
        self.he4_c12__o16 = np.nan
        self.n__p__weak__wc12 = np.nan
        self.he4_he4_he4__c12 = np.nan
        self.na23__ne23 = np.nan
        self.ne23__na23 = np.nan

# load data for na23 --> ne23
na23__ne23_table_path = _find_rate_file('23Na-23Ne_electroncapture.dat')
t_data2d = []
with open(na23__ne23_table_path) as tabular_file:
    for i, line in enumerate(tabular_file):
        if i < 6:
            continue
        line = line.strip()
        if not line:
            continue
        t_data2d.append(line.split())
na23__ne23_data = np.array(t_data2d, dtype=float)

# load data for ne23 --> na23
ne23__na23_table_path = _find_rate_file('23Ne-23Na_betadecay.dat')
t_data2d = []
with open(ne23__na23_table_path) as tabular_file:
    for i, line in enumerate(tabular_file):
        if i < 5:
            continue
        line = line.strip()
        if not line:
            continue
        t_data2d.append(line.split())
ne23__na23_data = np.array(t_data2d, dtype=float)

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def c12_c12__he4_ne20(rate_eval, tf):
    # c12 + c12 --> he4 + ne20
    rate = 0.0

    # cf88r
    rate += np.exp(  61.2863 + -84.165*tf.T913i + -1.56627*tf.T913
                  + -0.0736084*tf.T9 + -0.072797*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.c12_c12__he4_ne20 = rate

@numba.njit()
def c12_c12__n_mg23(rate_eval, tf):
    # c12 + c12 --> n + mg23
    rate = 0.0

    # cf88r
    rate += np.exp(  -12.8056 + -30.1485*tf.T9i + 11.4826*tf.T913
                  + 1.82849*tf.T9 + -0.34844*tf.T953)

    rate_eval.c12_c12__n_mg23 = rate

@numba.njit()
def c12_c12__p_na23(rate_eval, tf):
    # c12 + c12 --> p + na23
    rate = 0.0

    # cf88r
    rate += np.exp(  60.9649 + -84.165*tf.T913i + -1.4191*tf.T913
                  + -0.114619*tf.T9 + -0.070307*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.c12_c12__p_na23 = rate

@numba.njit()
def he4_c12__o16(rate_eval, tf):
    # c12 + he4 --> o16
    rate = 0.0

    # nac2 
    rate += np.exp(  69.6526 + -1.39254*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 70.3554*tf.lnT9)
    # nac2 
    rate += np.exp(  254.634 + -1.84097*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 137.303*tf.lnT9)

    rate_eval.he4_c12__o16 = rate

@numba.njit()
def n__p__weak__wc12(rate_eval, tf):
    # n --> p
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.78161)

    rate_eval.n__p__weak__wc12 = rate

@numba.njit()
def he4_he4_he4__c12(rate_eval, tf):
    # he4 + he4 + he4 --> c12
    rate = 0.0

    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)
    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)

    rate_eval.he4_he4_he4__c12 = rate

@numba.njit()
def na23__ne23(rate_eval, T, rhoY):
    # na23 --> ne23
    T_nearest = (na23__ne23_data[:, 1])[np.abs((na23__ne23_data[:, 1]) - T).argmin()]
    rhoY_nearest = (na23__ne23_data[:, 0])[np.abs((na23__ne23_data[:, 0]) - rhoY).argmin()]
    inde = np.where((na23__ne23_data[:, 1] == T_nearest) & (na23__ne23_data[:, 0] == rhoY_nearest))[0][0]
    rate_eval.na23__ne23 = na23__ne23_data[inde][5]

@numba.njit()
def ne23__na23(rate_eval, T, rhoY):
    # ne23 --> na23
    T_nearest = (ne23__na23_data[:, 1])[np.abs((ne23__na23_data[:, 1]) - T).argmin()]
    rhoY_nearest = (ne23__na23_data[:, 0])[np.abs((ne23__na23_data[:, 0]) - rhoY).argmin()]
    inde = np.where((ne23__na23_data[:, 1] == T_nearest) & (ne23__na23_data[:, 0] == rhoY_nearest))[0][0]
    rate_eval.ne23__na23 = ne23__na23_data[inde][5]

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    c12_c12__he4_ne20(rate_eval, tf)
    c12_c12__n_mg23(rate_eval, tf)
    c12_c12__p_na23(rate_eval, tf)
    he4_c12__o16(rate_eval, tf)
    n__p__weak__wc12(rate_eval, tf)
    he4_he4_he4__c12(rate_eval, tf)

    # tabular rates
    na23__ne23(rate_eval, T, rho*ye(Y))
    ne23__na23(rate_eval, T, rho*ye(Y))

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.c12_c12__he4_ne20 *= scor
        rate_eval.c12_c12__n_mg23 *= scor
        rate_eval.c12_c12__p_na23 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.he4_c12__o16 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.he4_he4_he4__c12 *= scor * scor2

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jn] = (
       -Y[jn]*rate_eval.n__p__weak__wc12
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__n_mg23
       )

    dYdt[jp] = (
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__p_na23
       +Y[jn]*rate_eval.n__p__weak__wc12
       )

    dYdt[jhe4] = (
       -rho*Y[jhe4]*Y[jc12]*rate_eval.he4_c12__o16
       -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.he4_he4_he4__c12
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__he4_ne20
       )

    dYdt[jc12] = (
       -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__he4_ne20
       -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__n_mg23
       -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__p_na23
       -rho*Y[jhe4]*Y[jc12]*rate_eval.he4_c12__o16
       +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.he4_he4_he4__c12
       )

    dYdt[jo16] = (
       +rho*Y[jhe4]*Y[jc12]*rate_eval.he4_c12__o16
       )

    dYdt[jne20] = (
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__he4_ne20
       )

    dYdt[jne23] = (
       -Y[jne23]*rate_eval.ne23__na23
       +Y[jna23]*rate_eval.na23__ne23
       )

    dYdt[jna23] = (
       -Y[jna23]*rate_eval.na23__ne23
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__p_na23
       +Y[jne23]*rate_eval.ne23__na23
       )

    dYdt[jmg23] = (
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.c12_c12__n_mg23
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    c12_c12__he4_ne20(rate_eval, tf)
    c12_c12__n_mg23(rate_eval, tf)
    c12_c12__p_na23(rate_eval, tf)
    he4_c12__o16(rate_eval, tf)
    n__p__weak__wc12(rate_eval, tf)
    he4_he4_he4__c12(rate_eval, tf)

    # tabular rates
    na23__ne23(rate_eval, T, rho*ye(Y))
    ne23__na23(rate_eval, T, rho*ye(Y))

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.c12_c12__he4_ne20 *= scor
        rate_eval.c12_c12__n_mg23 *= scor
        rate_eval.c12_c12__p_na23 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.he4_c12__o16 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.he4_he4_he4__c12 *= scor * scor2

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jn, jn] = (
       -rate_eval.n__p__weak__wc12
       )

    jac[jn, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__n_mg23
       )

    jac[jp, jn] = (
       +rate_eval.n__p__weak__wc12
       )

    jac[jp, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__p_na23
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jc12]*rate_eval.he4_c12__o16
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.he4_he4_he4__c12
       )

    jac[jhe4, jc12] = (
       -rho*Y[jhe4]*rate_eval.he4_c12__o16
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__he4_ne20
       )

    jac[jc12, jhe4] = (
       -rho*Y[jc12]*rate_eval.he4_c12__o16
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.he4_he4_he4__c12
       )

    jac[jc12, jc12] = (
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__he4_ne20
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__n_mg23
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__p_na23
       -rho*Y[jhe4]*rate_eval.he4_c12__o16
       )

    jac[jo16, jhe4] = (
       +rho*Y[jc12]*rate_eval.he4_c12__o16
       )

    jac[jo16, jc12] = (
       +rho*Y[jhe4]*rate_eval.he4_c12__o16
       )

    jac[jne20, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__he4_ne20
       )

    jac[jne23, jne23] = (
       -rate_eval.ne23__na23
       )

    jac[jne23, jna23] = (
       +rate_eval.na23__ne23
       )

    jac[jna23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__p_na23
       )

    jac[jna23, jne23] = (
       +rate_eval.ne23__na23
       )

    jac[jna23, jna23] = (
       -rate_eval.na23__ne23
       )

    jac[jmg23, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.c12_c12__n_mg23
       )

    return jac
