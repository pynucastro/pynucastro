import numba
import numpy as np
from scipy import constants
from numba.experimental import jitclass

from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors
from pynucastro.screening import PlasmaState, ScreenFactors

jp = 0
jhe4 = 1
jc12 = 2
jc13 = 3
jn13 = 4
jn14 = 5
jn15 = 6
jo14 = 7
jo15 = 8
jo16 = 9
jo17 = 10
jo18 = 11
jf17 = 12
jf18 = 13
jf19 = 14
jne18 = 15
jne19 = 16
jne20 = 17
jmg22 = 18
jmg24 = 19
jfe56 = 20
nnuc = 21

A = np.zeros((nnuc), dtype=np.int32)

A[jp] = 1
A[jhe4] = 4
A[jc12] = 12
A[jc13] = 13
A[jn13] = 13
A[jn14] = 14
A[jn15] = 15
A[jo14] = 14
A[jo15] = 15
A[jo16] = 16
A[jo17] = 17
A[jo18] = 18
A[jf17] = 17
A[jf18] = 18
A[jf19] = 19
A[jne18] = 18
A[jne19] = 19
A[jne20] = 20
A[jmg22] = 22
A[jmg24] = 24
A[jfe56] = 56

Z = np.zeros((nnuc), dtype=np.int32)

Z[jp] = 1
Z[jhe4] = 2
Z[jc12] = 6
Z[jc13] = 6
Z[jn13] = 7
Z[jn14] = 7
Z[jn15] = 7
Z[jo14] = 8
Z[jo15] = 8
Z[jo16] = 8
Z[jo17] = 8
Z[jo18] = 8
Z[jf17] = 9
Z[jf18] = 9
Z[jf19] = 9
Z[jne18] = 10
Z[jne19] = 10
Z[jne20] = 10
Z[jmg22] = 12
Z[jmg24] = 12
Z[jfe56] = 26

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jp] = 0.0015040963030260536
mass[jhe4] = 0.0059735574925878256
mass[jc12] = 0.017909017027273523
mass[jc13] = 0.019406441930882663
mass[jn13] = 0.01940999951603316
mass[jn14] = 0.020898440903103103
mass[jn15] = 0.0223864338056853
mass[jo14] = 0.020906683076491985
mass[jo15] = 0.02239084645968795
mass[jo16] = 0.023871099858982767
mass[jo17] = 0.02536981167252093
mass[jo18] = 0.02686227133140636
mass[jf17] = 0.025374234423440733
mass[jf18] = 0.026864924401329426
mass[jf19] = 0.028353560468882166
mass[jne18] = 0.026872045275379234
mass[jne19] = 0.02835875072008801
mass[jne20] = 0.02983707929641827
mass[jmg22] = 0.032832557028702955
mass[jmg24] = 0.03579570996619953
mass[jfe56] = 0.08347830935425092

names = []
names.append("H1")
names.append("He4")
names.append("C12")
names.append("C13")
names.append("N13")
names.append("N14")
names.append("N15")
names.append("O14")
names.append("O15")
names.append("O16")
names.append("O17")
names.append("O18")
names.append("F17")
names.append("F18")
names.append("F19")
names.append("Ne18")
names.append("Ne19")
names.append("Ne20")
names.append("Mg22")
names.append("Mg24")
names.append("Fe56")

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
    enuc *= -1*constants.Avogadro
    return enuc

@jitclass([
    ("N13__C13__weak__wc12", numba.float64),
    ("O14__N14__weak__wc12", numba.float64),
    ("O15__N15__weak__wc12", numba.float64),
    ("F17__O17__weak__wc12", numba.float64),
    ("F18__O18__weak__wc12", numba.float64),
    ("Ne18__F18__weak__wc12", numba.float64),
    ("Ne19__F19__weak__wc12", numba.float64),
    ("p_C12__N13", numba.float64),
    ("He4_C12__O16", numba.float64),
    ("p_C13__N14", numba.float64),
    ("p_N13__O14", numba.float64),
    ("p_N14__O15", numba.float64),
    ("He4_N14__F18", numba.float64),
    ("p_N15__O16", numba.float64),
    ("He4_N15__F19", numba.float64),
    ("He4_O14__Ne18", numba.float64),
    ("He4_O15__Ne19", numba.float64),
    ("p_O16__F17", numba.float64),
    ("He4_O16__Ne20", numba.float64),
    ("p_O17__F18", numba.float64),
    ("p_O18__F19", numba.float64),
    ("p_F17__Ne18", numba.float64),
    ("p_F18__Ne19", numba.float64),
    ("p_F19__Ne20", numba.float64),
    ("He4_Ne18__Mg22", numba.float64),
    ("He4_Ne20__Mg24", numba.float64),
    ("C12_C12__He4_Ne20", numba.float64),
    ("He4_N13__p_O16", numba.float64),
    ("p_N15__He4_C12", numba.float64),
    ("He4_O14__p_F17", numba.float64),
    ("C12_O16__He4_Mg24", numba.float64),
    ("p_O17__He4_N14", numba.float64),
    ("p_O18__He4_N15", numba.float64),
    ("p_F18__He4_O15", numba.float64),
    ("p_F19__He4_O16", numba.float64),
    ("p_Ne20__He4_F17", numba.float64),
    ("He4_He4_He4__C12", numba.float64),
])
class RateEval:
    def __init__(self):
        self.N13__C13__weak__wc12 = np.nan
        self.O14__N14__weak__wc12 = np.nan
        self.O15__N15__weak__wc12 = np.nan
        self.F17__O17__weak__wc12 = np.nan
        self.F18__O18__weak__wc12 = np.nan
        self.Ne18__F18__weak__wc12 = np.nan
        self.Ne19__F19__weak__wc12 = np.nan
        self.p_C12__N13 = np.nan
        self.He4_C12__O16 = np.nan
        self.p_C13__N14 = np.nan
        self.p_N13__O14 = np.nan
        self.p_N14__O15 = np.nan
        self.He4_N14__F18 = np.nan
        self.p_N15__O16 = np.nan
        self.He4_N15__F19 = np.nan
        self.He4_O14__Ne18 = np.nan
        self.He4_O15__Ne19 = np.nan
        self.p_O16__F17 = np.nan
        self.He4_O16__Ne20 = np.nan
        self.p_O17__F18 = np.nan
        self.p_O18__F19 = np.nan
        self.p_F17__Ne18 = np.nan
        self.p_F18__Ne19 = np.nan
        self.p_F19__Ne20 = np.nan
        self.He4_Ne18__Mg22 = np.nan
        self.He4_Ne20__Mg24 = np.nan
        self.C12_C12__He4_Ne20 = np.nan
        self.He4_N13__p_O16 = np.nan
        self.p_N15__He4_C12 = np.nan
        self.He4_O14__p_F17 = np.nan
        self.C12_O16__He4_Mg24 = np.nan
        self.p_O17__He4_N14 = np.nan
        self.p_O18__He4_N15 = np.nan
        self.p_F18__He4_O15 = np.nan
        self.p_F19__He4_O16 = np.nan
        self.p_Ne20__He4_F17 = np.nan
        self.He4_He4_He4__C12 = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def N13__C13__weak__wc12(rate_eval, tf):
    # N13 --> C13
    rate = 0.0

    # wc12w
    rate += np.exp(  -6.7601)

    rate_eval.N13__C13__weak__wc12 = rate

@numba.njit()
def O14__N14__weak__wc12(rate_eval, tf):
    # O14 --> N14
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.62354)

    rate_eval.O14__N14__weak__wc12 = rate

@numba.njit()
def O15__N15__weak__wc12(rate_eval, tf):
    # O15 --> N15
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.17053)

    rate_eval.O15__N15__weak__wc12 = rate

@numba.njit()
def F17__O17__weak__wc12(rate_eval, tf):
    # F17 --> O17
    rate = 0.0

    # wc12w
    rate += np.exp(  -4.53318)

    rate_eval.F17__O17__weak__wc12 = rate

@numba.njit()
def F18__O18__weak__wc12(rate_eval, tf):
    # F18 --> O18
    rate = 0.0

    # wc12w
    rate += np.exp(  -9.15982)

    rate_eval.F18__O18__weak__wc12 = rate

@numba.njit()
def Ne18__F18__weak__wc12(rate_eval, tf):
    # Ne18 --> F18
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.879336)

    rate_eval.Ne18__F18__weak__wc12 = rate

@numba.njit()
def Ne19__F19__weak__wc12(rate_eval, tf):
    # Ne19 --> F19
    rate = 0.0

    # wc12w
    rate += np.exp(  -3.21142)

    rate_eval.Ne19__F19__weak__wc12 = rate

@numba.njit()
def p_C12__N13(rate_eval, tf):
    # C12 + p --> N13
    rate = 0.0

    # ls09n
    rate += np.exp(  17.1482 + -13.692*tf.T913i + -0.230881*tf.T913
                  + 4.44362*tf.T9 + -3.15898*tf.T953 + -0.666667*tf.lnT9)
    # ls09r
    rate += np.exp(  17.5428 + -3.77849*tf.T9i + -5.10735*tf.T913i + -2.24111*tf.T913
                  + 0.148883*tf.T9 + -1.5*tf.lnT9)

    rate_eval.p_C12__N13 = rate

@numba.njit()
def He4_C12__O16(rate_eval, tf):
    # C12 + He4 --> O16
    rate = 0.0

    # nac2 
    rate += np.exp(  254.634 + -1.84097*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 137.303*tf.lnT9)
    # nac2 
    rate += np.exp(  69.6526 + -1.39254*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 70.3554*tf.lnT9)

    rate_eval.He4_C12__O16 = rate

@numba.njit()
def p_C13__N14(rate_eval, tf):
    # C13 + p --> N14
    rate = 0.0

    # nacrr
    rate += np.exp(  15.1825 + -13.5543*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.5155 + -13.72*tf.T913i + -0.450018*tf.T913
                  + 3.70823*tf.T9 + -1.70545*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9637 + -5.78147*tf.T9i + -0.196703*tf.T913
                  + 0.142126*tf.T9 + -0.0238912*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_C13__N14 = rate

@numba.njit()
def p_N13__O14(rate_eval, tf):
    # N13 + p --> O14
    rate = 0.0

    # lg06r
    rate += np.exp(  10.9971 + -6.12602*tf.T9i + 1.57122*tf.T913i
                  + -1.5*tf.lnT9)
    # lg06n
    rate += np.exp(  18.1356 + -15.1676*tf.T913i + 0.0955166*tf.T913
                  + 3.0659*tf.T9 + -0.507339*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_N13__O14 = rate

@numba.njit()
def p_N14__O15(rate_eval, tf):
    # N14 + p --> O15
    rate = 0.0

    # im05n
    rate += np.exp(  17.01 + -15.193*tf.T913i + -0.161954*tf.T913
                  + -7.52123*tf.T9 + -0.987565*tf.T953 + -0.666667*tf.lnT9)
    # im05r
    rate += np.exp(  6.73578 + -4.891*tf.T9i
                  + 0.0682*tf.lnT9)
    # im05r
    rate += np.exp(  7.65444 + -2.998*tf.T9i
                  + -1.5*tf.lnT9)
    # im05n
    rate += np.exp(  20.1169 + -15.193*tf.T913i + -4.63975*tf.T913
                  + 9.73458*tf.T9 + -9.55051*tf.T953 + 0.333333*tf.lnT9)

    rate_eval.p_N14__O15 = rate

@numba.njit()
def He4_N14__F18(rate_eval, tf):
    # N14 + He4 --> F18
    rate = 0.0

    # il10n
    rate += np.exp(  21.5339 + -36.2504*tf.T913i
                  + -5.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  13.8995 + -10.9656*tf.T9i + -5.6227*tf.T913i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  0.196838 + -5.16034*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_N14__F18 = rate

@numba.njit()
def p_N15__O16(rate_eval, tf):
    # N15 + p --> O16
    rate = 0.0

    # li10n
    rate += np.exp(  20.0176 + -15.24*tf.T913i + 0.334926*tf.T913
                  + 4.59088*tf.T9 + -4.78468*tf.T953 + -0.666667*tf.lnT9)
    # li10r
    rate += np.exp(  14.5444 + -10.2295*tf.T9i
                  + 0.0459037*tf.T9 + -1.5*tf.lnT9)
    # li10r
    rate += np.exp(  6.59056 + -2.92315*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_N15__O16 = rate

@numba.njit()
def He4_N15__F19(rate_eval, tf):
    # N15 + He4 --> F19
    rate = 0.0

    # il10r
    rate += np.exp(  -9.41892 + -4.17795*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  25.3916 + -36.2324*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -28.7989 + -4.19986*tf.T9i + 35.4292*tf.T913
                  + -5.5767*tf.T9 + 0.441293*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  3.5342 + -6.98462*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_N15__F19 = rate

@numba.njit()
def He4_O14__Ne18(rate_eval, tf):
    # O14 + He4 --> Ne18
    rate = 0.0

    # wh87r
    rate += np.exp(  -4.69948 + -12.159*tf.T9i
                  + 5.0*tf.lnT9)
    # wh87r
    rate += np.exp(  3.52636 + -22.61*tf.T9i
                  + -1.5*tf.lnT9)
    # wh87r
    rate += np.exp(  -2.15417 + -11.73*tf.T9i
                  + -1.5*tf.lnT9)
    # wh87n
    rate += np.exp(  26.4429 + -39.38*tf.T913i + -0.0772187*tf.T913
                  + -0.635361*tf.T9 + 0.106236*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O14__Ne18 = rate

@numba.njit()
def He4_O15__Ne19(rate_eval, tf):
    # O15 + He4 --> Ne19
    rate = 0.0

    # dc11r
    rate += np.exp(  -32.2496 + -4.20439*tf.T9i + -3.24609*tf.T913i + 44.4647*tf.T913
                  + -9.79962*tf.T9 + 0.841782*tf.T953 + -1.5*tf.lnT9)
    # dc11r
    rate += np.exp(  -0.0452465 + -5.88439*tf.T9i
                  + -1.5*tf.lnT9)
    # dc11n
    rate += np.exp(  26.2914 + -39.578*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O15__Ne19 = rate

@numba.njit()
def p_O16__F17(rate_eval, tf):
    # O16 + p --> F17
    rate = 0.0

    # ia08n
    rate += np.exp(  19.0904 + -16.696*tf.T913i + -1.16252*tf.T913
                  + 0.267703*tf.T9 + -0.0338411*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_O16__F17 = rate

@numba.njit()
def He4_O16__Ne20(rate_eval, tf):
    # O16 + He4 --> Ne20
    rate = 0.0

    # co10r
    rate += np.exp(  9.50848 + -12.7643*tf.T9i + -3.65925*tf.T913
                  + 0.714224*tf.T9 + -0.00107508*tf.T953 + -1.5*tf.lnT9)
    # co10r
    rate += np.exp(  3.88571 + -10.3585*tf.T9i
                  + -1.5*tf.lnT9)
    # co10n
    rate += np.exp(  23.903 + -39.7262*tf.T913i + -0.210799*tf.T913
                  + 0.442879*tf.T9 + -0.0797753*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O16__Ne20 = rate

@numba.njit()
def p_O17__F18(rate_eval, tf):
    # O17 + p --> F18
    rate = 0.0

    # il10n
    rate += np.exp(  15.8929 + -16.4035*tf.T913i + 4.31885*tf.T913
                  + -0.709921*tf.T9 + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  9.39048 + -6.22828*tf.T9i + 2.31435*tf.T913
                  + -0.302835*tf.T9 + 0.020133*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -13.077 + -0.746296*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O17__F18 = rate

@numba.njit()
def p_O18__F19(rate_eval, tf):
    # O18 + p --> F19
    rate = 0.0

    # il10r
    rate += np.exp(  -35.0079 + -0.244743*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.917 + -16.7246*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.26876 + -6.7253*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  5.07648 + -1.65681*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_O18__F19 = rate

@numba.njit()
def p_F17__Ne18(rate_eval, tf):
    # F17 + p --> Ne18
    rate = 0.0

    # cb09 
    rate += np.exp(  -7.84708 + -0.0323504*tf.T9i + -14.2191*tf.T913i + 34.0647*tf.T913
                  + -16.5698*tf.T9 + 2.48116*tf.T953 + -2.13376*tf.lnT9)
    # cb09 
    rate += np.exp(  27.5778 + -4.95969*tf.T9i + -21.3249*tf.T913i + -0.230774*tf.T913
                  + 0.917931*tf.T9 + -0.0440377*tf.T953 + -7.36014*tf.lnT9)

    rate_eval.p_F17__Ne18 = rate

@numba.njit()
def p_F18__Ne19(rate_eval, tf):
    # F18 + p --> Ne19
    rate = 0.0

    # il10n
    rate += np.exp(  57.4084 + -21.4023*tf.T913i + -93.766*tf.T913
                  + 179.258*tf.T9 + -202.561*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -5.85727 + -2.89147*tf.T9i + 13.1683*tf.T913
                  + -1.92023*tf.T9 + 0.16901*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -29.449 + -0.39895*tf.T9i + 22.4903*tf.T913
                  + 0.307872*tf.T9 + -0.296226*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_F18__Ne19 = rate

@numba.njit()
def p_F19__Ne20(rate_eval, tf):
    # F19 + p --> Ne20
    rate = 0.0

    # nacrr
    rate += np.exp(  -5.63093 + -7.74414*tf.T9i + 31.6442*tf.T913i + -58.6563*tf.T913
                  + 67.7365*tf.T9 + -22.9721*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  12.3816 + -1.71383*tf.T9i + -11.3832*tf.T913i + 5.47872*tf.T913
                  + -1.07203*tf.T9 + 0.11196*tf.T953 + -1.5*tf.lnT9)
    # nacrn
    rate += np.exp(  18.2807 + -18.116*tf.T913i + -1.4622*tf.T913
                  + 6.95113*tf.T9 + -2.90366*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F19__Ne20 = rate

@numba.njit()
def He4_Ne18__Mg22(rate_eval, tf):
    # Ne18 + He4 --> Mg22
    rate = 0.0

    # ths8r
    rate += np.exp(  32.8865 + -46.4859*tf.T913i + 0.956741*tf.T913
                  + -0.914402*tf.T9 + 0.0722478*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ne18__Mg22 = rate

@numba.njit()
def He4_Ne20__Mg24(rate_eval, tf):
    # Ne20 + He4 --> Mg24
    rate = 0.0

    # il10r
    rate += np.exp(  -38.7055 + -2.50605*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  24.5058 + -46.2525*tf.T913i + 5.58901*tf.T913
                  + 7.61843*tf.T9 + -3.683*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  -8.79827 + -12.7809*tf.T9i + 16.9229*tf.T913
                  + -2.57325*tf.T9 + 0.208997*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  1.98307 + -9.22026*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_Ne20__Mg24 = rate

@numba.njit()
def C12_C12__He4_Ne20(rate_eval, tf):
    # C12 + C12 --> He4 + Ne20
    rate = 0.0

    # cf88r
    rate += np.exp(  61.2863 + -84.165*tf.T913i + -1.56627*tf.T913
                  + -0.0736084*tf.T9 + -0.072797*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.C12_C12__He4_Ne20 = rate

@numba.njit()
def He4_N13__p_O16(rate_eval, tf):
    # N13 + He4 --> p + O16
    rate = 0.0

    # cf88n
    rate += np.exp(  40.4644 + -35.829*tf.T913i + -0.530275*tf.T913
                  + -0.982462*tf.T9 + 0.0808059*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_N13__p_O16 = rate

@numba.njit()
def p_N15__He4_C12(rate_eval, tf):
    # N15 + p --> He4 + C12
    rate = 0.0

    # nacrn
    rate += np.exp(  27.4764 + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.57522 + -1.1638*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.8972 + -7.406*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -4.87347 + -2.02117*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_N15__He4_C12 = rate

@numba.njit()
def He4_O14__p_F17(rate_eval, tf):
    # O14 + He4 --> p + F17
    rate = 0.0

    # Ha96n
    rate += np.exp(  40.8358 + -39.388*tf.T913i + -17.4673*tf.T913
                  + 35.3029*tf.T9 + -24.8162*tf.T953 + -0.666667*tf.lnT9)
    # Ha96r
    rate += np.exp(  16.3087 + -22.51*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  11.1184 + -13.6*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  -106.091 + -0.453036*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  12.1289 + -12.0223*tf.T9i
                  + -1.5*tf.lnT9)
    # Ha96r
    rate += np.exp(  18.6518 + -26.0*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O14__p_F17 = rate

@numba.njit()
def C12_O16__He4_Mg24(rate_eval, tf):
    # O16 + C12 --> He4 + Mg24
    rate = 0.0

    # cf88r
    rate += np.exp(  48.5341 + 0.37204*tf.T9i + -133.413*tf.T913i + 50.1572*tf.T913
                  + -3.15987*tf.T9 + 0.0178251*tf.T953 + -23.7027*tf.lnT9)

    rate_eval.C12_O16__He4_Mg24 = rate

@numba.njit()
def p_O17__He4_N14(rate_eval, tf):
    # O17 + p --> He4 + N14
    rate = 0.0

    # il10r
    rate += np.exp(  5.5336 + -2.11477*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -7.20763 + -0.753395*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  19.579 + -16.9078*tf.T913i
                  + -2.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  10.174 + -4.95865*tf.T9i + 5.10182*tf.T913
                  + 0.379373*tf.T9 + -0.0672515*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_O17__He4_N14 = rate

@numba.njit()
def p_O18__He4_N15(rate_eval, tf):
    # O18 + p --> He4 + N15
    rate = 0.0

    # il10r
    rate += np.exp(  10.2725 + -1.663*tf.T9i
                  + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -27.9044 + -0.245884*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  26.9671 + -16.6979*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  8.94352 + -5.32335*tf.T9i + 11.6568*tf.T913
                  + -2.16303*tf.T9 + 0.209965*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_O18__He4_N15 = rate

@numba.njit()
def p_F18__He4_O15(rate_eval, tf):
    # F18 + p --> He4 + O15
    rate = 0.0

    # il10n
    rate += np.exp(  62.0058 + -21.4023*tf.T913i + -80.8891*tf.T913
                  + 134.6*tf.T9 + -126.504*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  1.75704 + -3.01675*tf.T9i + 13.3223*tf.T913
                  + -1.36696*tf.T9 + 0.0757363*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  -31.7388 + -0.376432*tf.T9i + 61.738*tf.T913
                  + -108.29*tf.T9 + -34.2365*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_F18__He4_O15 = rate

@numba.njit()
def p_F19__He4_O16(rate_eval, tf):
    # F19 + p --> He4 + O16
    rate = 0.0

    # nacr 
    rate += np.exp(  8.239 + -2.46828*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  -52.7043 + -0.12765*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  26.2916 + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  14.3586 + -3.286*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  15.1955 + -3.75185*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_F19__He4_O16 = rate

@numba.njit()
def p_Ne20__He4_F17(rate_eval, tf):
    # Ne20 + p --> He4 + F17
    rate = 0.0

    # nacr 
    rate += np.exp(  41.563 + -47.9266*tf.T9i + -43.18*tf.T913i + 4.46827*tf.T913
                  + -1.63915*tf.T9 + 0.123483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ne20__He4_F17 = rate

@numba.njit()
def He4_He4_He4__C12(rate_eval, tf):
    # He4 + He4 + He4 --> C12
    rate = 0.0

    # fy05r
    rate += np.exp(  -24.3505 + -4.12656*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -13.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  -11.7884 + -1.02446*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + -2.16667*tf.lnT9)
    # fy05n
    rate += np.exp(  -0.971052 + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + -1.33333*tf.lnT9)

    rate_eval.He4_He4_He4__C12 = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    N13__C13__weak__wc12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    Ne19__F19__weak__wc12(rate_eval, tf)
    p_C12__N13(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    He4_N14__F18(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    He4_O14__Ne18(rate_eval, tf)
    He4_O15__Ne19(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    p_F18__Ne19(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_Ne18__Mg22(rate_eval, tf)
    He4_Ne20__Mg24(rate_eval, tf)
    C12_C12__He4_Ne20(rate_eval, tf)
    He4_N13__p_O16(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    C12_O16__He4_Mg24(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N14__F18 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__Ne18 *= scor
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O15__Ne19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor
        rate_eval.p_O17__He4_N14 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__Ne19 *= scor
        rate_eval.p_F18__He4_O15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne18__Mg22 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__Mg24 *= scor

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12__He4_Ne20 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N13__p_O16 *= scor

        scn_fac = ScreenFactors(6, 12, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_O16__He4_Mg24 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jp] = (
       -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    dYdt[jhe4] = (
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    dYdt[jc12] = (
       -rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       -2*5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20
       -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12
       )

    dYdt[jc13] = (
       -rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       +Y[jn13]*rate_eval.N13__C13__weak__wc12
       )

    dYdt[jn13] = (
       -Y[jn13]*rate_eval.N13__C13__weak__wc12
       -rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jp]*Y[jc12]*rate_eval.p_C12__N13
       )

    dYdt[jn14] = (
       -rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18
       +Y[jo14]*rate_eval.O14__N14__weak__wc12
       +rho*Y[jp]*Y[jc13]*rate_eval.p_C13__N14
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    dYdt[jn15] = (
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +Y[jo15]*rate_eval.O15__N15__weak__wc12
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    dYdt[jo14] = (
       -Y[jo14]*rate_eval.O14__N14__weak__wc12
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jp]*Y[jn13]*rate_eval.p_N13__O14
       )

    dYdt[jo15] = (
       -Y[jo15]*rate_eval.O15__N15__weak__wc12
       -rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19
       +rho*Y[jp]*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    dYdt[jo16] = (
       -rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jhe4]*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    dYdt[jo17] = (
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jp]*Y[jo17]*rate_eval.p_O17__He4_N14
       +Y[jf17]*rate_eval.F17__O17__weak__wc12
       )

    dYdt[jo18] = (
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       +Y[jf18]*rate_eval.F18__O18__weak__wc12
       )

    dYdt[jf17] = (
       -Y[jf17]*rate_eval.F17__O17__weak__wc12
       -rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       +rho*Y[jp]*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__p_F17
       +rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    dYdt[jf18] = (
       -Y[jf18]*rate_eval.F18__O18__weak__wc12
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jp]*Y[jf18]*rate_eval.p_F18__He4_O15
       +Y[jne18]*rate_eval.Ne18__F18__weak__wc12
       +rho*Y[jhe4]*Y[jn14]*rate_eval.He4_N14__F18
       +rho*Y[jp]*Y[jo17]*rate_eval.p_O17__F18
       )

    dYdt[jf19] = (
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +Y[jne19]*rate_eval.Ne19__F19__weak__wc12
       +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       )

    dYdt[jne18] = (
       -Y[jne18]*rate_eval.Ne18__F18__weak__wc12
       -rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22
       +rho*Y[jhe4]*Y[jo14]*rate_eval.He4_O14__Ne18
       +rho*Y[jp]*Y[jf17]*rate_eval.p_F17__Ne18
       )

    dYdt[jne19] = (
       -Y[jne19]*rate_eval.Ne19__F19__weak__wc12
       +rho*Y[jhe4]*Y[jo15]*rate_eval.He4_O15__Ne19
       +rho*Y[jp]*Y[jf18]*rate_eval.p_F18__Ne19
       )

    dYdt[jne20] = (
       -rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jp]*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__Ne20
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__Ne20
       +5.00000000000000e-01*rho*Y[jc12]**2*rate_eval.C12_C12__He4_Ne20
       )

    dYdt[jmg22] = (
       +rho*Y[jhe4]*Y[jne18]*rate_eval.He4_Ne18__Mg22
       )

    dYdt[jmg24] = (
       +rho*Y[jhe4]*Y[jne20]*rate_eval.He4_Ne20__Mg24
       +rho*Y[jc12]*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       )

    dYdt[jfe56] = 0.0

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    N13__C13__weak__wc12(rate_eval, tf)
    O14__N14__weak__wc12(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    F17__O17__weak__wc12(rate_eval, tf)
    F18__O18__weak__wc12(rate_eval, tf)
    Ne18__F18__weak__wc12(rate_eval, tf)
    Ne19__F19__weak__wc12(rate_eval, tf)
    p_C12__N13(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    p_C13__N14(rate_eval, tf)
    p_N13__O14(rate_eval, tf)
    p_N14__O15(rate_eval, tf)
    He4_N14__F18(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    He4_O14__Ne18(rate_eval, tf)
    He4_O15__Ne19(rate_eval, tf)
    p_O16__F17(rate_eval, tf)
    He4_O16__Ne20(rate_eval, tf)
    p_O17__F18(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    p_F17__Ne18(rate_eval, tf)
    p_F18__Ne19(rate_eval, tf)
    p_F19__Ne20(rate_eval, tf)
    He4_Ne18__Mg22(rate_eval, tf)
    He4_Ne20__Mg24(rate_eval, tf)
    C12_C12__He4_Ne20(rate_eval, tf)
    He4_N13__p_O16(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_O14__p_F17(rate_eval, tf)
    C12_O16__He4_Mg24(rate_eval, tf)
    p_O17__He4_N14(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_F18__He4_O15(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    p_Ne20__He4_F17(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C12__N13 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor

        scn_fac = ScreenFactors(1, 1, 6, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_C13__N14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N13__O14 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N14__O15 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N14__F18 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 14)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O14__Ne18 *= scor
        rate_eval.He4_O14__p_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O15__Ne19 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O16__F17 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__Ne20 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O17__F18 *= scor
        rate_eval.p_O17__He4_N14 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F17__Ne18 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F18__Ne19 *= scor
        rate_eval.p_F18__He4_O15 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__Ne20 *= scor
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne18__Mg22 *= scor

        scn_fac = ScreenFactors(2, 4, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ne20__Mg24 *= scor

        scn_fac = ScreenFactors(6, 12, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_C12__He4_Ne20 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 13)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N13__p_O16 *= scor

        scn_fac = ScreenFactors(6, 12, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.C12_O16__He4_Mg24 *= scor

        scn_fac = ScreenFactors(1, 1, 10, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ne20__He4_F17 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jp, jp] = (
       -rho*Y[jc12]*rate_eval.p_C12__N13
       -rho*Y[jc13]*rate_eval.p_C13__N14
       -rho*Y[jn13]*rate_eval.p_N13__O14
       -rho*Y[jn14]*rate_eval.p_N14__O15
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jo16]*rate_eval.p_O16__F17
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       -rho*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    jac[jp, jhe4] = (
       +rho*Y[jn13]*rate_eval.He4_N13__p_O16
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jp, jc12] = (
       -rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jp, jc13] = (
       -rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jp, jn13] = (
       -rho*Y[jp]*rate_eval.p_N13__O14
       +rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jp, jn14] = (
       -rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jp, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jp, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jp, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jp, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__F18
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jp, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jp, jf17] = (
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jp, jf18] = (
       -rho*Y[jp]*rate_eval.p_F18__Ne19
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jp, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jp, jne20] = (
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jhe4, jp] = (
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jn14]*rate_eval.He4_N14__F18
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jo15]*rate_eval.He4_O15__Ne19
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       -rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       -rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jn13]*rate_eval.He4_N13__p_O16
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       )

    jac[jhe4, jc12] = (
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       +rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jhe4, jn13] = (
       -rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jhe4, jn14] = (
       -rho*Y[jhe4]*rate_eval.He4_N14__F18
       )

    jac[jhe4, jn15] = (
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jhe4, jo14] = (
       -rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jhe4, jo15] = (
       -rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       )

    jac[jhe4, jo16] = (
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       +rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jhe4, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jhe4, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jhe4, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jhe4, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jhe4, jne18] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       )

    jac[jhe4, jne20] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jc12, jp] = (
       -rho*Y[jc12]*rate_eval.p_C12__N13
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       )

    jac[jc12, jc12] = (
       -rho*Y[jp]*rate_eval.p_C12__N13
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       -2*5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       -rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jc12, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jo16] = (
       -rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jc13, jp] = (
       -rho*Y[jc13]*rate_eval.p_C13__N14
       )

    jac[jc13, jc13] = (
       -rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jc13, jn13] = (
       +rate_eval.N13__C13__weak__wc12
       )

    jac[jn13, jp] = (
       -rho*Y[jn13]*rate_eval.p_N13__O14
       +rho*Y[jc12]*rate_eval.p_C12__N13
       )

    jac[jn13, jhe4] = (
       -rho*Y[jn13]*rate_eval.He4_N13__p_O16
       )

    jac[jn13, jc12] = (
       +rho*Y[jp]*rate_eval.p_C12__N13
       )

    jac[jn13, jn13] = (
       -rate_eval.N13__C13__weak__wc12
       -rho*Y[jp]*rate_eval.p_N13__O14
       -rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jn14, jp] = (
       -rho*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jc13]*rate_eval.p_C13__N14
       +rho*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    jac[jn14, jhe4] = (
       -rho*Y[jn14]*rate_eval.He4_N14__F18
       )

    jac[jn14, jc13] = (
       +rho*Y[jp]*rate_eval.p_C13__N14
       )

    jac[jn14, jn14] = (
       -rho*Y[jp]*rate_eval.p_N14__O15
       -rho*Y[jhe4]*rate_eval.He4_N14__F18
       )

    jac[jn14, jo14] = (
       +rate_eval.O14__N14__weak__wc12
       )

    jac[jn14, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jn15, jp] = (
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jhe4] = (
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       )

    jac[jn15, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jn15, jo15] = (
       +rate_eval.O15__N15__weak__wc12
       )

    jac[jn15, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jo14, jp] = (
       +rho*Y[jn13]*rate_eval.p_N13__O14
       )

    jac[jo14, jhe4] = (
       -rho*Y[jo14]*rate_eval.He4_O14__Ne18
       -rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jo14, jn13] = (
       +rho*Y[jp]*rate_eval.p_N13__O14
       )

    jac[jo14, jo14] = (
       -rate_eval.O14__N14__weak__wc12
       -rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       -rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jo15, jp] = (
       +rho*Y[jn14]*rate_eval.p_N14__O15
       +rho*Y[jf18]*rate_eval.p_F18__He4_O15
       )

    jac[jo15, jhe4] = (
       -rho*Y[jo15]*rate_eval.He4_O15__Ne19
       )

    jac[jo15, jn14] = (
       +rho*Y[jp]*rate_eval.p_N14__O15
       )

    jac[jo15, jo15] = (
       -rate_eval.O15__N15__weak__wc12
       -rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       )

    jac[jo15, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jo16, jp] = (
       -rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jhe4] = (
       -rho*Y[jo16]*rate_eval.He4_O16__Ne20
       +rho*Y[jc12]*rate_eval.He4_C12__O16
       +rho*Y[jn13]*rate_eval.He4_N13__p_O16
       )

    jac[jo16, jc12] = (
       -rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       +rho*Y[jhe4]*rate_eval.He4_C12__O16
       )

    jac[jo16, jn13] = (
       +rho*Y[jhe4]*rate_eval.He4_N13__p_O16
       )

    jac[jo16, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jo16, jo16] = (
       -rho*Y[jp]*rate_eval.p_O16__F17
       -rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       -rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jo16, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jo17, jp] = (
       -rho*Y[jo17]*rate_eval.p_O17__F18
       -rho*Y[jo17]*rate_eval.p_O17__He4_N14
       )

    jac[jo17, jo17] = (
       -rho*Y[jp]*rate_eval.p_O17__F18
       -rho*Y[jp]*rate_eval.p_O17__He4_N14
       )

    jac[jo17, jf17] = (
       +rate_eval.F17__O17__weak__wc12
       )

    jac[jo18, jp] = (
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jo18, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jo18, jf18] = (
       +rate_eval.F18__O18__weak__wc12
       )

    jac[jf17, jp] = (
       -rho*Y[jf17]*rate_eval.p_F17__Ne18
       +rho*Y[jo16]*rate_eval.p_O16__F17
       +rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf17, jhe4] = (
       +rho*Y[jo14]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__p_F17
       )

    jac[jf17, jo16] = (
       +rho*Y[jp]*rate_eval.p_O16__F17
       )

    jac[jf17, jf17] = (
       -rate_eval.F17__O17__weak__wc12
       -rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jf17, jne20] = (
       +rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jf18, jp] = (
       -rho*Y[jf18]*rate_eval.p_F18__Ne19
       -rho*Y[jf18]*rate_eval.p_F18__He4_O15
       +rho*Y[jo17]*rate_eval.p_O17__F18
       )

    jac[jf18, jhe4] = (
       +rho*Y[jn14]*rate_eval.He4_N14__F18
       )

    jac[jf18, jn14] = (
       +rho*Y[jhe4]*rate_eval.He4_N14__F18
       )

    jac[jf18, jo17] = (
       +rho*Y[jp]*rate_eval.p_O17__F18
       )

    jac[jf18, jf18] = (
       -rate_eval.F18__O18__weak__wc12
       -rho*Y[jp]*rate_eval.p_F18__Ne19
       -rho*Y[jp]*rate_eval.p_F18__He4_O15
       )

    jac[jf18, jne18] = (
       +rate_eval.Ne18__F18__weak__wc12
       )

    jac[jf19, jp] = (
       -rho*Y[jf19]*rate_eval.p_F19__Ne20
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jf19, jhe4] = (
       +rho*Y[jn15]*rate_eval.He4_N15__F19
       )

    jac[jf19, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__F19
       )

    jac[jf19, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jf19, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__Ne20
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jf19, jne19] = (
       +rate_eval.Ne19__F19__weak__wc12
       )

    jac[jne18, jp] = (
       +rho*Y[jf17]*rate_eval.p_F17__Ne18
       )

    jac[jne18, jhe4] = (
       -rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       +rho*Y[jo14]*rate_eval.He4_O14__Ne18
       )

    jac[jne18, jo14] = (
       +rho*Y[jhe4]*rate_eval.He4_O14__Ne18
       )

    jac[jne18, jf17] = (
       +rho*Y[jp]*rate_eval.p_F17__Ne18
       )

    jac[jne18, jne18] = (
       -rate_eval.Ne18__F18__weak__wc12
       -rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       )

    jac[jne19, jp] = (
       +rho*Y[jf18]*rate_eval.p_F18__Ne19
       )

    jac[jne19, jhe4] = (
       +rho*Y[jo15]*rate_eval.He4_O15__Ne19
       )

    jac[jne19, jo15] = (
       +rho*Y[jhe4]*rate_eval.He4_O15__Ne19
       )

    jac[jne19, jf18] = (
       +rho*Y[jp]*rate_eval.p_F18__Ne19
       )

    jac[jne19, jne19] = (
       -rate_eval.Ne19__F19__weak__wc12
       )

    jac[jne20, jp] = (
       -rho*Y[jne20]*rate_eval.p_Ne20__He4_F17
       +rho*Y[jf19]*rate_eval.p_F19__Ne20
       )

    jac[jne20, jhe4] = (
       -rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       +rho*Y[jo16]*rate_eval.He4_O16__Ne20
       )

    jac[jne20, jc12] = (
       +5.00000000000000e-01*rho*2*Y[jc12]*rate_eval.C12_C12__He4_Ne20
       )

    jac[jne20, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__Ne20
       )

    jac[jne20, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__Ne20
       )

    jac[jne20, jne20] = (
       -rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       -rho*Y[jp]*rate_eval.p_Ne20__He4_F17
       )

    jac[jmg22, jhe4] = (
       +rho*Y[jne18]*rate_eval.He4_Ne18__Mg22
       )

    jac[jmg22, jne18] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne18__Mg22
       )

    jac[jmg24, jhe4] = (
       +rho*Y[jne20]*rate_eval.He4_Ne20__Mg24
       )

    jac[jmg24, jc12] = (
       +rho*Y[jo16]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jmg24, jo16] = (
       +rho*Y[jc12]*rate_eval.C12_O16__He4_Mg24
       )

    jac[jmg24, jne20] = (
       +rho*Y[jhe4]*rate_eval.He4_Ne20__Mg24
       )

    jac[jfe56, jp] = 0.0

    jac[jfe56, jhe4] = 0.0

    jac[jfe56, jc12] = 0.0

    jac[jfe56, jc13] = 0.0

    jac[jfe56, jn13] = 0.0

    jac[jfe56, jn14] = 0.0

    jac[jfe56, jn15] = 0.0

    jac[jfe56, jo14] = 0.0

    jac[jfe56, jo15] = 0.0

    jac[jfe56, jo16] = 0.0

    jac[jfe56, jo17] = 0.0

    jac[jfe56, jo18] = 0.0

    jac[jfe56, jf17] = 0.0

    jac[jfe56, jf18] = 0.0

    jac[jfe56, jf19] = 0.0

    jac[jfe56, jne18] = 0.0

    jac[jfe56, jne19] = 0.0

    jac[jfe56, jne20] = 0.0

    jac[jfe56, jmg22] = 0.0

    jac[jfe56, jmg24] = 0.0

    jac[jfe56, jfe56] = 0.0

    return jac
