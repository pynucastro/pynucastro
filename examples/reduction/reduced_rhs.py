import numba
import numpy as np
from scipy import constants
from numba.experimental import jitclass

from pynucastro.rates import TableIndex, TableInterpolator, TabularRate, Tfactors
from pynucastro.screening import PlasmaState, ScreenFactors

jp = 0
jd = 1
jhe3 = 2
jhe4 = 3
jli7 = 4
jbe7 = 5
jb8 = 6
jc12 = 7
jn15 = 8
jo15 = 9
jo16 = 10
jo17 = 11
jo18 = 12
jf16 = 13
jf19 = 14
jf20 = 15
js28 = 16
js29 = 17
jcl29 = 18
jcl30 = 19
jar32 = 20
jar33 = 21
jk33 = 22
jk34 = 23
jca35 = 24
jca36 = 25
jca37 = 26
jca38 = 27
jsc36 = 28
jsc37 = 29
jsc38 = 30
jsc39 = 31
jti39 = 32
jti40 = 33
jv40 = 34
jv41 = 35
jcr43 = 36
jcr44 = 37
jmn44 = 38
jmn45 = 39
jfe45 = 40
jfe46 = 41
jfe47 = 42
jfe48 = 43
jco47 = 44
jco48 = 45
jco49 = 46
jni48 = 47
nnuc = 48

A = np.zeros((nnuc), dtype=np.int32)

A[jp] = 1
A[jd] = 2
A[jhe3] = 3
A[jhe4] = 4
A[jli7] = 7
A[jbe7] = 7
A[jb8] = 8
A[jc12] = 12
A[jn15] = 15
A[jo15] = 15
A[jo16] = 16
A[jo17] = 17
A[jo18] = 18
A[jf16] = 16
A[jf19] = 19
A[jf20] = 20
A[js28] = 28
A[js29] = 29
A[jcl29] = 29
A[jcl30] = 30
A[jar32] = 32
A[jar33] = 33
A[jk33] = 33
A[jk34] = 34
A[jca35] = 35
A[jca36] = 36
A[jca37] = 37
A[jca38] = 38
A[jsc36] = 36
A[jsc37] = 37
A[jsc38] = 38
A[jsc39] = 39
A[jti39] = 39
A[jti40] = 40
A[jv40] = 40
A[jv41] = 41
A[jcr43] = 43
A[jcr44] = 44
A[jmn44] = 44
A[jmn45] = 45
A[jfe45] = 45
A[jfe46] = 46
A[jfe47] = 47
A[jfe48] = 48
A[jco47] = 47
A[jco48] = 48
A[jco49] = 49
A[jni48] = 48

Z = np.zeros((nnuc), dtype=np.int32)

Z[jp] = 1
Z[jd] = 1
Z[jhe3] = 2
Z[jhe4] = 2
Z[jli7] = 3
Z[jbe7] = 4
Z[jb8] = 5
Z[jc12] = 6
Z[jn15] = 7
Z[jo15] = 8
Z[jo16] = 8
Z[jo17] = 8
Z[jo18] = 8
Z[jf16] = 9
Z[jf19] = 9
Z[jf20] = 9
Z[js28] = 16
Z[js29] = 16
Z[jcl29] = 17
Z[jcl30] = 17
Z[jar32] = 18
Z[jar33] = 18
Z[jk33] = 19
Z[jk34] = 19
Z[jca35] = 20
Z[jca36] = 20
Z[jca37] = 20
Z[jca38] = 20
Z[jsc36] = 21
Z[jsc37] = 21
Z[jsc38] = 21
Z[jsc39] = 21
Z[jti39] = 22
Z[jti40] = 22
Z[jv40] = 23
Z[jv41] = 23
Z[jcr43] = 24
Z[jcr44] = 24
Z[jmn44] = 25
Z[jmn45] = 25
Z[jfe45] = 26
Z[jfe46] = 26
Z[jfe47] = 26
Z[jfe48] = 26
Z[jco47] = 27
Z[jco48] = 27
Z[jco49] = 27
Z[jni48] = 28

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jp] = 0.0015040963030260536
mass[jd] = 0.003005881882823487
mass[jhe3] = 0.004501176676575961
mass[jhe4] = 0.0059735574925878256
mass[jli7] = 0.010470810406543588
mass[jbe7] = 0.010472191322584432
mass[jb8] = 0.011976069136782909
mass[jc12] = 0.017909017027273523
mass[jn15] = 0.0223864338056853
mass[jo15] = 0.02239084645968795
mass[jo16] = 0.023871099858982767
mass[jo17] = 0.02536981167252093
mass[jo18] = 0.02686227133140636
mass[jf16] = 0.023895792605265975
mass[jf19] = 0.028353560468882166
mass[jf20] = 0.02984833373331198
mass[js28] = 0.04179422725587193
mass[js29] = 0.043275167348072074
mass[jcl29] = 0.04330258699898635
mass[jcl30] = 0.04478003274394775
mass[jar32] = 0.0477538533099306
mass[jar33] = 0.04923476151881573
mass[jk33] = 0.049261877236822536
mass[jk34] = 0.05074026025511483
mass[jca35] = 0.0522429482929449
mass[jca36] = 0.05371671704253126
mass[jca37] = 0.05519842281494481
mass[jca38] = 0.05667654563975173
mass[jsc36] = 0.05375292623445966
mass[jsc37] = 0.055225525395103205
mass[jsc38] = 0.056705078002338316
mass[jsc39] = 0.058181597689205264
mass[jti39] = 0.05820831078022393
mass[jti40] = 0.059682319856305406
mass[jv40] = 0.05971670256687105
mass[jv41] = 0.06118963818460774
mass[jcr43] = 0.06417082139309448
mass[jcr44] = 0.0656448945562413
mass[jmn44] = 0.0656783480043592
mass[jmn45] = 0.06715083501263837
mass[jfe45] = 0.06718190121757164
mass[jfe46] = 0.06865317057160897
mass[jfe47] = 0.07013222650408754
mass[jfe48] = 0.07160721290791573
mass[jco47] = 0.07016066513934104
mass[jco48] = 0.0716388398746709
mass[jco49] = 0.07311281690721969
mass[jni48] = 0.07166519568030019

names = []
names.append("H1")
names.append("H2")
names.append("He3")
names.append("He4")
names.append("Li7")
names.append("Be7")
names.append("B8")
names.append("C12")
names.append("N15")
names.append("O15")
names.append("O16")
names.append("O17")
names.append("O18")
names.append("F16")
names.append("F19")
names.append("F20")
names.append("S28")
names.append("S29")
names.append("Cl29")
names.append("Cl30")
names.append("Ar32")
names.append("Ar33")
names.append("K33")
names.append("K34")
names.append("Ca35")
names.append("Ca36")
names.append("Ca37")
names.append("Ca38")
names.append("Sc36")
names.append("Sc37")
names.append("Sc38")
names.append("Sc39")
names.append("Ti39")
names.append("Ti40")
names.append("V40")
names.append("V41")
names.append("Cr43")
names.append("Cr44")
names.append("Mn44")
names.append("Mn45")
names.append("Fe45")
names.append("Fe46")
names.append("Fe47")
names.append("Fe48")
names.append("Co47")
names.append("Co48")
names.append("Co49")
names.append("Ni48")

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
    ("Be7__Li7__weak__electron_capture", numba.float64),
    ("O15__N15__weak__wc12", numba.float64),
    ("He3__p_d", numba.float64),
    ("He4__d_d", numba.float64),
    ("Be7__He4_He3", numba.float64),
    ("B8__p_Be7", numba.float64),
    ("B8__He4_He4__weak__wc12", numba.float64),
    ("O16__p_N15", numba.float64),
    ("O16__He4_C12", numba.float64),
    ("F19__p_O18", numba.float64),
    ("F19__He4_N15", numba.float64),
    ("C12__He4_He4_He4", numba.float64),
    ("p_p__d__weak__bet_pos_", numba.float64),
    ("p_p__d__weak__electron_capture", numba.float64),
    ("p_d__He3", numba.float64),
    ("d_d__He4", numba.float64),
    ("p_He3__He4__weak__bet_pos_", numba.float64),
    ("He4_He3__Be7", numba.float64),
    ("p_Be7__B8", numba.float64),
    ("He4_C12__O16", numba.float64),
    ("p_N15__O16", numba.float64),
    ("He4_N15__F19", numba.float64),
    ("p_O18__F19", numba.float64),
    ("d_He3__p_He4", numba.float64),
    ("p_He4__d_He3", numba.float64),
    ("He4_He4__p_Li7", numba.float64),
    ("p_Li7__He4_He4", numba.float64),
    ("He4_C12__p_N15", numba.float64),
    ("p_N15__He4_C12", numba.float64),
    ("He4_N15__p_O18", numba.float64),
    ("He4_O16__p_F19", numba.float64),
    ("He4_O17__p_F20", numba.float64),
    ("p_O18__He4_N15", numba.float64),
    ("p_F19__He4_O16", numba.float64),
    ("p_F20__He4_O17", numba.float64),
    ("He3_He3__p_p_He4", numba.float64),
    ("d_Be7__p_He4_He4", numba.float64),
    ("He3_Be7__p_p_He4_He4", numba.float64),
    ("He4_He4_He4__C12", numba.float64),
    ("p_p_He4__He3_He3", numba.float64),
    ("p_He4_He4__d_Be7", numba.float64),
    ("p_p_He4_He4__He3_Be7", numba.float64),
    ("p_O15__F16", numba.float64),
    ("F16__p_O15", numba.float64),
    ("p_S28__Cl29", numba.float64),
    ("He4_S28__Ar32", numba.float64),
    ("p_S29__Cl30", numba.float64),
    ("He4_S29__Ar33", numba.float64),
    ("Cl29__S29__weak__bqa_pos_", numba.float64),
    ("Cl29__p_S28", numba.float64),
    ("He4_Cl29__K33", numba.float64),
    ("He4_Cl29__p_Ar32", numba.float64),
    ("Ar32__He4_S28", numba.float64),
    ("p_Ar32__K33", numba.float64),
    ("He4_Ar32__Ca36", numba.float64),
    ("p_Ar32__He4_Cl29", numba.float64),
    ("p_Ca38__Sc39", numba.float64),
    ("p_p_Ca38__He4_Ca36", numba.float64),
    ("Cl30__p_S29", numba.float64),
    ("He4_Cl30__K34", numba.float64),
    ("He4_Cl30__p_Ar33", numba.float64),
    ("Ar33__He4_S29", numba.float64),
    ("p_Ar33__K34", numba.float64),
    ("He4_Ar33__Ca37", numba.float64),
    ("p_Ar33__He4_Cl30", numba.float64),
    ("K33__Ar33__weak__bqa_pos_", numba.float64),
    ("K33__p_Ar32", numba.float64),
    ("K33__He4_Cl29", numba.float64),
    ("He4_K33__Sc37", numba.float64),
    ("He4_K33__p_Ca36", numba.float64),
    ("Ca36__He4_Ar32", numba.float64),
    ("p_Ca36__Sc37", numba.float64),
    ("He4_Ca36__Ti40", numba.float64),
    ("p_Ca36__He4_K33", numba.float64),
    ("He4_Ca36__p_Sc39", numba.float64),
    ("He4_Ca36__p_p_Ca38", numba.float64),
    ("Sc39__p_Ca38", numba.float64),
    ("p_Sc39__Ti40", numba.float64),
    ("p_Sc39__He4_Ca36", numba.float64),
    ("K34__p_Ar33", numba.float64),
    ("K34__He4_Cl30", numba.float64),
    ("p_K34__Ca35", numba.float64),
    ("He4_K34__Sc38", numba.float64),
    ("He4_K34__p_Ca37", numba.float64),
    ("Ca37__He4_Ar33", numba.float64),
    ("p_Ca37__Sc38", numba.float64),
    ("p_Ca37__He4_K34", numba.float64),
    ("Sc37__Ca37__weak__bqa_pos_", numba.float64),
    ("Sc37__p_Ca36", numba.float64),
    ("Sc37__He4_K33", numba.float64),
    ("He4_Sc37__V41", numba.float64),
    ("He4_Sc37__p_Ti40", numba.float64),
    ("Ti40__p_Sc39", numba.float64),
    ("Ti40__He4_Ca36", numba.float64),
    ("p_Ti40__V41", numba.float64),
    ("He4_Ti40__Cr44", numba.float64),
    ("p_Ti40__He4_Sc37", numba.float64),
    ("Ca35__p_K34", numba.float64),
    ("p_Ca35__Sc36", numba.float64),
    ("He4_Ca35__Ti39", numba.float64),
    ("He4_Ca35__p_Sc38", numba.float64),
    ("Sc38__Ca38__weak__mo97", numba.float64),
    ("Sc38__p_Ca37", numba.float64),
    ("Sc38__He4_K34", numba.float64),
    ("p_Sc38__Ti39", numba.float64),
    ("p_Sc38__He4_Ca35", numba.float64),
    ("V41__p_Ti40", numba.float64),
    ("V41__He4_Sc37", numba.float64),
    ("He4_V41__Mn45", numba.float64),
    ("He4_V41__p_Cr44", numba.float64),
    ("Cr44__He4_Ti40", numba.float64),
    ("p_Cr44__Mn45", numba.float64),
    ("He4_Cr44__Fe48", numba.float64),
    ("p_Cr44__He4_V41", numba.float64),
    ("Sc36__Ca36__weak__bqa_pos_", numba.float64),
    ("Sc36__p_Ca35", numba.float64),
    ("He4_Sc36__V40", numba.float64),
    ("He4_Sc36__p_Ti39", numba.float64),
    ("Ti39__Sc39__weak__wc17", numba.float64),
    ("Ti39__p_Sc38", numba.float64),
    ("Ti39__p_Ca38__weak__wc12", numba.float64),
    ("Ti39__He4_Ca35", numba.float64),
    ("p_Ti39__V40", numba.float64),
    ("He4_Ti39__Cr43", numba.float64),
    ("p_Ti39__He4_Sc36", numba.float64),
    ("Mn45__p_Cr44", numba.float64),
    ("Mn45__He4_V41", numba.float64),
    ("p_Mn45__Fe46", numba.float64),
    ("He4_Mn45__Co49", numba.float64),
    ("He4_Mn45__p_Fe48", numba.float64),
    ("Fe48__He4_Cr44", numba.float64),
    ("p_Fe48__Co49", numba.float64),
    ("p_Fe48__He4_Mn45", numba.float64),
    ("V40__Ti40__weak__bqa_pos_", numba.float64),
    ("V40__p_Ti39", numba.float64),
    ("V40__He4_Sc36", numba.float64),
    ("He4_V40__Mn44", numba.float64),
    ("He4_V40__p_Cr43", numba.float64),
    ("Cr43__He4_Ti39", numba.float64),
    ("p_Cr43__Mn44", numba.float64),
    ("He4_Cr43__Fe47", numba.float64),
    ("p_Cr43__He4_V40", numba.float64),
    ("Fe46__p_Mn45", numba.float64),
    ("p_Fe46__Co47", numba.float64),
    ("He4_Fe46__p_Co49", numba.float64),
    ("Co49__p_Fe48", numba.float64),
    ("Co49__He4_Mn45", numba.float64),
    ("p_Co49__He4_Fe46", numba.float64),
    ("Mn44__Cr44__weak__bqa_pos_", numba.float64),
    ("Mn44__p_Cr43", numba.float64),
    ("Mn44__He4_V40", numba.float64),
    ("p_Mn44__Fe45", numba.float64),
    ("He4_Mn44__Co48", numba.float64),
    ("He4_Mn44__p_Fe47", numba.float64),
    ("Fe47__He4_Cr43", numba.float64),
    ("p_Fe47__Co48", numba.float64),
    ("p_Fe47__He4_Mn44", numba.float64),
    ("Co47__Fe47__weak__bqa_pos_", numba.float64),
    ("Co47__p_Fe46", numba.float64),
    ("p_Co47__Ni48", numba.float64),
    ("Fe45__Mn45__weak__wc17", numba.float64),
    ("Fe45__p_Mn44", numba.float64),
    ("Fe45__p_Cr44__weak__wc17", numba.float64),
    ("He4_Fe45__p_Co48", numba.float64),
    ("Co48__Fe48__weak__bqa_pos_", numba.float64),
    ("Co48__p_Fe47", numba.float64),
    ("Co48__He4_Mn44", numba.float64),
    ("p_Co48__He4_Fe45", numba.float64),
    ("Ni48__Co48__weak__wc17", numba.float64),
    ("Ni48__p_Co47", numba.float64),
])
class RateEval:
    def __init__(self):
        self.Be7__Li7__weak__electron_capture = np.nan
        self.O15__N15__weak__wc12 = np.nan
        self.He3__p_d = np.nan
        self.He4__d_d = np.nan
        self.Be7__He4_He3 = np.nan
        self.B8__p_Be7 = np.nan
        self.B8__He4_He4__weak__wc12 = np.nan
        self.O16__p_N15 = np.nan
        self.O16__He4_C12 = np.nan
        self.F19__p_O18 = np.nan
        self.F19__He4_N15 = np.nan
        self.C12__He4_He4_He4 = np.nan
        self.p_p__d__weak__bet_pos_ = np.nan
        self.p_p__d__weak__electron_capture = np.nan
        self.p_d__He3 = np.nan
        self.d_d__He4 = np.nan
        self.p_He3__He4__weak__bet_pos_ = np.nan
        self.He4_He3__Be7 = np.nan
        self.p_Be7__B8 = np.nan
        self.He4_C12__O16 = np.nan
        self.p_N15__O16 = np.nan
        self.He4_N15__F19 = np.nan
        self.p_O18__F19 = np.nan
        self.d_He3__p_He4 = np.nan
        self.p_He4__d_He3 = np.nan
        self.He4_He4__p_Li7 = np.nan
        self.p_Li7__He4_He4 = np.nan
        self.He4_C12__p_N15 = np.nan
        self.p_N15__He4_C12 = np.nan
        self.He4_N15__p_O18 = np.nan
        self.He4_O16__p_F19 = np.nan
        self.He4_O17__p_F20 = np.nan
        self.p_O18__He4_N15 = np.nan
        self.p_F19__He4_O16 = np.nan
        self.p_F20__He4_O17 = np.nan
        self.He3_He3__p_p_He4 = np.nan
        self.d_Be7__p_He4_He4 = np.nan
        self.He3_Be7__p_p_He4_He4 = np.nan
        self.He4_He4_He4__C12 = np.nan
        self.p_p_He4__He3_He3 = np.nan
        self.p_He4_He4__d_Be7 = np.nan
        self.p_p_He4_He4__He3_Be7 = np.nan
        self.p_O15__F16 = np.nan
        self.F16__p_O15 = np.nan
        self.p_S28__Cl29 = np.nan
        self.He4_S28__Ar32 = np.nan
        self.p_S29__Cl30 = np.nan
        self.He4_S29__Ar33 = np.nan
        self.Cl29__S29__weak__bqa_pos_ = np.nan
        self.Cl29__p_S28 = np.nan
        self.He4_Cl29__K33 = np.nan
        self.He4_Cl29__p_Ar32 = np.nan
        self.Ar32__He4_S28 = np.nan
        self.p_Ar32__K33 = np.nan
        self.He4_Ar32__Ca36 = np.nan
        self.p_Ar32__He4_Cl29 = np.nan
        self.p_Ca38__Sc39 = np.nan
        self.p_p_Ca38__He4_Ca36 = np.nan
        self.Cl30__p_S29 = np.nan
        self.He4_Cl30__K34 = np.nan
        self.He4_Cl30__p_Ar33 = np.nan
        self.Ar33__He4_S29 = np.nan
        self.p_Ar33__K34 = np.nan
        self.He4_Ar33__Ca37 = np.nan
        self.p_Ar33__He4_Cl30 = np.nan
        self.K33__Ar33__weak__bqa_pos_ = np.nan
        self.K33__p_Ar32 = np.nan
        self.K33__He4_Cl29 = np.nan
        self.He4_K33__Sc37 = np.nan
        self.He4_K33__p_Ca36 = np.nan
        self.Ca36__He4_Ar32 = np.nan
        self.p_Ca36__Sc37 = np.nan
        self.He4_Ca36__Ti40 = np.nan
        self.p_Ca36__He4_K33 = np.nan
        self.He4_Ca36__p_Sc39 = np.nan
        self.He4_Ca36__p_p_Ca38 = np.nan
        self.Sc39__p_Ca38 = np.nan
        self.p_Sc39__Ti40 = np.nan
        self.p_Sc39__He4_Ca36 = np.nan
        self.K34__p_Ar33 = np.nan
        self.K34__He4_Cl30 = np.nan
        self.p_K34__Ca35 = np.nan
        self.He4_K34__Sc38 = np.nan
        self.He4_K34__p_Ca37 = np.nan
        self.Ca37__He4_Ar33 = np.nan
        self.p_Ca37__Sc38 = np.nan
        self.p_Ca37__He4_K34 = np.nan
        self.Sc37__Ca37__weak__bqa_pos_ = np.nan
        self.Sc37__p_Ca36 = np.nan
        self.Sc37__He4_K33 = np.nan
        self.He4_Sc37__V41 = np.nan
        self.He4_Sc37__p_Ti40 = np.nan
        self.Ti40__p_Sc39 = np.nan
        self.Ti40__He4_Ca36 = np.nan
        self.p_Ti40__V41 = np.nan
        self.He4_Ti40__Cr44 = np.nan
        self.p_Ti40__He4_Sc37 = np.nan
        self.Ca35__p_K34 = np.nan
        self.p_Ca35__Sc36 = np.nan
        self.He4_Ca35__Ti39 = np.nan
        self.He4_Ca35__p_Sc38 = np.nan
        self.Sc38__Ca38__weak__mo97 = np.nan
        self.Sc38__p_Ca37 = np.nan
        self.Sc38__He4_K34 = np.nan
        self.p_Sc38__Ti39 = np.nan
        self.p_Sc38__He4_Ca35 = np.nan
        self.V41__p_Ti40 = np.nan
        self.V41__He4_Sc37 = np.nan
        self.He4_V41__Mn45 = np.nan
        self.He4_V41__p_Cr44 = np.nan
        self.Cr44__He4_Ti40 = np.nan
        self.p_Cr44__Mn45 = np.nan
        self.He4_Cr44__Fe48 = np.nan
        self.p_Cr44__He4_V41 = np.nan
        self.Sc36__Ca36__weak__bqa_pos_ = np.nan
        self.Sc36__p_Ca35 = np.nan
        self.He4_Sc36__V40 = np.nan
        self.He4_Sc36__p_Ti39 = np.nan
        self.Ti39__Sc39__weak__wc17 = np.nan
        self.Ti39__p_Sc38 = np.nan
        self.Ti39__p_Ca38__weak__wc12 = np.nan
        self.Ti39__He4_Ca35 = np.nan
        self.p_Ti39__V40 = np.nan
        self.He4_Ti39__Cr43 = np.nan
        self.p_Ti39__He4_Sc36 = np.nan
        self.Mn45__p_Cr44 = np.nan
        self.Mn45__He4_V41 = np.nan
        self.p_Mn45__Fe46 = np.nan
        self.He4_Mn45__Co49 = np.nan
        self.He4_Mn45__p_Fe48 = np.nan
        self.Fe48__He4_Cr44 = np.nan
        self.p_Fe48__Co49 = np.nan
        self.p_Fe48__He4_Mn45 = np.nan
        self.V40__Ti40__weak__bqa_pos_ = np.nan
        self.V40__p_Ti39 = np.nan
        self.V40__He4_Sc36 = np.nan
        self.He4_V40__Mn44 = np.nan
        self.He4_V40__p_Cr43 = np.nan
        self.Cr43__He4_Ti39 = np.nan
        self.p_Cr43__Mn44 = np.nan
        self.He4_Cr43__Fe47 = np.nan
        self.p_Cr43__He4_V40 = np.nan
        self.Fe46__p_Mn45 = np.nan
        self.p_Fe46__Co47 = np.nan
        self.He4_Fe46__p_Co49 = np.nan
        self.Co49__p_Fe48 = np.nan
        self.Co49__He4_Mn45 = np.nan
        self.p_Co49__He4_Fe46 = np.nan
        self.Mn44__Cr44__weak__bqa_pos_ = np.nan
        self.Mn44__p_Cr43 = np.nan
        self.Mn44__He4_V40 = np.nan
        self.p_Mn44__Fe45 = np.nan
        self.He4_Mn44__Co48 = np.nan
        self.He4_Mn44__p_Fe47 = np.nan
        self.Fe47__He4_Cr43 = np.nan
        self.p_Fe47__Co48 = np.nan
        self.p_Fe47__He4_Mn44 = np.nan
        self.Co47__Fe47__weak__bqa_pos_ = np.nan
        self.Co47__p_Fe46 = np.nan
        self.p_Co47__Ni48 = np.nan
        self.Fe45__Mn45__weak__wc17 = np.nan
        self.Fe45__p_Mn44 = np.nan
        self.Fe45__p_Cr44__weak__wc17 = np.nan
        self.He4_Fe45__p_Co48 = np.nan
        self.Co48__Fe48__weak__bqa_pos_ = np.nan
        self.Co48__p_Fe47 = np.nan
        self.Co48__He4_Mn44 = np.nan
        self.p_Co48__He4_Fe45 = np.nan
        self.Ni48__Co48__weak__wc17 = np.nan
        self.Ni48__p_Co47 = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def Be7__Li7__weak__electron_capture(rate_eval, tf):
    # Be7 --> Li7
    rate = 0.0

    #   ecw
    rate += np.exp(  -23.8328 + 3.02033*tf.T913
                  + -0.0742132*tf.T9 + -0.00792386*tf.T953 + -0.650113*tf.lnT9)

    rate_eval.Be7__Li7__weak__electron_capture = rate

@numba.njit()
def O15__N15__weak__wc12(rate_eval, tf):
    # O15 --> N15
    rate = 0.0

    # wc12w
    rate += np.exp(  -5.17053)

    rate_eval.O15__N15__weak__wc12 = rate

@numba.njit()
def He3__p_d(rate_eval, tf):
    # He3 --> p + d
    rate = 0.0

    # de04 
    rate += np.exp(  32.4383 + -63.7435*tf.T9i + -3.7208*tf.T913i + 0.198654*tf.T913
                  + 1.83333*tf.lnT9)
    # de04n
    rate += np.exp(  31.032 + -63.7435*tf.T9i + -3.7208*tf.T913i + 0.871782*tf.T913
                  + 0.833333*tf.lnT9)

    rate_eval.He3__p_d = rate

@numba.njit()
def He4__d_d(rate_eval, tf):
    # He4 --> d + d
    rate = 0.0

    # nacrn
    rate += np.exp(  28.2984 + -276.744*tf.T9i + -4.26166*tf.T913i + -0.119233*tf.T913
                  + 0.778829*tf.T9 + -0.0925203*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.He4__d_d = rate

@numba.njit()
def Be7__He4_He3(rate_eval, tf):
    # Be7 --> He4 + He3
    rate = 0.0

    # cd08n
    rate += np.exp(  38.7379 + -18.4059*tf.T9i + -12.8271*tf.T913i + -0.0308225*tf.T913
                  + -0.654685*tf.T9 + 0.0896331*tf.T953 + 0.833333*tf.lnT9)
    # cd08n
    rate += np.exp(  40.8355 + -18.4059*tf.T9i + -12.8271*tf.T913i + -3.8126*tf.T913
                  + 0.0942285*tf.T9 + -0.00301018*tf.T953 + 2.83333*tf.lnT9)

    rate_eval.Be7__He4_He3 = rate

@numba.njit()
def B8__p_Be7(rate_eval, tf):
    # B8 --> p + Be7
    rate = 0.0

    # nacrn
    rate += np.exp(  35.8138 + -1.58982*tf.T9i + -10.264*tf.T913i + -0.203472*tf.T913
                  + 0.121083*tf.T9 + -0.00700063*tf.T953 + 0.833333*tf.lnT9)
    # nacrr
    rate += np.exp(  31.0163 + -8.93482*tf.T9i)

    rate_eval.B8__p_Be7 = rate

@numba.njit()
def B8__He4_He4__weak__wc12(rate_eval, tf):
    # B8 --> He4 + He4
    rate = 0.0

    # wc12w
    rate += np.exp(  -0.105148)

    rate_eval.B8__He4_He4__weak__wc12 = rate

@numba.njit()
def O16__p_N15(rate_eval, tf):
    # O16 --> p + N15
    rate = 0.0

    # li10r
    rate += np.exp(  38.8465 + -150.962*tf.T9i
                  + 0.0459037*tf.T9)
    # li10r
    rate += np.exp(  30.8927 + -143.656*tf.T9i)
    # li10n
    rate += np.exp(  44.3197 + -140.732*tf.T9i + -15.24*tf.T913i + 0.334926*tf.T913
                  + 4.59088*tf.T9 + -4.78468*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.O16__p_N15 = rate

@numba.njit()
def O16__He4_C12(rate_eval, tf):
    # O16 --> He4 + C12
    rate = 0.0

    # nac2 
    rate += np.exp(  279.295 + -84.9515*tf.T9i + 103.411*tf.T913i + -420.567*tf.T913
                  + 64.0874*tf.T9 + -12.4624*tf.T953 + 138.803*tf.lnT9)
    # nac2 
    rate += np.exp(  94.3131 + -84.503*tf.T9i + 58.9128*tf.T913i + -148.273*tf.T913
                  + 9.08324*tf.T9 + -0.541041*tf.T953 + 71.8554*tf.lnT9)

    rate_eval.O16__He4_C12 = rate

@numba.njit()
def F19__p_O18(rate_eval, tf):
    # F19 --> p + O18
    rate = 0.0

    # il10n
    rate += np.exp(  42.8485 + -92.7757*tf.T9i + -16.7246*tf.T913i
                  + -3.0*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  30.2003 + -99.501*tf.T9i + 3.99059*tf.T913
                  + -0.593127*tf.T9 + 0.0877534*tf.T953)
    # il10r
    rate += np.exp(  28.008 + -94.4325*tf.T9i)
    # il10r
    rate += np.exp(  -12.0764 + -93.0204*tf.T9i)

    rate_eval.F19__p_O18 = rate

@numba.njit()
def F19__He4_N15(rate_eval, tf):
    # F19 --> He4 + N15
    rate = 0.0

    # il10r
    rate += np.exp(  15.3186 + -50.7554*tf.T9i)
    # il10n
    rate += np.exp(  50.1291 + -46.5774*tf.T9i + -36.2324*tf.T913i
                  + -2.0*tf.T953 + 0.833333*tf.lnT9)
    # il10r
    rate += np.exp(  -4.06142 + -50.7773*tf.T9i + 35.4292*tf.T913
                  + -5.5767*tf.T9 + 0.441293*tf.T953)
    # il10r
    rate += np.exp(  28.2717 + -53.5621*tf.T9i)

    rate_eval.F19__He4_N15 = rate

@numba.njit()
def C12__He4_He4_He4(rate_eval, tf):
    # C12 --> He4 + He4 + He4
    rate = 0.0

    # fy05n
    rate += np.exp(  45.7734 + -84.4227*tf.T9i + -37.06*tf.T913i + 29.3493*tf.T913
                  + -115.507*tf.T9 + -10.0*tf.T953 + 1.66667*tf.lnT9)
    # fy05r
    rate += np.exp(  22.394 + -88.5493*tf.T9i + -13.49*tf.T913i + 21.4259*tf.T913
                  + -1.34769*tf.T9 + 0.0879816*tf.T953 + -10.1653*tf.lnT9)
    # fy05r
    rate += np.exp(  34.9561 + -85.4472*tf.T9i + -23.57*tf.T913i + 20.4886*tf.T913
                  + -12.9882*tf.T9 + -20.0*tf.T953 + 0.83333*tf.lnT9)

    rate_eval.C12__He4_He4_He4 = rate

@numba.njit()
def p_p__d__weak__bet_pos_(rate_eval, tf):
    # p + p --> d
    rate = 0.0

    # bet+w
    rate += np.exp(  -34.7863 + -3.51193*tf.T913i + 3.10086*tf.T913
                  + -0.198314*tf.T9 + 0.0126251*tf.T953 + -1.02517*tf.lnT9)

    rate_eval.p_p__d__weak__bet_pos_ = rate

@numba.njit()
def p_p__d__weak__electron_capture(rate_eval, tf):
    # p + p --> d
    rate = 0.0

    #   ecw
    rate += np.exp(  -43.6499 + -0.00246064*tf.T9i + -2.7507*tf.T913i + -0.424877*tf.T913
                  + 0.015987*tf.T9 + -0.000690875*tf.T953 + -0.207625*tf.lnT9)

    rate_eval.p_p__d__weak__electron_capture = rate

@numba.njit()
def p_d__He3(rate_eval, tf):
    # d + p --> He3
    rate = 0.0

    # de04 
    rate += np.exp(  8.93525 + -3.7208*tf.T913i + 0.198654*tf.T913
                  + 0.333333*tf.lnT9)
    # de04n
    rate += np.exp(  7.52898 + -3.7208*tf.T913i + 0.871782*tf.T913
                  + -0.666667*tf.lnT9)

    rate_eval.p_d__He3 = rate

@numba.njit()
def d_d__He4(rate_eval, tf):
    # d + d --> He4
    rate = 0.0

    # nacrn
    rate += np.exp(  3.78177 + -4.26166*tf.T913i + -0.119233*tf.T913
                  + 0.778829*tf.T9 + -0.0925203*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.d_d__He4 = rate

@numba.njit()
def p_He3__He4__weak__bet_pos_(rate_eval, tf):
    # He3 + p --> He4
    rate = 0.0

    # bet+w
    rate += np.exp(  -27.7611 + -4.30107e-12*tf.T9i + -6.141*tf.T913i + -1.93473e-09*tf.T913
                  + 2.04145e-10*tf.T9 + -1.80372e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_He3__He4__weak__bet_pos_ = rate

@numba.njit()
def He4_He3__Be7(rate_eval, tf):
    # He3 + He4 --> Be7
    rate = 0.0

    # cd08n
    rate += np.exp(  17.7075 + -12.8271*tf.T913i + -3.8126*tf.T913
                  + 0.0942285*tf.T9 + -0.00301018*tf.T953 + 1.33333*tf.lnT9)
    # cd08n
    rate += np.exp(  15.6099 + -12.8271*tf.T913i + -0.0308225*tf.T913
                  + -0.654685*tf.T9 + 0.0896331*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_He3__Be7 = rate

@numba.njit()
def p_Be7__B8(rate_eval, tf):
    # Be7 + p --> B8
    rate = 0.0

    # nacrn
    rate += np.exp(  12.5315 + -10.264*tf.T913i + -0.203472*tf.T913
                  + 0.121083*tf.T9 + -0.00700063*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  7.73399 + -7.345*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Be7__B8 = rate

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
def d_He3__p_He4(rate_eval, tf):
    # He3 + d --> p + He4
    rate = 0.0

    # de04 
    rate += np.exp(  24.6839 + -7.182*tf.T913i + 0.473288*tf.T913
                  + 1.46847*tf.T9 + -27.9603*tf.T953 + -0.666667*tf.lnT9)
    # de04 
    rate += np.exp(  41.2969 + -7.182*tf.T913i + -17.1349*tf.T913
                  + 1.36908*tf.T9 + -0.0814423*tf.T953 + 3.35395*tf.lnT9)

    rate_eval.d_He3__p_He4 = rate

@numba.njit()
def p_He4__d_He3(rate_eval, tf):
    # p + He4 --> d + He3
    rate = 0.0

    # de04 
    rate += np.exp(  43.0037 + -212.977*tf.T9i + -7.182*tf.T913i + -17.1349*tf.T913
                  + 1.36908*tf.T9 + -0.0814423*tf.T953 + 3.35395*tf.lnT9)
    # de04 
    rate += np.exp(  26.3907 + -212.977*tf.T9i + -7.182*tf.T913i + 0.473288*tf.T913
                  + 1.46847*tf.T9 + -27.9603*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_He4__d_He3 = rate

@numba.njit()
def He4_He4__p_Li7(rate_eval, tf):
    # He4 + He4 --> p + Li7
    rate = 0.0

    # de04r
    rate += np.exp(  23.4325 + -227.465*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  21.9764 + -201.312*tf.T9i + -8.4727*tf.T913i + 0.297934*tf.T913
                  + 0.0582335*tf.T9 + -0.00413383*tf.T953 + -0.666667*tf.lnT9)
    # de04r
    rate += np.exp(  15.7864 + -205.79*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  13.4902 + -201.312*tf.T9i + -8.4727*tf.T913i + 0.417943*tf.T913
                  + 5.34565*tf.T9 + -4.8684*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_He4__p_Li7 = rate

@numba.njit()
def p_Li7__He4_He4(rate_eval, tf):
    # Li7 + p --> He4 + He4
    rate = 0.0

    # de04 
    rate += np.exp(  11.9576 + -8.4727*tf.T913i + 0.417943*tf.T913
                  + 5.34565*tf.T9 + -4.8684*tf.T953 + -0.666667*tf.lnT9)
    # de04r
    rate += np.exp(  21.8999 + -26.1527*tf.T9i
                  + -1.5*tf.lnT9)
    # de04 
    rate += np.exp(  20.4438 + -8.4727*tf.T913i + 0.297934*tf.T913
                  + 0.0582335*tf.T9 + -0.00413383*tf.T953 + -0.666667*tf.lnT9)
    # de04r
    rate += np.exp(  14.2538 + -4.478*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.p_Li7__He4_He4 = rate

@numba.njit()
def He4_C12__p_N15(rate_eval, tf):
    # C12 + He4 --> p + N15
    rate = 0.0

    # nacrn
    rate += np.exp(  27.118 + -57.6279*tf.T9i + -15.253*tf.T913i + 1.59318*tf.T913
                  + 2.4479*tf.T9 + -2.19708*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  -6.93365 + -58.7917*tf.T9i + 22.7105*tf.T913
                  + -2.90707*tf.T9 + 0.205754*tf.T953 + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  20.5388 + -65.034*tf.T9i
                  + -1.5*tf.lnT9)
    # nacrr
    rate += np.exp(  -5.2319 + -59.6491*tf.T9i + 30.8497*tf.T913
                  + -8.50433*tf.T9 + -1.54426*tf.T953 + -1.5*tf.lnT9)

    rate_eval.He4_C12__p_N15 = rate

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
def He4_N15__p_O18(rate_eval, tf):
    # N15 + He4 --> p + O18
    rate = 0.0

    # il10r
    rate += np.exp(  -29.7104 + -46.4444*tf.T9i
                  + -1.5*tf.lnT9)
    # il10n
    rate += np.exp(  25.1611 + -46.1986*tf.T9i + -16.6979*tf.T913i
                  + -3.0*tf.T953 + -0.666667*tf.lnT9)
    # il10r
    rate += np.exp(  7.13756 + -51.5219*tf.T9i + 11.6568*tf.T913
                  + -2.16303*tf.T9 + 0.209965*tf.T953 + -1.5*tf.lnT9)
    # il10r
    rate += np.exp(  8.46654 + -47.8616*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_N15__p_O18 = rate

@numba.njit()
def He4_O16__p_F19(rate_eval, tf):
    # O16 + He4 --> p + F19
    rate = 0.0

    # nacr 
    rate += np.exp(  -53.1397 + -94.2866*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  25.8562 + -94.1589*tf.T9i + -18.116*tf.T913i
                  + 1.86674*tf.T9 + -7.5666*tf.T953 + -0.666667*tf.lnT9)
    # nacrr
    rate += np.exp(  13.9232 + -97.4449*tf.T9i
                  + -0.21103*tf.T9 + 2.87702*tf.lnT9)
    # nacr 
    rate += np.exp(  14.7601 + -97.9108*tf.T9i
                  + -1.5*tf.lnT9)
    # nacr 
    rate += np.exp(  7.80363 + -96.6272*tf.T9i
                  + -1.5*tf.lnT9)

    rate_eval.He4_O16__p_F19 = rate

@numba.njit()
def He4_O17__p_F20(rate_eval, tf):
    # O17 + He4 --> p + F20
    rate = 0.0

    #  wagn
    rate += np.exp(  27.9452 + -65.6392*tf.T9i + -18.13*tf.T913i + -2.20482e-09*tf.T913
                  + 2.4033e-10*tf.T9 + -2.17403e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_O17__p_F20 = rate

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
def p_F20__He4_O17(rate_eval, tf):
    # F20 + p --> He4 + O17
    rate = 0.0

    #  wagn
    rate += np.exp(  29.27 + -4.58659e-12*tf.T9i + -18.13*tf.T913i + -2.20482e-09*tf.T913
                  + 2.4033e-10*tf.T9 + -2.17403e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_F20__He4_O17 = rate

@numba.njit()
def He3_He3__p_p_He4(rate_eval, tf):
    # He3 + He3 --> p + p + He4
    rate = 0.0

    # nacrn
    rate += np.exp(  24.7788 + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He3_He3__p_p_He4 = rate

@numba.njit()
def d_Be7__p_He4_He4(rate_eval, tf):
    # Be7 + d --> p + He4 + He4
    rate = 0.0

    # cf88n
    rate += np.exp(  27.6987 + -12.428*tf.T913i
                  + -0.666667*tf.lnT9)

    rate_eval.d_Be7__p_He4_He4 = rate

@numba.njit()
def He3_Be7__p_p_He4_He4(rate_eval, tf):
    # Be7 + He3 --> p + p + He4 + He4
    rate = 0.0

    # mafon
    rate += np.exp(  31.7435 + -5.45213e-12*tf.T9i + -21.793*tf.T913i + -1.98126e-09*tf.T913
                  + 1.84204e-10*tf.T9 + -1.46403e-11*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He3_Be7__p_p_He4_He4 = rate

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

@numba.njit()
def p_p_He4__He3_He3(rate_eval, tf):
    # p + p + He4 --> He3 + He3
    rate = 0.0

    # nacrn
    rate += np.exp(  2.98257 + -149.222*tf.T9i + -12.277*tf.T913i + -0.103699*tf.T913
                  + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -2.16667*tf.lnT9)

    rate_eval.p_p_He4__He3_He3 = rate

@numba.njit()
def p_He4_He4__d_Be7(rate_eval, tf):
    # p + He4 + He4 --> d + Be7
    rate = 0.0

    # cf88n
    rate += np.exp(  6.97069 + -194.561*tf.T9i + -12.428*tf.T913i
                  + -2.16667*tf.lnT9)

    rate_eval.p_He4_He4__d_Be7 = rate

@numba.njit()
def p_p_He4_He4__He3_Be7(rate_eval, tf):
    # p + p + He4 + He4 --> He3 + Be7
    rate = 0.0

    # mafon
    rate += np.exp(  -13.1807 + -130.807*tf.T9i + -21.793*tf.T913i + -1.98126e-09*tf.T913
                  + 1.84204e-10*tf.T9 + -1.46403e-11*tf.T953 + -3.66667*tf.lnT9)

    rate_eval.p_p_He4_He4__He3_Be7 = rate

@numba.njit()
def p_O15__F16(rate_eval, tf):
    # O15 + p --> F16
    rate = 0.0

    # rpsmr
    rate += np.exp(  -3.01222 + -8.51539*tf.T9i + 74.5955*tf.T913i + -73.794*tf.T913
                  + 2.87936*tf.T9 + -0.117988*tf.T953 + 47.5366*tf.lnT9)

    rate_eval.p_O15__F16 = rate

@numba.njit()
def F16__p_O15(rate_eval, tf):
    # F16 --> p + O15
    rate = 0.0

    # rpsmr
    rate += np.exp(  21.2899 + -2.29538*tf.T9i + 74.5955*tf.T913i + -73.794*tf.T913
                  + 2.87936*tf.T9 + -0.117988*tf.T953 + 49.0366*tf.lnT9)

    rate_eval.F16__p_O15 = rate

@numba.njit()
def p_S28__Cl29(rate_eval, tf):
    # S28 + p --> Cl29
    rate = 0.0

    # ths8r
    rate += np.exp(  -16.9432 + -20.6656*tf.T9i + 19.1618*tf.T913
                  + -2.43847*tf.T9 + 0.194115*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_S28__Cl29 = rate

@numba.njit()
def He4_S28__Ar32(rate_eval, tf):
    # S28 + He4 --> Ar32
    rate = 0.0

    # ths8r
    rate += np.exp(  40.0908 + -65.031*tf.T913i + 0.269344*tf.T913
                  + -0.818256*tf.T9 + 0.0515944*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_S28__Ar32 = rate

@numba.njit()
def p_S29__Cl30(rate_eval, tf):
    # S29 + p --> Cl30
    rate = 0.0

    # ths8r
    rate += np.exp(  -19.6182 + -3.6393*tf.T9i + 21.3406*tf.T913
                  + -3.25229*tf.T9 + 0.307741*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_S29__Cl30 = rate

@numba.njit()
def He4_S29__Ar33(rate_eval, tf):
    # S29 + He4 --> Ar33
    rate = 0.0

    # ths8r
    rate += np.exp(  41.2275 + -65.124*tf.T913i + -0.021112*tf.T913
                  + -0.799069*tf.T9 + 0.0579505*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_S29__Ar33 = rate

@numba.njit()
def Cl29__S29__weak__bqa_pos_(rate_eval, tf):
    # Cl29 --> S29
    rate = 0.0

    # bqa+w
    rate += np.exp(  4.1343)

    rate_eval.Cl29__S29__weak__bqa_pos_ = rate

@numba.njit()
def Cl29__p_S28(rate_eval, tf):
    # Cl29 --> p + S28
    rate = 0.0

    # ths8r
    rate += np.exp(  5.32364 + 19.1618*tf.T913
                  + -2.43847*tf.T9 + 0.194115*tf.T953)

    rate_eval.Cl29__p_S28 = rate

@numba.njit()
def He4_Cl29__K33(rate_eval, tf):
    # Cl29 + He4 --> K33
    rate = 0.0

    # ths8r
    rate += np.exp(  40.6053 + -67.8117*tf.T913i + -0.294831*tf.T913
                  + -0.394384*tf.T9 + 0.0004138*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cl29__K33 = rate

@numba.njit()
def He4_Cl29__p_Ar32(rate_eval, tf):
    # Cl29 + He4 --> p + Ar32
    rate = 0.0

    # ths8r
    rate += np.exp(  51.5091 + -67.8117*tf.T913i + 0.267424*tf.T913
                  + -0.726288*tf.T9 + 0.0481218*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cl29__p_Ar32 = rate

@numba.njit()
def Ar32__He4_S28(rate_eval, tf):
    # Ar32 --> He4 + S28
    rate = 0.0

    # ths8r
    rate += np.exp(  64.9826 + -100.94*tf.T9i + -65.031*tf.T913i + 0.269344*tf.T913
                  + -0.818256*tf.T9 + 0.0515944*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ar32__He4_S28 = rate

@numba.njit()
def p_Ar32__K33(rate_eval, tf):
    # Ar32 + p --> K33
    rate = 0.0

    # ths8r
    rate += np.exp(  -17.7268 + -19.4286*tf.T9i + 20.2635*tf.T913
                  + -2.76523*tf.T9 + 0.235117*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ar32__K33 = rate

@numba.njit()
def He4_Ar32__Ca36(rate_eval, tf):
    # Ar32 + He4 --> Ca36
    rate = 0.0

    # ths8r
    rate += np.exp(  43.8241 + -70.713*tf.T913i + -0.0778489*tf.T913
                  + -0.466075*tf.T9 + -0.00979269*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ar32__Ca36 = rate

@numba.njit()
def p_Ar32__He4_Cl29(rate_eval, tf):
    # Ar32 + p --> He4 + Cl29
    rate = 0.0

    # ths8r
    rate += np.exp(  54.134 + -121.605*tf.T9i + -67.8117*tf.T913i + 0.267424*tf.T913
                  + -0.726288*tf.T9 + 0.0481218*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ar32__He4_Cl29 = rate

@numba.njit()
def p_Ca38__Sc39(rate_eval, tf):
    # Ca38 + p --> Sc39
    rate = 0.0

    # ths8r
    rate += np.exp(  -17.7938 + -6.98055*tf.T9i + 20.4307*tf.T913
                  + -2.70468*tf.T9 + 0.242514*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ca38__Sc39 = rate

@numba.njit()
def p_p_Ca38__He4_Ca36(rate_eval, tf):
    # Ca38 + p + p --> He4 + Ca36
    rate = 0.0

    # fkthr
    rate += np.exp(  39.1952 + -41.7127*tf.T9i + -27.3832*tf.T913i + -55.6339*tf.T913
                  + 2.1009*tf.T9 + -0.113521*tf.T953 + 28.5739*tf.lnT9)

    rate_eval.p_p_Ca38__He4_Ca36 = rate

@numba.njit()
def Cl30__p_S29(rate_eval, tf):
    # Cl30 --> p + S29
    rate = 0.0

    # ths8r
    rate += np.exp(  3.88248 + 0.000155511*tf.T9i + 21.3406*tf.T913
                  + -3.25229*tf.T9 + 0.307741*tf.T953)

    rate_eval.Cl30__p_S29 = rate

@numba.njit()
def He4_Cl30__K34(rate_eval, tf):
    # Cl30 + He4 --> K34
    rate = 0.0

    # ths8r
    rate += np.exp(  41.9433 + -67.9024*tf.T913i + -0.618087*tf.T913
                  + -0.498925*tf.T9 + 0.0172983*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cl30__K34 = rate

@numba.njit()
def He4_Cl30__p_Ar33(rate_eval, tf):
    # Cl30 + He4 --> p + Ar33
    rate = 0.0

    # ths8r
    rate += np.exp(  52.3051 + -67.9024*tf.T913i + 0.19456*tf.T913
                  + -0.781433*tf.T9 + 0.0533255*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cl30__p_Ar33 = rate

@numba.njit()
def Ar33__He4_S29(rate_eval, tf):
    # Ar33 --> He4 + S29
    rate = 0.0

    # ths8r
    rate += np.exp(  67.2244 + -100.373*tf.T9i + -65.124*tf.T913i + -0.021112*tf.T913
                  + -0.799069*tf.T9 + 0.0579505*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ar33__He4_S29 = rate

@numba.njit()
def p_Ar33__K34(rate_eval, tf):
    # Ar33 + p --> K34
    rate = 0.0

    # ths8r
    rate += np.exp(  -18.9682 + -7.12714*tf.T9i + 20.7703*tf.T913
                  + -2.95588*tf.T9 + 0.257828*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ar33__K34 = rate

@numba.njit()
def He4_Ar33__Ca37(rate_eval, tf):
    # Ar33 + He4 --> Ca37
    rate = 0.0

    # ths8r
    rate += np.exp(  45.9124 + -70.792*tf.T913i + -3.32023*tf.T913
                  + 0.874694*tf.T9 + -0.191591*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ar33__Ca37 = rate

@numba.njit()
def p_Ar33__He4_Cl30(rate_eval, tf):
    # Ar33 + p --> He4 + Cl30
    rate = 0.0

    # ths8r
    rate += np.exp(  54.8012 + -104.017*tf.T9i + -67.9024*tf.T913i + 0.19456*tf.T913
                  + -0.781433*tf.T9 + 0.0533255*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ar33__He4_Cl30 = rate

@numba.njit()
def K33__Ar33__weak__bqa_pos_(rate_eval, tf):
    # K33 --> Ar33
    rate = 0.0

    # bqa+w
    rate += np.exp(  2.81555)

    rate_eval.K33__Ar33__weak__bqa_pos_ = rate

@numba.njit()
def K33__p_Ar32(rate_eval, tf):
    # K33 --> p + Ar32
    rate = 0.0

    # ths8r
    rate += np.exp(  5.23957 + 0.000838384*tf.T9i + 20.2635*tf.T913
                  + -2.76523*tf.T9 + 0.235117*tf.T953)

    rate_eval.K33__p_Ar32 = rate

@numba.njit()
def K33__He4_Cl29(rate_eval, tf):
    # K33 --> He4 + Cl29
    rate = 0.0

    # ths8r
    rate += np.exp(  66.1966 + -102.181*tf.T9i + -67.8117*tf.T913i + -0.294831*tf.T913
                  + -0.394384*tf.T9 + 0.0004138*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.K33__He4_Cl29 = rate

@numba.njit()
def He4_K33__Sc37(rate_eval, tf):
    # K33 + He4 --> Sc37
    rate = 0.0

    # ths8r
    rate += np.exp(  43.4914 + -73.3916*tf.T913i + -0.362928*tf.T913
                  + -0.500382*tf.T9 + 0.00594182*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_K33__Sc37 = rate

@numba.njit()
def He4_K33__p_Ca36(rate_eval, tf):
    # K33 + He4 --> p + Ca36
    rate = 0.0

    # ths8r
    rate += np.exp(  53.6861 + -73.3916*tf.T913i + 0.478953*tf.T913
                  + -0.846479*tf.T9 + 0.0605207*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_K33__p_Ca36 = rate

@numba.njit()
def Ca36__He4_Ar32(rate_eval, tf):
    # Ca36 --> He4 + Ar32
    rate = 0.0

    # ths8r
    rate += np.exp(  68.7395 + -77.3332*tf.T9i + -70.713*tf.T913i + -0.0778489*tf.T913
                  + -0.466075*tf.T9 + -0.00979269*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ca36__He4_Ar32 = rate

@numba.njit()
def p_Ca36__Sc37(rate_eval, tf):
    # Ca36 + p --> Sc37
    rate = 0.0

    # ths8r
    rate += np.exp(  -14.5883 + -23.1091*tf.T9i + 16.7352*tf.T913
                  + -1.67848*tf.T9 + 0.120144*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ca36__Sc37 = rate

@numba.njit()
def He4_Ca36__Ti40(rate_eval, tf):
    # Ca36 + He4 --> Ti40
    rate = 0.0

    # ths8r
    rate += np.exp(  49.6465 + -76.1732*tf.T913i + -4.87682*tf.T913
                  + 1.31056*tf.T9 + -0.240342*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ca36__Ti40 = rate

@numba.njit()
def p_Ca36__He4_K33(rate_eval, tf):
    # Ca36 + p --> He4 + K33
    rate = 0.0

    # ths8r
    rate += np.exp(  55.635 + -96.766*tf.T9i + -73.3916*tf.T913i + 0.478953*tf.T913
                  + -0.846479*tf.T9 + 0.0605207*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ca36__He4_K33 = rate

@numba.njit()
def He4_Ca36__p_Sc39(rate_eval, tf):
    # Ca36 + He4 --> p + Sc39
    rate = 0.0

    # ths8r
    rate += np.exp(  55.7563 + -76.1732*tf.T913i + -0.623677*tf.T913
                  + -0.571836*tf.T9 + 0.0260149*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ca36__p_Sc39 = rate

@numba.njit()
def He4_Ca36__p_p_Ca38(rate_eval, tf):
    # Ca36 + He4 --> p + p + Ca38
    rate = 0.0

    # fkthr
    rate += np.exp(  60.9026 + -1.47986*tf.T9i + -27.3832*tf.T913i + -55.6339*tf.T913
                  + 2.1009*tf.T9 + -0.113521*tf.T953 + 30.0739*tf.lnT9)

    rate_eval.He4_Ca36__p_p_Ca38 = rate

@numba.njit()
def Sc39__p_Ca38(rate_eval, tf):
    # Sc39 --> p + Ca38
    rate = 0.0

    # ths8r
    rate += np.exp(  3.79347 + 0.000298215*tf.T9i + 20.4307*tf.T913
                  + -2.70468*tf.T9 + 0.242514*tf.T953)

    rate_eval.Sc39__p_Ca38 = rate

@numba.njit()
def p_Sc39__Ti40(rate_eval, tf):
    # Sc39 + p --> Ti40
    rate = 0.0

    # ths8r
    rate += np.exp(  42.3709 + -32.1484*tf.T913i + -15.0165*tf.T913
                  + 1.54165*tf.T9 + -0.0793546*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Sc39__Ti40 = rate

@numba.njit()
def p_Sc39__He4_Ca36(rate_eval, tf):
    # Sc39 + p --> He4 + Ca36
    rate = 0.0

    # ths8r
    rate += np.exp(  54.9431 + -33.2438*tf.T9i + -76.1732*tf.T913i + -0.623677*tf.T913
                  + -0.571836*tf.T9 + 0.0260149*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Sc39__He4_Ca36 = rate

@numba.njit()
def K34__p_Ar33(rate_eval, tf):
    # K34 --> p + Ar33
    rate = 0.0

    # ths8r
    rate += np.exp(  4.28733 + 20.7703*tf.T913
                  + -2.95588*tf.T9 + 0.257828*tf.T953)

    rate_eval.K34__p_Ar33 = rate

@numba.njit()
def K34__He4_Cl30(rate_eval, tf):
    # K34 --> He4 + Cl30
    rate = 0.0

    # ths8r
    rate += np.exp(  67.6948 + -96.8892*tf.T9i + -67.9024*tf.T913i + -0.618087*tf.T913
                  + -0.498925*tf.T9 + 0.0172983*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.K34__He4_Cl30 = rate

@numba.njit()
def p_K34__Ca35(rate_eval, tf):
    # K34 + p --> Ca35
    rate = 0.0

    # ths8r
    rate += np.exp(  36.6033 + -30.0364*tf.T913i + -5.6016*tf.T913
                  + -1.67829*tf.T9 + 0.284814*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_K34__Ca35 = rate

@numba.njit()
def He4_K34__Sc38(rate_eval, tf):
    # K34 + He4 --> Sc38
    rate = 0.0

    # ths8r
    rate += np.exp(  44.1795 + -73.4689*tf.T913i + -1.11837*tf.T913
                  + -0.423848*tf.T9 + 0.00242462*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_K34__Sc38 = rate

@numba.njit()
def He4_K34__p_Ca37(rate_eval, tf):
    # K34 + He4 --> p + Ca37
    rate = 0.0

    # ths8r
    rate += np.exp(  54.6796 + -73.4689*tf.T913i + -0.350491*tf.T913
                  + -0.675654*tf.T9 + 0.0423497*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_K34__p_Ca37 = rate

@numba.njit()
def Ca37__He4_Ar33(rate_eval, tf):
    # Ca37 --> He4 + Ar33
    rate = 0.0

    # ths8r
    rate += np.exp(  70.1397 + -71.6748*tf.T9i + -70.792*tf.T913i + -3.32023*tf.T913
                  + 0.874694*tf.T9 + -0.191591*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ca37__He4_Ar33 = rate

@numba.njit()
def p_Ca37__Sc38(rate_eval, tf):
    # Ca37 + p --> Sc38
    rate = 0.0

    # ths8r
    rate += np.exp(  -18.4754 + -10.557*tf.T9i + 19.6599*tf.T913
                  + -2.5569*tf.T9 + 0.213755*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ca37__Sc38 = rate

@numba.njit()
def p_Ca37__He4_K34(rate_eval, tf):
    # Ca37 + p --> He4 + K34
    rate = 0.0

    # ths8r
    rate += np.exp(  55.6514 + -78.802*tf.T9i + -73.4689*tf.T913i + -0.350491*tf.T913
                  + -0.675654*tf.T9 + 0.0423497*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ca37__He4_K34 = rate

@numba.njit()
def Sc37__Ca37__weak__bqa_pos_(rate_eval, tf):
    # Sc37 --> Ca37
    rate = 0.0

    # bqa+w
    rate += np.exp(  2.51589)

    rate_eval.Sc37__Ca37__weak__bqa_pos_ = rate

@numba.njit()
def Sc37__p_Ca36(rate_eval, tf):
    # Sc37 --> p + Ca36
    rate = 0.0

    # ths8r
    rate += np.exp(  6.99693 + 16.7352*tf.T913
                  + -1.67848*tf.T9 + 0.120144*tf.T953)

    rate_eval.Sc37__p_Ca36 = rate

@numba.njit()
def Sc37__He4_K33(rate_eval, tf):
    # Sc37 --> He4 + K33
    rate = 0.0

    # ths8r
    rate += np.exp(  67.0255 + -73.656*tf.T9i + -73.3916*tf.T913i + -0.362928*tf.T913
                  + -0.500382*tf.T9 + 0.00594182*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Sc37__He4_K33 = rate

@numba.njit()
def He4_Sc37__V41(rate_eval, tf):
    # Sc37 + He4 --> V41
    rate = 0.0

    # ths8r
    rate += np.exp(  46.1396 + -78.7634*tf.T913i + -0.810841*tf.T913
                  + -0.485009*tf.T9 + 0.00342062*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Sc37__V41 = rate

@numba.njit()
def He4_Sc37__p_Ti40(rate_eval, tf):
    # Sc37 + He4 --> p + Ti40
    rate = 0.0

    # ths8r
    rate += np.exp(  56.9345 + -78.7634*tf.T913i + -1.03282*tf.T913
                  + -0.424405*tf.T9 + 0.00880205*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Sc37__p_Ti40 = rate

@numba.njit()
def Ti40__p_Sc39(rate_eval, tf):
    # Ti40 --> p + Sc39
    rate = 0.0

    # ths8r
    rate += np.exp(  68.118 + -22.8761*tf.T9i + -32.1484*tf.T913i + -15.0165*tf.T913
                  + 1.54165*tf.T9 + -0.0793546*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ti40__p_Sc39 = rate

@numba.njit()
def Ti40__He4_Ca36(rate_eval, tf):
    # Ti40 --> He4 + Ca36
    rate = 0.0

    # ths8r
    rate += np.exp(  74.5805 + -56.1174*tf.T9i + -76.1732*tf.T913i + -4.87682*tf.T913
                  + 1.31056*tf.T9 + -0.240342*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ti40__He4_Ca36 = rate

@numba.njit()
def p_Ti40__V41(rate_eval, tf):
    # Ti40 + p --> V41
    rate = 0.0

    # ths8r
    rate += np.exp(  -16.7267 + -15.7392*tf.T9i + 17.7727*tf.T913
                  + -1.89337*tf.T9 + 0.124478*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ti40__V41 = rate

@numba.njit()
def He4_Ti40__Cr44(rate_eval, tf):
    # Ti40 + He4 --> Cr44
    rate = 0.0

    # ths8r
    rate += np.exp(  50.2229 + -81.4428*tf.T913i + -3.13395*tf.T913
                  + 0.0095614*tf.T9 + -0.0629252*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ti40__Cr44 = rate

@numba.njit()
def p_Ti40__He4_Sc37(rate_eval, tf):
    # Ti40 + p --> He4 + Sc37
    rate = 0.0

    # ths8r
    rate += np.exp(  60.2833 + -79.2265*tf.T9i + -78.7634*tf.T913i + -1.03282*tf.T913
                  + -0.424405*tf.T9 + 0.00880205*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ti40__He4_Sc37 = rate

@numba.njit()
def Ca35__p_K34(rate_eval, tf):
    # Ca35 --> p + K34
    rate = 0.0

    # ths8r
    rate += np.exp(  60.6711 + -13.9947*tf.T9i + -30.0364*tf.T913i + -5.6016*tf.T913
                  + -1.67829*tf.T9 + 0.284814*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ca35__p_K34 = rate

@numba.njit()
def p_Ca35__Sc36(rate_eval, tf):
    # Ca35 + p --> Sc36
    rate = 0.0

    # ths8r
    rate += np.exp(  -17.041 + -23.2906*tf.T9i + 19.5625*tf.T913
                  + -2.69513*tf.T9 + 0.226635*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ca35__Sc36 = rate

@numba.njit()
def He4_Ca35__Ti39(rate_eval, tf):
    # Ca35 + He4 --> Ti39
    rate = 0.0

    # ths8r
    rate += np.exp(  45.1043 + -76.1016*tf.T913i + -0.913501*tf.T913
                  + -0.502372*tf.T9 + 0.00757111*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ca35__Ti39 = rate

@numba.njit()
def He4_Ca35__p_Sc38(rate_eval, tf):
    # Ca35 + He4 --> p + Sc38
    rate = 0.0

    # ths8r
    rate += np.exp(  55.7921 + -76.1016*tf.T913i + -0.977608*tf.T913
                  + -0.399735*tf.T9 + 0.00566881*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ca35__p_Sc38 = rate

@numba.njit()
def Sc38__Ca38__weak__mo97(rate_eval, tf):
    # Sc38 --> Ca38
    rate = 0.0

    # mo97w
    rate += np.exp(  2.79645)

    rate_eval.Sc38__Ca38__weak__mo97 = rate

@numba.njit()
def Sc38__p_Ca37(rate_eval, tf):
    # Sc38 --> p + Ca37
    rate = 0.0

    # ths8r
    rate += np.exp(  4.96722 + 19.6599*tf.T913
                  + -2.5569*tf.T9 + 0.213755*tf.T953)

    rate_eval.Sc38__p_Ca37 = rate

@numba.njit()
def Sc38__He4_K34(rate_eval, tf):
    # Sc38 --> He4 + K34
    rate = 0.0

    # ths8r
    rate += np.exp(  68.5939 + -68.2451*tf.T9i + -73.4689*tf.T913i + -1.11837*tf.T913
                  + -0.423848*tf.T9 + 0.00242462*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Sc38__He4_K34 = rate

@numba.njit()
def p_Sc38__Ti39(rate_eval, tf):
    # Sc38 + p --> Ti39
    rate = 0.0

    # ths8r
    rate += np.exp(  37.621 + -32.1413*tf.T913i + -4.55276*tf.T913
                  + -2.24702*tf.T9 + 0.357394*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Sc38__Ti39 = rate

@numba.njit()
def p_Sc38__He4_Ca35(rate_eval, tf):
    # Sc38 + p --> He4 + Ca35
    rate = 0.0

    # ths8r
    rate += np.exp(  56.1387 + -54.2505*tf.T9i + -76.1016*tf.T913i + -0.977608*tf.T913
                  + -0.399735*tf.T9 + 0.00566881*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Sc38__He4_Ca35 = rate

@numba.njit()
def V41__p_Ti40(rate_eval, tf):
    # V41 --> p + Ti40
    rate = 0.0

    # ths8r
    rate += np.exp(  6.24879 + 0.000660879*tf.T9i + 17.7727*tf.T913
                  + -1.89337*tf.T9 + 0.124478*tf.T953)

    rate_eval.V41__p_Ti40 = rate

@numba.njit()
def V41__He4_Sc37(rate_eval, tf):
    # V41 --> He4 + Sc37
    rate = 0.0

    # ths8r
    rate += np.exp(  72.4639 + -63.49*tf.T9i + -78.7634*tf.T913i + -0.810841*tf.T913
                  + -0.485009*tf.T9 + 0.00342062*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.V41__He4_Sc37 = rate

@numba.njit()
def He4_V41__Mn45(rate_eval, tf):
    # V41 + He4 --> Mn45
    rate = 0.0

    # ths8r
    rate += np.exp(  49.1422 + -83.9551*tf.T913i + -1.22789*tf.T913
                  + -0.494458*tf.T9 + 0.00588069*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_V41__Mn45 = rate

@numba.njit()
def He4_V41__p_Cr44(rate_eval, tf):
    # V41 + He4 --> p + Cr44
    rate = 0.0

    # ths8r
    rate += np.exp(  59.3909 + -83.9551*tf.T913i + -1.57038*tf.T913
                  + -0.345223*tf.T9 + 0.0032652*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_V41__p_Cr44 = rate

@numba.njit()
def Cr44__He4_Ti40(rate_eval, tf):
    # Cr44 --> He4 + Ti40
    rate = 0.0

    # ths8r
    rate += np.exp(  75.172 + -81.6451*tf.T9i + -81.4428*tf.T913i + -3.13395*tf.T913
                  + 0.0095614*tf.T9 + -0.0629252*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Cr44__He4_Ti40 = rate

@numba.njit()
def p_Cr44__Mn45(rate_eval, tf):
    # Cr44 + p --> Mn45
    rate = 0.0

    # ths8r
    rate += np.exp(  -16.0243 + -12.3249*tf.T9i + 17.9363*tf.T913
                  + -1.84483*tf.T9 + 0.126945*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Cr44__Mn45 = rate

@numba.njit()
def He4_Cr44__Fe48(rate_eval, tf):
    # Cr44 + He4 --> Fe48
    rate = 0.0

    # ths8r
    rate += np.exp(  52.0335 + -86.5458*tf.T913i + 1.97255*tf.T913
                  + -1.56115*tf.T9 + 0.065487*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cr44__Fe48 = rate

@numba.njit()
def p_Cr44__He4_V41(rate_eval, tf):
    # Cr44 + p --> He4 + V41
    rate = 0.0

    # ths8r
    rate += np.exp(  61.3644 + -97.3885*tf.T9i + -83.9551*tf.T913i + -1.57038*tf.T913
                  + -0.345223*tf.T9 + 0.0032652*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Cr44__He4_V41 = rate

@numba.njit()
def Sc36__Ca36__weak__bqa_pos_(rate_eval, tf):
    # Sc36 --> Ca36
    rate = 0.0

    # bqa+w
    rate += np.exp(  4.36649)

    rate_eval.Sc36__Ca36__weak__bqa_pos_ = rate

@numba.njit()
def Sc36__p_Ca35(rate_eval, tf):
    # Sc36 --> p + Ca35
    rate = 0.0

    # ths8r
    rate += np.exp(  5.11844 + 19.5625*tf.T913
                  + -2.69513*tf.T9 + 0.226635*tf.T953)

    rate_eval.Sc36__p_Ca35 = rate

@numba.njit()
def He4_Sc36__V40(rate_eval, tf):
    # Sc36 + He4 --> V40
    rate = 0.0

    # ths8r
    rate += np.exp(  44.4873 + -78.6932*tf.T913i + -0.0530652*tf.T913
                  + -0.629411*tf.T9 + 0.0237465*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Sc36__V40 = rate

@numba.njit()
def He4_Sc36__p_Ti39(rate_eval, tf):
    # Sc36 + He4 --> p + Ti39
    rate = 0.0

    # ths8r
    rate += np.exp(  56.9288 + -78.6932*tf.T913i + -1.07636*tf.T913
                  + -0.414773*tf.T9 + 0.00798428*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Sc36__p_Ti39 = rate

@numba.njit()
def Ti39__Sc39__weak__wc17(rate_eval, tf):
    # Ti39 --> Sc39
    rate = 0.0

    # wc17w
    rate += np.exp(  -112.022)

    rate_eval.Ti39__Sc39__weak__wc17 = rate

@numba.njit()
def Ti39__p_Sc38(rate_eval, tf):
    # Ti39 --> p + Sc38
    rate = 0.0

    # ths8r
    rate += np.exp(  61.5109 + -9.88712*tf.T9i + -32.1413*tf.T913i + -4.55276*tf.T913
                  + -2.24702*tf.T9 + 0.357394*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ti39__p_Sc38 = rate

@numba.njit()
def Ti39__p_Ca38__weak__wc12(rate_eval, tf):
    # Ti39 --> p + Ca38
    rate = 0.0

    # wc12w
    rate += np.exp(  3.10725)

    rate_eval.Ti39__p_Ca38__weak__wc12 = rate

@numba.njit()
def Ti39__He4_Ca35(rate_eval, tf):
    # Ti39 --> He4 + Ca35
    rate = 0.0

    # ths8r
    rate += np.exp(  69.3408 + -64.1399*tf.T9i + -76.1016*tf.T913i + -0.913501*tf.T913
                  + -0.502372*tf.T9 + 0.00757111*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Ti39__He4_Ca35 = rate

@numba.njit()
def p_Ti39__V40(rate_eval, tf):
    # Ti39 + p --> V40
    rate = 0.0

    # ths8r
    rate += np.exp(  -17.7509 + -17.8829*tf.T9i + 18.7676*tf.T913
                  + -2.45087*tf.T9 + 0.201364*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Ti39__V40 = rate

@numba.njit()
def He4_Ti39__Cr43(rate_eval, tf):
    # Ti39 + He4 --> Cr43
    rate = 0.0

    # ths8r
    rate += np.exp(  48.236 + -81.3803*tf.T913i + -1.18001*tf.T913
                  + -0.50423*tf.T9 + 0.00423975*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Ti39__Cr43 = rate

@numba.njit()
def p_Ti39__He4_Sc36(rate_eval, tf):
    # Ti39 + p --> He4 + Sc36
    rate = 0.0

    # ths8r
    rate += np.exp(  59.006 + -87.4314*tf.T9i + -78.6932*tf.T913i + -1.07636*tf.T913
                  + -0.414773*tf.T9 + 0.00798428*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Ti39__He4_Sc36 = rate

@numba.njit()
def Mn45__p_Cr44(rate_eval, tf):
    # Mn45 --> p + Cr44
    rate = 0.0

    # ths8r
    rate += np.exp(  5.85591 + 0.000537914*tf.T9i + 17.9363*tf.T913
                  + -1.84483*tf.T9 + 0.126945*tf.T953)

    rate_eval.Mn45__p_Cr44 = rate

@numba.njit()
def Mn45__He4_V41(rate_eval, tf):
    # Mn45 --> He4 + V41
    rate = 0.0

    # ths8r
    rate += np.exp(  72.9959 + -85.0631*tf.T9i + -83.9551*tf.T913i + -1.22789*tf.T913
                  + -0.494458*tf.T9 + 0.00588069*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mn45__He4_V41 = rate

@numba.njit()
def p_Mn45__Fe46(rate_eval, tf):
    # Mn45 + p --> Fe46
    rate = 0.0

    # ths8r
    rate += np.exp(  41.2369 + -36.1516*tf.T913i + -11.3353*tf.T913
                  + 0.290321*tf.T9 + 0.0486483*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mn45__Fe46 = rate

@numba.njit()
def He4_Mn45__Co49(rate_eval, tf):
    # Mn45 + He4 --> Co49
    rate = 0.0

    # ths8r
    rate += np.exp(  51.4799 + -88.989*tf.T913i + -1.46275*tf.T913
                  + -0.528983*tf.T9 + 0.0100936*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mn45__Co49 = rate

@numba.njit()
def He4_Mn45__p_Fe48(rate_eval, tf):
    # Mn45 + He4 --> p + Fe48
    rate = 0.0

    # ths8r
    rate += np.exp(  61.3147 + -88.989*tf.T913i + -1.30257*tf.T913
                  + -0.417112*tf.T9 + 0.00356355*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mn45__p_Fe48 = rate

@numba.njit()
def Fe48__He4_Cr44(rate_eval, tf):
    # Fe48 --> He4 + Cr44
    rate = 0.0

    # ths8r
    rate += np.exp(  76.995 + -82.6695*tf.T9i + -86.5458*tf.T913i + 1.97255*tf.T913
                  + -1.56115*tf.T9 + 0.065487*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Fe48__He4_Cr44 = rate

@numba.njit()
def p_Fe48__Co49(rate_eval, tf):
    # Fe48 + p --> Co49
    rate = 0.0

    # ths8r
    rate += np.exp(  -17.0528 + -10.5047*tf.T9i + 19.0735*tf.T913
                  + -2.23228*tf.T9 + 0.159842*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Fe48__Co49 = rate

@numba.njit()
def p_Fe48__He4_Mn45(rate_eval, tf):
    # Fe48 + p --> He4 + Mn45
    rate = 0.0

    # ths8r
    rate += np.exp(  64.396 + -94.9985*tf.T9i + -88.989*tf.T913i + -1.30257*tf.T913
                  + -0.417112*tf.T9 + 0.00356355*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Fe48__He4_Mn45 = rate

@numba.njit()
def V40__Ti40__weak__bqa_pos_(rate_eval, tf):
    # V40 --> Ti40
    rate = 0.0

    # bqa+w
    rate += np.exp(  4.50002)

    rate_eval.V40__Ti40__weak__bqa_pos_ = rate

@numba.njit()
def V40__p_Ti39(rate_eval, tf):
    # V40 --> p + Ti39
    rate = 0.0

    # ths8r
    rate += np.exp(  5.69366 + 0.000752445*tf.T9i + 18.7676*tf.T913
                  + -2.45087*tf.T9 + 0.201364*tf.T953)

    rate_eval.V40__p_Ti39 = rate

@numba.njit()
def V40__He4_Sc36(rate_eval, tf):
    # V40 --> He4 + Sc36
    rate = 0.0

    # ths8r
    rate += np.exp(  70.0091 + -69.5448*tf.T9i + -78.6932*tf.T913i + -0.0530652*tf.T913
                  + -0.629411*tf.T9 + 0.0237465*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.V40__He4_Sc36 = rate

@numba.njit()
def He4_V40__Mn44(rate_eval, tf):
    # V40 + He4 --> Mn44
    rate = 0.0

    # ths8r
    rate += np.exp(  48.0646 + -83.8937*tf.T913i + -0.548514*tf.T913
                  + -0.638906*tf.T9 + 0.0171963*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_V40__Mn44 = rate

@numba.njit()
def He4_V40__p_Cr43(rate_eval, tf):
    # V40 + He4 --> p + Cr43
    rate = 0.0

    # ths8r
    rate += np.exp(  59.9001 + -83.8937*tf.T913i + -1.36185*tf.T913
                  + -0.476387*tf.T9 + 0.0148535*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_V40__p_Cr43 = rate

@numba.njit()
def Cr43__He4_Ti39(rate_eval, tf):
    # Cr43 --> He4 + Ti39
    rate = 0.0

    # ths8r
    rate += np.exp(  73.1815 + -70.3022*tf.T9i + -81.3803*tf.T913i + -1.18001*tf.T913
                  + -0.50423*tf.T9 + 0.00423975*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Cr43__He4_Ti39 = rate

@numba.njit()
def p_Cr43__Mn44(rate_eval, tf):
    # Cr43 + p --> Mn44
    rate = 0.0

    # ths8r
    rate += np.exp(  -19.4919 + -14.4247*tf.T9i + 19.201*tf.T913
                  + -2.67118*tf.T9 + 0.226275*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Cr43__Mn44 = rate

@numba.njit()
def He4_Cr43__Fe47(rate_eval, tf):
    # Cr43 + He4 --> Fe47
    rate = 0.0

    # ths8r
    rate += np.exp(  50.9622 + -86.4906*tf.T913i + -1.66143*tf.T913
                  + -0.451078*tf.T9 + 6.28872e-05*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Cr43__Fe47 = rate

@numba.njit()
def p_Cr43__He4_V40(rate_eval, tf):
    # Cr43 + p --> He4 + V40
    rate = 0.0

    # ths8r
    rate += np.exp(  61.4011 + -88.182*tf.T9i + -83.8937*tf.T913i + -1.36185*tf.T913
                  + -0.476387*tf.T9 + 0.0148535*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Cr43__He4_V40 = rate

@numba.njit()
def Fe46__p_Mn45(rate_eval, tf):
    # Fe46 --> p + Mn45
    rate = 0.0

    # ths8r
    rate += np.exp(  66.7014 + -12.8479*tf.T9i + -36.1516*tf.T913i + -11.3353*tf.T913
                  + 0.290321*tf.T9 + 0.0486483*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Fe46__p_Mn45 = rate

@numba.njit()
def p_Fe46__Co47(rate_eval, tf):
    # Fe46 + p --> Co47
    rate = 0.0

    # ths8r
    rate += np.exp(  -18.5553 + -17.6735*tf.T9i + 20.309*tf.T913
                  + -2.64768*tf.T9 + 0.201647*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Fe46__Co47 = rate

@numba.njit()
def He4_Fe46__p_Co49(rate_eval, tf):
    # Fe46 + He4 --> p + Co49
    rate = 0.0

    # ths8r
    rate += np.exp(  62.2446 + -91.401*tf.T913i + -1.38716*tf.T913
                  + -0.407203*tf.T9 + 0.00261659*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Fe46__p_Co49 = rate

@numba.njit()
def Co49__p_Fe48(rate_eval, tf):
    # Co49 --> p + Fe48
    rate = 0.0

    # ths8r
    rate += np.exp(  5.23566 + 0.000493736*tf.T9i + 19.0735*tf.T913
                  + -2.23228*tf.T9 + 0.159842*tf.T953)

    rate_eval.Co49__p_Fe48 = rate

@numba.njit()
def Co49__He4_Mn45(rate_eval, tf):
    # Co49 --> He4 + Mn45
    rate = 0.0

    # ths8r
    rate += np.exp(  76.8496 + -84.4933*tf.T9i + -88.989*tf.T913i + -1.46275*tf.T913
                  + -0.528983*tf.T9 + 0.0100936*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Co49__He4_Mn45 = rate

@numba.njit()
def p_Co49__He4_Fe46(rate_eval, tf):
    # Co49 + p --> He4 + Fe46
    rate = 0.0

    # ths8r
    rate += np.exp(  62.1498 + -71.6454*tf.T9i + -91.401*tf.T913i + -1.38716*tf.T913
                  + -0.407203*tf.T9 + 0.00261659*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Co49__He4_Fe46 = rate

@numba.njit()
def Mn44__Cr44__weak__bqa_pos_(rate_eval, tf):
    # Mn44 --> Cr44
    rate = 0.0

    # bqa+w
    rate += np.exp(  3.95348)

    rate_eval.Mn44__Cr44__weak__bqa_pos_ = rate

@numba.njit()
def Mn44__p_Cr43(rate_eval, tf):
    # Mn44 --> p + Cr43
    rate = 0.0

    # ths8r
    rate += np.exp(  3.95624 + 19.201*tf.T913
                  + -2.67118*tf.T9 + 0.226275*tf.T953)

    rate_eval.Mn44__p_Cr43 = rate

@numba.njit()
def Mn44__He4_V40(rate_eval, tf):
    # Mn44 --> He4 + V40
    rate = 0.0

    # ths8r
    rate += np.exp(  73.0137 + -73.7573*tf.T9i + -83.8937*tf.T913i + -0.548514*tf.T913
                  + -0.638906*tf.T9 + 0.0171963*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Mn44__He4_V40 = rate

@numba.njit()
def p_Mn44__Fe45(rate_eval, tf):
    # Mn44 + p --> Fe45
    rate = 0.0

    # ths8r
    rate += np.exp(  40.6028 + -36.1457*tf.T913i + -6.20399*tf.T913
                  + -2.04399*tf.T9 + 0.328581*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Mn44__Fe45 = rate

@numba.njit()
def He4_Mn44__Co48(rate_eval, tf):
    # Mn44 + He4 --> Co48
    rate = 0.0

    # ths8r
    rate += np.exp(  52.4663 + -88.9347*tf.T913i + -2.41015*tf.T913
                  + -0.458817*tf.T9 + 0.0164617*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mn44__Co48 = rate

@numba.njit()
def He4_Mn44__p_Fe47(rate_eval, tf):
    # Mn44 + He4 --> p + Fe47
    rate = 0.0

    # ths8r
    rate += np.exp(  63.3625 + -88.9347*tf.T913i + -2.86881*tf.T913
                  + -0.327748*tf.T9 + 0.012567*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Mn44__p_Fe47 = rate

@numba.njit()
def Fe47__He4_Cr43(rate_eval, tf):
    # Fe47 --> He4 + Cr43
    rate = 0.0

    # ths8r
    rate += np.exp(  75.2276 + -86.7594*tf.T9i + -86.4906*tf.T913i + -1.66143*tf.T913
                  + -0.451078*tf.T9 + 6.28872e-05*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Fe47__He4_Cr43 = rate

@numba.njit()
def p_Fe47__Co48(rate_eval, tf):
    # Fe47 + p --> Co48
    rate = 0.0

    # ths8r
    rate += np.exp(  -19.3572 + -9.377*tf.T9i + 20.4773*tf.T913
                  + -2.90303*tf.T9 + 0.254241*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Fe47__Co48 = rate

@numba.njit()
def p_Fe47__He4_Mn44(rate_eval, tf):
    # Fe47 + p --> He4 + Mn44
    rate = 0.0

    # ths8r
    rate += np.exp(  64.1798 + -101.185*tf.T9i + -88.9347*tf.T913i + -2.86881*tf.T913
                  + -0.327748*tf.T9 + 0.012567*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Fe47__He4_Mn44 = rate

@numba.njit()
def Co47__Fe47__weak__bqa_pos_(rate_eval, tf):
    # Co47 --> Fe47
    rate = 0.0

    # bqa+w
    rate += np.exp(  4.2091)

    rate_eval.Co47__Fe47__weak__bqa_pos_ = rate

@numba.njit()
def Co47__p_Fe46(rate_eval, tf):
    # Co47 --> p + Fe46
    rate = 0.0

    # ths8r
    rate += np.exp(  4.42497 + 0.000682203*tf.T9i + 20.309*tf.T913
                  + -2.64768*tf.T9 + 0.201647*tf.T953)

    rate_eval.Co47__p_Fe46 = rate

@numba.njit()
def p_Co47__Ni48(rate_eval, tf):
    # Co47 + p --> Ni48
    rate = 0.0

    # ths8r
    rate += np.exp(  -18.8912 + -14.2057*tf.T9i + 18.9366*tf.T913
                  + -2.28963*tf.T9 + 0.162055*tf.T953 + -1.5*tf.lnT9)

    rate_eval.p_Co47__Ni48 = rate

@numba.njit()
def Fe45__Mn45__weak__wc17(rate_eval, tf):
    # Fe45 --> Mn45
    rate = 0.0

    # wc17w
    rate += np.exp(  5.59269)

    rate_eval.Fe45__Mn45__weak__wc17 = rate

@numba.njit()
def Fe45__p_Mn44(rate_eval, tf):
    # Fe45 --> p + Mn44
    rate = 0.0

    # ths8r
    rate += np.exp(  64.498 + -1.26455*tf.T9i + -36.1457*tf.T913i + -6.20399*tf.T913
                  + -2.04399*tf.T9 + 0.328581*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Fe45__p_Mn44 = rate

@numba.njit()
def Fe45__p_Cr44__weak__wc17(rate_eval, tf):
    # Fe45 --> p + Cr44
    rate = 0.0

    # wc17w
    rate += np.exp(  4.24393)

    rate_eval.Fe45__p_Cr44__weak__wc17 = rate

@numba.njit()
def He4_Fe45__p_Co48(rate_eval, tf):
    # Fe45 + He4 --> p + Co48
    rate = 0.0

    # ths8r
    rate += np.exp(  62.1224 + -91.3476*tf.T913i + -1.40215*tf.T913
                  + -0.406715*tf.T9 + 0.00387095*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.He4_Fe45__p_Co48 = rate

@numba.njit()
def Co48__Fe48__weak__bqa_pos_(rate_eval, tf):
    # Co48 --> Fe48
    rate = 0.0

    # bqa+w
    rate += np.exp(  3.94599)

    rate_eval.Co48__Fe48__weak__bqa_pos_ = rate

@numba.njit()
def Co48__p_Fe47(rate_eval, tf):
    # Co48 --> p + Fe47
    rate = 0.0

    # ths8r
    rate += np.exp(  4.45043 + 0.000408275*tf.T9i + 20.4773*tf.T913
                  + -2.90303*tf.T9 + 0.254241*tf.T953)

    rate_eval.Co48__p_Fe47 = rate

@numba.njit()
def Co48__He4_Mn44(rate_eval, tf):
    # Co48 --> He4 + Mn44
    rate = 0.0

    # ths8r
    rate += np.exp(  77.0913 + -91.8074*tf.T9i + -88.9347*tf.T913i + -2.41015*tf.T913
                  + -0.458817*tf.T9 + 0.0164617*tf.T953 + 0.833333*tf.lnT9)

    rate_eval.Co48__He4_Mn44 = rate

@numba.njit()
def p_Co48__He4_Fe45(rate_eval, tf):
    # Co48 + p --> He4 + Fe45
    rate = 0.0

    # ths8r
    rate += np.exp(  62.8523 + -90.5428*tf.T9i + -91.3476*tf.T913i + -1.40215*tf.T913
                  + -0.406715*tf.T9 + 0.00387095*tf.T953 + -0.666667*tf.lnT9)

    rate_eval.p_Co48__He4_Fe45 = rate

@numba.njit()
def Ni48__Co48__weak__wc17(rate_eval, tf):
    # Ni48 --> Co48
    rate = 0.0

    # wc17w
    rate += np.exp(  5.7993)

    rate_eval.Ni48__Co48__weak__wc17 = rate

@numba.njit()
def Ni48__p_Co47(rate_eval, tf):
    # Ni48 --> p + Co47
    rate = 0.0

    # ths8r
    rate += np.exp(  5.47605 + 0.00056019*tf.T9i + 18.9366*tf.T913
                  + -2.28963*tf.T9 + 0.162055*tf.T953)

    rate_eval.Ni48__p_Co47 = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    Be7__Li7__weak__electron_capture(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    He3__p_d(rate_eval, tf)
    He4__d_d(rate_eval, tf)
    Be7__He4_He3(rate_eval, tf)
    B8__p_Be7(rate_eval, tf)
    B8__He4_He4__weak__wc12(rate_eval, tf)
    O16__p_N15(rate_eval, tf)
    O16__He4_C12(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    F19__He4_N15(rate_eval, tf)
    C12__He4_He4_He4(rate_eval, tf)
    p_p__d__weak__bet_pos_(rate_eval, tf)
    p_p__d__weak__electron_capture(rate_eval, tf)
    p_d__He3(rate_eval, tf)
    d_d__He4(rate_eval, tf)
    p_He3__He4__weak__bet_pos_(rate_eval, tf)
    He4_He3__Be7(rate_eval, tf)
    p_Be7__B8(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    d_He3__p_He4(rate_eval, tf)
    p_He4__d_He3(rate_eval, tf)
    He4_He4__p_Li7(rate_eval, tf)
    p_Li7__He4_He4(rate_eval, tf)
    He4_C12__p_N15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_N15__p_O18(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    He4_O17__p_F20(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    p_F20__He4_O17(rate_eval, tf)
    He3_He3__p_p_He4(rate_eval, tf)
    d_Be7__p_He4_He4(rate_eval, tf)
    He3_Be7__p_p_He4_He4(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)
    p_p_He4__He3_He3(rate_eval, tf)
    p_He4_He4__d_Be7(rate_eval, tf)
    p_p_He4_He4__He3_Be7(rate_eval, tf)
    p_O15__F16(rate_eval, tf)
    F16__p_O15(rate_eval, tf)
    p_S28__Cl29(rate_eval, tf)
    He4_S28__Ar32(rate_eval, tf)
    p_S29__Cl30(rate_eval, tf)
    He4_S29__Ar33(rate_eval, tf)
    Cl29__S29__weak__bqa_pos_(rate_eval, tf)
    Cl29__p_S28(rate_eval, tf)
    He4_Cl29__K33(rate_eval, tf)
    He4_Cl29__p_Ar32(rate_eval, tf)
    Ar32__He4_S28(rate_eval, tf)
    p_Ar32__K33(rate_eval, tf)
    He4_Ar32__Ca36(rate_eval, tf)
    p_Ar32__He4_Cl29(rate_eval, tf)
    p_Ca38__Sc39(rate_eval, tf)
    p_p_Ca38__He4_Ca36(rate_eval, tf)
    Cl30__p_S29(rate_eval, tf)
    He4_Cl30__K34(rate_eval, tf)
    He4_Cl30__p_Ar33(rate_eval, tf)
    Ar33__He4_S29(rate_eval, tf)
    p_Ar33__K34(rate_eval, tf)
    He4_Ar33__Ca37(rate_eval, tf)
    p_Ar33__He4_Cl30(rate_eval, tf)
    K33__Ar33__weak__bqa_pos_(rate_eval, tf)
    K33__p_Ar32(rate_eval, tf)
    K33__He4_Cl29(rate_eval, tf)
    He4_K33__Sc37(rate_eval, tf)
    He4_K33__p_Ca36(rate_eval, tf)
    Ca36__He4_Ar32(rate_eval, tf)
    p_Ca36__Sc37(rate_eval, tf)
    He4_Ca36__Ti40(rate_eval, tf)
    p_Ca36__He4_K33(rate_eval, tf)
    He4_Ca36__p_Sc39(rate_eval, tf)
    He4_Ca36__p_p_Ca38(rate_eval, tf)
    Sc39__p_Ca38(rate_eval, tf)
    p_Sc39__Ti40(rate_eval, tf)
    p_Sc39__He4_Ca36(rate_eval, tf)
    K34__p_Ar33(rate_eval, tf)
    K34__He4_Cl30(rate_eval, tf)
    p_K34__Ca35(rate_eval, tf)
    He4_K34__Sc38(rate_eval, tf)
    He4_K34__p_Ca37(rate_eval, tf)
    Ca37__He4_Ar33(rate_eval, tf)
    p_Ca37__Sc38(rate_eval, tf)
    p_Ca37__He4_K34(rate_eval, tf)
    Sc37__Ca37__weak__bqa_pos_(rate_eval, tf)
    Sc37__p_Ca36(rate_eval, tf)
    Sc37__He4_K33(rate_eval, tf)
    He4_Sc37__V41(rate_eval, tf)
    He4_Sc37__p_Ti40(rate_eval, tf)
    Ti40__p_Sc39(rate_eval, tf)
    Ti40__He4_Ca36(rate_eval, tf)
    p_Ti40__V41(rate_eval, tf)
    He4_Ti40__Cr44(rate_eval, tf)
    p_Ti40__He4_Sc37(rate_eval, tf)
    Ca35__p_K34(rate_eval, tf)
    p_Ca35__Sc36(rate_eval, tf)
    He4_Ca35__Ti39(rate_eval, tf)
    He4_Ca35__p_Sc38(rate_eval, tf)
    Sc38__Ca38__weak__mo97(rate_eval, tf)
    Sc38__p_Ca37(rate_eval, tf)
    Sc38__He4_K34(rate_eval, tf)
    p_Sc38__Ti39(rate_eval, tf)
    p_Sc38__He4_Ca35(rate_eval, tf)
    V41__p_Ti40(rate_eval, tf)
    V41__He4_Sc37(rate_eval, tf)
    He4_V41__Mn45(rate_eval, tf)
    He4_V41__p_Cr44(rate_eval, tf)
    Cr44__He4_Ti40(rate_eval, tf)
    p_Cr44__Mn45(rate_eval, tf)
    He4_Cr44__Fe48(rate_eval, tf)
    p_Cr44__He4_V41(rate_eval, tf)
    Sc36__Ca36__weak__bqa_pos_(rate_eval, tf)
    Sc36__p_Ca35(rate_eval, tf)
    He4_Sc36__V40(rate_eval, tf)
    He4_Sc36__p_Ti39(rate_eval, tf)
    Ti39__Sc39__weak__wc17(rate_eval, tf)
    Ti39__p_Sc38(rate_eval, tf)
    Ti39__p_Ca38__weak__wc12(rate_eval, tf)
    Ti39__He4_Ca35(rate_eval, tf)
    p_Ti39__V40(rate_eval, tf)
    He4_Ti39__Cr43(rate_eval, tf)
    p_Ti39__He4_Sc36(rate_eval, tf)
    Mn45__p_Cr44(rate_eval, tf)
    Mn45__He4_V41(rate_eval, tf)
    p_Mn45__Fe46(rate_eval, tf)
    He4_Mn45__Co49(rate_eval, tf)
    He4_Mn45__p_Fe48(rate_eval, tf)
    Fe48__He4_Cr44(rate_eval, tf)
    p_Fe48__Co49(rate_eval, tf)
    p_Fe48__He4_Mn45(rate_eval, tf)
    V40__Ti40__weak__bqa_pos_(rate_eval, tf)
    V40__p_Ti39(rate_eval, tf)
    V40__He4_Sc36(rate_eval, tf)
    He4_V40__Mn44(rate_eval, tf)
    He4_V40__p_Cr43(rate_eval, tf)
    Cr43__He4_Ti39(rate_eval, tf)
    p_Cr43__Mn44(rate_eval, tf)
    He4_Cr43__Fe47(rate_eval, tf)
    p_Cr43__He4_V40(rate_eval, tf)
    Fe46__p_Mn45(rate_eval, tf)
    p_Fe46__Co47(rate_eval, tf)
    He4_Fe46__p_Co49(rate_eval, tf)
    Co49__p_Fe48(rate_eval, tf)
    Co49__He4_Mn45(rate_eval, tf)
    p_Co49__He4_Fe46(rate_eval, tf)
    Mn44__Cr44__weak__bqa_pos_(rate_eval, tf)
    Mn44__p_Cr43(rate_eval, tf)
    Mn44__He4_V40(rate_eval, tf)
    p_Mn44__Fe45(rate_eval, tf)
    He4_Mn44__Co48(rate_eval, tf)
    He4_Mn44__p_Fe47(rate_eval, tf)
    Fe47__He4_Cr43(rate_eval, tf)
    p_Fe47__Co48(rate_eval, tf)
    p_Fe47__He4_Mn44(rate_eval, tf)
    Co47__Fe47__weak__bqa_pos_(rate_eval, tf)
    Co47__p_Fe46(rate_eval, tf)
    p_Co47__Ni48(rate_eval, tf)
    Fe45__Mn45__weak__wc17(rate_eval, tf)
    Fe45__p_Mn44(rate_eval, tf)
    Fe45__p_Cr44__weak__wc17(rate_eval, tf)
    He4_Fe45__p_Co48(rate_eval, tf)
    Co48__Fe48__weak__bqa_pos_(rate_eval, tf)
    Co48__p_Fe47(rate_eval, tf)
    Co48__He4_Mn44(rate_eval, tf)
    p_Co48__He4_Fe45(rate_eval, tf)
    Ni48__Co48__weak__wc17(rate_eval, tf)
    Ni48__p_Co47(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p__d__weak__bet_pos_ *= scor
        rate_eval.p_p__d__weak__electron_capture *= scor
        rate_eval.p_p_He4_He4__He3_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_d__He3 *= scor

        scn_fac = ScreenFactors(1, 2, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_d__He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He3__He4__weak__bet_pos_ *= scor

        scn_fac = ScreenFactors(2, 4, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He3__Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be7__B8 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor
        rate_eval.He4_C12__p_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor
        rate_eval.He4_N15__p_O18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He3__p_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4__d_He3 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He4__p_Li7 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li7__He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__p_F20 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F20__He4_O17 *= scor

        scn_fac = ScreenFactors(2, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_He3__p_p_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Be7__p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Be7__p_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_He4__He3_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4_He4__d_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O15__F16 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S28__Cl29 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S28__Ar32 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S29__Cl30 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S29__Ar33 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl29__K33 *= scor
        rate_eval.He4_Cl29__p_Ar32 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar32__K33 *= scor
        rate_eval.p_Ar32__He4_Cl29 *= scor

        scn_fac = ScreenFactors(2, 4, 18, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ar32__Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 38)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca38__Sc39 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_Ca38__He4_Ca36 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl30__K34 *= scor
        rate_eval.He4_Cl30__p_Ar33 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar33__K34 *= scor
        rate_eval.p_Ar33__He4_Cl30 *= scor

        scn_fac = ScreenFactors(2, 4, 18, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ar33__Ca37 *= scor

        scn_fac = ScreenFactors(2, 4, 19, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_K33__Sc37 *= scor
        rate_eval.He4_K33__p_Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca36__Sc37 *= scor
        rate_eval.p_Ca36__He4_K33 *= scor

        scn_fac = ScreenFactors(2, 4, 20, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ca36__Ti40 *= scor
        rate_eval.He4_Ca36__p_Sc39 *= scor
        rate_eval.He4_Ca36__p_p_Ca38 *= scor

        scn_fac = ScreenFactors(1, 1, 21, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Sc39__Ti40 *= scor
        rate_eval.p_Sc39__He4_Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 19, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_K34__Ca35 *= scor

        scn_fac = ScreenFactors(2, 4, 19, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_K34__Sc38 *= scor
        rate_eval.He4_K34__p_Ca37 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 37)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca37__Sc38 *= scor
        rate_eval.p_Ca37__He4_K34 *= scor

        scn_fac = ScreenFactors(2, 4, 21, 37)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Sc37__V41 *= scor
        rate_eval.He4_Sc37__p_Ti40 *= scor

        scn_fac = ScreenFactors(1, 1, 22, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ti40__V41 *= scor
        rate_eval.p_Ti40__He4_Sc37 *= scor

        scn_fac = ScreenFactors(2, 4, 22, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ti40__Cr44 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca35__Sc36 *= scor

        scn_fac = ScreenFactors(2, 4, 20, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ca35__Ti39 *= scor
        rate_eval.He4_Ca35__p_Sc38 *= scor

        scn_fac = ScreenFactors(1, 1, 21, 38)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Sc38__Ti39 *= scor
        rate_eval.p_Sc38__He4_Ca35 *= scor

        scn_fac = ScreenFactors(2, 4, 23, 41)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_V41__Mn45 *= scor
        rate_eval.He4_V41__p_Cr44 *= scor

        scn_fac = ScreenFactors(1, 1, 24, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cr44__Mn45 *= scor
        rate_eval.p_Cr44__He4_V41 *= scor

        scn_fac = ScreenFactors(2, 4, 24, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cr44__Fe48 *= scor

        scn_fac = ScreenFactors(2, 4, 21, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Sc36__V40 *= scor
        rate_eval.He4_Sc36__p_Ti39 *= scor

        scn_fac = ScreenFactors(1, 1, 22, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ti39__V40 *= scor
        rate_eval.p_Ti39__He4_Sc36 *= scor

        scn_fac = ScreenFactors(2, 4, 22, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ti39__Cr43 *= scor

        scn_fac = ScreenFactors(1, 1, 25, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mn45__Fe46 *= scor

        scn_fac = ScreenFactors(2, 4, 25, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mn45__Co49 *= scor
        rate_eval.He4_Mn45__p_Fe48 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 48)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe48__Co49 *= scor
        rate_eval.p_Fe48__He4_Mn45 *= scor

        scn_fac = ScreenFactors(2, 4, 23, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_V40__Mn44 *= scor
        rate_eval.He4_V40__p_Cr43 *= scor

        scn_fac = ScreenFactors(1, 1, 24, 43)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cr43__Mn44 *= scor
        rate_eval.p_Cr43__He4_V40 *= scor

        scn_fac = ScreenFactors(2, 4, 24, 43)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cr43__Fe47 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 46)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe46__Co47 *= scor

        scn_fac = ScreenFactors(2, 4, 26, 46)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Fe46__p_Co49 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 49)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co49__He4_Fe46 *= scor

        scn_fac = ScreenFactors(1, 1, 25, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mn44__Fe45 *= scor

        scn_fac = ScreenFactors(2, 4, 25, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mn44__Co48 *= scor
        rate_eval.He4_Mn44__p_Fe47 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 47)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe47__Co48 *= scor
        rate_eval.p_Fe47__He4_Mn44 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 47)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co47__Ni48 *= scor

        scn_fac = ScreenFactors(2, 4, 26, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Fe45__p_Co48 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 48)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co48__He4_Fe45 *= scor

    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jp] = (
       -2*5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p__d__weak__bet_pos_
       -2*5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p__d__weak__electron_capture
       -rho*Y[jp]*Y[jd]*rate_eval.p_d__He3
       -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3
       -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jp]*Y[jf20]*rate_eval.p_F20__He4_O17
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16
       -rho*Y[jp]*Y[js28]*rate_eval.p_S28__Cl29
       -rho*Y[jp]*Y[js29]*rate_eval.p_S29__Cl30
       -rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__K33
       -rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       -rho*Y[jp]*Y[jca38]*rate_eval.p_Ca38__Sc39
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       -rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__K34
       -rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       -rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__Sc37
       -rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__He4_K33
       -rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__Ti40
       -rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       -rho*Y[jp]*Y[jk34]*rate_eval.p_K34__Ca35
       -rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__Sc38
       -rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__He4_K34
       -rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__V41
       -rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       -rho*Y[jp]*Y[jca35]*rate_eval.p_Ca35__Sc36
       -rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__Ti39
       -rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       -rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__Mn45
       -rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       -rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__V40
       -rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       -rho*Y[jp]*Y[jmn45]*rate_eval.p_Mn45__Fe46
       -rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__Co49
       -rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       -rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__Mn44
       -rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       -rho*Y[jp]*Y[jfe46]*rate_eval.p_Fe46__Co47
       -rho*Y[jp]*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       -rho*Y[jp]*Y[jmn44]*rate_eval.p_Mn44__Fe45
       -rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__Co48
       -rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       -rho*Y[jp]*Y[jco47]*rate_eval.p_Co47__Ni48
       -rho*Y[jp]*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       +Y[jhe3]*rate_eval.He3__p_d
       +Y[jb8]*rate_eval.B8__p_Be7
       +Y[jo16]*rate_eval.O16__p_N15
       +Y[jf19]*rate_eval.F19__p_O18
       +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4
       +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7
       +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15
       +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__p_F20
       +2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4
       +rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       +2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +Y[jf16]*rate_eval.F16__p_O15
       +Y[jcl29]*rate_eval.Cl29__p_S28
       +rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       +Y[jcl30]*rate_eval.Cl30__p_S29
       +rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       +Y[jk33]*rate_eval.K33__p_Ar32
       +rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__p_Ca36
       +rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       +2*rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       +Y[jsc39]*rate_eval.Sc39__p_Ca38
       +Y[jk34]*rate_eval.K34__p_Ar33
       +rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__p_Ca37
       +Y[jsc37]*rate_eval.Sc37__p_Ca36
       +rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       +Y[jti40]*rate_eval.Ti40__p_Sc39
       +Y[jca35]*rate_eval.Ca35__p_K34
       +rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       +Y[jsc38]*rate_eval.Sc38__p_Ca37
       +Y[jv41]*rate_eval.V41__p_Ti40
       +rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__p_Cr44
       +Y[jsc36]*rate_eval.Sc36__p_Ca35
       +rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       +Y[jti39]*rate_eval.Ti39__p_Sc38
       +Y[jti39]*rate_eval.Ti39__p_Ca38__weak__wc12
       +Y[jmn45]*rate_eval.Mn45__p_Cr44
       +rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       +Y[jv40]*rate_eval.V40__p_Ti39
       +rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__p_Cr43
       +Y[jfe46]*rate_eval.Fe46__p_Mn45
       +rho*Y[jhe4]*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       +Y[jco49]*rate_eval.Co49__p_Fe48
       +Y[jmn44]*rate_eval.Mn44__p_Cr43
       +rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       +Y[jco47]*rate_eval.Co47__p_Fe46
       +Y[jfe45]*rate_eval.Fe45__p_Mn44
       +Y[jfe45]*rate_eval.Fe45__p_Cr44__weak__wc17
       +rho*Y[jhe4]*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       +Y[jco48]*rate_eval.Co48__p_Fe47
       +Y[jni48]*rate_eval.Ni48__p_Co47
       )

    dYdt[jd] = (
       -rho*Y[jp]*Y[jd]*rate_eval.p_d__He3
       -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__He4
       -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4
       -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       +Y[jhe3]*rate_eval.He3__p_d
       +2*Y[jhe4]*rate_eval.He4__d_d
       +5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p__d__weak__bet_pos_
       +5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p__d__weak__electron_capture
       +rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3
       +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       )

    dYdt[jhe3] = (
       -Y[jhe3]*rate_eval.He3__p_d
       -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7
       -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4
       -2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4
       -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +Y[jbe7]*rate_eval.Be7__He4_He3
       +rho*Y[jp]*Y[jd]*rate_eval.p_d__He3
       +rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3
       +2*5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       +2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    dYdt[jhe4] = (
       -Y[jhe4]*rate_eval.He4__d_d
       -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jp]*Y[jhe4]*rate_eval.p_He4__d_He3
       -2*5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__p_F20
       -3*1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12
       -5.00000000000000e-01*rho**2*Y[jp]**2*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -2*5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       -rho*Y[jhe4]*Y[js28]*rate_eval.He4_S28__Ar32
       -rho*Y[jhe4]*Y[js29]*rate_eval.He4_S29__Ar33
       -rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__K33
       -rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       -rho*Y[jhe4]*Y[jar32]*rate_eval.He4_Ar32__Ca36
       -rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__K34
       -rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       -rho*Y[jhe4]*Y[jar33]*rate_eval.He4_Ar33__Ca37
       -rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__Sc37
       -rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__p_Ca36
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       -rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__Sc38
       -rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__p_Ca37
       -rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__V41
       -rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       -rho*Y[jhe4]*Y[jti40]*rate_eval.He4_Ti40__Cr44
       -rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       -rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__Mn45
       -rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__p_Cr44
       -rho*Y[jhe4]*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       -rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__V40
       -rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       -rho*Y[jhe4]*Y[jti39]*rate_eval.He4_Ti39__Cr43
       -rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__Co49
       -rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       -rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__Mn44
       -rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__p_Cr43
       -rho*Y[jhe4]*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       -rho*Y[jhe4]*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       -rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__Co48
       -rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       -rho*Y[jhe4]*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       +Y[jbe7]*rate_eval.Be7__He4_He3
       +2*Y[jb8]*rate_eval.B8__He4_He4__weak__wc12
       +Y[jo16]*rate_eval.O16__He4_C12
       +Y[jf19]*rate_eval.F19__He4_N15
       +3*Y[jc12]*rate_eval.C12__He4_He4_He4
       +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d__He4
       +rho*Y[jp]*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3__p_He4
       +2*rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jp]*Y[jf20]*rate_eval.p_F20__He4_O17
       +5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3__p_p_He4
       +2*rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       +2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +Y[jar32]*rate_eval.Ar32__He4_S28
       +rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       +5.00000000000000e-01*rho**2*Y[jp]**2*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       +Y[jar33]*rate_eval.Ar33__He4_S29
       +rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       +Y[jk33]*rate_eval.K33__He4_Cl29
       +Y[jca36]*rate_eval.Ca36__He4_Ar32
       +rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__He4_K33
       +rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       +Y[jk34]*rate_eval.K34__He4_Cl30
       +Y[jca37]*rate_eval.Ca37__He4_Ar33
       +rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__He4_K34
       +Y[jsc37]*rate_eval.Sc37__He4_K33
       +Y[jti40]*rate_eval.Ti40__He4_Ca36
       +rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       +Y[jsc38]*rate_eval.Sc38__He4_K34
       +rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       +Y[jv41]*rate_eval.V41__He4_Sc37
       +Y[jcr44]*rate_eval.Cr44__He4_Ti40
       +rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       +Y[jti39]*rate_eval.Ti39__He4_Ca35
       +rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       +Y[jmn45]*rate_eval.Mn45__He4_V41
       +Y[jfe48]*rate_eval.Fe48__He4_Cr44
       +rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       +Y[jv40]*rate_eval.V40__He4_Sc36
       +Y[jcr43]*rate_eval.Cr43__He4_Ti39
       +rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       +Y[jco49]*rate_eval.Co49__He4_Mn45
       +rho*Y[jp]*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       +Y[jmn44]*rate_eval.Mn44__He4_V40
       +Y[jfe47]*rate_eval.Fe47__He4_Cr43
       +rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       +Y[jco48]*rate_eval.Co48__He4_Mn44
       +rho*Y[jp]*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       )

    dYdt[jli7] = (
       -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7__He4_He4
       +rho*ye(Y)*Y[jbe7]*rate_eval.Be7__Li7__weak__electron_capture
       +5.00000000000000e-01*rho*Y[jhe4]**2*rate_eval.He4_He4__p_Li7
       )

    dYdt[jbe7] = (
       -rho*ye(Y)*Y[jbe7]*rate_eval.Be7__Li7__weak__electron_capture
       -Y[jbe7]*rate_eval.Be7__He4_He3
       -rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8
       -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +Y[jb8]*rate_eval.B8__p_Be7
       +rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3__Be7
       +5.00000000000000e-01*rho**2*Y[jp]*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       +2.50000000000000e-01*rho**3*Y[jp]**2*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    dYdt[jb8] = (
       -Y[jb8]*rate_eval.B8__p_Be7
       -Y[jb8]*rate_eval.B8__He4_He4__weak__wc12
       +rho*Y[jp]*Y[jbe7]*rate_eval.p_Be7__B8
       )

    dYdt[jc12] = (
       -Y[jc12]*rate_eval.C12__He4_He4_He4
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15
       +Y[jo16]*rate_eval.O16__He4_C12
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       +1.66666666666667e-01*rho**2*Y[jhe4]**3*rate_eval.He4_He4_He4__C12
       )

    dYdt[jn15] = (
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jp]*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18
       +Y[jo15]*rate_eval.O15__N15__weak__wc12
       +Y[jo16]*rate_eval.O16__p_N15
       +Y[jf19]*rate_eval.F19__He4_N15
       +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__p_N15
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    dYdt[jo15] = (
       -Y[jo15]*rate_eval.O15__N15__weak__wc12
       -rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16
       +Y[jf16]*rate_eval.F16__p_O15
       )

    dYdt[jo16] = (
       -Y[jo16]*rate_eval.O16__p_N15
       -Y[jo16]*rate_eval.O16__He4_C12
       -rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jhe4]*Y[jc12]*rate_eval.He4_C12__O16
       +rho*Y[jp]*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    dYdt[jo17] = (
       -rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__p_F20
       +rho*Y[jp]*Y[jf20]*rate_eval.p_F20__He4_O17
       )

    dYdt[jo18] = (
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jp]*Y[jo18]*rate_eval.p_O18__He4_N15
       +Y[jf19]*rate_eval.F19__p_O18
       +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__p_O18
       )

    dYdt[jf16] = (
       -Y[jf16]*rate_eval.F16__p_O15
       +rho*Y[jp]*Y[jo15]*rate_eval.p_O15__F16
       )

    dYdt[jf19] = (
       -Y[jf19]*rate_eval.F19__p_O18
       -Y[jf19]*rate_eval.F19__He4_N15
       -rho*Y[jp]*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jhe4]*Y[jn15]*rate_eval.He4_N15__F19
       +rho*Y[jp]*Y[jo18]*rate_eval.p_O18__F19
       +rho*Y[jhe4]*Y[jo16]*rate_eval.He4_O16__p_F19
       )

    dYdt[jf20] = (
       -rho*Y[jp]*Y[jf20]*rate_eval.p_F20__He4_O17
       +rho*Y[jhe4]*Y[jo17]*rate_eval.He4_O17__p_F20
       )

    dYdt[js28] = (
       -rho*Y[jp]*Y[js28]*rate_eval.p_S28__Cl29
       -rho*Y[jhe4]*Y[js28]*rate_eval.He4_S28__Ar32
       +Y[jcl29]*rate_eval.Cl29__p_S28
       +Y[jar32]*rate_eval.Ar32__He4_S28
       )

    dYdt[js29] = (
       -rho*Y[jp]*Y[js29]*rate_eval.p_S29__Cl30
       -rho*Y[jhe4]*Y[js29]*rate_eval.He4_S29__Ar33
       +Y[jcl29]*rate_eval.Cl29__S29__weak__bqa_pos_
       +Y[jcl30]*rate_eval.Cl30__p_S29
       +Y[jar33]*rate_eval.Ar33__He4_S29
       )

    dYdt[jcl29] = (
       -Y[jcl29]*rate_eval.Cl29__S29__weak__bqa_pos_
       -Y[jcl29]*rate_eval.Cl29__p_S28
       -rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__K33
       -rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       +rho*Y[jp]*Y[js28]*rate_eval.p_S28__Cl29
       +rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       +Y[jk33]*rate_eval.K33__He4_Cl29
       )

    dYdt[jcl30] = (
       -Y[jcl30]*rate_eval.Cl30__p_S29
       -rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__K34
       -rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       +rho*Y[jp]*Y[js29]*rate_eval.p_S29__Cl30
       +rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       +Y[jk34]*rate_eval.K34__He4_Cl30
       )

    dYdt[jar32] = (
       -Y[jar32]*rate_eval.Ar32__He4_S28
       -rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__K33
       -rho*Y[jhe4]*Y[jar32]*rate_eval.He4_Ar32__Ca36
       -rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       +rho*Y[jhe4]*Y[js28]*rate_eval.He4_S28__Ar32
       +rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       +Y[jk33]*rate_eval.K33__p_Ar32
       +Y[jca36]*rate_eval.Ca36__He4_Ar32
       )

    dYdt[jar33] = (
       -Y[jar33]*rate_eval.Ar33__He4_S29
       -rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__K34
       -rho*Y[jhe4]*Y[jar33]*rate_eval.He4_Ar33__Ca37
       -rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       +rho*Y[jhe4]*Y[js29]*rate_eval.He4_S29__Ar33
       +rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       +Y[jk33]*rate_eval.K33__Ar33__weak__bqa_pos_
       +Y[jk34]*rate_eval.K34__p_Ar33
       +Y[jca37]*rate_eval.Ca37__He4_Ar33
       )

    dYdt[jk33] = (
       -Y[jk33]*rate_eval.K33__Ar33__weak__bqa_pos_
       -Y[jk33]*rate_eval.K33__p_Ar32
       -Y[jk33]*rate_eval.K33__He4_Cl29
       -rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__Sc37
       -rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__p_Ca36
       +rho*Y[jhe4]*Y[jcl29]*rate_eval.He4_Cl29__K33
       +rho*Y[jp]*Y[jar32]*rate_eval.p_Ar32__K33
       +rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__He4_K33
       +Y[jsc37]*rate_eval.Sc37__He4_K33
       )

    dYdt[jk34] = (
       -Y[jk34]*rate_eval.K34__p_Ar33
       -Y[jk34]*rate_eval.K34__He4_Cl30
       -rho*Y[jp]*Y[jk34]*rate_eval.p_K34__Ca35
       -rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__Sc38
       -rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__p_Ca37
       +rho*Y[jhe4]*Y[jcl30]*rate_eval.He4_Cl30__K34
       +rho*Y[jp]*Y[jar33]*rate_eval.p_Ar33__K34
       +rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__He4_K34
       +Y[jca35]*rate_eval.Ca35__p_K34
       +Y[jsc38]*rate_eval.Sc38__He4_K34
       )

    dYdt[jca35] = (
       -Y[jca35]*rate_eval.Ca35__p_K34
       -rho*Y[jp]*Y[jca35]*rate_eval.p_Ca35__Sc36
       -rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       +rho*Y[jp]*Y[jk34]*rate_eval.p_K34__Ca35
       +rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       +Y[jsc36]*rate_eval.Sc36__p_Ca35
       +Y[jti39]*rate_eval.Ti39__He4_Ca35
       )

    dYdt[jca36] = (
       -Y[jca36]*rate_eval.Ca36__He4_Ar32
       -rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__Sc37
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__He4_K33
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       +rho*Y[jhe4]*Y[jar32]*rate_eval.He4_Ar32__Ca36
       +5.00000000000000e-01*rho**2*Y[jp]**2*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       +rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__p_Ca36
       +rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       +Y[jsc37]*rate_eval.Sc37__p_Ca36
       +Y[jti40]*rate_eval.Ti40__He4_Ca36
       +Y[jsc36]*rate_eval.Sc36__Ca36__weak__bqa_pos_
       )

    dYdt[jca37] = (
       -Y[jca37]*rate_eval.Ca37__He4_Ar33
       -rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__Sc38
       -rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__He4_K34
       +rho*Y[jhe4]*Y[jar33]*rate_eval.He4_Ar33__Ca37
       +rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__p_Ca37
       +Y[jsc37]*rate_eval.Sc37__Ca37__weak__bqa_pos_
       +Y[jsc38]*rate_eval.Sc38__p_Ca37
       )

    dYdt[jca38] = (
       -rho*Y[jp]*Y[jca38]*rate_eval.p_Ca38__Sc39
       -5.00000000000000e-01*rho**2*Y[jp]**2*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       +rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       +Y[jsc39]*rate_eval.Sc39__p_Ca38
       +Y[jsc38]*rate_eval.Sc38__Ca38__weak__mo97
       +Y[jti39]*rate_eval.Ti39__p_Ca38__weak__wc12
       )

    dYdt[jsc36] = (
       -Y[jsc36]*rate_eval.Sc36__Ca36__weak__bqa_pos_
       -Y[jsc36]*rate_eval.Sc36__p_Ca35
       -rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__V40
       -rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       +rho*Y[jp]*Y[jca35]*rate_eval.p_Ca35__Sc36
       +rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       +Y[jv40]*rate_eval.V40__He4_Sc36
       )

    dYdt[jsc37] = (
       -Y[jsc37]*rate_eval.Sc37__Ca37__weak__bqa_pos_
       -Y[jsc37]*rate_eval.Sc37__p_Ca36
       -Y[jsc37]*rate_eval.Sc37__He4_K33
       -rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__V41
       -rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       +rho*Y[jhe4]*Y[jk33]*rate_eval.He4_K33__Sc37
       +rho*Y[jp]*Y[jca36]*rate_eval.p_Ca36__Sc37
       +rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       +Y[jv41]*rate_eval.V41__He4_Sc37
       )

    dYdt[jsc38] = (
       -Y[jsc38]*rate_eval.Sc38__Ca38__weak__mo97
       -Y[jsc38]*rate_eval.Sc38__p_Ca37
       -Y[jsc38]*rate_eval.Sc38__He4_K34
       -rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__Ti39
       -rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       +rho*Y[jhe4]*Y[jk34]*rate_eval.He4_K34__Sc38
       +rho*Y[jp]*Y[jca37]*rate_eval.p_Ca37__Sc38
       +rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       +Y[jti39]*rate_eval.Ti39__p_Sc38
       )

    dYdt[jsc39] = (
       -Y[jsc39]*rate_eval.Sc39__p_Ca38
       -rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__Ti40
       -rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       +rho*Y[jp]*Y[jca38]*rate_eval.p_Ca38__Sc39
       +rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       +Y[jti40]*rate_eval.Ti40__p_Sc39
       +Y[jti39]*rate_eval.Ti39__Sc39__weak__wc17
       )

    dYdt[jti39] = (
       -Y[jti39]*rate_eval.Ti39__Sc39__weak__wc17
       -Y[jti39]*rate_eval.Ti39__p_Sc38
       -Y[jti39]*rate_eval.Ti39__p_Ca38__weak__wc12
       -Y[jti39]*rate_eval.Ti39__He4_Ca35
       -rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__V40
       -rho*Y[jhe4]*Y[jti39]*rate_eval.He4_Ti39__Cr43
       -rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       +rho*Y[jhe4]*Y[jca35]*rate_eval.He4_Ca35__Ti39
       +rho*Y[jp]*Y[jsc38]*rate_eval.p_Sc38__Ti39
       +rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       +Y[jv40]*rate_eval.V40__p_Ti39
       +Y[jcr43]*rate_eval.Cr43__He4_Ti39
       )

    dYdt[jti40] = (
       -Y[jti40]*rate_eval.Ti40__p_Sc39
       -Y[jti40]*rate_eval.Ti40__He4_Ca36
       -rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__V41
       -rho*Y[jhe4]*Y[jti40]*rate_eval.He4_Ti40__Cr44
       -rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       +rho*Y[jhe4]*Y[jca36]*rate_eval.He4_Ca36__Ti40
       +rho*Y[jp]*Y[jsc39]*rate_eval.p_Sc39__Ti40
       +rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       +Y[jv41]*rate_eval.V41__p_Ti40
       +Y[jcr44]*rate_eval.Cr44__He4_Ti40
       +Y[jv40]*rate_eval.V40__Ti40__weak__bqa_pos_
       )

    dYdt[jv40] = (
       -Y[jv40]*rate_eval.V40__Ti40__weak__bqa_pos_
       -Y[jv40]*rate_eval.V40__p_Ti39
       -Y[jv40]*rate_eval.V40__He4_Sc36
       -rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__Mn44
       -rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__p_Cr43
       +rho*Y[jhe4]*Y[jsc36]*rate_eval.He4_Sc36__V40
       +rho*Y[jp]*Y[jti39]*rate_eval.p_Ti39__V40
       +rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       +Y[jmn44]*rate_eval.Mn44__He4_V40
       )

    dYdt[jv41] = (
       -Y[jv41]*rate_eval.V41__p_Ti40
       -Y[jv41]*rate_eval.V41__He4_Sc37
       -rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__Mn45
       -rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__p_Cr44
       +rho*Y[jhe4]*Y[jsc37]*rate_eval.He4_Sc37__V41
       +rho*Y[jp]*Y[jti40]*rate_eval.p_Ti40__V41
       +rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       +Y[jmn45]*rate_eval.Mn45__He4_V41
       )

    dYdt[jcr43] = (
       -Y[jcr43]*rate_eval.Cr43__He4_Ti39
       -rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__Mn44
       -rho*Y[jhe4]*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       -rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       +rho*Y[jhe4]*Y[jti39]*rate_eval.He4_Ti39__Cr43
       +rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__p_Cr43
       +Y[jmn44]*rate_eval.Mn44__p_Cr43
       +Y[jfe47]*rate_eval.Fe47__He4_Cr43
       )

    dYdt[jcr44] = (
       -Y[jcr44]*rate_eval.Cr44__He4_Ti40
       -rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__Mn45
       -rho*Y[jhe4]*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       -rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       +rho*Y[jhe4]*Y[jti40]*rate_eval.He4_Ti40__Cr44
       +rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__p_Cr44
       +Y[jmn45]*rate_eval.Mn45__p_Cr44
       +Y[jfe48]*rate_eval.Fe48__He4_Cr44
       +Y[jmn44]*rate_eval.Mn44__Cr44__weak__bqa_pos_
       +Y[jfe45]*rate_eval.Fe45__p_Cr44__weak__wc17
       )

    dYdt[jmn44] = (
       -Y[jmn44]*rate_eval.Mn44__Cr44__weak__bqa_pos_
       -Y[jmn44]*rate_eval.Mn44__p_Cr43
       -Y[jmn44]*rate_eval.Mn44__He4_V40
       -rho*Y[jp]*Y[jmn44]*rate_eval.p_Mn44__Fe45
       -rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__Co48
       -rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       +rho*Y[jhe4]*Y[jv40]*rate_eval.He4_V40__Mn44
       +rho*Y[jp]*Y[jcr43]*rate_eval.p_Cr43__Mn44
       +rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       +Y[jfe45]*rate_eval.Fe45__p_Mn44
       +Y[jco48]*rate_eval.Co48__He4_Mn44
       )

    dYdt[jmn45] = (
       -Y[jmn45]*rate_eval.Mn45__p_Cr44
       -Y[jmn45]*rate_eval.Mn45__He4_V41
       -rho*Y[jp]*Y[jmn45]*rate_eval.p_Mn45__Fe46
       -rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__Co49
       -rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       +rho*Y[jhe4]*Y[jv41]*rate_eval.He4_V41__Mn45
       +rho*Y[jp]*Y[jcr44]*rate_eval.p_Cr44__Mn45
       +rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       +Y[jfe46]*rate_eval.Fe46__p_Mn45
       +Y[jco49]*rate_eval.Co49__He4_Mn45
       +Y[jfe45]*rate_eval.Fe45__Mn45__weak__wc17
       )

    dYdt[jfe45] = (
       -Y[jfe45]*rate_eval.Fe45__Mn45__weak__wc17
       -Y[jfe45]*rate_eval.Fe45__p_Mn44
       -Y[jfe45]*rate_eval.Fe45__p_Cr44__weak__wc17
       -rho*Y[jhe4]*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       +rho*Y[jp]*Y[jmn44]*rate_eval.p_Mn44__Fe45
       +rho*Y[jp]*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       )

    dYdt[jfe46] = (
       -Y[jfe46]*rate_eval.Fe46__p_Mn45
       -rho*Y[jp]*Y[jfe46]*rate_eval.p_Fe46__Co47
       -rho*Y[jhe4]*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       +rho*Y[jp]*Y[jmn45]*rate_eval.p_Mn45__Fe46
       +rho*Y[jp]*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       +Y[jco47]*rate_eval.Co47__p_Fe46
       )

    dYdt[jfe47] = (
       -Y[jfe47]*rate_eval.Fe47__He4_Cr43
       -rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__Co48
       -rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       +rho*Y[jhe4]*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       +rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       +Y[jco47]*rate_eval.Co47__Fe47__weak__bqa_pos_
       +Y[jco48]*rate_eval.Co48__p_Fe47
       )

    dYdt[jfe48] = (
       -Y[jfe48]*rate_eval.Fe48__He4_Cr44
       -rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__Co49
       -rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       +rho*Y[jhe4]*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       +rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       +Y[jco49]*rate_eval.Co49__p_Fe48
       +Y[jco48]*rate_eval.Co48__Fe48__weak__bqa_pos_
       )

    dYdt[jco47] = (
       -Y[jco47]*rate_eval.Co47__Fe47__weak__bqa_pos_
       -Y[jco47]*rate_eval.Co47__p_Fe46
       -rho*Y[jp]*Y[jco47]*rate_eval.p_Co47__Ni48
       +rho*Y[jp]*Y[jfe46]*rate_eval.p_Fe46__Co47
       +Y[jni48]*rate_eval.Ni48__p_Co47
       )

    dYdt[jco48] = (
       -Y[jco48]*rate_eval.Co48__Fe48__weak__bqa_pos_
       -Y[jco48]*rate_eval.Co48__p_Fe47
       -Y[jco48]*rate_eval.Co48__He4_Mn44
       -rho*Y[jp]*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       +rho*Y[jhe4]*Y[jmn44]*rate_eval.He4_Mn44__Co48
       +rho*Y[jp]*Y[jfe47]*rate_eval.p_Fe47__Co48
       +rho*Y[jhe4]*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       +Y[jni48]*rate_eval.Ni48__Co48__weak__wc17
       )

    dYdt[jco49] = (
       -Y[jco49]*rate_eval.Co49__p_Fe48
       -Y[jco49]*rate_eval.Co49__He4_Mn45
       -rho*Y[jp]*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       +rho*Y[jhe4]*Y[jmn45]*rate_eval.He4_Mn45__Co49
       +rho*Y[jp]*Y[jfe48]*rate_eval.p_Fe48__Co49
       +rho*Y[jhe4]*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       )

    dYdt[jni48] = (
       -Y[jni48]*rate_eval.Ni48__Co48__weak__wc17
       -Y[jni48]*rate_eval.Ni48__p_Co47
       +rho*Y[jp]*Y[jco47]*rate_eval.p_Co47__Ni48
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    # reaclib rates
    Be7__Li7__weak__electron_capture(rate_eval, tf)
    O15__N15__weak__wc12(rate_eval, tf)
    He3__p_d(rate_eval, tf)
    He4__d_d(rate_eval, tf)
    Be7__He4_He3(rate_eval, tf)
    B8__p_Be7(rate_eval, tf)
    B8__He4_He4__weak__wc12(rate_eval, tf)
    O16__p_N15(rate_eval, tf)
    O16__He4_C12(rate_eval, tf)
    F19__p_O18(rate_eval, tf)
    F19__He4_N15(rate_eval, tf)
    C12__He4_He4_He4(rate_eval, tf)
    p_p__d__weak__bet_pos_(rate_eval, tf)
    p_p__d__weak__electron_capture(rate_eval, tf)
    p_d__He3(rate_eval, tf)
    d_d__He4(rate_eval, tf)
    p_He3__He4__weak__bet_pos_(rate_eval, tf)
    He4_He3__Be7(rate_eval, tf)
    p_Be7__B8(rate_eval, tf)
    He4_C12__O16(rate_eval, tf)
    p_N15__O16(rate_eval, tf)
    He4_N15__F19(rate_eval, tf)
    p_O18__F19(rate_eval, tf)
    d_He3__p_He4(rate_eval, tf)
    p_He4__d_He3(rate_eval, tf)
    He4_He4__p_Li7(rate_eval, tf)
    p_Li7__He4_He4(rate_eval, tf)
    He4_C12__p_N15(rate_eval, tf)
    p_N15__He4_C12(rate_eval, tf)
    He4_N15__p_O18(rate_eval, tf)
    He4_O16__p_F19(rate_eval, tf)
    He4_O17__p_F20(rate_eval, tf)
    p_O18__He4_N15(rate_eval, tf)
    p_F19__He4_O16(rate_eval, tf)
    p_F20__He4_O17(rate_eval, tf)
    He3_He3__p_p_He4(rate_eval, tf)
    d_Be7__p_He4_He4(rate_eval, tf)
    He3_Be7__p_p_He4_He4(rate_eval, tf)
    He4_He4_He4__C12(rate_eval, tf)
    p_p_He4__He3_He3(rate_eval, tf)
    p_He4_He4__d_Be7(rate_eval, tf)
    p_p_He4_He4__He3_Be7(rate_eval, tf)
    p_O15__F16(rate_eval, tf)
    F16__p_O15(rate_eval, tf)
    p_S28__Cl29(rate_eval, tf)
    He4_S28__Ar32(rate_eval, tf)
    p_S29__Cl30(rate_eval, tf)
    He4_S29__Ar33(rate_eval, tf)
    Cl29__S29__weak__bqa_pos_(rate_eval, tf)
    Cl29__p_S28(rate_eval, tf)
    He4_Cl29__K33(rate_eval, tf)
    He4_Cl29__p_Ar32(rate_eval, tf)
    Ar32__He4_S28(rate_eval, tf)
    p_Ar32__K33(rate_eval, tf)
    He4_Ar32__Ca36(rate_eval, tf)
    p_Ar32__He4_Cl29(rate_eval, tf)
    p_Ca38__Sc39(rate_eval, tf)
    p_p_Ca38__He4_Ca36(rate_eval, tf)
    Cl30__p_S29(rate_eval, tf)
    He4_Cl30__K34(rate_eval, tf)
    He4_Cl30__p_Ar33(rate_eval, tf)
    Ar33__He4_S29(rate_eval, tf)
    p_Ar33__K34(rate_eval, tf)
    He4_Ar33__Ca37(rate_eval, tf)
    p_Ar33__He4_Cl30(rate_eval, tf)
    K33__Ar33__weak__bqa_pos_(rate_eval, tf)
    K33__p_Ar32(rate_eval, tf)
    K33__He4_Cl29(rate_eval, tf)
    He4_K33__Sc37(rate_eval, tf)
    He4_K33__p_Ca36(rate_eval, tf)
    Ca36__He4_Ar32(rate_eval, tf)
    p_Ca36__Sc37(rate_eval, tf)
    He4_Ca36__Ti40(rate_eval, tf)
    p_Ca36__He4_K33(rate_eval, tf)
    He4_Ca36__p_Sc39(rate_eval, tf)
    He4_Ca36__p_p_Ca38(rate_eval, tf)
    Sc39__p_Ca38(rate_eval, tf)
    p_Sc39__Ti40(rate_eval, tf)
    p_Sc39__He4_Ca36(rate_eval, tf)
    K34__p_Ar33(rate_eval, tf)
    K34__He4_Cl30(rate_eval, tf)
    p_K34__Ca35(rate_eval, tf)
    He4_K34__Sc38(rate_eval, tf)
    He4_K34__p_Ca37(rate_eval, tf)
    Ca37__He4_Ar33(rate_eval, tf)
    p_Ca37__Sc38(rate_eval, tf)
    p_Ca37__He4_K34(rate_eval, tf)
    Sc37__Ca37__weak__bqa_pos_(rate_eval, tf)
    Sc37__p_Ca36(rate_eval, tf)
    Sc37__He4_K33(rate_eval, tf)
    He4_Sc37__V41(rate_eval, tf)
    He4_Sc37__p_Ti40(rate_eval, tf)
    Ti40__p_Sc39(rate_eval, tf)
    Ti40__He4_Ca36(rate_eval, tf)
    p_Ti40__V41(rate_eval, tf)
    He4_Ti40__Cr44(rate_eval, tf)
    p_Ti40__He4_Sc37(rate_eval, tf)
    Ca35__p_K34(rate_eval, tf)
    p_Ca35__Sc36(rate_eval, tf)
    He4_Ca35__Ti39(rate_eval, tf)
    He4_Ca35__p_Sc38(rate_eval, tf)
    Sc38__Ca38__weak__mo97(rate_eval, tf)
    Sc38__p_Ca37(rate_eval, tf)
    Sc38__He4_K34(rate_eval, tf)
    p_Sc38__Ti39(rate_eval, tf)
    p_Sc38__He4_Ca35(rate_eval, tf)
    V41__p_Ti40(rate_eval, tf)
    V41__He4_Sc37(rate_eval, tf)
    He4_V41__Mn45(rate_eval, tf)
    He4_V41__p_Cr44(rate_eval, tf)
    Cr44__He4_Ti40(rate_eval, tf)
    p_Cr44__Mn45(rate_eval, tf)
    He4_Cr44__Fe48(rate_eval, tf)
    p_Cr44__He4_V41(rate_eval, tf)
    Sc36__Ca36__weak__bqa_pos_(rate_eval, tf)
    Sc36__p_Ca35(rate_eval, tf)
    He4_Sc36__V40(rate_eval, tf)
    He4_Sc36__p_Ti39(rate_eval, tf)
    Ti39__Sc39__weak__wc17(rate_eval, tf)
    Ti39__p_Sc38(rate_eval, tf)
    Ti39__p_Ca38__weak__wc12(rate_eval, tf)
    Ti39__He4_Ca35(rate_eval, tf)
    p_Ti39__V40(rate_eval, tf)
    He4_Ti39__Cr43(rate_eval, tf)
    p_Ti39__He4_Sc36(rate_eval, tf)
    Mn45__p_Cr44(rate_eval, tf)
    Mn45__He4_V41(rate_eval, tf)
    p_Mn45__Fe46(rate_eval, tf)
    He4_Mn45__Co49(rate_eval, tf)
    He4_Mn45__p_Fe48(rate_eval, tf)
    Fe48__He4_Cr44(rate_eval, tf)
    p_Fe48__Co49(rate_eval, tf)
    p_Fe48__He4_Mn45(rate_eval, tf)
    V40__Ti40__weak__bqa_pos_(rate_eval, tf)
    V40__p_Ti39(rate_eval, tf)
    V40__He4_Sc36(rate_eval, tf)
    He4_V40__Mn44(rate_eval, tf)
    He4_V40__p_Cr43(rate_eval, tf)
    Cr43__He4_Ti39(rate_eval, tf)
    p_Cr43__Mn44(rate_eval, tf)
    He4_Cr43__Fe47(rate_eval, tf)
    p_Cr43__He4_V40(rate_eval, tf)
    Fe46__p_Mn45(rate_eval, tf)
    p_Fe46__Co47(rate_eval, tf)
    He4_Fe46__p_Co49(rate_eval, tf)
    Co49__p_Fe48(rate_eval, tf)
    Co49__He4_Mn45(rate_eval, tf)
    p_Co49__He4_Fe46(rate_eval, tf)
    Mn44__Cr44__weak__bqa_pos_(rate_eval, tf)
    Mn44__p_Cr43(rate_eval, tf)
    Mn44__He4_V40(rate_eval, tf)
    p_Mn44__Fe45(rate_eval, tf)
    He4_Mn44__Co48(rate_eval, tf)
    He4_Mn44__p_Fe47(rate_eval, tf)
    Fe47__He4_Cr43(rate_eval, tf)
    p_Fe47__Co48(rate_eval, tf)
    p_Fe47__He4_Mn44(rate_eval, tf)
    Co47__Fe47__weak__bqa_pos_(rate_eval, tf)
    Co47__p_Fe46(rate_eval, tf)
    p_Co47__Ni48(rate_eval, tf)
    Fe45__Mn45__weak__wc17(rate_eval, tf)
    Fe45__p_Mn44(rate_eval, tf)
    Fe45__p_Cr44__weak__wc17(rate_eval, tf)
    He4_Fe45__p_Co48(rate_eval, tf)
    Co48__Fe48__weak__bqa_pos_(rate_eval, tf)
    Co48__p_Fe47(rate_eval, tf)
    Co48__He4_Mn44(rate_eval, tf)
    p_Co48__He4_Fe45(rate_eval, tf)
    Ni48__Co48__weak__wc17(rate_eval, tf)
    Ni48__p_Co47(rate_eval, tf)

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p__d__weak__bet_pos_ *= scor
        rate_eval.p_p__d__weak__electron_capture *= scor
        rate_eval.p_p_He4_He4__He3_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_d__He3 *= scor

        scn_fac = ScreenFactors(1, 2, 1, 2)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_d__He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He3__He4__weak__bet_pos_ *= scor

        scn_fac = ScreenFactors(2, 4, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He3__Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Be7__B8 *= scor

        scn_fac = ScreenFactors(2, 4, 6, 12)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_C12__O16 *= scor
        rate_eval.He4_C12__p_N15 *= scor

        scn_fac = ScreenFactors(1, 1, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_N15__O16 *= scor
        rate_eval.p_N15__He4_C12 *= scor

        scn_fac = ScreenFactors(2, 4, 7, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_N15__F19 *= scor
        rate_eval.He4_N15__p_O18 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 18)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O18__F19 *= scor
        rate_eval.p_O18__He4_N15 *= scor

        scn_fac = ScreenFactors(1, 2, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_He3__p_He4 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4__d_He3 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_He4__p_Li7 *= scor

        scn_fac = ScreenFactors(1, 1, 3, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Li7__He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 16)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O16__p_F19 *= scor

        scn_fac = ScreenFactors(2, 4, 8, 17)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_O17__p_F20 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 19)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F19__He4_O16 *= scor

        scn_fac = ScreenFactors(1, 1, 9, 20)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_F20__He4_O17 *= scor

        scn_fac = ScreenFactors(2, 3, 2, 3)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_He3__p_p_He4 *= scor

        scn_fac = ScreenFactors(1, 2, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.d_Be7__p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 3, 4, 7)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He3_Be7__p_p_He4_He4 *= scor

        scn_fac = ScreenFactors(2, 4, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        scn_fac2 = ScreenFactors(2, 4, 4, 8)
        scor2 = screen_func(plasma_state, scn_fac2)
        rate_eval.He4_He4_He4__C12 *= scor * scor2

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_He4__He3_He3 *= scor

        scn_fac = ScreenFactors(1, 1, 2, 4)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_He4_He4__d_Be7 *= scor

        scn_fac = ScreenFactors(1, 1, 8, 15)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_O15__F16 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S28__Cl29 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 28)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S28__Ar32 *= scor

        scn_fac = ScreenFactors(1, 1, 16, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_S29__Cl30 *= scor

        scn_fac = ScreenFactors(2, 4, 16, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_S29__Ar33 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 29)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl29__K33 *= scor
        rate_eval.He4_Cl29__p_Ar32 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar32__K33 *= scor
        rate_eval.p_Ar32__He4_Cl29 *= scor

        scn_fac = ScreenFactors(2, 4, 18, 32)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ar32__Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 38)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca38__Sc39 *= scor

        scn_fac = ScreenFactors(1, 1, 1, 1)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_p_Ca38__He4_Ca36 *= scor

        scn_fac = ScreenFactors(2, 4, 17, 30)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cl30__K34 *= scor
        rate_eval.He4_Cl30__p_Ar33 *= scor

        scn_fac = ScreenFactors(1, 1, 18, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ar33__K34 *= scor
        rate_eval.p_Ar33__He4_Cl30 *= scor

        scn_fac = ScreenFactors(2, 4, 18, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ar33__Ca37 *= scor

        scn_fac = ScreenFactors(2, 4, 19, 33)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_K33__Sc37 *= scor
        rate_eval.He4_K33__p_Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca36__Sc37 *= scor
        rate_eval.p_Ca36__He4_K33 *= scor

        scn_fac = ScreenFactors(2, 4, 20, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ca36__Ti40 *= scor
        rate_eval.He4_Ca36__p_Sc39 *= scor
        rate_eval.He4_Ca36__p_p_Ca38 *= scor

        scn_fac = ScreenFactors(1, 1, 21, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Sc39__Ti40 *= scor
        rate_eval.p_Sc39__He4_Ca36 *= scor

        scn_fac = ScreenFactors(1, 1, 19, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_K34__Ca35 *= scor

        scn_fac = ScreenFactors(2, 4, 19, 34)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_K34__Sc38 *= scor
        rate_eval.He4_K34__p_Ca37 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 37)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca37__Sc38 *= scor
        rate_eval.p_Ca37__He4_K34 *= scor

        scn_fac = ScreenFactors(2, 4, 21, 37)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Sc37__V41 *= scor
        rate_eval.He4_Sc37__p_Ti40 *= scor

        scn_fac = ScreenFactors(1, 1, 22, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ti40__V41 *= scor
        rate_eval.p_Ti40__He4_Sc37 *= scor

        scn_fac = ScreenFactors(2, 4, 22, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ti40__Cr44 *= scor

        scn_fac = ScreenFactors(1, 1, 20, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ca35__Sc36 *= scor

        scn_fac = ScreenFactors(2, 4, 20, 35)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ca35__Ti39 *= scor
        rate_eval.He4_Ca35__p_Sc38 *= scor

        scn_fac = ScreenFactors(1, 1, 21, 38)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Sc38__Ti39 *= scor
        rate_eval.p_Sc38__He4_Ca35 *= scor

        scn_fac = ScreenFactors(2, 4, 23, 41)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_V41__Mn45 *= scor
        rate_eval.He4_V41__p_Cr44 *= scor

        scn_fac = ScreenFactors(1, 1, 24, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cr44__Mn45 *= scor
        rate_eval.p_Cr44__He4_V41 *= scor

        scn_fac = ScreenFactors(2, 4, 24, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cr44__Fe48 *= scor

        scn_fac = ScreenFactors(2, 4, 21, 36)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Sc36__V40 *= scor
        rate_eval.He4_Sc36__p_Ti39 *= scor

        scn_fac = ScreenFactors(1, 1, 22, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Ti39__V40 *= scor
        rate_eval.p_Ti39__He4_Sc36 *= scor

        scn_fac = ScreenFactors(2, 4, 22, 39)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Ti39__Cr43 *= scor

        scn_fac = ScreenFactors(1, 1, 25, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mn45__Fe46 *= scor

        scn_fac = ScreenFactors(2, 4, 25, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mn45__Co49 *= scor
        rate_eval.He4_Mn45__p_Fe48 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 48)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe48__Co49 *= scor
        rate_eval.p_Fe48__He4_Mn45 *= scor

        scn_fac = ScreenFactors(2, 4, 23, 40)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_V40__Mn44 *= scor
        rate_eval.He4_V40__p_Cr43 *= scor

        scn_fac = ScreenFactors(1, 1, 24, 43)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Cr43__Mn44 *= scor
        rate_eval.p_Cr43__He4_V40 *= scor

        scn_fac = ScreenFactors(2, 4, 24, 43)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Cr43__Fe47 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 46)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe46__Co47 *= scor

        scn_fac = ScreenFactors(2, 4, 26, 46)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Fe46__p_Co49 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 49)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co49__He4_Fe46 *= scor

        scn_fac = ScreenFactors(1, 1, 25, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Mn44__Fe45 *= scor

        scn_fac = ScreenFactors(2, 4, 25, 44)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Mn44__Co48 *= scor
        rate_eval.He4_Mn44__p_Fe47 *= scor

        scn_fac = ScreenFactors(1, 1, 26, 47)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Fe47__Co48 *= scor
        rate_eval.p_Fe47__He4_Mn44 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 47)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co47__Ni48 *= scor

        scn_fac = ScreenFactors(2, 4, 26, 45)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.He4_Fe45__p_Co48 *= scor

        scn_fac = ScreenFactors(1, 1, 27, 48)
        scor = screen_func(plasma_state, scn_fac)
        rate_eval.p_Co48__He4_Fe45 *= scor

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jp, jp] = (
       -2*5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p__d__weak__bet_pos_
       -2*5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p__d__weak__electron_capture
       -rho*Y[jd]*rate_eval.p_d__He3
       -rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jbe7]*rate_eval.p_Be7__B8
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jhe4]*rate_eval.p_He4__d_He3
       -rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       -rho*Y[jf20]*rate_eval.p_F20__He4_O17
       -2*5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       -rho*Y[jo15]*rate_eval.p_O15__F16
       -rho*Y[js28]*rate_eval.p_S28__Cl29
       -rho*Y[js29]*rate_eval.p_S29__Cl30
       -rho*Y[jar32]*rate_eval.p_Ar32__K33
       -rho*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       -rho*Y[jca38]*rate_eval.p_Ca38__Sc39
       -2*5.00000000000000e-01*rho**2*2*Y[jp]*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       -rho*Y[jar33]*rate_eval.p_Ar33__K34
       -rho*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       -rho*Y[jca36]*rate_eval.p_Ca36__Sc37
       -rho*Y[jca36]*rate_eval.p_Ca36__He4_K33
       -rho*Y[jsc39]*rate_eval.p_Sc39__Ti40
       -rho*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       -rho*Y[jk34]*rate_eval.p_K34__Ca35
       -rho*Y[jca37]*rate_eval.p_Ca37__Sc38
       -rho*Y[jca37]*rate_eval.p_Ca37__He4_K34
       -rho*Y[jti40]*rate_eval.p_Ti40__V41
       -rho*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       -rho*Y[jca35]*rate_eval.p_Ca35__Sc36
       -rho*Y[jsc38]*rate_eval.p_Sc38__Ti39
       -rho*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       -rho*Y[jcr44]*rate_eval.p_Cr44__Mn45
       -rho*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       -rho*Y[jti39]*rate_eval.p_Ti39__V40
       -rho*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       -rho*Y[jmn45]*rate_eval.p_Mn45__Fe46
       -rho*Y[jfe48]*rate_eval.p_Fe48__Co49
       -rho*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       -rho*Y[jcr43]*rate_eval.p_Cr43__Mn44
       -rho*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       -rho*Y[jfe46]*rate_eval.p_Fe46__Co47
       -rho*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       -rho*Y[jmn44]*rate_eval.p_Mn44__Fe45
       -rho*Y[jfe47]*rate_eval.p_Fe47__Co48
       -rho*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       -rho*Y[jco47]*rate_eval.p_Co47__Ni48
       -rho*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jp, jd] = (
       -rho*Y[jp]*rate_eval.p_d__He3
       +rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jp, jhe3] = (
       -rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       +rate_eval.He3__p_d
       +rho*Y[jd]*rate_eval.d_He3__p_He4
       +2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       +2*rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jp, jhe4] = (
       -rho*Y[jp]*rate_eval.p_He4__d_He3
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       -5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       +rho*Y[jc12]*rate_eval.He4_C12__p_N15
       +rho*Y[jn15]*rate_eval.He4_N15__p_O18
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jo17]*rate_eval.He4_O17__p_F20
       +rho*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       +rho*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       +rho*Y[jk33]*rate_eval.He4_K33__p_Ca36
       +rho*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       +2*rho*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       +rho*Y[jk34]*rate_eval.He4_K34__p_Ca37
       +rho*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       +rho*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       +rho*Y[jv41]*rate_eval.He4_V41__p_Cr44
       +rho*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       +rho*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       +rho*Y[jv40]*rate_eval.He4_V40__p_Cr43
       +rho*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       +rho*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       +rho*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jp, jli7] = (
       -rho*Y[jp]*rate_eval.p_Li7__He4_He4
       )

    jac[jp, jbe7] = (
       -rho*Y[jp]*rate_eval.p_Be7__B8
       +rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       +2*rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jp, jb8] = (
       +rate_eval.B8__p_Be7
       )

    jac[jp, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       )

    jac[jp, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       +rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jp, jo15] = (
       -rho*Y[jp]*rate_eval.p_O15__F16
       )

    jac[jp, jo16] = (
       +rate_eval.O16__p_N15
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jp, jo17] = (
       +rho*Y[jhe4]*rate_eval.He4_O17__p_F20
       )

    jac[jp, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jp, jf16] = (
       +rate_eval.F16__p_O15
       )

    jac[jp, jf19] = (
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       +rate_eval.F19__p_O18
       )

    jac[jp, jf20] = (
       -rho*Y[jp]*rate_eval.p_F20__He4_O17
       )

    jac[jp, js28] = (
       -rho*Y[jp]*rate_eval.p_S28__Cl29
       )

    jac[jp, js29] = (
       -rho*Y[jp]*rate_eval.p_S29__Cl30
       )

    jac[jp, jcl29] = (
       +rate_eval.Cl29__p_S28
       +rho*Y[jhe4]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jp, jcl30] = (
       +rate_eval.Cl30__p_S29
       +rho*Y[jhe4]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jp, jar32] = (
       -rho*Y[jp]*rate_eval.p_Ar32__K33
       -rho*Y[jp]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jp, jar33] = (
       -rho*Y[jp]*rate_eval.p_Ar33__K34
       -rho*Y[jp]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jp, jk33] = (
       +rate_eval.K33__p_Ar32
       +rho*Y[jhe4]*rate_eval.He4_K33__p_Ca36
       )

    jac[jp, jk34] = (
       -rho*Y[jp]*rate_eval.p_K34__Ca35
       +rate_eval.K34__p_Ar33
       +rho*Y[jhe4]*rate_eval.He4_K34__p_Ca37
       )

    jac[jp, jca35] = (
       -rho*Y[jp]*rate_eval.p_Ca35__Sc36
       +rate_eval.Ca35__p_K34
       +rho*Y[jhe4]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jp, jca36] = (
       -rho*Y[jp]*rate_eval.p_Ca36__Sc37
       -rho*Y[jp]*rate_eval.p_Ca36__He4_K33
       +rho*Y[jhe4]*rate_eval.He4_Ca36__p_Sc39
       +2*rho*Y[jhe4]*rate_eval.He4_Ca36__p_p_Ca38
       )

    jac[jp, jca37] = (
       -rho*Y[jp]*rate_eval.p_Ca37__Sc38
       -rho*Y[jp]*rate_eval.p_Ca37__He4_K34
       )

    jac[jp, jca38] = (
       -rho*Y[jp]*rate_eval.p_Ca38__Sc39
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_Ca38__He4_Ca36
       )

    jac[jp, jsc36] = (
       +rate_eval.Sc36__p_Ca35
       +rho*Y[jhe4]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jp, jsc37] = (
       +rate_eval.Sc37__p_Ca36
       +rho*Y[jhe4]*rate_eval.He4_Sc37__p_Ti40
       )

    jac[jp, jsc38] = (
       -rho*Y[jp]*rate_eval.p_Sc38__Ti39
       -rho*Y[jp]*rate_eval.p_Sc38__He4_Ca35
       +rate_eval.Sc38__p_Ca37
       )

    jac[jp, jsc39] = (
       -rho*Y[jp]*rate_eval.p_Sc39__Ti40
       -rho*Y[jp]*rate_eval.p_Sc39__He4_Ca36
       +rate_eval.Sc39__p_Ca38
       )

    jac[jp, jti39] = (
       -rho*Y[jp]*rate_eval.p_Ti39__V40
       -rho*Y[jp]*rate_eval.p_Ti39__He4_Sc36
       +rate_eval.Ti39__p_Sc38
       +rate_eval.Ti39__p_Ca38__weak__wc12
       )

    jac[jp, jti40] = (
       -rho*Y[jp]*rate_eval.p_Ti40__V41
       -rho*Y[jp]*rate_eval.p_Ti40__He4_Sc37
       +rate_eval.Ti40__p_Sc39
       )

    jac[jp, jv40] = (
       +rate_eval.V40__p_Ti39
       +rho*Y[jhe4]*rate_eval.He4_V40__p_Cr43
       )

    jac[jp, jv41] = (
       +rate_eval.V41__p_Ti40
       +rho*Y[jhe4]*rate_eval.He4_V41__p_Cr44
       )

    jac[jp, jcr43] = (
       -rho*Y[jp]*rate_eval.p_Cr43__Mn44
       -rho*Y[jp]*rate_eval.p_Cr43__He4_V40
       )

    jac[jp, jcr44] = (
       -rho*Y[jp]*rate_eval.p_Cr44__Mn45
       -rho*Y[jp]*rate_eval.p_Cr44__He4_V41
       )

    jac[jp, jmn44] = (
       -rho*Y[jp]*rate_eval.p_Mn44__Fe45
       +rate_eval.Mn44__p_Cr43
       +rho*Y[jhe4]*rate_eval.He4_Mn44__p_Fe47
       )

    jac[jp, jmn45] = (
       -rho*Y[jp]*rate_eval.p_Mn45__Fe46
       +rate_eval.Mn45__p_Cr44
       +rho*Y[jhe4]*rate_eval.He4_Mn45__p_Fe48
       )

    jac[jp, jfe45] = (
       +rate_eval.Fe45__p_Mn44
       +rate_eval.Fe45__p_Cr44__weak__wc17
       +rho*Y[jhe4]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jp, jfe46] = (
       -rho*Y[jp]*rate_eval.p_Fe46__Co47
       +rate_eval.Fe46__p_Mn45
       +rho*Y[jhe4]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jp, jfe47] = (
       -rho*Y[jp]*rate_eval.p_Fe47__Co48
       -rho*Y[jp]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jp, jfe48] = (
       -rho*Y[jp]*rate_eval.p_Fe48__Co49
       -rho*Y[jp]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jp, jco47] = (
       -rho*Y[jp]*rate_eval.p_Co47__Ni48
       +rate_eval.Co47__p_Fe46
       )

    jac[jp, jco48] = (
       -rho*Y[jp]*rate_eval.p_Co48__He4_Fe45
       +rate_eval.Co48__p_Fe47
       )

    jac[jp, jco49] = (
       -rho*Y[jp]*rate_eval.p_Co49__He4_Fe46
       +rate_eval.Co49__p_Fe48
       )

    jac[jp, jni48] = (
       +rate_eval.Ni48__p_Co47
       )

    jac[jd, jp] = (
       -rho*Y[jd]*rate_eval.p_d__He3
       +5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p__d__weak__bet_pos_
       +5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p__d__weak__electron_capture
       +rho*Y[jhe4]*rate_eval.p_He4__d_He3
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       )

    jac[jd, jd] = (
       -rho*Y[jp]*rate_eval.p_d__He3
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__He4
       -rho*Y[jhe3]*rate_eval.d_He3__p_He4
       -rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jd, jhe3] = (
       -rho*Y[jd]*rate_eval.d_He3__p_He4
       +rate_eval.He3__p_d
       )

    jac[jd, jhe4] = (
       +2*rate_eval.He4__d_d
       +rho*Y[jp]*rate_eval.p_He4__d_He3
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       )

    jac[jd, jbe7] = (
       -rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jhe3, jp] = (
       -rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jd]*rate_eval.p_d__He3
       +rho*Y[jhe4]*rate_eval.p_He4__d_He3
       +2*5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       +2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jhe3, jd] = (
       -rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +rho*Y[jp]*rate_eval.p_d__He3
       )

    jac[jhe3, jhe3] = (
       -rate_eval.He3__p_d
       -rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       -rho*Y[jhe4]*rate_eval.He4_He3__Be7
       -rho*Y[jd]*rate_eval.d_He3__p_He4
       -2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       -rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe3, jhe4] = (
       -rho*Y[jhe3]*rate_eval.He4_He3__Be7
       +rho*Y[jp]*rate_eval.p_He4__d_He3
       +2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       +2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jhe3, jbe7] = (
       -rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       +rate_eval.Be7__He4_He3
       )

    jac[jhe4, jp] = (
       -rho*Y[jhe4]*rate_eval.p_He4__d_He3
       -5.00000000000000e-01*rho**2*2*Y[jp]*Y[jhe4]*rate_eval.p_p_He4__He3_He3
       -2*5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       +rho*Y[jhe3]*rate_eval.p_He3__He4__weak__bet_pos_
       +2*rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jf20]*rate_eval.p_F20__He4_O17
       +rho*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       +5.00000000000000e-01*rho**2*2*Y[jp]*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       +rho*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       +rho*Y[jca36]*rate_eval.p_Ca36__He4_K33
       +rho*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       +rho*Y[jca37]*rate_eval.p_Ca37__He4_K34
       +rho*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       +rho*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       +rho*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       +rho*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       +rho*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       +rho*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       +rho*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       +rho*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       +rho*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jhe4, jd] = (
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d__He4
       +rho*Y[jhe3]*rate_eval.d_He3__p_He4
       +2*rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jhe4, jhe3] = (
       -rho*Y[jhe4]*rate_eval.He4_He3__Be7
       +rho*Y[jp]*rate_eval.p_He3__He4__weak__bet_pos_
       +rho*Y[jd]*rate_eval.d_He3__p_He4
       +5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3__p_p_He4
       +2*rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe4, jhe4] = (
       -rate_eval.He4__d_d
       -rho*Y[jhe3]*rate_eval.He4_He3__Be7
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jp]*rate_eval.p_He4__d_He3
       -2*5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       -rho*Y[jc12]*rate_eval.He4_C12__p_N15
       -rho*Y[jn15]*rate_eval.He4_N15__p_O18
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       -rho*Y[jo17]*rate_eval.He4_O17__p_F20
       -3*1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       -5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_He4__He3_He3
       -2*5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       -2*2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       -rho*Y[js28]*rate_eval.He4_S28__Ar32
       -rho*Y[js29]*rate_eval.He4_S29__Ar33
       -rho*Y[jcl29]*rate_eval.He4_Cl29__K33
       -rho*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       -rho*Y[jar32]*rate_eval.He4_Ar32__Ca36
       -rho*Y[jcl30]*rate_eval.He4_Cl30__K34
       -rho*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       -rho*Y[jar33]*rate_eval.He4_Ar33__Ca37
       -rho*Y[jk33]*rate_eval.He4_K33__Sc37
       -rho*Y[jk33]*rate_eval.He4_K33__p_Ca36
       -rho*Y[jca36]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       -rho*Y[jk34]*rate_eval.He4_K34__Sc38
       -rho*Y[jk34]*rate_eval.He4_K34__p_Ca37
       -rho*Y[jsc37]*rate_eval.He4_Sc37__V41
       -rho*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       -rho*Y[jti40]*rate_eval.He4_Ti40__Cr44
       -rho*Y[jca35]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       -rho*Y[jv41]*rate_eval.He4_V41__Mn45
       -rho*Y[jv41]*rate_eval.He4_V41__p_Cr44
       -rho*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       -rho*Y[jsc36]*rate_eval.He4_Sc36__V40
       -rho*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       -rho*Y[jti39]*rate_eval.He4_Ti39__Cr43
       -rho*Y[jmn45]*rate_eval.He4_Mn45__Co49
       -rho*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       -rho*Y[jv40]*rate_eval.He4_V40__Mn44
       -rho*Y[jv40]*rate_eval.He4_V40__p_Cr43
       -rho*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       -rho*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       -rho*Y[jmn44]*rate_eval.He4_Mn44__Co48
       -rho*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       -rho*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jhe4, jli7] = (
       +2*rho*Y[jp]*rate_eval.p_Li7__He4_He4
       )

    jac[jhe4, jbe7] = (
       +rate_eval.Be7__He4_He3
       +2*rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       +2*rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jhe4, jb8] = (
       +2*rate_eval.B8__He4_He4__weak__wc12
       )

    jac[jhe4, jc12] = (
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       +3*rate_eval.C12__He4_He4_He4
       )

    jac[jhe4, jn15] = (
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       -rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jhe4, jo16] = (
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       +rate_eval.O16__He4_C12
       )

    jac[jhe4, jo17] = (
       -rho*Y[jhe4]*rate_eval.He4_O17__p_F20
       )

    jac[jhe4, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jhe4, jf19] = (
       +rate_eval.F19__He4_N15
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jhe4, jf20] = (
       +rho*Y[jp]*rate_eval.p_F20__He4_O17
       )

    jac[jhe4, js28] = (
       -rho*Y[jhe4]*rate_eval.He4_S28__Ar32
       )

    jac[jhe4, js29] = (
       -rho*Y[jhe4]*rate_eval.He4_S29__Ar33
       )

    jac[jhe4, jcl29] = (
       -rho*Y[jhe4]*rate_eval.He4_Cl29__K33
       -rho*Y[jhe4]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jhe4, jcl30] = (
       -rho*Y[jhe4]*rate_eval.He4_Cl30__K34
       -rho*Y[jhe4]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jhe4, jar32] = (
       -rho*Y[jhe4]*rate_eval.He4_Ar32__Ca36
       +rate_eval.Ar32__He4_S28
       +rho*Y[jp]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jhe4, jar33] = (
       -rho*Y[jhe4]*rate_eval.He4_Ar33__Ca37
       +rate_eval.Ar33__He4_S29
       +rho*Y[jp]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jhe4, jk33] = (
       -rho*Y[jhe4]*rate_eval.He4_K33__Sc37
       -rho*Y[jhe4]*rate_eval.He4_K33__p_Ca36
       +rate_eval.K33__He4_Cl29
       )

    jac[jhe4, jk34] = (
       -rho*Y[jhe4]*rate_eval.He4_K34__Sc38
       -rho*Y[jhe4]*rate_eval.He4_K34__p_Ca37
       +rate_eval.K34__He4_Cl30
       )

    jac[jhe4, jca35] = (
       -rho*Y[jhe4]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jhe4]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jhe4, jca36] = (
       -rho*Y[jhe4]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jhe4]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jhe4]*rate_eval.He4_Ca36__p_p_Ca38
       +rate_eval.Ca36__He4_Ar32
       +rho*Y[jp]*rate_eval.p_Ca36__He4_K33
       )

    jac[jhe4, jca37] = (
       +rate_eval.Ca37__He4_Ar33
       +rho*Y[jp]*rate_eval.p_Ca37__He4_K34
       )

    jac[jhe4, jca38] = (
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_Ca38__He4_Ca36
       )

    jac[jhe4, jsc36] = (
       -rho*Y[jhe4]*rate_eval.He4_Sc36__V40
       -rho*Y[jhe4]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jhe4, jsc37] = (
       -rho*Y[jhe4]*rate_eval.He4_Sc37__V41
       -rho*Y[jhe4]*rate_eval.He4_Sc37__p_Ti40
       +rate_eval.Sc37__He4_K33
       )

    jac[jhe4, jsc38] = (
       +rate_eval.Sc38__He4_K34
       +rho*Y[jp]*rate_eval.p_Sc38__He4_Ca35
       )

    jac[jhe4, jsc39] = (
       +rho*Y[jp]*rate_eval.p_Sc39__He4_Ca36
       )

    jac[jhe4, jti39] = (
       -rho*Y[jhe4]*rate_eval.He4_Ti39__Cr43
       +rate_eval.Ti39__He4_Ca35
       +rho*Y[jp]*rate_eval.p_Ti39__He4_Sc36
       )

    jac[jhe4, jti40] = (
       -rho*Y[jhe4]*rate_eval.He4_Ti40__Cr44
       +rate_eval.Ti40__He4_Ca36
       +rho*Y[jp]*rate_eval.p_Ti40__He4_Sc37
       )

    jac[jhe4, jv40] = (
       -rho*Y[jhe4]*rate_eval.He4_V40__Mn44
       -rho*Y[jhe4]*rate_eval.He4_V40__p_Cr43
       +rate_eval.V40__He4_Sc36
       )

    jac[jhe4, jv41] = (
       -rho*Y[jhe4]*rate_eval.He4_V41__Mn45
       -rho*Y[jhe4]*rate_eval.He4_V41__p_Cr44
       +rate_eval.V41__He4_Sc37
       )

    jac[jhe4, jcr43] = (
       -rho*Y[jhe4]*rate_eval.He4_Cr43__Fe47
       +rate_eval.Cr43__He4_Ti39
       +rho*Y[jp]*rate_eval.p_Cr43__He4_V40
       )

    jac[jhe4, jcr44] = (
       -rho*Y[jhe4]*rate_eval.He4_Cr44__Fe48
       +rate_eval.Cr44__He4_Ti40
       +rho*Y[jp]*rate_eval.p_Cr44__He4_V41
       )

    jac[jhe4, jmn44] = (
       -rho*Y[jhe4]*rate_eval.He4_Mn44__Co48
       -rho*Y[jhe4]*rate_eval.He4_Mn44__p_Fe47
       +rate_eval.Mn44__He4_V40
       )

    jac[jhe4, jmn45] = (
       -rho*Y[jhe4]*rate_eval.He4_Mn45__Co49
       -rho*Y[jhe4]*rate_eval.He4_Mn45__p_Fe48
       +rate_eval.Mn45__He4_V41
       )

    jac[jhe4, jfe45] = (
       -rho*Y[jhe4]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jhe4, jfe46] = (
       -rho*Y[jhe4]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jhe4, jfe47] = (
       +rate_eval.Fe47__He4_Cr43
       +rho*Y[jp]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jhe4, jfe48] = (
       +rate_eval.Fe48__He4_Cr44
       +rho*Y[jp]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jhe4, jco48] = (
       +rate_eval.Co48__He4_Mn44
       +rho*Y[jp]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jhe4, jco49] = (
       +rate_eval.Co49__He4_Mn45
       +rho*Y[jp]*rate_eval.p_Co49__He4_Fe46
       )

    jac[jli7, jp] = (
       -rho*Y[jli7]*rate_eval.p_Li7__He4_He4
       )

    jac[jli7, jhe4] = (
       +5.00000000000000e-01*rho*2*Y[jhe4]*rate_eval.He4_He4__p_Li7
       )

    jac[jli7, jli7] = (
       -rho*Y[jp]*rate_eval.p_Li7__He4_He4
       )

    jac[jli7, jbe7] = (
       +rho*ye(Y)*rate_eval.Be7__Li7__weak__electron_capture
       )

    jac[jbe7, jp] = (
       -rho*Y[jbe7]*rate_eval.p_Be7__B8
       +5.00000000000000e-01*rho**2*Y[jhe4]**2*rate_eval.p_He4_He4__d_Be7
       +2.50000000000000e-01*rho**3*2*Y[jp]*Y[jhe4]**2*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jbe7, jd] = (
       -rho*Y[jbe7]*rate_eval.d_Be7__p_He4_He4
       )

    jac[jbe7, jhe3] = (
       -rho*Y[jbe7]*rate_eval.He3_Be7__p_p_He4_He4
       +rho*Y[jhe4]*rate_eval.He4_He3__Be7
       )

    jac[jbe7, jhe4] = (
       +rho*Y[jhe3]*rate_eval.He4_He3__Be7
       +5.00000000000000e-01*rho**2*Y[jp]*2*Y[jhe4]*rate_eval.p_He4_He4__d_Be7
       +2.50000000000000e-01*rho**3*Y[jp]**2*2*Y[jhe4]*rate_eval.p_p_He4_He4__He3_Be7
       )

    jac[jbe7, jbe7] = (
       -rho*ye(Y)*rate_eval.Be7__Li7__weak__electron_capture
       -rate_eval.Be7__He4_He3
       -rho*Y[jp]*rate_eval.p_Be7__B8
       -rho*Y[jd]*rate_eval.d_Be7__p_He4_He4
       -rho*Y[jhe3]*rate_eval.He3_Be7__p_p_He4_He4
       )

    jac[jbe7, jb8] = (
       +rate_eval.B8__p_Be7
       )

    jac[jb8, jp] = (
       +rho*Y[jbe7]*rate_eval.p_Be7__B8
       )

    jac[jb8, jbe7] = (
       +rho*Y[jp]*rate_eval.p_Be7__B8
       )

    jac[jb8, jb8] = (
       -rate_eval.B8__p_Be7
       -rate_eval.B8__He4_He4__weak__wc12
       )

    jac[jc12, jp] = (
       +rho*Y[jn15]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jhe4] = (
       -rho*Y[jc12]*rate_eval.He4_C12__O16
       -rho*Y[jc12]*rate_eval.He4_C12__p_N15
       +1.66666666666667e-01*rho**2*3*Y[jhe4]**2*rate_eval.He4_He4_He4__C12
       )

    jac[jc12, jc12] = (
       -rate_eval.C12__He4_He4_He4
       -rho*Y[jhe4]*rate_eval.He4_C12__O16
       -rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       )

    jac[jc12, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__He4_C12
       )

    jac[jc12, jo16] = (
       +rate_eval.O16__He4_C12
       )

    jac[jn15, jp] = (
       -rho*Y[jn15]*rate_eval.p_N15__O16
       -rho*Y[jn15]*rate_eval.p_N15__He4_C12
       +rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jhe4] = (
       -rho*Y[jn15]*rate_eval.He4_N15__F19
       -rho*Y[jn15]*rate_eval.He4_N15__p_O18
       +rho*Y[jc12]*rate_eval.He4_C12__p_N15
       )

    jac[jn15, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12__p_N15
       )

    jac[jn15, jn15] = (
       -rho*Y[jp]*rate_eval.p_N15__O16
       -rho*Y[jhe4]*rate_eval.He4_N15__F19
       -rho*Y[jp]*rate_eval.p_N15__He4_C12
       -rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jn15, jo15] = (
       +rate_eval.O15__N15__weak__wc12
       )

    jac[jn15, jo16] = (
       +rate_eval.O16__p_N15
       )

    jac[jn15, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jn15, jf19] = (
       +rate_eval.F19__He4_N15
       )

    jac[jo15, jp] = (
       -rho*Y[jo15]*rate_eval.p_O15__F16
       )

    jac[jo15, jo15] = (
       -rate_eval.O15__N15__weak__wc12
       -rho*Y[jp]*rate_eval.p_O15__F16
       )

    jac[jo15, jf16] = (
       +rate_eval.F16__p_O15
       )

    jac[jo16, jp] = (
       +rho*Y[jn15]*rate_eval.p_N15__O16
       +rho*Y[jf19]*rate_eval.p_F19__He4_O16
       )

    jac[jo16, jhe4] = (
       -rho*Y[jo16]*rate_eval.He4_O16__p_F19
       +rho*Y[jc12]*rate_eval.He4_C12__O16
       )

    jac[jo16, jc12] = (
       +rho*Y[jhe4]*rate_eval.He4_C12__O16
       )

    jac[jo16, jn15] = (
       +rho*Y[jp]*rate_eval.p_N15__O16
       )

    jac[jo16, jo16] = (
       -rate_eval.O16__p_N15
       -rate_eval.O16__He4_C12
       -rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jo16, jf19] = (
       +rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jo17, jp] = (
       +rho*Y[jf20]*rate_eval.p_F20__He4_O17
       )

    jac[jo17, jhe4] = (
       -rho*Y[jo17]*rate_eval.He4_O17__p_F20
       )

    jac[jo17, jo17] = (
       -rho*Y[jhe4]*rate_eval.He4_O17__p_F20
       )

    jac[jo17, jf20] = (
       +rho*Y[jp]*rate_eval.p_F20__He4_O17
       )

    jac[jo18, jp] = (
       -rho*Y[jo18]*rate_eval.p_O18__F19
       -rho*Y[jo18]*rate_eval.p_O18__He4_N15
       )

    jac[jo18, jhe4] = (
       +rho*Y[jn15]*rate_eval.He4_N15__p_O18
       )

    jac[jo18, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__p_O18
       )

    jac[jo18, jo18] = (
       -rho*Y[jp]*rate_eval.p_O18__F19
       -rho*Y[jp]*rate_eval.p_O18__He4_N15
       )

    jac[jo18, jf19] = (
       +rate_eval.F19__p_O18
       )

    jac[jf16, jp] = (
       +rho*Y[jo15]*rate_eval.p_O15__F16
       )

    jac[jf16, jo15] = (
       +rho*Y[jp]*rate_eval.p_O15__F16
       )

    jac[jf16, jf16] = (
       -rate_eval.F16__p_O15
       )

    jac[jf19, jp] = (
       -rho*Y[jf19]*rate_eval.p_F19__He4_O16
       +rho*Y[jo18]*rate_eval.p_O18__F19
       )

    jac[jf19, jhe4] = (
       +rho*Y[jn15]*rate_eval.He4_N15__F19
       +rho*Y[jo16]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jn15] = (
       +rho*Y[jhe4]*rate_eval.He4_N15__F19
       )

    jac[jf19, jo16] = (
       +rho*Y[jhe4]*rate_eval.He4_O16__p_F19
       )

    jac[jf19, jo18] = (
       +rho*Y[jp]*rate_eval.p_O18__F19
       )

    jac[jf19, jf19] = (
       -rate_eval.F19__p_O18
       -rate_eval.F19__He4_N15
       -rho*Y[jp]*rate_eval.p_F19__He4_O16
       )

    jac[jf20, jp] = (
       -rho*Y[jf20]*rate_eval.p_F20__He4_O17
       )

    jac[jf20, jhe4] = (
       +rho*Y[jo17]*rate_eval.He4_O17__p_F20
       )

    jac[jf20, jo17] = (
       +rho*Y[jhe4]*rate_eval.He4_O17__p_F20
       )

    jac[jf20, jf20] = (
       -rho*Y[jp]*rate_eval.p_F20__He4_O17
       )

    jac[js28, jp] = (
       -rho*Y[js28]*rate_eval.p_S28__Cl29
       )

    jac[js28, jhe4] = (
       -rho*Y[js28]*rate_eval.He4_S28__Ar32
       )

    jac[js28, js28] = (
       -rho*Y[jp]*rate_eval.p_S28__Cl29
       -rho*Y[jhe4]*rate_eval.He4_S28__Ar32
       )

    jac[js28, jcl29] = (
       +rate_eval.Cl29__p_S28
       )

    jac[js28, jar32] = (
       +rate_eval.Ar32__He4_S28
       )

    jac[js29, jp] = (
       -rho*Y[js29]*rate_eval.p_S29__Cl30
       )

    jac[js29, jhe4] = (
       -rho*Y[js29]*rate_eval.He4_S29__Ar33
       )

    jac[js29, js29] = (
       -rho*Y[jp]*rate_eval.p_S29__Cl30
       -rho*Y[jhe4]*rate_eval.He4_S29__Ar33
       )

    jac[js29, jcl29] = (
       +rate_eval.Cl29__S29__weak__bqa_pos_
       )

    jac[js29, jcl30] = (
       +rate_eval.Cl30__p_S29
       )

    jac[js29, jar33] = (
       +rate_eval.Ar33__He4_S29
       )

    jac[jcl29, jp] = (
       +rho*Y[js28]*rate_eval.p_S28__Cl29
       +rho*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jcl29, jhe4] = (
       -rho*Y[jcl29]*rate_eval.He4_Cl29__K33
       -rho*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jcl29, js28] = (
       +rho*Y[jp]*rate_eval.p_S28__Cl29
       )

    jac[jcl29, jcl29] = (
       -rate_eval.Cl29__S29__weak__bqa_pos_
       -rate_eval.Cl29__p_S28
       -rho*Y[jhe4]*rate_eval.He4_Cl29__K33
       -rho*Y[jhe4]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jcl29, jar32] = (
       +rho*Y[jp]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jcl29, jk33] = (
       +rate_eval.K33__He4_Cl29
       )

    jac[jcl30, jp] = (
       +rho*Y[js29]*rate_eval.p_S29__Cl30
       +rho*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jcl30, jhe4] = (
       -rho*Y[jcl30]*rate_eval.He4_Cl30__K34
       -rho*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jcl30, js29] = (
       +rho*Y[jp]*rate_eval.p_S29__Cl30
       )

    jac[jcl30, jcl30] = (
       -rate_eval.Cl30__p_S29
       -rho*Y[jhe4]*rate_eval.He4_Cl30__K34
       -rho*Y[jhe4]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jcl30, jar33] = (
       +rho*Y[jp]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jcl30, jk34] = (
       +rate_eval.K34__He4_Cl30
       )

    jac[jar32, jp] = (
       -rho*Y[jar32]*rate_eval.p_Ar32__K33
       -rho*Y[jar32]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jar32, jhe4] = (
       -rho*Y[jar32]*rate_eval.He4_Ar32__Ca36
       +rho*Y[js28]*rate_eval.He4_S28__Ar32
       +rho*Y[jcl29]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jar32, js28] = (
       +rho*Y[jhe4]*rate_eval.He4_S28__Ar32
       )

    jac[jar32, jcl29] = (
       +rho*Y[jhe4]*rate_eval.He4_Cl29__p_Ar32
       )

    jac[jar32, jar32] = (
       -rate_eval.Ar32__He4_S28
       -rho*Y[jp]*rate_eval.p_Ar32__K33
       -rho*Y[jhe4]*rate_eval.He4_Ar32__Ca36
       -rho*Y[jp]*rate_eval.p_Ar32__He4_Cl29
       )

    jac[jar32, jk33] = (
       +rate_eval.K33__p_Ar32
       )

    jac[jar32, jca36] = (
       +rate_eval.Ca36__He4_Ar32
       )

    jac[jar33, jp] = (
       -rho*Y[jar33]*rate_eval.p_Ar33__K34
       -rho*Y[jar33]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jar33, jhe4] = (
       -rho*Y[jar33]*rate_eval.He4_Ar33__Ca37
       +rho*Y[js29]*rate_eval.He4_S29__Ar33
       +rho*Y[jcl30]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jar33, js29] = (
       +rho*Y[jhe4]*rate_eval.He4_S29__Ar33
       )

    jac[jar33, jcl30] = (
       +rho*Y[jhe4]*rate_eval.He4_Cl30__p_Ar33
       )

    jac[jar33, jar33] = (
       -rate_eval.Ar33__He4_S29
       -rho*Y[jp]*rate_eval.p_Ar33__K34
       -rho*Y[jhe4]*rate_eval.He4_Ar33__Ca37
       -rho*Y[jp]*rate_eval.p_Ar33__He4_Cl30
       )

    jac[jar33, jk33] = (
       +rate_eval.K33__Ar33__weak__bqa_pos_
       )

    jac[jar33, jk34] = (
       +rate_eval.K34__p_Ar33
       )

    jac[jar33, jca37] = (
       +rate_eval.Ca37__He4_Ar33
       )

    jac[jk33, jp] = (
       +rho*Y[jar32]*rate_eval.p_Ar32__K33
       +rho*Y[jca36]*rate_eval.p_Ca36__He4_K33
       )

    jac[jk33, jhe4] = (
       -rho*Y[jk33]*rate_eval.He4_K33__Sc37
       -rho*Y[jk33]*rate_eval.He4_K33__p_Ca36
       +rho*Y[jcl29]*rate_eval.He4_Cl29__K33
       )

    jac[jk33, jcl29] = (
       +rho*Y[jhe4]*rate_eval.He4_Cl29__K33
       )

    jac[jk33, jar32] = (
       +rho*Y[jp]*rate_eval.p_Ar32__K33
       )

    jac[jk33, jk33] = (
       -rate_eval.K33__Ar33__weak__bqa_pos_
       -rate_eval.K33__p_Ar32
       -rate_eval.K33__He4_Cl29
       -rho*Y[jhe4]*rate_eval.He4_K33__Sc37
       -rho*Y[jhe4]*rate_eval.He4_K33__p_Ca36
       )

    jac[jk33, jca36] = (
       +rho*Y[jp]*rate_eval.p_Ca36__He4_K33
       )

    jac[jk33, jsc37] = (
       +rate_eval.Sc37__He4_K33
       )

    jac[jk34, jp] = (
       -rho*Y[jk34]*rate_eval.p_K34__Ca35
       +rho*Y[jar33]*rate_eval.p_Ar33__K34
       +rho*Y[jca37]*rate_eval.p_Ca37__He4_K34
       )

    jac[jk34, jhe4] = (
       -rho*Y[jk34]*rate_eval.He4_K34__Sc38
       -rho*Y[jk34]*rate_eval.He4_K34__p_Ca37
       +rho*Y[jcl30]*rate_eval.He4_Cl30__K34
       )

    jac[jk34, jcl30] = (
       +rho*Y[jhe4]*rate_eval.He4_Cl30__K34
       )

    jac[jk34, jar33] = (
       +rho*Y[jp]*rate_eval.p_Ar33__K34
       )

    jac[jk34, jk34] = (
       -rate_eval.K34__p_Ar33
       -rate_eval.K34__He4_Cl30
       -rho*Y[jp]*rate_eval.p_K34__Ca35
       -rho*Y[jhe4]*rate_eval.He4_K34__Sc38
       -rho*Y[jhe4]*rate_eval.He4_K34__p_Ca37
       )

    jac[jk34, jca35] = (
       +rate_eval.Ca35__p_K34
       )

    jac[jk34, jca37] = (
       +rho*Y[jp]*rate_eval.p_Ca37__He4_K34
       )

    jac[jk34, jsc38] = (
       +rate_eval.Sc38__He4_K34
       )

    jac[jca35, jp] = (
       -rho*Y[jca35]*rate_eval.p_Ca35__Sc36
       +rho*Y[jk34]*rate_eval.p_K34__Ca35
       +rho*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       )

    jac[jca35, jhe4] = (
       -rho*Y[jca35]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jca35, jk34] = (
       +rho*Y[jp]*rate_eval.p_K34__Ca35
       )

    jac[jca35, jca35] = (
       -rate_eval.Ca35__p_K34
       -rho*Y[jp]*rate_eval.p_Ca35__Sc36
       -rho*Y[jhe4]*rate_eval.He4_Ca35__Ti39
       -rho*Y[jhe4]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jca35, jsc36] = (
       +rate_eval.Sc36__p_Ca35
       )

    jac[jca35, jsc38] = (
       +rho*Y[jp]*rate_eval.p_Sc38__He4_Ca35
       )

    jac[jca35, jti39] = (
       +rate_eval.Ti39__He4_Ca35
       )

    jac[jca36, jp] = (
       -rho*Y[jca36]*rate_eval.p_Ca36__Sc37
       -rho*Y[jca36]*rate_eval.p_Ca36__He4_K33
       +5.00000000000000e-01*rho**2*2*Y[jp]*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       +rho*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       )

    jac[jca36, jhe4] = (
       -rho*Y[jca36]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       +rho*Y[jar32]*rate_eval.He4_Ar32__Ca36
       +rho*Y[jk33]*rate_eval.He4_K33__p_Ca36
       )

    jac[jca36, jar32] = (
       +rho*Y[jhe4]*rate_eval.He4_Ar32__Ca36
       )

    jac[jca36, jk33] = (
       +rho*Y[jhe4]*rate_eval.He4_K33__p_Ca36
       )

    jac[jca36, jca36] = (
       -rate_eval.Ca36__He4_Ar32
       -rho*Y[jp]*rate_eval.p_Ca36__Sc37
       -rho*Y[jhe4]*rate_eval.He4_Ca36__Ti40
       -rho*Y[jp]*rate_eval.p_Ca36__He4_K33
       -rho*Y[jhe4]*rate_eval.He4_Ca36__p_Sc39
       -rho*Y[jhe4]*rate_eval.He4_Ca36__p_p_Ca38
       )

    jac[jca36, jca38] = (
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_Ca38__He4_Ca36
       )

    jac[jca36, jsc36] = (
       +rate_eval.Sc36__Ca36__weak__bqa_pos_
       )

    jac[jca36, jsc37] = (
       +rate_eval.Sc37__p_Ca36
       )

    jac[jca36, jsc39] = (
       +rho*Y[jp]*rate_eval.p_Sc39__He4_Ca36
       )

    jac[jca36, jti40] = (
       +rate_eval.Ti40__He4_Ca36
       )

    jac[jca37, jp] = (
       -rho*Y[jca37]*rate_eval.p_Ca37__Sc38
       -rho*Y[jca37]*rate_eval.p_Ca37__He4_K34
       )

    jac[jca37, jhe4] = (
       +rho*Y[jar33]*rate_eval.He4_Ar33__Ca37
       +rho*Y[jk34]*rate_eval.He4_K34__p_Ca37
       )

    jac[jca37, jar33] = (
       +rho*Y[jhe4]*rate_eval.He4_Ar33__Ca37
       )

    jac[jca37, jk34] = (
       +rho*Y[jhe4]*rate_eval.He4_K34__p_Ca37
       )

    jac[jca37, jca37] = (
       -rate_eval.Ca37__He4_Ar33
       -rho*Y[jp]*rate_eval.p_Ca37__Sc38
       -rho*Y[jp]*rate_eval.p_Ca37__He4_K34
       )

    jac[jca37, jsc37] = (
       +rate_eval.Sc37__Ca37__weak__bqa_pos_
       )

    jac[jca37, jsc38] = (
       +rate_eval.Sc38__p_Ca37
       )

    jac[jca38, jp] = (
       -rho*Y[jca38]*rate_eval.p_Ca38__Sc39
       -5.00000000000000e-01*rho**2*2*Y[jp]*Y[jca38]*rate_eval.p_p_Ca38__He4_Ca36
       )

    jac[jca38, jhe4] = (
       +rho*Y[jca36]*rate_eval.He4_Ca36__p_p_Ca38
       )

    jac[jca38, jca36] = (
       +rho*Y[jhe4]*rate_eval.He4_Ca36__p_p_Ca38
       )

    jac[jca38, jca38] = (
       -rho*Y[jp]*rate_eval.p_Ca38__Sc39
       -5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.p_p_Ca38__He4_Ca36
       )

    jac[jca38, jsc38] = (
       +rate_eval.Sc38__Ca38__weak__mo97
       )

    jac[jca38, jsc39] = (
       +rate_eval.Sc39__p_Ca38
       )

    jac[jca38, jti39] = (
       +rate_eval.Ti39__p_Ca38__weak__wc12
       )

    jac[jsc36, jp] = (
       +rho*Y[jca35]*rate_eval.p_Ca35__Sc36
       +rho*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       )

    jac[jsc36, jhe4] = (
       -rho*Y[jsc36]*rate_eval.He4_Sc36__V40
       -rho*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jsc36, jca35] = (
       +rho*Y[jp]*rate_eval.p_Ca35__Sc36
       )

    jac[jsc36, jsc36] = (
       -rate_eval.Sc36__Ca36__weak__bqa_pos_
       -rate_eval.Sc36__p_Ca35
       -rho*Y[jhe4]*rate_eval.He4_Sc36__V40
       -rho*Y[jhe4]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jsc36, jti39] = (
       +rho*Y[jp]*rate_eval.p_Ti39__He4_Sc36
       )

    jac[jsc36, jv40] = (
       +rate_eval.V40__He4_Sc36
       )

    jac[jsc37, jp] = (
       +rho*Y[jca36]*rate_eval.p_Ca36__Sc37
       +rho*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       )

    jac[jsc37, jhe4] = (
       -rho*Y[jsc37]*rate_eval.He4_Sc37__V41
       -rho*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       +rho*Y[jk33]*rate_eval.He4_K33__Sc37
       )

    jac[jsc37, jk33] = (
       +rho*Y[jhe4]*rate_eval.He4_K33__Sc37
       )

    jac[jsc37, jca36] = (
       +rho*Y[jp]*rate_eval.p_Ca36__Sc37
       )

    jac[jsc37, jsc37] = (
       -rate_eval.Sc37__Ca37__weak__bqa_pos_
       -rate_eval.Sc37__p_Ca36
       -rate_eval.Sc37__He4_K33
       -rho*Y[jhe4]*rate_eval.He4_Sc37__V41
       -rho*Y[jhe4]*rate_eval.He4_Sc37__p_Ti40
       )

    jac[jsc37, jti40] = (
       +rho*Y[jp]*rate_eval.p_Ti40__He4_Sc37
       )

    jac[jsc37, jv41] = (
       +rate_eval.V41__He4_Sc37
       )

    jac[jsc38, jp] = (
       -rho*Y[jsc38]*rate_eval.p_Sc38__Ti39
       -rho*Y[jsc38]*rate_eval.p_Sc38__He4_Ca35
       +rho*Y[jca37]*rate_eval.p_Ca37__Sc38
       )

    jac[jsc38, jhe4] = (
       +rho*Y[jk34]*rate_eval.He4_K34__Sc38
       +rho*Y[jca35]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jsc38, jk34] = (
       +rho*Y[jhe4]*rate_eval.He4_K34__Sc38
       )

    jac[jsc38, jca35] = (
       +rho*Y[jhe4]*rate_eval.He4_Ca35__p_Sc38
       )

    jac[jsc38, jca37] = (
       +rho*Y[jp]*rate_eval.p_Ca37__Sc38
       )

    jac[jsc38, jsc38] = (
       -rate_eval.Sc38__Ca38__weak__mo97
       -rate_eval.Sc38__p_Ca37
       -rate_eval.Sc38__He4_K34
       -rho*Y[jp]*rate_eval.p_Sc38__Ti39
       -rho*Y[jp]*rate_eval.p_Sc38__He4_Ca35
       )

    jac[jsc38, jti39] = (
       +rate_eval.Ti39__p_Sc38
       )

    jac[jsc39, jp] = (
       -rho*Y[jsc39]*rate_eval.p_Sc39__Ti40
       -rho*Y[jsc39]*rate_eval.p_Sc39__He4_Ca36
       +rho*Y[jca38]*rate_eval.p_Ca38__Sc39
       )

    jac[jsc39, jhe4] = (
       +rho*Y[jca36]*rate_eval.He4_Ca36__p_Sc39
       )

    jac[jsc39, jca36] = (
       +rho*Y[jhe4]*rate_eval.He4_Ca36__p_Sc39
       )

    jac[jsc39, jca38] = (
       +rho*Y[jp]*rate_eval.p_Ca38__Sc39
       )

    jac[jsc39, jsc39] = (
       -rate_eval.Sc39__p_Ca38
       -rho*Y[jp]*rate_eval.p_Sc39__Ti40
       -rho*Y[jp]*rate_eval.p_Sc39__He4_Ca36
       )

    jac[jsc39, jti39] = (
       +rate_eval.Ti39__Sc39__weak__wc17
       )

    jac[jsc39, jti40] = (
       +rate_eval.Ti40__p_Sc39
       )

    jac[jti39, jp] = (
       -rho*Y[jti39]*rate_eval.p_Ti39__V40
       -rho*Y[jti39]*rate_eval.p_Ti39__He4_Sc36
       +rho*Y[jsc38]*rate_eval.p_Sc38__Ti39
       )

    jac[jti39, jhe4] = (
       -rho*Y[jti39]*rate_eval.He4_Ti39__Cr43
       +rho*Y[jca35]*rate_eval.He4_Ca35__Ti39
       +rho*Y[jsc36]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jti39, jca35] = (
       +rho*Y[jhe4]*rate_eval.He4_Ca35__Ti39
       )

    jac[jti39, jsc36] = (
       +rho*Y[jhe4]*rate_eval.He4_Sc36__p_Ti39
       )

    jac[jti39, jsc38] = (
       +rho*Y[jp]*rate_eval.p_Sc38__Ti39
       )

    jac[jti39, jti39] = (
       -rate_eval.Ti39__Sc39__weak__wc17
       -rate_eval.Ti39__p_Sc38
       -rate_eval.Ti39__p_Ca38__weak__wc12
       -rate_eval.Ti39__He4_Ca35
       -rho*Y[jp]*rate_eval.p_Ti39__V40
       -rho*Y[jhe4]*rate_eval.He4_Ti39__Cr43
       -rho*Y[jp]*rate_eval.p_Ti39__He4_Sc36
       )

    jac[jti39, jv40] = (
       +rate_eval.V40__p_Ti39
       )

    jac[jti39, jcr43] = (
       +rate_eval.Cr43__He4_Ti39
       )

    jac[jti40, jp] = (
       -rho*Y[jti40]*rate_eval.p_Ti40__V41
       -rho*Y[jti40]*rate_eval.p_Ti40__He4_Sc37
       +rho*Y[jsc39]*rate_eval.p_Sc39__Ti40
       )

    jac[jti40, jhe4] = (
       -rho*Y[jti40]*rate_eval.He4_Ti40__Cr44
       +rho*Y[jca36]*rate_eval.He4_Ca36__Ti40
       +rho*Y[jsc37]*rate_eval.He4_Sc37__p_Ti40
       )

    jac[jti40, jca36] = (
       +rho*Y[jhe4]*rate_eval.He4_Ca36__Ti40
       )

    jac[jti40, jsc37] = (
       +rho*Y[jhe4]*rate_eval.He4_Sc37__p_Ti40
       )

    jac[jti40, jsc39] = (
       +rho*Y[jp]*rate_eval.p_Sc39__Ti40
       )

    jac[jti40, jti40] = (
       -rate_eval.Ti40__p_Sc39
       -rate_eval.Ti40__He4_Ca36
       -rho*Y[jp]*rate_eval.p_Ti40__V41
       -rho*Y[jhe4]*rate_eval.He4_Ti40__Cr44
       -rho*Y[jp]*rate_eval.p_Ti40__He4_Sc37
       )

    jac[jti40, jv40] = (
       +rate_eval.V40__Ti40__weak__bqa_pos_
       )

    jac[jti40, jv41] = (
       +rate_eval.V41__p_Ti40
       )

    jac[jti40, jcr44] = (
       +rate_eval.Cr44__He4_Ti40
       )

    jac[jv40, jp] = (
       +rho*Y[jti39]*rate_eval.p_Ti39__V40
       +rho*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       )

    jac[jv40, jhe4] = (
       -rho*Y[jv40]*rate_eval.He4_V40__Mn44
       -rho*Y[jv40]*rate_eval.He4_V40__p_Cr43
       +rho*Y[jsc36]*rate_eval.He4_Sc36__V40
       )

    jac[jv40, jsc36] = (
       +rho*Y[jhe4]*rate_eval.He4_Sc36__V40
       )

    jac[jv40, jti39] = (
       +rho*Y[jp]*rate_eval.p_Ti39__V40
       )

    jac[jv40, jv40] = (
       -rate_eval.V40__Ti40__weak__bqa_pos_
       -rate_eval.V40__p_Ti39
       -rate_eval.V40__He4_Sc36
       -rho*Y[jhe4]*rate_eval.He4_V40__Mn44
       -rho*Y[jhe4]*rate_eval.He4_V40__p_Cr43
       )

    jac[jv40, jcr43] = (
       +rho*Y[jp]*rate_eval.p_Cr43__He4_V40
       )

    jac[jv40, jmn44] = (
       +rate_eval.Mn44__He4_V40
       )

    jac[jv41, jp] = (
       +rho*Y[jti40]*rate_eval.p_Ti40__V41
       +rho*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       )

    jac[jv41, jhe4] = (
       -rho*Y[jv41]*rate_eval.He4_V41__Mn45
       -rho*Y[jv41]*rate_eval.He4_V41__p_Cr44
       +rho*Y[jsc37]*rate_eval.He4_Sc37__V41
       )

    jac[jv41, jsc37] = (
       +rho*Y[jhe4]*rate_eval.He4_Sc37__V41
       )

    jac[jv41, jti40] = (
       +rho*Y[jp]*rate_eval.p_Ti40__V41
       )

    jac[jv41, jv41] = (
       -rate_eval.V41__p_Ti40
       -rate_eval.V41__He4_Sc37
       -rho*Y[jhe4]*rate_eval.He4_V41__Mn45
       -rho*Y[jhe4]*rate_eval.He4_V41__p_Cr44
       )

    jac[jv41, jcr44] = (
       +rho*Y[jp]*rate_eval.p_Cr44__He4_V41
       )

    jac[jv41, jmn45] = (
       +rate_eval.Mn45__He4_V41
       )

    jac[jcr43, jp] = (
       -rho*Y[jcr43]*rate_eval.p_Cr43__Mn44
       -rho*Y[jcr43]*rate_eval.p_Cr43__He4_V40
       )

    jac[jcr43, jhe4] = (
       -rho*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       +rho*Y[jti39]*rate_eval.He4_Ti39__Cr43
       +rho*Y[jv40]*rate_eval.He4_V40__p_Cr43
       )

    jac[jcr43, jti39] = (
       +rho*Y[jhe4]*rate_eval.He4_Ti39__Cr43
       )

    jac[jcr43, jv40] = (
       +rho*Y[jhe4]*rate_eval.He4_V40__p_Cr43
       )

    jac[jcr43, jcr43] = (
       -rate_eval.Cr43__He4_Ti39
       -rho*Y[jp]*rate_eval.p_Cr43__Mn44
       -rho*Y[jhe4]*rate_eval.He4_Cr43__Fe47
       -rho*Y[jp]*rate_eval.p_Cr43__He4_V40
       )

    jac[jcr43, jmn44] = (
       +rate_eval.Mn44__p_Cr43
       )

    jac[jcr43, jfe47] = (
       +rate_eval.Fe47__He4_Cr43
       )

    jac[jcr44, jp] = (
       -rho*Y[jcr44]*rate_eval.p_Cr44__Mn45
       -rho*Y[jcr44]*rate_eval.p_Cr44__He4_V41
       )

    jac[jcr44, jhe4] = (
       -rho*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       +rho*Y[jti40]*rate_eval.He4_Ti40__Cr44
       +rho*Y[jv41]*rate_eval.He4_V41__p_Cr44
       )

    jac[jcr44, jti40] = (
       +rho*Y[jhe4]*rate_eval.He4_Ti40__Cr44
       )

    jac[jcr44, jv41] = (
       +rho*Y[jhe4]*rate_eval.He4_V41__p_Cr44
       )

    jac[jcr44, jcr44] = (
       -rate_eval.Cr44__He4_Ti40
       -rho*Y[jp]*rate_eval.p_Cr44__Mn45
       -rho*Y[jhe4]*rate_eval.He4_Cr44__Fe48
       -rho*Y[jp]*rate_eval.p_Cr44__He4_V41
       )

    jac[jcr44, jmn44] = (
       +rate_eval.Mn44__Cr44__weak__bqa_pos_
       )

    jac[jcr44, jmn45] = (
       +rate_eval.Mn45__p_Cr44
       )

    jac[jcr44, jfe45] = (
       +rate_eval.Fe45__p_Cr44__weak__wc17
       )

    jac[jcr44, jfe48] = (
       +rate_eval.Fe48__He4_Cr44
       )

    jac[jmn44, jp] = (
       -rho*Y[jmn44]*rate_eval.p_Mn44__Fe45
       +rho*Y[jcr43]*rate_eval.p_Cr43__Mn44
       +rho*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jmn44, jhe4] = (
       -rho*Y[jmn44]*rate_eval.He4_Mn44__Co48
       -rho*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       +rho*Y[jv40]*rate_eval.He4_V40__Mn44
       )

    jac[jmn44, jv40] = (
       +rho*Y[jhe4]*rate_eval.He4_V40__Mn44
       )

    jac[jmn44, jcr43] = (
       +rho*Y[jp]*rate_eval.p_Cr43__Mn44
       )

    jac[jmn44, jmn44] = (
       -rate_eval.Mn44__Cr44__weak__bqa_pos_
       -rate_eval.Mn44__p_Cr43
       -rate_eval.Mn44__He4_V40
       -rho*Y[jp]*rate_eval.p_Mn44__Fe45
       -rho*Y[jhe4]*rate_eval.He4_Mn44__Co48
       -rho*Y[jhe4]*rate_eval.He4_Mn44__p_Fe47
       )

    jac[jmn44, jfe45] = (
       +rate_eval.Fe45__p_Mn44
       )

    jac[jmn44, jfe47] = (
       +rho*Y[jp]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jmn44, jco48] = (
       +rate_eval.Co48__He4_Mn44
       )

    jac[jmn45, jp] = (
       -rho*Y[jmn45]*rate_eval.p_Mn45__Fe46
       +rho*Y[jcr44]*rate_eval.p_Cr44__Mn45
       +rho*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jmn45, jhe4] = (
       -rho*Y[jmn45]*rate_eval.He4_Mn45__Co49
       -rho*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       +rho*Y[jv41]*rate_eval.He4_V41__Mn45
       )

    jac[jmn45, jv41] = (
       +rho*Y[jhe4]*rate_eval.He4_V41__Mn45
       )

    jac[jmn45, jcr44] = (
       +rho*Y[jp]*rate_eval.p_Cr44__Mn45
       )

    jac[jmn45, jmn45] = (
       -rate_eval.Mn45__p_Cr44
       -rate_eval.Mn45__He4_V41
       -rho*Y[jp]*rate_eval.p_Mn45__Fe46
       -rho*Y[jhe4]*rate_eval.He4_Mn45__Co49
       -rho*Y[jhe4]*rate_eval.He4_Mn45__p_Fe48
       )

    jac[jmn45, jfe45] = (
       +rate_eval.Fe45__Mn45__weak__wc17
       )

    jac[jmn45, jfe46] = (
       +rate_eval.Fe46__p_Mn45
       )

    jac[jmn45, jfe48] = (
       +rho*Y[jp]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jmn45, jco49] = (
       +rate_eval.Co49__He4_Mn45
       )

    jac[jfe45, jp] = (
       +rho*Y[jmn44]*rate_eval.p_Mn44__Fe45
       +rho*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jfe45, jhe4] = (
       -rho*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jfe45, jmn44] = (
       +rho*Y[jp]*rate_eval.p_Mn44__Fe45
       )

    jac[jfe45, jfe45] = (
       -rate_eval.Fe45__Mn45__weak__wc17
       -rate_eval.Fe45__p_Mn44
       -rate_eval.Fe45__p_Cr44__weak__wc17
       -rho*Y[jhe4]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jfe45, jco48] = (
       +rho*Y[jp]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jfe46, jp] = (
       -rho*Y[jfe46]*rate_eval.p_Fe46__Co47
       +rho*Y[jmn45]*rate_eval.p_Mn45__Fe46
       +rho*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       )

    jac[jfe46, jhe4] = (
       -rho*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jfe46, jmn45] = (
       +rho*Y[jp]*rate_eval.p_Mn45__Fe46
       )

    jac[jfe46, jfe46] = (
       -rate_eval.Fe46__p_Mn45
       -rho*Y[jp]*rate_eval.p_Fe46__Co47
       -rho*Y[jhe4]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jfe46, jco47] = (
       +rate_eval.Co47__p_Fe46
       )

    jac[jfe46, jco49] = (
       +rho*Y[jp]*rate_eval.p_Co49__He4_Fe46
       )

    jac[jfe47, jp] = (
       -rho*Y[jfe47]*rate_eval.p_Fe47__Co48
       -rho*Y[jfe47]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jfe47, jhe4] = (
       +rho*Y[jcr43]*rate_eval.He4_Cr43__Fe47
       +rho*Y[jmn44]*rate_eval.He4_Mn44__p_Fe47
       )

    jac[jfe47, jcr43] = (
       +rho*Y[jhe4]*rate_eval.He4_Cr43__Fe47
       )

    jac[jfe47, jmn44] = (
       +rho*Y[jhe4]*rate_eval.He4_Mn44__p_Fe47
       )

    jac[jfe47, jfe47] = (
       -rate_eval.Fe47__He4_Cr43
       -rho*Y[jp]*rate_eval.p_Fe47__Co48
       -rho*Y[jp]*rate_eval.p_Fe47__He4_Mn44
       )

    jac[jfe47, jco47] = (
       +rate_eval.Co47__Fe47__weak__bqa_pos_
       )

    jac[jfe47, jco48] = (
       +rate_eval.Co48__p_Fe47
       )

    jac[jfe48, jp] = (
       -rho*Y[jfe48]*rate_eval.p_Fe48__Co49
       -rho*Y[jfe48]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jfe48, jhe4] = (
       +rho*Y[jcr44]*rate_eval.He4_Cr44__Fe48
       +rho*Y[jmn45]*rate_eval.He4_Mn45__p_Fe48
       )

    jac[jfe48, jcr44] = (
       +rho*Y[jhe4]*rate_eval.He4_Cr44__Fe48
       )

    jac[jfe48, jmn45] = (
       +rho*Y[jhe4]*rate_eval.He4_Mn45__p_Fe48
       )

    jac[jfe48, jfe48] = (
       -rate_eval.Fe48__He4_Cr44
       -rho*Y[jp]*rate_eval.p_Fe48__Co49
       -rho*Y[jp]*rate_eval.p_Fe48__He4_Mn45
       )

    jac[jfe48, jco48] = (
       +rate_eval.Co48__Fe48__weak__bqa_pos_
       )

    jac[jfe48, jco49] = (
       +rate_eval.Co49__p_Fe48
       )

    jac[jco47, jp] = (
       -rho*Y[jco47]*rate_eval.p_Co47__Ni48
       +rho*Y[jfe46]*rate_eval.p_Fe46__Co47
       )

    jac[jco47, jfe46] = (
       +rho*Y[jp]*rate_eval.p_Fe46__Co47
       )

    jac[jco47, jco47] = (
       -rate_eval.Co47__Fe47__weak__bqa_pos_
       -rate_eval.Co47__p_Fe46
       -rho*Y[jp]*rate_eval.p_Co47__Ni48
       )

    jac[jco47, jni48] = (
       +rate_eval.Ni48__p_Co47
       )

    jac[jco48, jp] = (
       -rho*Y[jco48]*rate_eval.p_Co48__He4_Fe45
       +rho*Y[jfe47]*rate_eval.p_Fe47__Co48
       )

    jac[jco48, jhe4] = (
       +rho*Y[jmn44]*rate_eval.He4_Mn44__Co48
       +rho*Y[jfe45]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jco48, jmn44] = (
       +rho*Y[jhe4]*rate_eval.He4_Mn44__Co48
       )

    jac[jco48, jfe45] = (
       +rho*Y[jhe4]*rate_eval.He4_Fe45__p_Co48
       )

    jac[jco48, jfe47] = (
       +rho*Y[jp]*rate_eval.p_Fe47__Co48
       )

    jac[jco48, jco48] = (
       -rate_eval.Co48__Fe48__weak__bqa_pos_
       -rate_eval.Co48__p_Fe47
       -rate_eval.Co48__He4_Mn44
       -rho*Y[jp]*rate_eval.p_Co48__He4_Fe45
       )

    jac[jco48, jni48] = (
       +rate_eval.Ni48__Co48__weak__wc17
       )

    jac[jco49, jp] = (
       -rho*Y[jco49]*rate_eval.p_Co49__He4_Fe46
       +rho*Y[jfe48]*rate_eval.p_Fe48__Co49
       )

    jac[jco49, jhe4] = (
       +rho*Y[jmn45]*rate_eval.He4_Mn45__Co49
       +rho*Y[jfe46]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jco49, jmn45] = (
       +rho*Y[jhe4]*rate_eval.He4_Mn45__Co49
       )

    jac[jco49, jfe46] = (
       +rho*Y[jhe4]*rate_eval.He4_Fe46__p_Co49
       )

    jac[jco49, jfe48] = (
       +rho*Y[jp]*rate_eval.p_Fe48__Co49
       )

    jac[jco49, jco49] = (
       -rate_eval.Co49__p_Fe48
       -rate_eval.Co49__He4_Mn45
       -rho*Y[jp]*rate_eval.p_Co49__He4_Fe46
       )

    jac[jni48, jp] = (
       +rho*Y[jco47]*rate_eval.p_Co47__Ni48
       )

    jac[jni48, jco47] = (
       +rho*Y[jp]*rate_eval.p_Co47__Ni48
       )

    jac[jni48, jni48] = (
       -rate_eval.Ni48__Co48__weak__wc17
       -rate_eval.Ni48__p_Co47
       )

    return jac
