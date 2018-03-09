import numpy as np
from scipy.optimize import curve_fit
from pynucastro.nucdata import PartitionFunction, PartitionFunctionCollection

pfc = PartitionFunctionCollection()

pf_ni56 = pfc.get_partition_function('ni56')

t = pf_ni56.temperature
f = pf_ni56.partition_function

print(t)
print(f)

def fit(T, ashift, a0, a1, a2, a3, a4, a5, a6):
    T = T - ashift
    T9 = T/1.e9
    T9i = 1.0/T9
    T913i = np.power(T9i,1./3.)
    T913 = np.power(T9,1./3.)
    T953 = np.power(T9,5./3.)
    lnT9 = np.log(T9)
    f = np.exp(a0 +
               a1*T9i +
               a2*T913i +
               a3*T913 +
               a4*T9 +
               a5*T953 +
               a6*lnT9)
    f = f + 1
    return f

c = t[-1]*np.log(f[-1])
popt, pcov = curve_fit(fit, t, f, p0=(t[0], c, c, c, c, c, c, c))
pstd = np.sqrt(np.diag(pcov))

for v, s in zip(popt, pstd):
    print('{} +/- {}'.format(v, s))
