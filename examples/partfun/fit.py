import numpy as np
from scipy.optimize import curve_fit
from pynucastro.nucdata import PartitionFunction, PartitionFunctionCollection
import matplotlib.pyplot as plt

pfc = PartitionFunctionCollection()

pf_ni56 = pfc.get_partition_function('ni56')

t = np.log10(pf_ni56.temperature)
f = np.log10(pf_ni56.partition_function)

def fit(T, a0, a1):
    f = a0 * T**(3./2.) * np.exp(-a1/T)
    return f

c = 1.0
popt, pcov = curve_fit(fit, t, f, p0=(c, c))
pstd = np.sqrt(np.diag(pcov))

for v, s in zip(popt, pstd):
    print('{} +/- {}'.format(v, s))

xt = np.linspace(t[0], t[-1], num=1000)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(t, f, 'r.', label='ni56')
ax1.plot(xt, fit(xt, *popt), 'b-', label='fit')
ax1.set_ylabel('Log F')
ax2.plot(t, (np.power(10.0,fit(t, *popt))-np.power(10.0,f))/np.power(10.0, f), 'g-', label='fit-actual')
ax2.set_xlabel('Log T (K)')
ax2.set_ylabel('Fit Relative Error')
plt.savefig('fit_partition_function_{}.png'.format('ni56'), dpi=600)
plt.clf()
