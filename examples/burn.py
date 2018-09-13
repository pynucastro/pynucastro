from __future__ import print_function

import numpy as np
from scipy.integrate import ode
import network_module as network
import pynucastro
import matplotlib.pyplot as plt

def burn(Y0, rho, T, tmax, nsave):

    r = ode(network.rhs).set_integrator("vode", method="bdf",
                                    with_jacobian=False,
                                    atol=1.e-8, rtol=1.e-8,
                                    nsteps = 1500000, order=5) #, min_step=dt)

    t = 0.0

    r.set_initial_value(Y0, t)

    r.set_f_params(rho, T)

    dt = tmax/nsave

    t_out = []
    Y_out = [[] for _ in range(network.nnuc)]

    while r.successful() and r.t < tmax:
        print(r.t, r.y[:])
        r.integrate(r.t+dt)

        t_out.append(r.t)
        for i, yi in enumerate(r.y):
            Y_out[i].append(yi)

    return t_out, Y_out


if __name__ == "__main__":

    Y0 = np.zeros((network.nnuc), dtype=np.float64)

    # Initial conditions
    rho = 1.0e6
    T = 1.e9
    tmax = 1.0e-6
    nsteps = 10000
    Y0[network.ip] = 0.7
    Y0[network.ihe4] = 0.28
    Y0[network.ic12] = 0.02
    # End initial conditions
    
    Y0[:] = Y0[:]/network.A[:]

    t, Y_out = burn(Y0, rho, T, tmax, nsteps)

    for i in range(network.nnuc):
        plt.loglog(t, np.array(Y_out[i])*network.A[i], label=network.species_names[i])

    plt.xlabel("time [s]")
    plt.ylabel("mass fraction, X")

    plt.legend(frameon=False)

    plt.savefig("burn.png")
