import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode
from scipy import constants
from pynucastro import Composition, Nucleus

import rpNi56_rhs as rpni56
import CNO_rhs as cno
import reduced_rhs as reduced

m_u = constants.value('atomic mass constant energy equivalent in MeV')
MeV2erg = (constants.eV * constants.mega) / constants.erg

def evaluate_energy(comp, net):
    E = 0
    for i, nuc in enumerate([Nucleus.from_cache(name) for name in net.names]):
        E += nuc.A_nuc * comp[i]
    E *= m_u * MeV2erg * constants.Avogadro 
    return E

def burn(Y0, rho, T, tmax, nsave, net):

    r = ode(net.rhs).set_integrator("vode", method="bdf",
                                    with_jacobian=False,
                                    atol=1.e-8, rtol=1.e-8,
                                    nsteps = 1500000, order=5) #, min_step=dt)

    t = 0.0

    r.set_initial_value(Y0, t)

    r.set_f_params(rho, T)

    dt = tmax/nsave

    t_out = []
    E_out = []

    print(t, Y0[net.jp], Y0[net.jo14], evaluate_energy(r.y, net))

    istep = 1
    while r.successful() and istep <= nsave:
        r.integrate(t+dt*istep)

        if r.successful():
            print(r.t, r.y[net.jp], r.y[net.jo14], evaluate_energy(r.y, net))
            t_out.append(r.t)
            E_out.append(evaluate_energy(r.y, net))
            istep = istep + 1
        else:
            print(f"An integration error occurred at time {r.t}")

    return t_out, E_out 

def get_init_M(net):
    Y0 = np.zeros((net.nnuc), dtype=np.float64)

    Y0[net.jp] = 0.7
    Y0[net.jhe4] = 0.28
    Y0[net.jc12] = 0.01
    Y0[net.jo16] = 0.01

    return Y0[:]/net.A[:]

def Eval_Burn(net, rho, T, tmax=None):

    # initialize as mass fractions first
    Y0 = get_init_M(net)

    # estimate the H destruction time
    if tmax is None:
        Ydot = net.rhs(0.0, Y0, rho, T)

        tmax = 10.0*np.abs(Y0[net.jp]/Ydot[net.jp])
        print(f"tmax: {tmax}")

    nsteps = 100

    return burn(Y0, rho, T, tmax, nsteps, net), tmax

#    rets = burn(Y0, rho, T, tmax, nsteps, net)
#    return rets(0), rets(1), tmax

if __name__ == "__main__":

    start_time = time.time()

    nets = [cno, reduced, rpni56]
    labels = ["CNO extra", "Reduced Network", "rp-process to Ni56 network"]

    rho = 10000.0
    T = 1.e8
    tmax = None

    fig, ax = plt.subplots()
    ax.set_xlabel("time (s)")
    ax.set_ylabel("Energy (erg/g)")

    for label, net in zip(labels, nets):
        #if i == 0:
        try:
            print(f"Starting {label}")
            (cnot, cnoE), tmax = Eval_Burn(net, rho, T, tmax)
            ax.plot(cnot, cnoE[0] - cnoE, label = label)
        except:
            print(f"Failed integration of {label}")

    ax.legend(frameon=False)

    plt.show()
    fig.savefig("burn.png")
    print(f"Execution time is {time.time() - start_time} s")
