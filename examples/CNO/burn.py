import numpy as np
from scipy.integrate import ode
import cno_rhs as cno
import matplotlib.pyplot as plt

def burn(Y0, rho, T, tmax, nsave):

    r = ode(cno.rhs).set_integrator("vode", method="bdf",
                                    with_jacobian=False,
                                    atol=1.e-8, rtol=1.e-8,
                                    nsteps = 1500000, order=5) #, min_step=dt)

    t = 0.0

    r.set_initial_value(Y0, t)

    r.set_f_params(rho, T)

    dt = tmax/nsave

    t_out = []
    H_out = []
    He_out = []
    O14_out = []
    O15_out = []

    print(t, Y0[cno.ip], Y0[cno.io14])

    istep = 1
    while r.successful() and istep <= nsave:
        r.integrate(t+dt*istep)

        if r.successful():
            print(r.t, r.y[cno.ip], r.y[cno.io14])
            t_out.append(r.t)
            H_out.append(r.y[cno.ip])
            He_out.append(r.y[cno.ihe4])
            O14_out.append(r.y[cno.io14])
            O15_out.append(r.y[cno.io15])
            istep = istep + 1
        else:
            print("An integration error occurred at time {}".format(r.t))

    return t_out, H_out, He_out, O14_out, O15_out


if __name__ == "__main__":

    # initialize as mass fractions first
    Y0 = np.zeros((cno.nnuc), dtype=np.float64)

    Y0[cno.ip] = 0.7
    Y0[cno.ihe4] = 0.28
    Y0[cno.ic12] = 0.02

    Y0[:] = Y0[:]/cno.A[:]

    rho = 10000.0
    T = 1.e8

    # estimate the H destruction time
    Ydot = cno.rhs(0.0, Y0, rho, T)

    tmax = 10.0*np.abs(Y0[cno.ip]/Ydot[cno.ip])
    print("tmax: {}".format(tmax))

    nsteps = 100

    t, Y_H, Y_He, Y_O14, Y_O15 = burn(Y0, rho, T, tmax, nsteps)

    plt.loglog(t, np.array(Y_H)*cno.A[cno.ip], label="H1")
    plt.loglog(t, np.array(Y_He)*cno.A[cno.ihe4], label="He4")
    plt.loglog(t, np.array(Y_O14)*cno.A[cno.io14], label="O14")
    plt.loglog(t, np.array(Y_O15)*cno.A[cno.io15], label="O15")

    plt.xlabel("time [s]")
    plt.ylabel("mass fraction, X")

    plt.legend(frameon=False)

    plt.savefig("burn.png")
