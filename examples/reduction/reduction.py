from pynucastro import ReacLibLibrary, PythonNetwork
from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus
from pynucastro.reduction import mpi_importer, sens_analysis
from pynucastro.reduction.drgep import drgep
from pynucastro.reduction.generate_data import dataset
import numpy as np

def get_net_info(net, comp, rho, T):

    y_dict = comp.get_molar()
    # Can alternatively use the NumPy-based method (evaluate_ydots_arr)
    ydots_dict = net.evaluate_ydots(rho, T, comp)

    y = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    ydot = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    z = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    a = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    ebind = np.zeros(len(net.unique_nuclei), dtype=np.float64)
    m = np.zeros(len(net.unique_nuclei), dtype=np.float64)

    for i, n in enumerate(net.unique_nuclei):

        y[i] = y_dict[n]
        ydot[i] = ydots_dict[n]
        z[i] = n.Z
        a[i] = n.A
        ebind[i] = n.nucbind or 0.0
        m[i] = n.A_nuc * constants.m_u

    return (y, ydot, z, a, ebind, m)

def enuc_dot(net_info):
    """Calculate the nuclear energy generation rate."""

    return -np.sum(net_info[1] * net_info[5]) * constants.N_A * constants.c_light**2

def rel_err(x, x0):
    """Compute the relative error between two NumPy arrays."""

    return np.abs((x - x0) / x0)

def get_errfunc_enuc(net_old, conds):
    """Function for computing error in nuclear energy generation."""

    enucdot_list = []

    for comp, rho, T in conds:
        net_info_old = get_net_info(net_old, comp, rho, T)
        enucdot_list.append(enuc_dot(net_info_old))

    def erf(net_new):

        err = 0.0

        for cond, enucdot_old in zip(conds, enucdot_list):

            net_info_new = get_net_info(net_new, *cond)
            enucdot_new = enuc_dot(net_info_new)
            err = max(err, rel_err(enucdot_new, enucdot_old))

        return err

    return erf
