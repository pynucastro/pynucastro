# In this module file we construct N_{ele} and N_{pos} from Timmes & Arnett (1999) https://iopscience.iop.org/article/10.1086/313304,
# solve for eta from eq (7) using the reaction network "

from pynucastro.eos import fermi_dirac
import numpy as np
from pynucastro.constants import constants
from scipy import optimize

m_e = constants.m_e
m_u = constants.m_u
hbar = constants.hbar
c = constants.c_light
k = constants.k

def N_ele(eta, beta):
    integrand = 8*np.pi*np.sqrt(2)*(m_e*c*np.sqrt(beta))**3*(fermi_dirac(0.5, eta, beta) + fermi_dirac(1.5, eta, beta))/(hbar**3)
    return integrand

def N_pos(eta, beta):
    integrand = 8*np.pi*np.sqrt(2)*(m_e*c*np.sqrt(beta))**3*(fermi_dirac(0.5, -eta-2/beta, beta) + beta*fermi_dirac(1.5, -eta-2/beta, beta))/(hbar**3)
    return integrand

def newton_target_function(rho,T, eta):
    beta = k*T/(m_e*c**2)
    N_ele_matter = rho/m_u
    return N_ele(eta,beta) - N_pos(eta,beta) - N_ele_matter

def newton_eta(rho, T, initial_eta, max_iter):
    root = optimize.newton(func=lambda eta: newton_target_function(rho, T, eta), x0=initial_eta, tol=1.0e-15, maxiter=max_iter)
    return root

def brentq_eta(rho, T, low_eta, max_eta, max_iter):
    root = optimize.brentq(f=lambda eta: newton_target_function(rho, T, eta), a=low_eta, b=max_eta, xtol=1.0e-15, maxiter=max_iter)
    return root