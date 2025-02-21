# In this particular module, we write the routines from Aparicio (1998), described in
# https://iopscience.iop.org/article/10.1086/313121. The method consist on breaking down
# the Fermi-integrals (eta, theta) integration domain into four subintervals and apply
# Gauss-Legendre or the Gauss laguerre quadrature methods with 20 points over each domain
# to guarantee double-precision.

# In the first subinterval, a x=z^2 change of variable is applied to overcome the kernel
# singularity near the origin.

import numpy as np
from scipy.special import roots_legendre, roots_laguerre


def qfermi(k, eta, beta, x, first=False):

    if first:
        z = np.sqrt(x)
        f = 2*z**(2*k+1) * np.sqrt(1 + (z**2*beta/2))/(np.exp(z**2-eta) + 1) # There is a missing 2 factor in eq. 5
    else:
        f = x**k * np.sqrt(1 + (x*beta/2))/(np.exp(x-eta) + 1)

    return f

def qdfermi_deta(k, eta, beta, x, first=False):

    if first:
        z = np.sqrt(x)
        f = 2*z**(2*k+1) * np.sqrt(1 + (z**2*beta/2)) / ( (np.exp(z**2-eta) + 1) * (1+(np.exp(eta-z**2))) )
    else:
        f = x**k * np.sqrt(1 + (x*beta/2)) / ( (np.exp(x-eta) + 1) * (1 + np.exp(eta-x)) )

    return f

def qdfermi_dbeta(k, eta, beta, x, first=False):

    if first:
        z = np.sqrt(x)
        f = 2*z**(2*k+1) * np.sqrt(1 + (z**2*beta/2)) * z**2 / ( (np.exp(z**2-eta) + 1) * (4 + 2*beta*z**2) )
    else:
        f =  x**k * np.sqrt(1 + (x*beta/2)) * x / ( (np.exp(x-eta) + 1) * (4 + 2*beta*x) )

    return f


def fermi_points(eta):

    D =  3.3609
    sigma = 9.1186e-2

    xi = (1/sigma)*np.log(1+np.exp(sigma*(eta-D)))

    a1 = 6.7774
    b1 = 1.1418
    c1 = 2.9826

    a2 = 3.7601
    b2 = 9.3719e-2
    c2 = 2.1063e-2
    d2 = 3.1084e1
    e2 = 1.0056

    a3 = 7.5669
    b3 = 1.1695
    c3 = 7.5416e-1
    d3 = 6.6558
    e3 = -1.2819e-1

    X_a = (a1 + b1*xi + c1*xi**2)/(1 + c1*xi)
    X_b = (a2 + b2*xi + c2*d2*xi**2)/(1 + e2*xi + c2*xi**2)
    X_c = (a3 + b3*xi + c3*d3*xi**2)/(1 + e3*xi + c3*xi**2)

    S_1 = X_a - X_b
    S_2 = X_a
    S_3 = X_a + X_c


    points = [S_1, S_2, S_3]
    points = np.array(points, dtype=np.double)

    return points

def dfermi_points(eta):

    D = 4.99551e0
    sigma = 9.11856e-2

    xi = (1/sigma)*np.log(1+np.exp(sigma*(eta-D)))

    a1 = 6.77740e0
    b1 = 1.14180e0
    c1 = 2.98255e0

    a2 = 3.76010e0
    b2 = 9.37188e-2
    c2 = 2.10635e-2
    d2 = 3.95015e1
    e2 = 1.00557e0

    a3 = 7.56690e0
    b3 = 1.16953e0
    c3 = 7.54162e-1
    d3 = 7.64734e0
    e3 = -1.28190e-1

    X_a = (a1 + b1*xi + c1*xi**2)/(1 + c1*xi)
    X_b = (a2 + b2*xi + c2*d2*xi**2)/(1 + e2*xi + c2*xi**2)
    X_c = (a3 + b3*xi + c3*d3*xi**2)/(1 + e3*xi + c3*xi**2)

    S_1 = X_a - X_b
    S_2 = X_a
    S_3 = X_a + X_c

    points = [S_1, S_2, S_3]
    points = np.array(points, dtype=np.double)

    return points

def compute_fermi(k, eta, beta, n, function_fermi, der):

    point_per_interval = 50
    integral = 0

    if der:
        S_1, S_2, S_3 = dfermi_points(eta)
    else:
        S_1, S_2, S_3 = fermi_points(eta)

    if n == 0:
        a = 0.0
        b = np.sqrt(S_1)
        roots, weights, _ = roots_legendre(point_per_interval)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * function_fermi(k=k, eta=eta, beta=beta, x=scaled_roots, first=True).dot(weights)

    elif n == 1:
        a = S_1
        b = S_2
        roots, weights, _ = roots_legendre(point_per_interval)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * function_fermi(k=k, eta=eta, beta=beta, x=scaled_roots).dot(weights)

    elif n == 2:
        a = S_2
        b = S_3
        roots, weights, _ = roots_legendre(point_per_interval)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * function_fermi(k=k, eta=eta, beta=beta, x=scaled_roots).dot(weights)

    elif n == 3:
        roots, weights, _ = roots_laguerre(point_per_interval)
        scaled_roots = roots + S_3
        integral = (np.exp(scaled_roots) * function_fermi(k=k, eta=eta, beta=beta, x=scaled_roots)).dot(weights)

    return integral

def fermi(k, eta, beta):
    f = 0.0
    for n in range(4):
        f += compute_fermi(k, eta, beta, n, qfermi, der=False)
    return f

def dfermi_deta(k, eta, beta):
    f = 0.0
    for n in range(4):
        f += compute_fermi(k, eta, beta, n, qdfermi_deta, der=True)
    return f

def dfermi_dbeta(k, eta, beta):
    f = 0.0
    for n in range(4):
        f += compute_fermi(k, eta, beta, n, qdfermi_dbeta, der=True)

    return f
