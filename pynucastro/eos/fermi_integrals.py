# In this particular module, we write the routines from Aparicio (1998), described in
# https://iopscience.iop.org/article/10.1086/313121. The method consist on breaking down
# the Fermi-integrals (eta, theta) integration domain into four subintervals and apply
# Gauss-Legendre or the Gauss laguerre quadrature methods with 20 points over each domain
# to guarantee double-precision.

# In the first subinterval, a x=z^2 change of variable is applied to overcome the kernel
# singularity near the origin.

import numpy as np
from scipy.special import roots_legendre, roots_laguerre


def quadrature_function(k:float, eta:float, theta:float, x:float, first:bool)-> float:

    f:float = 0.0

    if(first):
        z:float = x**2
        f = z**(2*k+1) * np.sqrt(1 + (z**2*theta/2))/(np.exp(z**2-eta) + 1)
    else:
        f = x**k * np.sqrt(1 + (x*theta/2))/(np.exp(x-eta) + 1)

    return f

def interval_points(eta:float)-> "np.array64":

    D:float =  3.3609
    sigma:float = 9.1186e-2

    xi:float = (1/sigma)*np.log(1+np.exp(sigma*(eta-D)))

    a1 = 6.774
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
    X_b = (a2 + b2*xi + c2*d2*xi**2)/(1 + e2*xi + c2*d2*xi**2)
    X_c = (a3 + b3*xi + c3*d3*xi**2)/(1 + e3*xi + c3*xi)

    points = [X_a, X_b, X_c]
    points = np.array(points, dtype=np.array64)

    return points

def compute_quadrature(k: float, eta:float, theta:float, int n)->float

    X_a, X_b, X_c = interval_points(eta)

    if(n = 0):
        a = 0.0
        b = X_a
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=false).dot(weights)

    elif(n = 1):
        a = X_a
        b = X_b
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=false).dot(weights)

    elif(n==2):
        a = X_b
        b = X_c
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=false).dot(weights)

    elif(n=3):
        roots, weights = roots_laguerre(20)
        scaled_roots = roots + X_c
        integral = np.exp(x) * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=false).dot(weights)

    return integral

def fermi_dirac(k, eta, theta):

    sum = 0.0

    for n in range(4):
        sum += compute_quadrature(k, eta, theta, n)

    return sum
