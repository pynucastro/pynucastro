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

    if(first):
        z = np.sqrt(x)
        f = z**(2*k+1) * np.sqrt(1 + (z**2*theta/2))/(np.exp(z**2-eta) + 1)
    else:
        f = x**k * np.sqrt(1 + (x*theta/2))/(np.exp(x-eta) + 1)

    return f

def interval_points(eta:float)-> float:

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

def compute_quadrature(k: float, eta:float, theta:float, n:int)->float:

    S_1, S_2, S_3 = interval_points(eta)

    if(n == 0):
        a = 0.0
        b = np.sqrt(S_1)
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=True).dot(weights)

    elif(n == 1):
        a = S_1
        b = S_2
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=False).dot(weights)

    elif(n == 2):
        a = S_2
        b = S_3
        roots, weights = roots_legendre(20)
        scaled_roots = (b-a)/2 * roots + (a+b)/2
        integral = (b-a)/2 * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=False).dot(weights)

    elif(n == 3):
        roots, weights = roots_laguerre(20)
        scaled_roots = roots + S_3
        integral = (np.exp(scaled_roots) * quadrature_function(k=k, eta=eta, theta=theta, x=scaled_roots, first=False)).dot(weights)

    return integral

def fermi_dirac(k:float, eta:float, theta:float):

    sum = 0.0

    for n in range(4):
        sum += compute_quadrature(k, eta, theta, n)

    return sum
