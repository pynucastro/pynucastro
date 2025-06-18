# In this particular module, we write the routines from Aparicio (1998) and Gong (2001),
# described in https://iopscience.iop.org/article/10.1086/313121. and
# https://www.sciencedirect.com/science/article/abs/pii/S001046550100145X?via%3Dihub.
# The method consist on breaking down the Fermi-integrals (eta, theta) integration domain
# into four subintervals and apply Gauss-Legendre or the Gauss laguerre quadrature methods
# with 20 points over each domain to guarantee accuracy beyond double-precision.

# In the first subinterval, a x=z^2 change of variable is applied to overcome the kernel
# singularity near the origin.

import numpy as np
from scipy.special import roots_laguerre, roots_legendre

class FermiIntegrals:

    def __init__(self, k, eta, beta, point_per_interval=50):
        self.eta = eta
        self.beta = beta
        self.k = k
        self.point_per_interval = point_per_interval

    def _kernel(self, x, first=False):

        eta = self.eta
        beta = self.beta
        k = self.k

        if first:
            xsq = x*x
            if ((xsq-eta) < 1.0e2):
                xdk = x**(2*k + 1)
                factor = np.sqrt(1 + (xsq*beta/2))
                denom = np.exp(xsq - eta) + 1
                denomi = 1/denom
                kernel = 2*xdk*factor*denomi
            else:
                xdk =  x**(2*k + 1)
                factor = np.sqrt(1 + (xsq*beta/2))
                denomi = np.exp(eta - xsq)
                kernel = 2*xdk*factor*denomi
        else:
            if ((x-eta) < 1.0e2):
                xdk = x**k
                factor = np.sqrt(1 + (x*beta/2))
                denom = np.exp(x - eta) + 1
                denomi = 1/denom
                kernel = xdk * factor * denomi
            else:
                xdk = x**k
                factor = np.sqrt(1 + (x*beta/2))
                denomi = np.exp(eta - x)
                kernel = xdk * factor * denomi

        return kernel

    def qfermi(self, x, first=False):
        kernel = self._kernel(x, first)
        weight = 1
        f = kernel * weight
        return f


    def qdfermi_deta(self, x, first=False):

        eta = self.eta
        kernel = self._kernel(x, first)

        if first:
            xsq = x*x
            if ((xsq-eta) < 1.0e2):
                factor = np.exp(xsq-eta)
                denom =  1 + factor
                weight = factor/denom
            else:
                weight = 1
        else:
            if ((x-eta) < 1.0e2):
                factor = np.exp(x-eta)
                denom =  1 + factor
                weight = factor/denom
            else:
                weight = 1
        f = kernel * weight
        return f


    def qdfermi_dbeta(self, x, first=False):

        beta = self.beta

        kernel = self._kernel(x, first)

        if first:
            xsq = x*x
            weight = xsq / (4 + 2*beta*xsq)
        else:
            weight = x / (4 + 2*beta*x)
        f = kernel * weight
        return f


    def fermi_points(self):

        eta = self.eta

        D = 3.3609
        sigma = 9.1186e-2

        xi = (1/sigma)*np.log(1+np.exp(sigma*(eta-D)))

        a1 = 6.7774
        b1 = 1.1418
        c1 = 2.9826

        a2 = 3.7601
        b2 = 9.3719e-2
        c2 = 2.1064e-2
        d2 = 3.1084e1
        e2 = 1.0056

        a3 = 7.5669
        b3 = 1.1695
        c3 = 7.5416e-1
        d3 = 6.6559
        e3 = -1.2819

        X_a = (a1 + b1*xi + c1*xi**2)/(1 + c1*xi)
        X_b = (a2 + b2*xi + c2*d2*xi**2)/(1 + e2*xi + c2*xi**2)
        X_c = (a3 + b3*xi + c3*d3*xi**2)/(1 + e3*xi + c3*xi**2)

        S_1 = X_a - X_b
        S_2 = X_a
        S_3 = X_a + X_c

        points = [S_1, S_2, S_3]
        points = np.array(points, dtype=np.double)

        return points


    def dfermi_points(self):

        eta = self.eta

        D = 4.99551
        sigma = 9.11856e-2

        xi = (1/sigma)*np.log(1+np.exp(sigma*(eta-D)))

        a1 = 6.77740
        b1 = 1.14180
        c1 = 2.98255

        a2 = 3.76010
        b2 = 9.37188e-2
        c2 = 2.10635e-2
        d2 = 3.95015e1
        e2 = 1.00557

        a3 = 7.56690
        b3 = 1.16953
        c3 = 7.54162
        d3 = 7.64734
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

    def _compute_legendre(self, function_fermi, a, b, first=False):

        point_per_interval = self.point_per_interval/2
        roots, weights = roots_legendre(point_per_interval)
        integral = 0

        for i in range(len(roots)):
            x_1 =  (a+b)/2 + (b-a)/2 * roots[i]
            x_2 =  (a+b)/2 - (b-a)/2 * roots[i]
            fval_1 = function_fermi(x_1, first)
            fval_2 = function_fermi(x_2, first)
            integral += (fval_1 + fval_2) * weights[i]

        integral *= (b-a)/2
        return integral

    def _compute_laguerre(self, function_fermi, a):

        point_per_interval = self.point_per_interval
        roots, weights = roots_laguerre(point_per_interval)
        integral = 0

        scaled_roots = roots + a
        for i, x in enumerate(scaled_roots):
            integral += np.exp(roots[i]) * function_fermi(x, first=False) * weights[i]

        return integral

    def compute_fermi(self, n, function_fermi, der):

        integral = 0

        if der:
            S_1, S_2, S_3 = self.dfermi_points()
        else:
            S_1, S_2, S_3 = self.fermi_points()

        if n == 0:
            a = 0.0
            b = np.sqrt(S_1)
            integral = self._compute_legendre(function_fermi, a, b, first=True)
        elif n == 1:
            a = S_1
            b = S_2
            integral = self._compute_legendre(function_fermi, a, b)
        elif n == 2:
            a = S_2
            b = S_3
            integral = self._compute_legendre(function_fermi, a, b)
        elif n == 3:
            a = S_3
            integral = self._compute_laguerre(function_fermi, a)

        return integral

    def fd(self):
        f = 0.0
        for n in range(4):
            f += self.compute_fermi(n, self.qfermi, der=False)
        return f

    def dfdeta(self):
        f = 0.0
        for n in range(4):
            f += self.compute_fermi(n, self.qdfermi_deta, der=True)
        return f

    def dfdbeta(self):
        f = 0.0
        for n in range(4):
            f += self.compute_fermi(n, self.qdfermi_dbeta, der=False)
        return f
