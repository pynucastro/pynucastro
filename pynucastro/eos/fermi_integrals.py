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

    def _qderivatives(self, x, eta_der, beta_der, first=False):

        eta = self.eta
        beta = self.beta
        k = self.k

        if first:
            xsq = x*x
            xdk = x**(2*k + 1)
            if ((xsq-eta) < 2.0e2):
                sqrt_factor = np.sqrt(1 + (xsq*beta/2))
                exp_factor = np.exp(xsq - eta)
                denom = exp_factor + 1
                denomi = 1/denom
                kernel = 2*xdk*sqrt_factor*denomi
                denom2i =  1/(4 + 2*beta*xsq)

                fval = kernel
                df_deta = fval*exp_factor*denomi
                df_deta2 = (2*exp_factor*denomi - 1)*df_deta
                df_dbeta = fval*xsq*denom2i
                df_dbeta2 = -df_dbeta*xsq*denom2i
                df_deta_dbeta = df_dbeta*exp_factor*denomi
            else:
                sqrt_factor = np.sqrt(1 + (xsq*beta/2))
                denomi = np.exp(eta - xsq)
                kernel = 2*xdk*sqrt_factor*denomi
                denom2i =  1/(4 + 2*beta*xsq)

                fval = kernel
                df_deta = fval
                df_deta2 = fval
                df_dbeta = fval*xsq*denom2i
                df_dbeta2 = -df_dbeta*xsq*denom2i
                df_deta_dbeta = df_dbeta
        else:
            if ((x-eta) < 2.0e2):
                xdk = x**k
                sqrt_factor = np.sqrt(1 + (x*beta/2))
                exp_factor = np.exp(x - eta)
                denom = exp_factor + 1
                denomi = 1/denom
                kernel = xdk*sqrt_factor*denomi
                denom2i =  1/(4 + 2*beta*x)

                fval = kernel
                df_deta = fval*exp_factor*denomi
                df_deta2 = (2*exp_factor*denomi - 1)*df_deta
                df_dbeta = fval*x*denom2i
                df_dbeta2 = -df_dbeta*x*denom2i
                df_deta_dbeta = df_dbeta*exp_factor*denomi
            else:
                xdk = x**k
                sqrt_factor = np.sqrt(1 + (x*beta/2))
                denomi = np.exp(eta - x)
                kernel = xdk*sqrt_factor*denomi
                denom2i =  1/(4 + 2*beta*x)

                fval = kernel
                df_deta = fval
                df_deta2 = fval
                df_dbeta = fval*x*denom2i
                df_dbeta2 = -df_dbeta*x*denom2i
                df_deta_dbeta = df_dbeta

        if(eta_der == 0 and beta_der == 0):
            return fval
        elif(eta_der == 1 and beta_der == 0):
            return df_deta
        elif(eta_der == 2 and beta_der == 0):
            return df_deta2
        elif(eta_der == 0 and beta_der == 1):
            return df_dbeta
        elif(eta_der == 0 and beta_der == 2):
            return df_dbeta2
        elif(eta_der == 1 and beta_der == 1):
            return df_deta_dbeta
        else:
            raise ValueError("Invalid derivative order")

    def _fermi_points(self, eta_der=0):

        eta = self.eta

        D = [3.3609, 4.99551, 3.93830, 4.17444]
        sigma = [9.1186e-2, 9.11856e-2, 9.11856e-2, 9.11856e-2]

        xi = [(1/s)*np.log(1+np.exp(s*(eta-d))) for d, s in zip(D,sigma)]

        a1 = np.array([6.7774, 6.77740, 6.77740, 6.77740])
        b1 = np.array([1.1418, 1.14180, 1.14180, 1.14180])
        c1 = np.array([2.9826, 2.98255, 2.98255, 2.98255])

        a2 = np.array([3.7601, 3.76010, 3.76010, 3.76010])
        b2 = np.array([9.3719e-2, 9.37188e-2, 9.37188e-2, 9.37188e-2])
        c2 = np.array([2.1064e-2, 2.10635e-2, 2.10635e-2, 2.10635e-2])
        d2 = np.array([3.1084e1, 3.95015e1, 3.14499e1, 3.05412e1])
        e2 = np.array([1.0056, 1.00557, 1.00557, 1.00557])

        a3 = np.array([7.5669, 7.56690, 7.56690, 7.56690])
        b3 = np.array([1.1695, 1.16953, 1.16953, 1.16953])
        c3 = np.array([7.5416e-1, 7.54162, 7.54162, 7.54162])
        d3 = np.array([6.6559, 7.64734, 6.86346, 7.88030])
        e3 = np.array([-1.2819, -1.28190e-1, -1.28190e-1, -1.28190e-1])

        X_a = (a1 + b1*xi + c1*xi**2)/(1 + c1*xi)
        X_b = (a2 + b2*xi + c2*d2*xi**2)/(1 + e2*xi + c2*xi**2)
        X_c = (a3 + b3*xi + c3*d3*xi**2)/(1 + e3*xi + c3*xi**2)

        S_1 = X_a - X_b
        S_2 = X_a
        S_3 = X_a + X_c

        points = [S_1[eta_der], S_2[eta_der], S_3[eta_der]]
        points = np.array(points, dtype=np.double)

        return points

    def _compute_legendre(self, a, b, eta_der, beta_der, first=False):

        point_per_interval = self.point_per_interval//2
        roots, weights = roots_legendre(point_per_interval)
        integral = 0

        # We set the correspondence (a+b)/2 -> 0 and map the (-1,0) and (0,1)
        # intervals separately.

        for weight, root in zip(weights,roots):
            x_1 = (a+b)/2 + (b-a)/2 * root
            x_2 = (a+b)/2 - (b-a)/2 * root
            qder_1 = self._qderivatives(x_1, eta_der, beta_der, first)
            qder_2 = self._qderivatives(x_2, eta_der, beta_der, first)
            integral += (qder_1 + qder_2) * weight

        integral *= (b-a)/2
        return integral

    def _compute_laguerre(self, a, eta_der, beta_der):

        point_per_interval = self.point_per_interval
        roots, weights = roots_laguerre(point_per_interval)
        integral = 0

        scaled_roots = roots + a
        for sroot, root, weight in zip(scaled_roots, roots, weights):
            integral += np.exp(root) * self._qderivatives(sroot, eta_der, beta_der, first=False) * weight

        return integral

    def _compute_fermi(self, eta_der, beta_der):

        integral = 0
        S_1, S_2, S_3 = self._fermi_points(eta_der)

        a1 = 0.0
        b2 = np.sqrt(S_1)
        integral += self._compute_legendre(a1, b2, eta_der, beta_der, first=True)

        a3 = S_1
        b4 = S_2
        integral += self._compute_legendre(a3, b4, eta_der, beta_der)

        a5 = S_2
        b6 = S_3
        integral += self._compute_legendre(a5, b6, eta_der, beta_der)

        a7 = S_3
        integral += self._compute_laguerre(a7, eta_der, beta_der)

        return integral

    def evaluate(self):

        fval = self._compute_fermi(eta_der=0, beta_der=0)
        df_deta = self._compute_fermi(eta_der=1, beta_der=0)
        df_deta2 = self._compute_fermi(eta_der=2, beta_der=0)
        df_dbeta = self._compute_fermi(eta_der=0, beta_der=1)
        df_dbeta2 = self._compute_fermi(eta_der=0, beta_der=2)
        df_deta_dbeta = self._compute_fermi(eta_der=1, beta_der=1)

        return np.array([fval, df_deta, df_deta2, df_dbeta, df_dbeta2, df_deta_dbeta])
