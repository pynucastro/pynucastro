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


class BreakPoints:

    def __init__(self, *, itype="F"):

        if itype == "F":
            self.D = 3.36091
            self.sigma = 9.11856e-2
            self.a = np.array([6.77740, 3.76010, 7.56690])
            self.b = np.array([1.14180, 9.37188e-2, 1.16953])
            self.c = np.array([2.98255, 2.10635e-2, 7.54162e-1])
            self.d = np.array([0.0, 3.10836e1, 6.65585])
            self.e = np.array([0.0, 1.00557, -1.28190e-1])

        elif itype == "dF/deta":
            self.D = 4.99551
            self.sigma = 9.11856e-2
            self.a = np.array([6.77740, 3.76010, 7.56690])
            self.b = np.array([1.14180, 9.37188e-2, 1.16953])
            self.c = np.array([2.98255, 2.10635e-2, 7.54162])
            self.d = np.array([0.0, 3.95015e1, 7.64734])
            self.e = np.array([0.0, 1.00557, -1.28190e-1])

        elif itype == "d2F/deta2":
            self.D = 3.93830
            self.sigma = 9.11856e-2
            self.a = np.array([6.77740, 3.76010, 7.56690])
            self.b = np.array([1.14180, 9.37188e-2, 1.16953])
            self.c = np.array([2.98255, 2.10635e-2, 7.54162])
            self.d = np.array([0.0, 3.14499e1, 6.86346])
            self.e = np.array([0.0, 1.00557, -1.28190e-1])

        elif itype == "d3F/deta3":
            self.D = 4.17444
            self.sigma = 9.11856e-2
            self.a = np.array([6.77740, 3.76010, 7.56690])
            self.b = np.array([1.14180, 9.37188e-2, 1.16953])
            self.c = np.array([2.98255, 2.10635e-2, 7.54162])
            self.d = np.array([0.0, 3.05412e1, 7.88030])
            self.e = np.array([0.0, 1.00557, -1.28190e-1])

    def get_points(self, eta):

        xi = np.log(1.0 + np.exp(self.sigma * (eta - self.D))) / self.sigma

        X_a = ((self.a[0] + self.b[0] * xi + self.c[0] * xi**2) /
               (1.0 + self.c[0] * xi))

        X_b = ((self.a[1] + self.b[1] * xi + self.c[1] * self.d[1] * xi**2) /
               (1.0 + self.e[1] * xi + self.c[1] * xi**2))

        X_c = ((self.a[2] + self.b[2] * xi + self.c[2] * self.d[2] * xi**2) /
               (1.0 + self.e[2] * xi + self.c[2] * xi**2))

        return X_a - X_b, X_a, X_a + X_c


class FermiIntegrals:
    r"""Construct the integral

    .. math::

       F_k(\eta, \beta) = \int_0^\infty
           \frac{x^k [1 + (x\beta/2)]^{1/2}}{e^{x-\eta} + 1} dx

    using the method from Gong et al. 2001.  This splits the
    integration into 4 intervals, and uses Legendre quadrature for the
    first 3 parts and Laguerre quadrature for the last interval.  For
    the first interval, a change of variables is done to avoid
    singularities (effectively integrating in terms of momentum
    instead of energy).
    """

    def __init__(self, k, eta, beta, point_per_interval=100):
        self.eta = eta
        self.beta = beta
        self.k = k

        # for the Legendre quadrature, we use 1/2 the weights, so
        # we require the number of quadrature points to be even
        assert point_per_interval % 2 == 0

        self.point_per_interval = point_per_interval

        self.F = None
        self.dF_deta = None
        self.dF_dbeta = None
        self.d2F_deta2 = None
        self.d2F_detadbeta = None
        self.d2F_dbeta2 = None

    def __str__(self):
        fstr = ""
        fstr += f"F        = {self.F}\n"
        fstr += f"dF/dη    = {self.dF_deta}\n"
        fstr += f"dF/dβ    = {self.dF_dbeta}\n"
        fstr += f"d²F/dη²  = {self.d2F_deta2}\n"
        fstr += f"d²F/dηdβ = {self.d2F_detadbeta}\n"
        fstr += f"d²F/dβ²  = {self.d2F_dbeta2}\n"
        return fstr

    def _integrand(self, x, *,
                   eta_der=0, beta_der=0,
                   interval=0):

        assert 0 <= eta_der <= 2
        assert 0 <= beta_der <= 2

        eta = self.eta
        beta = self.beta
        k = self.k

        if interval == 0:
            # we want to work in terms of x**2
            # see Aparicio 1998 (but note they are missing a factor of 2
            # in the conversion from x to z).

            xsq = x * x
            sqrt_term = np.sqrt(1.0 + (0.5*xsq*beta))
            num = 2.0 * x**(2*k + 1.0) * sqrt_term
            denomi = 1.0 / (np.exp(xsq - eta) + 1.0)
            test = xsq - eta < -700.0

            # now construct the integrand for what we are actual computing
            if eta_der == 0 and beta_der == 0:
                if test:
                    return num
                return num * denomi

            if eta_der == 1 and beta_der == 0:
                # this is IB = 1 from Gong et al.
                # this corresponds to eq A.1 in terms of x**2
                if test:
                    return 0.0
                return num / (np.exp(xsq - eta) + 2.0 + np.exp(eta - xsq))

            if eta_der == 0 and beta_der == 1:
                # this is IB = 2 from Gong et al.
                # this corresponds to eq A.2 in terms of x**2
                if test:
                    return 0.5 * x**(2.0*k + 3.0) / sqrt_term
                return 0.5 * x**(2.0*k + 3.0) * denomi / sqrt_term

            if eta_der == 2 and beta_der == 0:
                # this is IB = 3 from Gong et al.
                # this corresponds to eq A.3 in terms of x**2
                if test:
                    return 0.0
                return num / (np.exp(xsq - eta) + 2.0 + np.exp(eta - xsq)) * \
                    ((np.exp(xsq - eta) - 1.0) / (np.exp(xsq - eta) + 1.0))

            if eta_der == 1 and beta_der == 1:
                # this is IB = 4 from Gong et al.
                # this corresponds to eq A.4 in terms of x**2
                if test:
                    return 0.0
                return 0.5 * x**(2.0*k + 3.0) / \
                    (np.exp(xsq - eta) + 2.0 + np.exp(eta - xsq)) / sqrt_term

            if eta_der == 0 and beta_der == 2:
                # this is IB = 5 from Gong et al.
                # this corresponds to eq A.5 in terms of x**2
                if test:
                    return -0.125 * x**(2.0*k + 5.0) / sqrt_term**3
                return -0.125 * x**(2.0*k + 5.0) * denomi / sqrt_term**3

        else:
            # we will work in terms of x

            sqrt_term = np.sqrt(1.0 + (0.5*x*beta))
            num = x**k * sqrt_term
            denomi = 1.0 / (np.exp(x - eta) + 1.0)
            test = x - eta < 700

            # now construct the integrand for what we are actual computing
            if eta_der == 0 and beta_der == 0:
                if test:
                    return num * denomi
                return 0.0

            if eta_der == 1 and beta_der == 0:
                # this is IB = 1 from Gong et al.
                # this corresponds to eq A.1
                if test:
                    return num / (np.exp(x - eta) + 2.0 + np.exp(eta - x))
                return 0.0

            if eta_der == 0 and beta_der == 1:
                # this is IB = 2 from Gong et al.
                # this corresponds to eq A.2
                if test:
                    return 0.25 * x**(k + 1.0) * denomi / sqrt_term
                return 0.0

            if eta_der == 2 and beta_der == 0:
                # this is IB = 3 from Gong et al.
                # this corresponds to eq A.3
                if test:
                    return num / (np.exp(x - eta) + 2.0 + np.exp(eta - x)) * \
                    ((np.exp(x - eta) - 1.0) / (np.exp(x - eta) + 1.0))
                return 0.0

            if eta_der == 1 and beta_der == 1:
                # this is IB = 4 from Gong et al.
                # this corresponds to eq A.4
                if test:
                    return 0.25 * x**(k + 1.0) / \
                        (np.exp(x - eta) + 2.0 + np.exp(eta - x)) / sqrt_term
                return 0.0

            if eta_der == 0 and beta_der == 2:
                # this is IB = 5 from Gong et al.
                # this corresponds to eq A.5
                if test:
                    return -0.0625 * x**(k + 2.0) * denomi / sqrt_term**3
                return 0.0

    def _compute_legendre(self, a, b, eta_der, beta_der, interval=None):

        roots, weights = roots_legendre(self.point_per_interval)
        N = self.point_per_interval//2

        integral = 0

        # We set the correspondence (a+b)/2 -> 0 and map the (-1,0) and (0,1)
        # intervals separately.

        for weight, root in zip(weights[N:], roots[N:]):
            x_1 = (a+b)/2 + (b-a)/2 * root
            x_2 = (a+b)/2 - (b-a)/2 * root
            qder_1 = self._integrand(x_1, eta_der=eta_der, beta_der=beta_der,
                                     interval=interval)
            qder_2 = self._integrand(x_2, eta_der=eta_der, beta_der=beta_der,
                                     interval=interval)
            integral += (qder_1 + qder_2) * weight

        integral *= (b-a)/2
        return integral

    def _compute_laguerre(self, a, eta_der, beta_der, interval=100):

        point_per_interval = self.point_per_interval
        roots, weights = roots_laguerre(point_per_interval)
        integral = 0

        scaled_roots = roots + a
        for sroot, root, weight in zip(scaled_roots, roots, weights):
            integral += np.exp(root) * self._integrand(sroot, eta_der=eta_der, beta_der=beta_der, interval=interval) * weight

        return integral

    def _compute_fermi(self, eta_der, beta_der):

        if eta_der == 0:
            bp = BreakPoints(itype="F")
        elif eta_der == 1:
            bp = BreakPoints(itype="dF/deta")
        elif eta_der == 2:
            bp = BreakPoints(itype="d2F/deta2")
        elif eta_der == 3:
            bp = BreakPoints(itype="d3F/deta3")

        S_1, S_2, S_3 = bp.get_points(self.eta)
        print(S_1, S_2, S_3)

        I0 = self._compute_legendre(0.0, np.sqrt(S_1),
                                    eta_der, beta_der, interval=0)

        I1 = self._compute_legendre(S_1, S_2,
                                    eta_der, beta_der, interval=1)

        I2 = self._compute_legendre(S_2, S_3,
                                    eta_der, beta_der, interval=2)

        I3 = self._compute_laguerre(S_3, eta_der, beta_der)

        print("intervals = ", I0, I1, I2, I3)

        return I0 + I1 + I2 + I3

    def evaluate(self):

        self.F = self._compute_fermi(eta_der=0, beta_der=0)
        self.dF_deta = self._compute_fermi(eta_der=1, beta_der=0)
        self.d2F_deta2 = self._compute_fermi(eta_der=2, beta_der=0)
        self.dF_dbeta = self._compute_fermi(eta_der=0, beta_der=1)
        self.d2F_detadbeta = self._compute_fermi(eta_der=1, beta_der=1)
        self.d2F_dbeta2 = self._compute_fermi(eta_der=0, beta_der=2)
