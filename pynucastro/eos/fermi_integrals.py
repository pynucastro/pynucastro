"""Classes that enable the computation of Fermi-Dirac integrals to
high precision.  This uses the methods of Aparicio (1998) and Gong
(2001), described in https://iopscience.iop.org/article/10.1086/313121
and
https://www.sciencedirect.com/science/article/abs/pii/S001046550100145X?via%3Dihub.

The method consist on breaking down the Fermi-integrals (eta, theta)
integration domain into four subintervals and apply Gauss-Legendre or
the Gauss laguerre quadrature methods with 200 points over each domain
to guarantee accuracy beyond double-precision.

In the first subinterval, a x=z^2 change of variable is applied to
overcome the kernel singularity near the origin.

"""

import numba
import numpy as np

from .quadrature_weights import w_lag, w_leg, x_lag, x_leg


class BreakPoints:
    """The break points used in splitting the integral from [0, inf] into
    4 separate integrals over smaller domains.  These are described in
    Gong et al. 2001.  Different choices of the break points are made
    depending on whether we are computing the integral itself or one
    of its derivatives with respect to eta (the degeneracy parameter).

    Parameters
    ----------
    itype : str
        the type of integral we are doing---mainly this concerns the
        eta derivatives.  Valid choices are:

        * "F" : the Fermi-Dirac integral

        * "dF/deta" : the first derivative of the F-D integral with
          respect to eta.

        * "d2F/deta2" : the second derivative of the F-D integral with
          respect to eta.

        * "d3F/deta3" : the third derivative of the F-D integral with
          respect to eta.

    """

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

        else:
            raise ValueError("invalid itype")

    def get_points(self, eta):
        """Return the 3 break points that define the 4 sub-intervals
        of integration

        Parameters
        ----------
        eta : float
            The degeneracy parameter.

        Returns
        -------
        S_1, S_2, S_3 : float
            The break points.

        """

        term = self.sigma * (eta - self.D)
        if term <= 50:
            xi = np.log(1.0 + np.exp(term)) / self.sigma
        else:
            xi = eta - self.D

        X_a = ((self.a[0] + self.b[0] * xi + self.c[0] * xi**2) /
               (1.0 + self.c[0] * xi))

        X_b = ((self.a[1] + self.b[1] * xi + self.c[1] * self.d[1] * xi**2) /
               (1.0 + self.e[1] * xi + self.c[1] * xi**2))

        X_c = ((self.a[2] + self.b[2] * xi + self.c[2] * self.d[2] * xi**2) /
               (1.0 + self.e[2] * xi + self.c[2] * xi**2))

        return X_a - X_b, X_a, X_a + X_c


@numba.njit()
def _kernel_p(x, k, eta, beta,
              eta_der=0, beta_der=0):

    result = 0.0

    # we want to work in terms of x**2
    # see Aparicio 1998 (but note they are missing a factor of 2
    # in the conversion from x to z).

    xsq = x * x
    sqrt_term = np.sqrt(1.0 + 0.5 * x * x * beta)
    num = 2.0 * x**(2*k + 1.0) * sqrt_term
    # this is what we are usually exponentiating
    delta = xsq - eta
    cosh_delta = np.cosh(delta)  # this is 0.5 * (exp(xsq - eta) + exp(eta - xsq))
    testm = delta < -700.0
    if testm:
        denomi = 1.0
        exp_delta = 0.0
    else:
        exp_delta = np.exp(delta)
        inv_exp_delta = 1.0 / exp_delta
        # 1 / (exp(x**2 - eta) + 1) rewritten
        denomi = inv_exp_delta / (1.0 + inv_exp_delta)

    # now construct the integrand for what we are actual computing
    if eta_der == 0 and beta_der == 0:
        result = num
        if not testm:
            result *= denomi

    elif eta_der == 1 and beta_der == 0:
        # this is IB = 1 from Gong et al.
        # this corresponds to eq A.1 in terms of x**2
        if not testm:
            result = num / (2.0 * (1.0 + cosh_delta))

    elif eta_der == 0 and beta_der == 1:
        # this is IB = 2 from Gong et al.
        # this corresponds to eq A.2 in terms of x**2
        result = 0.5 * x**(2.0*k + 3.0) / sqrt_term
        if not testm:
            result *= denomi

    elif eta_der == 2 and beta_der == 0:
        # this is IB = 3 from Gong et al.
        # this corresponds to eq A.3 in terms of x**2
        if not testm:
            result = num / (2.0 * (1.0 + cosh_delta)) * \
                ((np.exp(xsq - eta) - 1.0) / (np.exp(xsq - eta) + 1.0))

    elif eta_der == 1 and beta_der == 1:
        # this is IB = 4 from Gong et al.
        # this corresponds to eq A.4 in terms of x**2
        if not testm:
            result = 0.5 * x**(2.0*k + 3.0) / (2.0 * (1.0 + cosh_delta)) / sqrt_term

    elif eta_der == 0 and beta_der == 2:
        # this is IB = 5 from Gong et al.
        # this corresponds to eq A.5 in terms of x**2
        result = -0.125 * x**(2.0*k + 5.0) / sqrt_term**3
        if not testm:
            result *= denomi

    return result


@numba.njit()
def _kernel_E(x, k, eta, beta,
              eta_der=0, beta_der=0):

    result = 0.0

    # we will work in terms of x

    sqrt_term = np.sqrt(1.0 + 0.5 * x * beta)
    num = x**k * sqrt_term
    # this is what we are usually exponentiating
    delta = x - eta
    cosh_delta = np.cosh(delta)
    testm = x - eta < -700
    if testm:
        denomi = 1.0
        exp_delta = 0.0
    else:
        exp_delta = np.exp(delta)
        inv_exp_delta = 1.0 / exp_delta
        # 1 / (exp(x - eta) + 1) rewritten
        denomi = inv_exp_delta / (1.0 + inv_exp_delta)

    # now construct the integrand for what we are actual computing
    if eta_der == 0 and beta_der == 0:
        result = num
        if not testm:
            result *= denomi

    elif eta_der == 1 and beta_der == 0:
        # this is IB = 1 from Gong et al.
        # this corresponds to eq A.1
        if not testm:
            result = num / (2.0 * (1.0 + cosh_delta))

    elif eta_der == 0 and beta_der == 1:
        # this is IB = 2 from Gong et al.
        # this corresponds to eq A.2
        result = 0.25 * x**(k + 1.0)  / sqrt_term
        if not testm:
            result *= denomi

    elif eta_der == 2 and beta_der == 0:
        # this is IB = 3 from Gong et al.
        # this corresponds to eq A.3
        if not testm:
            result = num / (2.0 * (1.0 + cosh_delta)) * \
                ((1.0 - np.exp(eta - x)) / (1.0 + np.exp(eta - x)))

    elif eta_der == 1 and beta_der == 1:
        # this is IB = 4 from Gong et al.
        # this corresponds to eq A.4
        if not testm:
            result = 0.25 * x**(k + 1.0) / (2.0 * (1.0 + cosh_delta)) / sqrt_term

    elif eta_der == 0 and beta_der == 2:
        # this is IB = 5 from Gong et al.
        # this corresponds to eq A.5
        result = -0.0625 * x**(k + 2.0) / sqrt_term**3
        if not testm:
            result *= denomi

    return result


class FermiIntegral:
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

    First and second derivatives with respect to $\eta$ and $\beta$
    are supported.

    Parameters
    ----------
    k : float
        The index of the Fermi-Dirac integral
    eta : float
        The degeneracy parameter ($\eta = \mu / k_B T$)
    beta : float
        The dimenionless temperature (relativity parameter)
        ($\beta = k_B T / m_e c^2$)

    """

    def __init__(self, k, eta, beta):
        self.eta = eta
        self.beta = beta
        self.k = k

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

    def _compute_legendre(self, kernel, a, b,
                          k, eta, beta,
                          eta_der, beta_der):

        N = len(x_leg)//2

        integral = 0

        # We set the correspondence (a+b)/2 -> 0 and map the (-1,0) and (0,1)
        # intervals separately.

        fac1 = (a + b)/2
        fac2 = (b - a)/2

        for root, weight in zip(x_leg[N:], w_leg[N:]):
            x_1 = fac1 + fac2 * root
            x_2 = fac1 - fac2 * root
            qder_1 = kernel(x_1, k, eta, beta,
                            eta_der=eta_der, beta_der=beta_der)
            qder_2 = kernel(x_2, k, eta, beta,
                            eta_der=eta_der, beta_der=beta_der)
            integral += (qder_1 + qder_2) * weight

        integral *= (b - a) / 2
        return integral

    def _compute_laguerre(self, a, eta_der, beta_der):

        # Laguerre quadrature solves and integral of the form:
        #
        #   ∞
        # ∫  f(x) exp(-x) dx ~ ∑ f(x_i) w_i
        #  0
        #
        # where x_i are the nodes (roots of the Laguerre polynomial)
        # and w_i are the weights
        #
        # We want to integrate just f(x) from a instead of 0, then we
        # change variables x = z + a
        #
        #   ∞            ∞
        # ∫  f(x) dx = ∫  f(z + a) dz
        #  a            0
        #
        # Since Laguerre quadrature includes an exp(-x) kernel, we
        # need to scale our function by exp(x), giving
        #
        #   ∞             ∞
        # ∫  f(x) dx -> ∫  [f(x + a) exp(x)] exp(-x) dx
        #  a             0
        #
        #             ~ ∑ f(x_i) w_i exp(x) f(x + a)

        # note: the w_lag already have the exp(x) term included
        # for root, weight in zip(x_lag, w_lag):
        #     I = _kernel_E(root + a, self.k, self.eta, self.beta,
        #                   eta_der=eta_der, beta_der=beta_der)
        #     integral += I * weight

        integral = sum(_kernel_E(root + a, self.k, self.eta, self.beta,
                                 eta_der=eta_der, beta_der=beta_der) * weight for
                       root, weight in zip(x_lag, w_lag))

        return integral

    def _compute_fermi(self, eta_der, beta_der):

        bp = None
        if eta_der == 0:
            bp = BreakPoints(itype="F")
        elif eta_der == 1:
            bp = BreakPoints(itype="dF/deta")
        elif eta_der == 2:
            bp = BreakPoints(itype="d2F/deta2")
        elif eta_der == 3:
            bp = BreakPoints(itype="d3F/deta3")

        S_1, S_2, S_3 = bp.get_points(self.eta)

        I0 = self._compute_legendre(_kernel_p, 0.0, np.sqrt(S_1),
                                    self.k, self.eta, self.beta,
                                    eta_der, beta_der)

        I1 = self._compute_legendre(_kernel_E, S_1, S_2,
                                    self.k, self.eta, self.beta,
                                    eta_der, beta_der)

        I2 = self._compute_legendre(_kernel_E, S_2, S_3,
                                    self.k, self.eta, self.beta,
                                    eta_der, beta_der)

        I3 = self._compute_laguerre(S_3, eta_der, beta_der)

        return I0 + I1 + I2 + I3

    def evaluate(self, *, do_first_derivs=True, do_second_derivs=True):
        """Perform the integration for the Fermi-Dirac function
        and its derivatives (if desired)

        Parameters
        ----------
        do_first_derivs : bool
            do we compute the first derivatives?
        do_second_derivs: bool
            do we compute the second derivatives?

        """

        self.F = self._compute_fermi(eta_der=0, beta_der=0)

        if do_first_derivs:
            self.dF_deta = self._compute_fermi(eta_der=1, beta_der=0)
            self.dF_dbeta = self._compute_fermi(eta_der=0, beta_der=1)

        if do_second_derivs:
            self.d2F_deta2 = self._compute_fermi(eta_der=2, beta_der=0)
            self.d2F_detadbeta = self._compute_fermi(eta_der=1, beta_der=1)
            self.d2F_dbeta2 = self._compute_fermi(eta_der=0, beta_der=2)
