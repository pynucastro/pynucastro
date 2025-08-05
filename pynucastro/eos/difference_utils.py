"""Some high-order finite-difference approximations for the EOS."""

import numpy as np

# the coefficients for these can be found in many sources, including:
# https://en.wikipedia.org/wiki/Finite_difference_coefficient


def fourth_order_diff(func, x0, delta, *, component=None):
    """Compute a 4th order accurate centered difference approximation
    of a function, and allow us to specify the component of the object
    that is returned (if applicable)

    Parameters
    ----------
    func : Callable
        the function to difference, assumed to be of the form `func(x)`
    x0 : float
        the point at which to approximate the derivative
    delta : float
        the step-size to use
    component : str
        if func returns an object, use this component for the derivative.

    Returns
    -------
    float

    """

    fvals = []
    for i in [-2, -1, 0, 1, 2]:
        if i == 0:
            fvals.append(None)
            continue
        _f = func(x0 + i*delta)
        if component:
            fvals.append(getattr(_f, component))
        else:
            fvals.append(_f)

    deriv = (fvals[0] - 8.0 * fvals[1] + 8.0 * fvals[3] - fvals[4]) / (12 * delta)
    return deriv


def sixth_order_diff(func, x0, delta, *, component=None):
    """Compute a 6th order accurate centered difference approximation
    of a function, and allow us to specify the component of the object
    that is returned (if applicable)

    Parameters
    ----------
    func : Callable
        the function to difference, assumed to be of the form `func(x)`
    x0 : float
        the point at which to approximate the derivative
    delta : float
        the step-size to use
    component : str
        if func returns an object, use this component for the derivative.

    Returns
    -------
    float

    """

    fvals = []
    for i in [-3, -2, -1, 0, 1, 2, 3]:
        if i == 0:
            fvals.append(None)
            continue
        _f = func(x0 + i*delta)
        if component:
            fvals.append(getattr(_f, component))
        else:
            fvals.append(_f)

    deriv = (-fvals[0] + 9.0 * fvals[1] - 45.0 * fvals[2] + 45.0 * fvals[4] - 9.0 * fvals[5] + fvals[6]) / (60 * delta)
    return deriv


def eighth_order_diff(func, x0, delta, *, component=None):
    """Compute an 8th order accurate centered difference approximation
    of a function, and allow us to specify the component of the object
    that is returned (if applicable)

    Parameters
    ----------
    func : Callable
        the function to difference, assumed to be of the form `func(x)`
    x0 : float
        the point at which to approximate the derivative
    delta : float
        the step-size to use
    component : str
        if func returns an object, use this component for the derivative.

    Returns
    -------
    float

    """

    fvals = []
    for i in [-4, -3, -2, -1, 0, 1, 2, 3, 4]:
        if i == 0:
            fvals.append(None)
            continue
        _f = func(x0 + i*delta)
        if component:
            fvals.append(getattr(_f, component))
        else:
            fvals.append(_f)

    deriv = ((1.0 / 280.0) * (fvals[0] - fvals[8]) +
             (4.0 / 105.0) * (-fvals[1] + fvals[7]) +
             0.2 * (fvals[2] - fvals[6]) +
             0.8 * (-fvals[3] + fvals[5])) / delta
    return deriv


def adaptive_diff(func, x0, h, *, component=None, max_levels=10):
    """Perform an adaptive centered-difference estimate to the first
    derivative using Richardson-extrapolation / Ridders' method.


    Parameters
    ----------
    func : Callable
        the function to difference, assumed to be of the form `func(x)`
    x0 : float
        the point at which to approximate the derivative
    delta : float
        the step-size to use
    component : str
        if func returns an object, use this component for the derivative.
    max_levels : int
        the number of levels / iterations in the tableau used in
        extrapolating the deriviative to higher order

    Returns
    -------
    deriv : float
        an estimate of the derivative at x0
    err : float
        the estimated error in the difference approximation

    """

    # see C.J.F. Ridders, "Accurate computation of F'x(x) and F''(x)
    # Advances in Engineering Software, 4, 2, 1978

    # we'll use a second-order centered difference, starting with the
    # initial h.  this has the form:
    #
    #  D_h(f) = 0.5 * (f(x0 + h) - f(x0 - h)) / h + c (h**2) + O(h**4)
    #
    # the key point is that C depends only on the function f (or its
    # derivatives) at point x0.  That means we can compute D_{h/2}(f)
    # and then combine the 2 for a better estimate of the derivative:
    #
    # D' = (4 D_{h/2}(f) - D_h(f)) / 3 + O(h**4)

    # In Ridders' method, we build a tableau of the form:
    #
    #
    #   A(1, 1)  A(1, 2)  A(1, 3)  A(1, 4) ...
    #            A(2, 2)  A(2, 3)  A(2, 4) ...
    #                     A(3, 3)  A(3, 3) ...
    #                              A(4, 4) ...
    #
    # where the entries in each row come via extrapolation of the data
    # in the previous row.
    #
    # In this notation, A(1, 1) = D_h(f),
    #                   A(1, 2) = D_{h/2}(f)
    #                   A(2, 2) = D'
    #
    # we then have:
    #
    # A(1, n) = 0.5 * (f(x + h/2**(n-1)) - f(x - h/2**(n-1))) / (h/2**(n-1))
    #
    # A(m, n) = (4**(m-1) A(m-1, n) - A(m-1, n-1)) / (4**(m-1) - 1)

    # Note: we use 0-based indexing below

    A = np.zeros((max_levels, max_levels), dtype=np.float64)

    assert h > 0, "step size, h, must be positive"

    _h = h

    jump = 2

    _fp = func(x0 + _h)
    _fm = func(x0 - _h)
    if component:
        _fp = getattr(_fp, component)
        _fm = getattr(_fm, component)

    A[0, 0] = 0.5 * (_fp - _fm) / _h
    deriv = A[0, 0]

    err = np.finfo(np.float64).max

    for n in range(1, max_levels):
        _h /= jump
        _fp = func(x0 + _h)
        _fm = func(x0 - _h)
        if component:
            _fp = getattr(_fp, component)
            _fm = getattr(_fm, component)

        A[0, n] = 0.5 * (_fp - _fm) / _h

        # for the current column n, fill in all the rows m > 0
        for m in range(1, n+1):
            A[m, n] = ((jump * jump)**m * A[m-1, n] - A[m-1, n-1]) / ((jump * jump)**m - 1)

            # estimate the error so far, using the current row
            # and the estimates in the row above us
            err_m = max(np.abs(A[m, n] - A[m-1, n]), np.abs(A[m, n] - A[m-1, n-1]))

            if err_m < err:
                err = err_m
                deriv = A[m, n]

        # at this point, A(n, n) should be the most accurate estimate
        # for the derivative, but it's possible that things have
        # become worse -- let's check
        if np.abs(A[n, n] - A[n-1, n-1]) > 4.0 * err:
            break

    return float(deriv), float(err)
