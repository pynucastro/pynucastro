"""Some high-order finite-difference approximations for the EOS."""


def fourth_order_diff(func, x0, delta, component=None):
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
