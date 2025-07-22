"""Some high-order finite-difference approximations for the EOS."""


def fourth_order_rho(eos, state, component, delta):
    """Compute a 4th order accurate centered difference approximation
    to the density derivative of an EOS quantity.

    Parameters
    ----------
    eos : ElectronEOS
        the EOS object
    state : tuple
        (density, temperature, Composition)
    component : str
        the EOS quantity to differentiate
    delta : float
        the step-size to use

    Returns
    -------
    float

    """

    rho, T, comp = state

    fvals = []
    for i in [-2, -1, 0, 1, 2]:
        if i == 0:
            fvals.append(None)
            continue
        _es = eos.pe_state(rho + i*delta, T, comp)
        fvals.append(getattr(_es, component))

    deriv = (fvals[0] - 8.0 * fvals[1] + 8.0 * fvals[3] - fvals[4]) / (12 * delta)
    return deriv


def fourth_order_temp(eos, state, component, delta):
    """Compute a 4th order accurate centered difference approximation
    to the temperature derivative of an EOS quantity.

    Parameters
    ----------
    eos : ElectronEOS
        the EOS object
    state : tuple
        (density, temperature, Composition)
    component : str
        the EOS quantity to differentiate
    delta : float
        the step-size to use

    Returns
    -------
    float

    """

    rho, T, comp = state

    fvals = []
    for i in [-2, -1, 0, 1, 2]:
        if i == 0:
            fvals.append(None)
            continue
        _es = eos.pe_state(rho, T + i*delta, comp)
        fvals.append(getattr(_es, component))

    deriv = (fvals[0] - 8.0 * fvals[1] + 8.0 * fvals[3] - fvals[4]) / (12 * delta)
    return deriv
