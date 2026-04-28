"""A collection of methods for accessing data from a MESA models
as read by py_mesa_reader, and putting them into a form that
pynucastro can use.

"""

import re

from pynucastro.nucdata import Composition, Nucleus


def get_nuclei(model):
    """Return a list of nuclei contained in the MESA model.

    Parameters
    ----------
    model : MesaData
        The MESA model as read bt mesa_reader.MesaData

    Returns
    -------
    list(Nuclei)

    """

    # the species are written as log_he4, etc. in the
    # "bulk names" section of the header
    species_re = re.compile(r"log_([a-z]+[\d]+)")

    nuclei = []
    for field in model.bulk_names:
        if g := species_re.match(field):
            nuclei.append(Nucleus(g.group(1)))

    return nuclei


def get_zone_data(model, i, *, nuclei=None):
    """Return a tuple of (rho, T, composition) for a single zone from
    the MESA model.  These are the thermodynamic conditions needed to
    evaluate a pynucastro network's rates.

    Parameters
    ----------
    model : MesaData
        The MESA model as read bt mesa_reader.MesaData
    i : int
        The MESA zone index
    nuclei : list(Nucleus)
        The list of nuclei in the MESA model.  This will be
        read from the model if not provided, but if you are
        looping over multiple zones in succession, passing
        a cached list in can speed up the reading.
    Returns
    -------
    tuple(float, float, Composition)

    """

    if nuclei is None:
        nuclei = get_nuclei(model)

    rho = 10.0**model.bulk_data[i]["logRho"]
    T = 10.0**model.bulk_data[i]["logT"]
    comp = Composition(nuclei)

    return (rho, T, comp)


def get_all_data(model):
    """Return a dictionary keyed by MESA zone index containing the
    thermodynamic information (rho, T, composition) for all the
    zones in the model.

    Parameters
    ----------
    model : MesaData
        The MESA model as read bt mesa_reader.MesaData

    Returns
    -------
    dict[int, tuple(float, float, Composition)]

    """

    mesa_zones = {}

    # get the list of nuclei objects once, so we don't
    # need to repeat the construction.
    nuclei = get_nuclei(model)

    for i in range(model.header_data["num_zones"]):
        rho, T, comp = get_zone_data(model, i, nuclei=nuclei)
        mesa_zones[i] = (rho, T, comp)
