"""A collection of methods for accessing data from a MESA models
as read by py_mesa_reader, and putting them into a form that
pynucastro can use.

"""

import re

from pynucastro.nucdata import Composition, Nucleus, UnidentifiedElement


def get_nuclei(model):
    """Return a list of nuclei contained in the MESA model.

    Parameters
    ----------
    model : ``mesa_reader.MesaData``
        The MESA model as read by ``mesa_reader.MesaData``

    Returns
    -------
    list(Nucleus)

    """

    # the species may be written as log_he4, etc. in the
    # "bulk names" section of the header
    species_re = re.compile(r"log_([a-z]+[\d]+)")

    # sometimes they may be without the log_.  In this case,
    # restrict the element name to 1-2 characters
    species_nolog_re = re.compile(r"^([a-z]{1,2}[\d]+)$")

    nuclei = []
    for field in model.bulk_names:
        if g := species_re.match(field):
            nuclei.append(Nucleus(g.group(1)))
        elif g := species_nolog_re.match(field):
            try:
                nuc = Nucleus(g.group(1))
            except UnidentifiedElement:
                continue
            nuclei.append(nuc)

    return nuclei


def get_zone_data(model, i, *, nuclei=None):
    """Return a tuple of (rho, T, composition) for a single zone from
    the MESA model.  These are the thermodynamic conditions needed to
    evaluate a pynucastro network's rates.

    Parameters
    ----------
    model : ``mesa_reader.MesaData``
        The MESA model as read by ``mesa_reader.MesaData``
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
    for n in nuclei:
        # sometimes the species data is in terms of log
        try:
            val = 10.0**model.bulk_data[i][f"log_{str(n).lower()}"]
        except ValueError:
            val = model.bulk_data[i][f"{str(n).lower()}"]
        comp[n] = val
    return (rho, T, comp)


def get_all_data(model):
    """Return a dictionary keyed by MESA zone index containing the
    thermodynamic information (rho, T, composition) for all the
    zones in the model.

    Parameters
    ----------
    model : ``mesa_reader.MesaData``
        The MESA model as read by ``mesa_reader.MesaData``

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

    return mesa_zones
