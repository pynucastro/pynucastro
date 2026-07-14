"""A collection of methods for accessing data from a MESA models
as read by py_mesa_reader, and putting them into a form that
pynucastro can use.

"""

import re

import numpy as np

from pynucastro.nucdata import Composition, Nucleus, UnidentifiedElement
from pynucastro.constants.constants import R_sun, M_sun


class MesaZoneState:
    """A container for a single MESA zone's data

    Attributes
    ----------
    r : float
        radius (cm)
    m : float
        mass shell (g)
    rho : float
        density (CGS)
    T : float
        temperature (K)
    comp : Composition
        the composition of the zone

    """

    def __init__(self, r, m, rho, T, comp):
        self.r = r
        self.m = m
        self.rho = rho
        self.T = T
        self.comp = comp


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
    MesaZoneState

    """

    if nuclei is None:
        nuclei = get_nuclei(model)

    r = 10.0**model.bulk_data[i]["logR"] * R_sun
    m = model.bulk_data[i]["mass"] * M_sun
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

    return MesaZoneState(r, m, rho, T, comp)


class MesaModel:
    """Read data from a MESA model and store the thermodynamic
    state for each zone in a pynucastro-compatible format.

    Parameters
    ----------
    model : ``mesa_reader.MesaData``
        The MESA model as read by ``mesa_reader.MesaData``

    """

    def __init__(self, model):

        self.mesa_zones = {}

        # get the list of nuclei objects once, so we don't
        # need to repeat the construction.
        nuclei = get_nuclei(model)

        for i in range(model.header_data["num_zones"]):
            self.mesa_zones[i] = get_zone_data(model, i, nuclei=nuclei)

    def get_data_array(self, var="T"):
        """Return an array of a single variable indexed by zone

        Parameters
        ----------
        var : str
           The name of the variable to access.  This must be a member
           of :py:obj:`MesaZoneState`.

        Returns
        -------
        numpy.ndarray

        """

        return np.asarray([getattr(m, var) for _, m in self.mesa_zones.items()])

    def get_peak_index(self, var):
        """Return the index in the MESA model with the highest
        value of var.

        Parameters
        ----------
        var : str
           The name of the variable to access.  This must be a member
           of :py:obj:`MesaZoneState`.

        Returns
        -------
        int

        """

        vs = self.get_data_array(var)
        return np.argmax(vs)

    def get_zone_data(self, i):
        """Return the thermodynamic state for a single zone in the
        MESA model

        Parameters
        ----------
        i : int
            The index into the MESA model

        Returns
        -------
        MesaZoneState

        """

        return self.mesa_zones[i]
