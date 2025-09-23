"""Methods to ease the creation of networks."""

from pynucastro.rates import (DerivedRate, ReacLibLibrary, ReacLibRate,
                              TabularLibrary)

from .amrexastro_cxx_network import AmrexAstroCxxNetwork
from .fortran_network import FortranNetwork
from .python_network import PythonNetwork
from .simple_cxx_network import SimpleCxxNetwork


def network_helper(nuclei, *,
                   network_type="python",
                   use_detailed_balance=True,
                   use_tabular_rates=True,
                   tabular_ordering=None):
    """Generate a basic network connecting all of the input nuclei
    using all of the rates known to pynucastro.

    Parameters
    ----------
    nuclei : Iterable(Nucleus)
        the nuclei to use for the network
    network_type : str
        the type of network to create.  Allowed values are:

        * "python" : create a :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`
        * "cxx" : create a :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`
        * "fortran" : create a :py:obj:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>`
        * "amrex" : create a :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`
    use_detailed_balanace : bool
        do we rederive inverse rates using detailed balance?
    use_tabular_rates : bool
        do we include tabulated weak rates?
    tabular_ordering : Iterable(str)
        if we are including tabular rates, a list of sources
        can be provided to specify which rate sources are used,
        and the priority that each source should have.

    Returns
    -------
    PythonNetwork, SimpleCxxNetwork, AmrexAstroCxxNetwork, FortranNetwork

    """

    rl = ReacLibLibrary()
    lib = rl.linking_nuclei(nuclei)

    if use_tabular_rates:
        if tabular_ordering:
            tl = TabularLibrary(ordering=tabular_ordering)
        else:
            tl = TabularLibrary()

        lib += tl.linking_nuclei(nuclei)

        # if we have both a tabular and ReacLib rate,
        # remove the ReacLib version

        lib.eliminate_duplicates(rate_type_preference="tabular")

    if use_detailed_balance:
        rates_to_derive = lib.backward().get_rates()

        # now for each of those derived rates, look to see if the pair exists
        for r in rates_to_derive:
            fr = lib.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                print(f"modifying {r} from {fr}")
                lib.remove_rate(r)
                d = DerivedRate(rate=fr, compute_Q=False, use_pf=True)
                lib.add_rate(d)

    if network_type == "python":
        return PythonNetwork(libraries=[lib])

    if network_type == "cxx":
        return SimpleCxxNetwork(libraries=[lib])

    if network_type == "fortran":
        return FortranNetwork(libraries=[lib])

    if network_type == "amrex":
        return AmrexAstroCxxNetwork(libraries=[lib])

    raise ValueError("invalid network_type")
