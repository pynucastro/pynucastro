from pynucastro.networks.helper import get_net_class, network_helper
from pynucastro.nucdata import Nucleus
from pynucastro.rates import ApproximateRate, BranchedRate, Library, ModifiedRate, ReacLibLibrary, make_ap_pg_rates


def aprox13(*, network_type="python"):
    """Create the aprox13 network.

    aprox13 contains 13 nuclei and is suitable for explosive He and C
    burning.  It uses the (α,p)(p,γ) approximation for the α-chain
    and approximates C+C, C+O, and O+O burning.

    Parameters
    ----------
    network_type : str
        The type of network to create.  Allowed values are:

        * "python" : create a :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`
        * "cxx" : create a :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`
        * "fortran" : create a :py:obj:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>`
        * "amrex" : create a :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`

    Returns
    -------
    PythonNetwork, SimpleCxxNetwork, AmrexAstroCxxNetwork, FortranNetwork

    """

    nuc = ["p", "he4",
           "c12", "o16", "ne20", "na23",
           "mg24", "al27", "si28", "p31", "s32",
           "cl35", "ar36", "k39", "ca40",
           "sc43", "ti44", "v47", "cr48",
           "mn51", "fe52", "co55", "ni56"]

    net = network_helper(nuc, network_type=network_type)

    net.make_CO_burning_approx(root_nuclei="C")
    net.make_CO_burning_approx(root_nuclei="O")
    net.make_CO_burning_approx(root_nuclei="CO")

    net.remove_nuclei(["na23", "al27", "p31"])

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43",
                                               "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43",
                       "v47", "mn51", "co55"])

    return net


def mesa_basic(*, network_type="python"):
    """Create the MESA basic.net network.

    MESA's basic.net contains 8 nuclei and approximates pp-I,II,III,
    CNO (but not hot-CNO), and He-burning.

    Parameters
    ----------
    network_type : str
        The type of network to create.  Allowed values are:

        * "python" : create a :py:obj:`PythonNetwork <pynucastro.networks.python_network.PythonNetwork>`
        * "cxx" : create a :py:obj:`SimpleCxxNetwork <pynucastro.networks.simple_cxx_network.SimpleCxxNetwork>`
        * "fortran" : create a :py:obj:`FortranNetwork <pynucastro.networks.fortran_network.FortranNetwork>`
        * "amrex" : create a :py:obj:`AmrexAstroCxxNetwork <pynucastro.networks.amrexastro_cxx_network.AmrexAstroCxxNetwork>`

    Returns
    -------
    PythonNetwork, SimpleCxxNetwork, AmrexAstroCxxNetwork, FortranNetwork

    """

    rl = ReacLibLibrary()

    # pp

    rpp, rpep = rl.get_rate_by_name("p(p,)d")

    rppp_he3 = ModifiedRate(rpp,
                            new_products=[Nucleus("he3")],
                            stoichiometry={Nucleus("p"): 3},
                            description="p(p,e⁺ν)d(p,γ)He3",
                            rate_source=rpp.src)

    rpepp_he3 = ModifiedRate(rpep,
                             new_products=[Nucleus("he3")],
                             stoichiometry={Nucleus("p"): 3},
                             description="p(pe⁻,ν)d(p,γ)He3",
                             rate_source=rpep.src)

    rhe3he3 = rl.get_rate_by_name("he3(he3,pp)he4")

    rhe3he4 = rl.get_rate_by_name("he4(he3,g)be7")
    rhe3he4p_2he4 = ModifiedRate(rhe3he4,
                                 new_reactants=[Nucleus("he4"),
                                                Nucleus("he3"),
                                                Nucleus("p")],
                                 new_products=[Nucleus("he4"),
                                               Nucleus("he4")],
                                 not_in_ydot_term=[Nucleus("p")],
                                 description="He4(He3,γ)Be7(e⁻,ν)Li7(p,α)He4")

    lib_pp = Library(rates=[rppp_he3, rpepp_he3,
                            rhe3he3, rhe3he4p_2he4])

    # CNO

    rc12pg = rl.get_rate_by_name("c12(p,g)n13")
    rn14pg = rl.get_rate_by_name("n14(p,g)o15")
    rn15pa = rl.get_rate_by_name("n15(p,a)c12")
    rn15pg = rl.get_rate_by_name("n15(p,g)o16")
    ro16pg = rl.get_rate_by_name("o16(p,g)f17")

    rc12_2p_n14 = ModifiedRate(rc12pg,
                               new_products=[Nucleus("n14")],
                               stoichiometry={Nucleus("p"): 2},
                               description="C12(p,γ)N13(,e⁺ν)C13(p,γ)N14")

    rn14_2p_c12 = BranchedRate(rn14pg,
                               primary_branch=rn15pa,
                               other_branch=rn15pg,
                               stoichiometry={Nucleus("p"): 2},
                               description="N14(p,γ)O15(,e⁺ν)N15(p,α)C12")

    rn14_2p_o16 = BranchedRate(rn14pg,
                               primary_branch=rn15pg,
                               other_branch=rn15pa,
                               stoichiometry={Nucleus("p"): 2},
                               description="N14(p,γ)O15(,e⁺ν)N15(p,γ)O16")

    ro16_2p_n14_a = ModifiedRate(ro16pg,
                                 new_products=[Nucleus("n14"),
                                               Nucleus("he4")],
                                 stoichiometry={Nucleus("p"): 2},
                                 description="O16(p,γ)F17(,e⁺ν)O17(p,α)N14")

    lib_cno = Library(rates=[rc12_2p_n14, rn14_2p_c12,
                             rn14_2p_o16, ro16_2p_n14_a])

    # He burning

    r3a = rl.get_rate_by_name("a(aa,g)c12")

    rc12ag, _ = make_ap_pg_rates(rl, "c12", "o16",
                                 use_detailed_balance=True)

    rn14ag = rl.get_rate_by_name("n14(a,g)f18")
    rn14ag_lite = ModifiedRate(rn14ag,
                               new_products=["ne20"],
                               stoichiometry={Nucleus("he4"): 1.5})

    ro16ag, _ = make_ap_pg_rates(rl, "o16", "ne20",
                                 use_detailed_balance=True)
    rne20ag, _ = make_ap_pg_rates(rl, "ne20", "mg24",
                                  use_detailed_balance=True)

    lib_he = Library(rates=[r3a, rc12ag,
                            rn14ag_lite, ro16ag, rne20ag])

    net_class = get_net_class(network_type=network_type)

    return net_class(libraries=[lib_pp, lib_cno, lib_he],
                     allowed_duplicates=[(rppp_he3, rpepp_he3)])
