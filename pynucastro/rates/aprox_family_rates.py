"""Functions to help construct rate approximations used in the
aprox13, aprox19, aprox21 family of reaction networks

"""

from pynucastro.nucdata import Nucleus
from pynucastro.rates.approximate_rates import ApproximateRate
from pynucastro.rates.library import _rate_name_to_nuc


def make_CO_approx_rates(all_rates, root_nuclei,
                         return_obsolete_rate_names=False):
    """We want to model a sequence like:

    ::

       N - Z
       ^
       |                      X
       |                 ... . .
       |            ...     .   .
       |       ...         .     .
       |    S ............E ..... F
       |
       +------------------------------> Z


    This comes up in a few cases:

    * C + C: S = C12, E = Ne20, F = Ng24, X = Na23
    * C + O: S = C12 + O16, E = Mg24, F = Si28, X = Al27
    * O + O: S = O16, E = Si28, F = S32, X = P31

    We start with 5 rates (+ their inverses):

    * λ1 = S(S,p)X
    * λ2 = S(S,α)E
    * λ3 = X(p,α)E
    * λ4 = X(p,γ)F
    * λ5 = E(α,γ)F

    and we want to combine these into 3 effective rates (and their inverses)
    by eliminating X assuming equilibrium.

    The resulting effective rates are:

    * S(S,p)X(p,γ)F
    * S(S,α)E + S(S,p)X(p,α)E
    * E(α,γ)F + E(α,p)X(p,γ)F

    Parameters
    ----------
    all_rates : list(Rate)
         the collection of rates that we will use for finding λ1 through λ5 and
         their inverses.  This could be via
         :py:func:`Library.get_rates <pynucastro.rates.library.Library.get_rates>`
         or :py:func:`RateCollection.get_rates <pynucastro.networks.rate_collection.RateCollection.get_rates>`.
    root_nuclei : str
         the identity of nucleus A in this approximation.  This can be
         "C", "CO" or "O"
    return_obsolete_rate_names : bool
         Do we return a list of rates names that can be removed if we use
         these approximate rates?

    Returns
    -------
    approximate_rates : list(ApproximateRate)
        A list of the approximated rates
    obsolete_rates : list(str), optional
        A list of the rate names, in the form "A(x,y)B" that are
        obsoleted by the new approximate rates.

    """

    assert root_nuclei in ["C", "CO", "O"]

    # create the nuclei
    if root_nuclei == "C":
        S1 = Nucleus("c12")
        S2 = Nucleus("c12")
    elif root_nuclei == "CO":
        S1 = Nucleus("c12")
        S2 = Nucleus("o16")
    else:
        S1 = Nucleus("o16")
        S2 = Nucleus("o16")

    # F is the heaviest nucleus we get, from direct S + S fusion
    F = S1 + S2

    # E is an alpha lighter
    E = F - Nucleus("he4")

    # X is a proton lighter than C
    X = F - Nucleus("p")

    print(S1, S2, E, F, X)

    # now get the rates (including inverses)
    rates = {}

    lambda1 = f"{S1}({S2},p){X}"
    lambda1_r = f"{X}(p,{S2}){S1}"
    reactants, products = _rate_name_to_nuc(lambda1)
    rates["S(S,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["X(p,S)S"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda2 = f"{S1}({S2},a){E}"
    lambda2_r = f"{E}(a,{S2}){S1}"
    reactants, products = _rate_name_to_nuc(lambda2)
    rates["S(S,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["E(a,S)S"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda3 = f"{X}(p,a){E}"
    lambda3_r = f"{E}(a,p){X}"
    reactants, products = _rate_name_to_nuc(lambda3)
    rates["X(p,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["E(a,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda4 = f"{X}(p,g){F}"
    lambda4_r = f"{F}(g,p){X}"
    reactants, products = _rate_name_to_nuc(lambda4)
    rates["X(p,g)F"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["F(g,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda5 = f"{E}(a,g){F}"
    lambda5_r = f"{F}(g,a){E}"
    reactants, products = _rate_name_to_nuc(lambda5)
    rates["E(a,g)F"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["F(g,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    new_rates = []

    # S(S,p)X(p,γ)F (e.g., C12(C12,p)Na23(p,γ)Mg24)
    # this also needs X(p,α)E (e.g. Na23(p,α)Ne20) as an alternate
    # branching for normalization

    # for ApproximateRate, we write these in the form A <--> B
    # where A = S, B = F, C = E
    child_rates = {}

    child_rates["A(Y,p)X"] = rates["S(S,p)X"]
    child_rates["X(p,Y)A"] = rates["X(p,S)S"]

    child_rates["X(p,g)B"] = rates["X(p,g)F"]
    child_rates["B(g,p)X"] = rates["F(g,p)X"]

    child_rates["X(p,a)C"] = rates["X(p,a)E"]

    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=False, approx_type="Yp_pg"))
    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=True, approx_type="Yp_pg"))

    # S(S,α)E + S(S,p)X(p,α)E (e.g. C12(C12,α)Ne20 + C12(C12,p)Na23(p,α)Ne20)
    # this also needs X(p,g)F (e.g. Na23(p,γ)Mg24) as an alternate
    # branching for normalization

    # for ApproximateRate, we write these in the form A <--> B
    # where A = S, B = E, C = F
    child_rates = {}

    child_rates["A(Y,a)B"] = rates["S(S,a)E"]
    child_rates["B(a,Y)A"] = rates["E(a,S)S"]

    child_rates["A(Y,p)X"] = rates["S(S,p)X"]
    child_rates["X(p,Y)A"] = rates["X(p,S)S"]

    child_rates["X(p,a)B"] = rates["X(p,a)E"]
    child_rates["B(a,p)X"] = rates["E(a,p)X"]

    child_rates["X(p,g)C"] = rates["X(p,g)F"]

    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=False, approx_type="Yp_pa"))
    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=True, approx_type="Yp_pa"))

    # E(α,γ)F + E(α,p)X(p,γ)F (e.g. Ne20(α,γ)Mg24 + Ne20(α,p)Na23(p,γ)Mg24)
    # with X(p,S)S (e.g. Na23(p,C12)C12) as an alternate branch for
    # normalization.

    # for ApproximateRate, we write these in the form A <--> B
    # where A = E, B = F, and C = S
    child_rates = {}

    child_rates["A(a,g)B"] = rates["E(a,g)F"]
    child_rates["B(g,a)A"] = rates["F(g,a)E"]

    child_rates["A(a,p)X"] = rates["E(a,p)X"]
    child_rates["X(p,a)A"] = rates["X(p,a)E"]

    child_rates["X(p,g)B"] = rates["X(p,g)F"]
    child_rates["B(g,p)X"] = rates["F(g,p)X"]

    child_rates["X(p,Y)C"] = rates["X(p,S)S"]

    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=False, approx_type="ap_pg"))
    new_rates.append(ApproximateRate(child_rates,
                                     is_reverse=True, approx_type="ap_pg"))

    obsolete_rates = [lambda1, lambda2, lambda3, lambda4, lambda5,
                      lambda1_r, lambda2_r, lambda3_r, lambda4_r, lambda5_r]

    if return_obsolete_rate_names:
        return new_rates, obsolete_rates
    return new_rates
