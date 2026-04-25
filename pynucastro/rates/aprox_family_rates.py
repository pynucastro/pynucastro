from pynucastro.nucdata import Nucleus
from pynucastro.rates.library import _rate_name_to_nuc
from pynucastro.rates.approximate_rates import ApproximateRate


def make_CO_approximation(all_rates, root_nuclei):
    """We want to model a sequence like:

                     X
              ....  .
          ...      .   .
      ...          .     .
    S ............E ..... F

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
         or :py:func:`RateCollection.get_rates <pynucastro.rates.rate_collection.RateCollection.get_rates>`.
    root_nuclei : str
         the identity of nucleus A in this approximation.  This can be
         "C", "CO" or "O"

    Returns
    -------
    list(ApproximateRate)

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
    reactants, products = _rate_name_to_nuc(lambda1)
    rates["S(S,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["X(p,S)S"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda2 = f"{S1}({S2},a){E}"
    reactants, products = _rate_name_to_nuc(lambda2)
    rates["S(S,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["E(a,S)S"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda3 = f"{X}(p,a){E}"
    reactants, products = _rate_name_to_nuc(lambda3)
    rates["X(p,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["E(a,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda4 = f"{X}(p,g){F}"
    reactants, products = _rate_name_to_nuc(lambda4)
    rates["X(p,g)F"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["F(g,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    lambda5 = f"{E}(a,g){F}"
    reactants, products = _rate_name_to_nuc(lambda5)
    rates["E(a,g)F"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]
    rates["F(g,a)E"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(products) and
                           sorted(r.products) == sorted(reactants)][0]

    new_rates = []

    # S(S,p)X(p,γ)F
    child_rates = [rates["S(S,p)X"], rates["X(p,g)F"],
                   rates["X(p,S)S"], rates["F(g,p)X"]]
    new_rates.append(ApproximateRate(child_rates, is_reverse=False, approx_type="Yp_pg"))
    new_rates.append(ApproximateRate(child_rates, is_reverse=True, approx_type="Yp_pg"))

    for k, v in rates.items():
        print(k, v)
