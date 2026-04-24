from pynucastro.nucdata import Nucleus
from pynucastro.rates.library import _rate_name_to_nuc


def make_CO_approximation(all_rates, root_nuclei):
    """We want to model a sequence like:

                     X
              ....  .
          ...      .   .
      ...          .     .
    A ............ B ..... C

    This comes up in a few cases:

    * C + C: A = C12, B = Ne20, C = Ng24, X = Na23
    * C + O: A = C12 + O16, B = Mg24, C = Si28, X = Al27
    * O + O: A = O16, B = Si28, C = S32, X = P31

    We start with 5 rates (+ their inverses):

    * λ1 = A(A,p)X
    * λ2 = A(A,α)B
    * λ3 = X(p,α)B
    * λ4 = X(p,γ)C
    * λ5 = B(α,γ)C

    and we want to combine these into 3 effective rates (and their inverses)
    by eliminating X assuming equilibrium.

    The resulting effective rates are:

    * A(A,p)X(p,γ)C
    * A(A,α)B + A(A,p)X(p,α)B
    * B(α,γ) + B(α,p)X(p,γ)C

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
        A1 = Nucleus("c12")
        A2 = Nucleus("c12")
    elif root_nuclei == "CO":
        A1 = Nucleus("c12")
        A2 = Nucleus("o16")
    else:
        A1 = Nucleus("o16")
        A2 = Nucleus("o16")

    # C is the heaviest nucleus we get, from direct A + A fusion
    C = A1 + A2

    # B is an alpha lighter
    B = C - Nucleus("he4")

    # X is a proton lighter than C
    X = C - Nucleus("p")

    print(A1, A2, B, C, X)

    # now get the rates
    rates = {}

    lambda1 = f"{A1}({A2},p){X}"
    reactants, products = _rate_name_to_nuc(lambda1)
    rates["A(A,p)X"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]

    lambda2 = f"{A1}({A2},a){B}"
    reactants, products = _rate_name_to_nuc(lambda2)
    rates["A(A,a)B"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]

    lambda3 = f"{X}(p,a){B}"
    reactants, products = _rate_name_to_nuc(lambda3)
    rates["X(p,a)B"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]

    lambda4 = f"{X}(p,g){C}"
    reactants, products = _rate_name_to_nuc(lambda4)
    rates["X(p,g)C"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]

    lambda5 = f"{B}(a,g){C}"
    reactants, products = _rate_name_to_nuc(lambda5)
    rates["B(a,g)C"] = [r for r in all_rates
                        if sorted(r.reactants) == sorted(reactants) and
                           sorted(r.products) == sorted(products)][0]

    for k, v in rates.items():
        print(k, v)
