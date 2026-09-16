"""Some helper functions for determining which rates need screening"""


def get_screening_pair_set(rates):
    """Create a set of unique screening pairs for a list of rates.

    Parameters
    ----------
    rates : Iterable(Rate)
        A list of the rates in our network.  We expect this to be the
        complete list of rates (including any child rates of composite
        rate types), e.g., ``RateCollection.all_rates``

    Returns
    -------
    set(tuple(Nucleus, Nucleus))

    """

    # Create a full set of unique screening pairs across all rates
    unique_pairs = set()
    for r in rates:
        for scn_pair in r.screening_pairs:
            unique_pairs.add(scn_pair)

    return unique_pairs
