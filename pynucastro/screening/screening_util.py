"""Some helper functions for determining which rates need screening"""

def get_screening_pair_set(rates):
    """Create a set of unique screening pairs for a list of rates.

    Parameters
    ----------
    rates : Iterable(Rate)
        A list of the rates in our network.

    Returns
    -------
    set(tuple(Nucleus, Nucleus))

    """

    # we need to consider the child rates that come with ApproximateRate
    all_rates = []
    for r in rates:
        # check if rate is an ApproximateRate
        # by checking if it has attribute approx_type
        # Don't do isinstance(r, ApproximateRate) to avoid circular dependency
        if hasattr(r, "approx_type"):
            all_rates += r.get_child_rates()
        else:
            all_rates.append(r)

    # Create a full set of unique screening pairs across all rates
    unique_pairs = set()
    for r in all_rates:
        for scn_pair in r.screening_pairs:
            unique_pairs.add(scn_pair)

    return unique_pairs
