"""Tools for detecting if a list of rates contains duplicates."""

import collections

# there are some exceptions to the no-duplicate rates restriction.  We
# list them here by class name and then fname
ALLOWED_DUPLICATES = [
    {"ReacLibRate: p_p_to_d_weak_bet_pos_",
     "ReacLibRate: p_p_to_d_weak_electron_capture"}
]


def find_duplicate_rates(rate_list):
    """Given a list of rates, return a list of groups of duplicate
    rates

    Parameters
    ----------
    rate_list : list(Rate)
        the input list of rates

    Returns
    -------
    list(Rate)

    """

    # Group the rates into lists of potential duplicates, keyed by their
    # reactants and products.
    grouped_rates = collections.defaultdict(list)
    for rate in rate_list:
        grouped_rates[tuple(sorted(rate.reactants)),
                      tuple(sorted(rate.products))].append(rate)

    # any entry in grouped_rates containing more than one rate is a duplicate
    duplicates = [entry for entry in grouped_rates.values() if len(entry) > 1]

    return duplicates


def is_allowed_dupe(rate_list):
    """Check if any of the duplicates are allowed.  Some duplicates
    may be allowed since they represent distinct processes between the
    same endpoints (those rates are listed in ``ALLOWED_DUPLICATES``).
    Return `True` is the input set of rates is an allowed duplicate.

    Parameters
    ----------
    rate_list : list(Rate)
        a group of rates represented the same sequence that may be
        allowed duplicates.

    Returns
    -------
    bool

    """

    # make rate_list into a set of strings in the same format as
    # ALLOWED_DUPLICATES, then check if it matches any of the allowed sets
    key_set = {f"{r.__class__.__name__}: {r.fname}" for r in rate_list}
    return key_set in ALLOWED_DUPLICATES
