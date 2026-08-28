"""Tools for detecting if a list of rates contains duplicates."""

import collections

# put any exceptions to the no-duplicate rates restriction.
# For example, a set like
#  {"ReacLibRate: p + p --> d <reaclib_bet+>",
#   "ReacLibRate: p + p --> d <reaclib_ec>"},
#
# Different weak rate types are already allowed, so
# the above does not need to be listed here.
ALLOWED_DUPLICATES = []


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

    # now check to see if the duplicates are allowed weak rate pairs
    remove = []
    for dupe in duplicates:
        if len(dupe) != 2:
            continue
        if ((dupe[0].weak_type == "electron_capture" and dupe[1].weak_type == "beta_pos") or
            (dupe[1].weak_type == "electron_capture" and dupe[0].weak_type == "beta_pos")):
            # e⁻-capture and β⁺-decay have the same reactants and products, but
            # represent 2 distinct processes.  If they are separated out, then we
            # assume that they are designed to be used in tandem and of the same type
            if type(dupe[0]) == type(dupe[1]):
                remove.append(dupe)

    for dupe in remove:
        duplicates.remove(dupe)

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
    key_set = {f"{r.__class__.__name__}: {r.id}" for r in rate_list}
    return key_set in ALLOWED_DUPLICATES
