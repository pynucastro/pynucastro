import collections

# there are some exceptions to the no-duplicate rates restriction.  We
# list them here by class name and then fname
ALLOWED_DUPLICATES = [
    {"ReacLibRate: p_p__d__weak__bet_pos_",
     "ReacLibRate: p_p__d__weak__electron_capture"}
]


def find_duplicate_rates(rate_list):
    """given a list of rates, return a list of groups of duplicate rates"""

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
    """rate_list is a list of rates that provide the same connection
    in a network.  Return True if this is an allowed duplicate"""

    # make rate_list into a set of strings in the same format as
    # ALLOWED_DUPLICATES, then check if it matches any of the allowed sets
    key_set = {f"{r.__class__.__name__}: {r.fname}" for r in rate_list}
    return key_set in ALLOWED_DUPLICATES
