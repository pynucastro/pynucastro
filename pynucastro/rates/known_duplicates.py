# there are some exceptions to the no-duplicate rates restriction.  We
# list them here by class name and then fname
ALLOWED_DUPLICATES = [
    (("ReacLibRate: p_p__d__weak__bet_pos_"),
     ("ReacLibRate: p_p__d__weak__electron_capture"))
]


def find_duplicate_rates(rate_list):
    """given a list of rates, return a list of groups of duplicate
    rates"""

    duplicates = []
    for rate in rate_list:
        same_links = [q for q in rate_list
                      if q != rate and
                      sorted(q.reactants) == sorted(rate.reactants) and
                      sorted(q.products) == sorted(rate.products)]

        if same_links:
            new_entry = [rate] + same_links
            already_found = False
            # we may have already found this pair
            for dupe in duplicates:
                if new_entry[0] in dupe:
                    already_found = True
                    break
            if not already_found:
                duplicates.append(new_entry)

    return duplicates


def is_allowed_dupe(rate_list):
    """rate_list is a list of rates that provide the same connection
    in a network.  Return True if this is an allowed duplicate"""

    for allowed_dupe in ALLOWED_DUPLICATES:
        found = 0
        if len(rate_list) == len(allowed_dupe):
            found = 1
            for r in rate_list:
                rate_key = f"{r.__class__.__name__}: {r.fname}"
                if rate_key in allowed_dupe:
                    found *= 1
                else:
                    found *= 0
            if found:
                return True

    return False


class RateDuplicationError(Exception):
    """An error of multiple rates linking the same nuclei occurred"""
