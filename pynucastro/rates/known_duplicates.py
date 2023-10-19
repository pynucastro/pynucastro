# there are some exceptions to the no-duplicate rates restriction.  We
# list them here by class name and then fname
ALLOWED_DUPLICATES = [
    (("ReacLibRate: p_p__d__weak__bet_pos_"),
     ("ReacLibRate: p_p__d__weak__electron_capture"))
]


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
