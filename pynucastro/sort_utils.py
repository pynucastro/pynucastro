"""Functions used for sorting rates"""


class CircularRateDependency(Exception):
    """A circular dependency is detected in a sequence of rates"""


def topo_sort(rates):
    """Perform a topological sort on a list of rates.

    ``topo_sort`` takes a list of ``Rate`` objects and does a
    depth-first search using the child rates as dependencies to
    produce a topologically-sorted ordering.

    Parameters
    ----------
    rates : Iterable(Rate)
        The collection of rates to sort

    Returns
    -------
    list(Rate)

    """

    visited = set()

    # visiting will keep track of the rates we are currently
    # exploring to detect circular dependencies
    visiting = set()

    order = []

    def visit(rate):
        if rate in visiting:
            raise CircularRateDependency(f"circular dependency with {rate}")

        if rate in visited:
            return

        visiting.add(rate)

        if crates := rate.get_child_rates():
            for child_rate in crates:
                visit(child_rate)

        visiting.remove(rate)
        visited.add(rate)
        order.append(rate)

    for rate in rates:
        visit(rate)

    return order
