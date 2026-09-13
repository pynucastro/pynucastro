"""Functions used for sorting rates"""

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
    order = []

    def visit(rate):
        if rate in visited:
            return

        visited.add(rate.name)

        if rate.get_child_rates():
            for child_rate in rate.get_child_rates():
                visit(child_rate)

        order.append(rate)

    for rate in rates:
        visit(rate)

    return order
