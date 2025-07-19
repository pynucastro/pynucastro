"""Some helper functions for determining which rates need screening"""

from pynucastro.nucdata import Nucleus
from pynucastro.rates import ApproximateRate


class ScreeningPair:
    """A pair of nuclei that will have rate screening applied.  We
    store a list of all rates that match this pair of nuclei

    """

    def __init__(self, name, nuc1, nuc2, rate=None):
        self.name = name
        self.n1 = nuc1
        self.n2 = nuc2

        if rate is None:
            self.rates = []
        else:
            self.rates = [rate]

    def add_rate(self, rate):
        """Add a new rate to the screening pair.

        Parameters
        ----------
        rate : Rate
            The rate to add.

        """

        if rate not in self.rates:
            self.rates.append(rate)

    def __str__(self):
        ostr = f"screening for {self.n1} + {self.n2}\n"
        ostr += "rates:\n"
        for r in self.rates:
            ostr += f"  {r}\n"
        return ostr

    def __eq__(self, other):
        """Test for equality.  All we care about is whether the names
        are the same -- that conveys what the reaction is

        """

        return self.name == other.name


def get_screening_map(rates, *, symmetric_screening=False):
    """Create a screening map---this is just a list of ScreeningPair
    objects containing the information about nuclei pairs for
    screening If symmetric_screening=True, then for reverse rates, we
    screen using the forward rate nuclei (assuming that we got here
    via detailed balance).

    Parameters
    ----------
    rates : Iterable(Rate)
        A list of the rates in our network.
    symmetric_screening : bool
        Do we use the same screening factor for the forward and reverse
        rate of a single reaction process?

    Returns
    -------
    list(ScreeningPair)

    """
    screening_map = []

    # we need to consider the child rates that come with ApproximateRate
    all_rates = []
    for r in rates:
        if isinstance(r, ApproximateRate):
            all_rates += r.get_child_rates()
        else:
            all_rates.append(r)

    for r in all_rates:
        screen_nuclei = r.ion_screen
        if symmetric_screening:
            screen_nuclei = r.symmetric_screen

        # screen_nuclei may be [] if it is a decay, gamma-capture, or
        # neutron-capture
        if not screen_nuclei:
            continue

        nucs = "_".join([str(q) for q in screen_nuclei])

        scr = [q for q in screening_map if q.name == nucs]

        assert len(scr) <= 1

        if scr:
            # we already have the reactants in our map, so we
            # will already be doing the screening factors.
            # Just append this new rate to the list we are
            # keeping of the rates where this screening is
            # needed -- if the rate is already in the list, then
            # this is a no-op

            scr[0].add_rate(r)

            # if we got here because nuc == "He4_He4_He4",
            # then we also have to add to "He4_He4_He4_dummy"

            if nucs == "He4_He4_He4":
                scr2 = [q for q in screening_map if q.name == nucs + "_dummy"]
                assert len(scr2) == 1

                scr2[0].add_rate(r)

        else:

            # we handle 3-alpha specially -- we actually need
            # 2 screening factors for it

            if nucs == "He4_He4_He4":
                # he4 + he4
                scr1 = ScreeningPair(nucs, screen_nuclei[0], screen_nuclei[1], r)

                # he4 + be8
                be8 = Nucleus("Be8", dummy=True)
                scr2 = ScreeningPair(nucs + "_dummy", screen_nuclei[2], be8, r)

                screening_map.append(scr1)
                screening_map.append(scr2)

            else:
                scr1 = ScreeningPair(nucs, screen_nuclei[0], screen_nuclei[1], r)
                screening_map.append(scr1)

    return screening_map
