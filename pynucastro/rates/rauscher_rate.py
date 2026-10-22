""" Classes and methods to describe reaction rates using Rasucher's 2025 SMARAGD calculations."""

from pynucastro.rates.reaclib_rate import ReacLibRate, SingleSet


class RauscherRate(ReacLibRate):
    """A forward reaction rate using Rauscher fits in Reaclib format.

    Parameters
    ----------
    reactants : list(str), list(Nucleus)
        the reactants for the reactions
    products : list(str), list(Nucleus)
        the products for the reactions
    sets : list(SingleSet)
        the sets that make up the reaction rate
    Q : float
        energy release in MeV
    """

    def __init__(self, reactants, products, sets, Q):
        super().__init__(
            reactants=reactants,
            products=products,
            sets=sets,
            Q=Q,
            label="rauscher",
            rate_source="ra25",
        )

    @classmethod
    def from_record(cls, header, forward_line):
        """Construct a forward rate from Rauscher data.

        Parameters
        ----------
        header, forward_line : str
            the first line of every reaction is header,
            second line carries forward rate coefficients.

        Returns
        -------
        RauscherRate
            this gives the forward reaction rate

        """

        columns = header.split()
        if len(columns) != 15:
            raise ValueError("Expected 15 columns in the header")

        target, projectile, ejectile, final_nucleus = columns[:4]
        Q = float(columns[-1])

        coefficients = [float(value) for value in forward_line.split()]
        if len(coefficients) != 7:
            raise ValueError("Expected 7 coefficients")

        products = [final_nucleus]
        if ejectile != "gam":
            products.append(ejectile)

        fit = SingleSet(coefficients, labelprops="ra25  ")
        return cls(
            reactants=[target, projectile],
            products=products,
            sets=[fit],
            Q=Q,
        )
