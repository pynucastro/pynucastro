"""A collection of classes and methods to deal with collections of
rates that together make up a network.

"""

import collections
import functools
import math
import warnings
from operator import mul
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from ipywidgets import interact
from matplotlib.colors import SymLogNorm
from matplotlib.patches import ConnectionPatch
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.linalg import eigvals

# Import Rate
from pynucastro.constants import constants
from pynucastro.nucdata import Nucleus
from pynucastro.rates import (ApproximateRate, DerivedRate, Library,
                              ModifiedRate, Rate, RateFileError, RatePair,
                              TabularRate, find_duplicate_rates,
                              is_allowed_dupe, load_rate)
from pynucastro.rates.library import _rate_name_to_nuc, capitalize_rid
from pynucastro.screening import (get_screening_map, make_plasma_state,
                                  make_screen_factors)

mpl.rcParams['figure.dpi'] = 100

# for plotting a legend on the network plot
# the tuple is (dZ, dN)
RATE_LINES = {r"$(\alpha, p)$": (1, 2),
              r"$(\alpha, \gamma)$": (2, 2),
              r"$(\alpha, n)$": (2, 1),
              r"$(p, \gamma)$": (1, 0),
              r"$(n, \gamma)$": (0, 1),
              r"$\beta^-$": (1, -1),
              r"$\beta^+$": (-1, 1)}


class RateDuplicationError(Exception):
    """An error of multiple rates linking the same nuclei occurred"""


def _skip_xalpha(n, p, r):
    """Check if we should show an (a, x) or (x, a) rate.  Here, p is
    the product we want to link to

    """

    # first check if alpha is the heaviest nucleus on the RHS
    rhs_heavy = max(r.products)
    if not (rhs_heavy.Z == 2 and rhs_heavy.A == 4):

        # for rates that are A (x, alpha) B, where A and B are heavy nuclei,
        # don't show the connection of the nucleus to alpha, only show it to B
        if p.Z == 2 and p.A == 4:
            return True

        # likewise, hide A (alpha, x) B, unless A itself is an alpha
        c = r.reactants
        n_alpha = 0
        for nuc in c:
            if nuc.Z == 2 and nuc.A == 4:
                n_alpha += 1
        # if there is only 1 alpha and we are working on the alpha node,
        # then skip
        if n_alpha == 1 and n.Z == 2 and n.A == 4:
            return True

    return False


def _skip_xp(n, p, r):
    """Check if we should show an (p, x) or (x, p) rate.  Here, p is
    the product we want to link to

    """

    # for rates that are A (x, p) B, where A and B are heavy nuclei,
    # don't show the connection of the nucleus to p, only show it to B
    if p.Z == 1 and p.A == 1:
        return True

    # likewise, hide A (p, x) B, unless A itself is an p
    c = r.reactants
    n_p = 0
    for nuc in c:
        if nuc.Z == 1 and nuc.A == 1:
            n_p += 1
    # if there is only 1 p and we are working on the p node,
    # then skip
    if n_p == 1 and n.Z == 1 and n.A == 1:
        return True

    return False


class Composition(collections.UserDict):
    """A composition holds the mass fractions of the nuclei in a
    network.

    Parameters
    ----------
    nuclei : list, tuple
        an iterable of Nucleus objects
    small : float
        a floor for nuclei mass fractions, used as the default value

    """

    def __init__(self, nuclei, small=1.e-16):
        try:
            super().__init__({Nucleus.cast(k): small for k in nuclei})
        except TypeError:
            raise ValueError("must supply an iterable of Nucleus objects or strings") from None

    @property
    def X(self):
        """backwards-compatible getter for self.X"""
        return self.data

    @X.setter
    def X(self, new_value):
        """backwards-compatible setter for self.X"""
        self.data = new_value

    def __delitem__(self, key):
        super().__delitem__(Nucleus.cast(key))

    def __getitem__(self, key):
        return super().__getitem__(Nucleus.cast(key))

    def __setitem__(self, key, value):
        super().__setitem__(Nucleus.cast(key), value)

    def __repr__(self):
        return "Composition(" + super().__repr__() + ")"

    def __str__(self):
        return "".join(f"  X({k}) : {v}\n" for k, v in self.items())

    @property
    def A(self):
        """Nucleus molar masses

        Returns
        -------
        A : dict
            {Nucleus : A} pairs
        """
        return {n: n.A for n in self}

    @property
    def Z(self):
        """Nucleus charge

        Returns
        -------
        Z : dict
            {Nucleus : Z} pairs
        """
        return {n: n.Z for n in self}

    def get_nuclei(self):
        """Return a list of Nuclei objects that make up this
        composition.

        Returns
        -------
        list

        """
        return list(self)

    def get_molar(self):
        """Return a dictionary of molar fractions, Y = X/A.

        Returns
        -------
        molar : dict
            {Nucleus : Y}
        """
        return {k: v/k.A for k, v in self.items()}

    def get_sum_X(self):
        """Return the sum of the mass fractions.

        Returns
        -------
        float
        """
        return math.fsum(self.values())

    def set_solar_like(self, *, Z=0.02, half_life_thresh=None):
        """Approximate a solar abundance, setting p to 0.7, He4 to 0.3
        - Z and the remainder evenly distributed with Z.

        Parameters
        ----------
        Z : float
            The desired metalicity
        half_life_thresh : float
            The half life value below which to zero the mass fraction
            of a nucleus.  This prevents us from making a composition
            that is not really stable.

        """

        rem = Z/(len(self)-2)
        for k in self:
            if k == Nucleus("p"):
                self[k] = 0.7
            elif k.raw == "he4":
                self[k] = 0.3 - Z
            else:
                self[k] = rem

        self.normalize(half_life_thresh=half_life_thresh)

    def set_array(self, arr):
        """Set the mass fractions of all species to the values
        in arr, `get_nuclei()`

        Parameters
        ----------
        arr : list, tuple, numpy.ndarray
            input values of mass fractions
        """
        for i, k in enumerate(self):
            self[k] = arr[i]

    def set_all(self, xval: float):
        """Set all species to the same scalar value.

        Parameters
        ----------
        xval : float
            mass fraction value for all species
        """
        for k in self:
            self[k] = xval

    def set_equal(self):
        """Set all species to be equal"""
        self.set_all(1.0 / len(self))

    def set_random(self, alpha=None, seed=None):
        """Set all species using a Dirichlet distribution with
        parameters alpha and specified rng seed.

        Parameters
        ----------
        alpha : list, tuple, numpy.ndarray
            distribution length for the Dirichlet distribution
        seed : float
            seed for the random number generator
        """

        # initializes random seed
        rng = np.random.default_rng(seed)

        # default is a flat Dirichlet distribution
        if alpha is None:
            alpha = np.ones(len(self))

        fracs = rng.dirichlet(alpha)
        self.set_array(fracs)

        # ensures exact normalization
        self.normalize()

    def set_nuc(self, name, xval: float):
        """Set nuclei name to the mass fraction xval.

        Parameters
        ----------
        name : Nucleus
            the nucleus to set
        xval: float
        """
        self[name] = xval

    def normalize(self, *, half_life_thresh=None):
        """Normalize the mass fractions to sum to 1.

        Parameters
        ----------
        half_life_thresh : float
            The half life value below which to zero the mass fraction
            of a nucleus.  This prevents us from making a composition
            that is not really stable.

        """

        if half_life_thresh is not None:
            for k in self:
                if k.tau != "stable" and k.tau is not None:
                    if k.tau < half_life_thresh:
                        self[k] = 0.0

        X_sum = self.get_sum_X()

        for k in self:
            self[k] /= X_sum

    @property
    def ye(self):
        """Return the electron fraction of the composition

        Returns
        -------
        float
        """
        electron_frac = math.fsum(self[n] * n.Z / n.A for n in self) / self.get_sum_X()
        return electron_frac

    @property
    def abar(self):
        """Return the mean molecular weight

        Returns
        -------
        float
        """
        abar = math.fsum(self[n] / n.A for n in self)
        return 1. / abar

    @property
    def zbar(self):
        """Return the mean charge, Zbar

        Returns
        -------
        float
        """
        return self.abar * self.ye

    def bin_as(self, nuclei, *, verbose=False, exclude=None):
        """Given a list of nuclei, return a new Composition object
        with the current composition mass fractions binned into the
        new nuclei.

        Parameters
        ----------
        nuclei : list
            Input nuclei (either as string names or
            Nucleus objects) defining the new composition.
        verbose : bool
            Output more information
        exclude : bool
            List of nuclei in `nuclei` that only
            exact matches from the original composition can
            map into

        Returns
        -------
        new_composition : Composition
            The new binned composition
        """

        nuclei = Nucleus.cast_list(nuclei)
        exclude = Nucleus.cast_list(exclude, allow_None=True)

        # sort the input nuclei by A, then Z
        nuclei.sort(key=lambda n: (n.A, n.Z))

        # create the new composition
        new_comp = Composition(nuclei)

        # first do any exact matches if we provided an exclude list
        if exclude is None:
            exclude = []

        for ex_nuc in exclude:
            # if the exclude nucleus is in both our original
            # composition and the reduced composition, then set
            # the abundance in the new, reduced composition and
            # remove the nucleus from consideration for the other
            # original nuclei
            if ex_nuc in nuclei and ex_nuc in self:
                nuclei.remove(ex_nuc)
                new_comp[ex_nuc] = self[ex_nuc]
                if verbose:
                    print(f"storing {ex_nuc} as {ex_nuc}")

            else:
                raise ValueError("cannot use exclude if nucleus is not present in both the original and new compostion")

        # loop over our original nuclei.  Find the new nucleus such
        # that n_orig.A >= n_new.A.  If there are multiple, then do
        # the same for Z
        for old_n, v in self.items():

            if old_n in exclude:
                # we should have already dealt with this above
                continue

            candidates = [q for q in nuclei if old_n.A >= q.A]
            # if candidates is empty, then all of the nuclei are heavier than
            # old_n, so just put its composition in the first new nucleus
            # (which will be the lightest)
            if not candidates:
                match_nuc = nuclei[0]
            else:
                max_A = max(q.A for q in candidates)
                match_A = [q for q in candidates if q.A == max_A]
                if len(match_A) > 1:
                    match_Z = [q for q in sorted(match_A, key=lambda p: p.Z) if old_n.Z >= q.Z]
                    if not match_Z:
                        # our nucleus has a Z less than any of the Z's in match_A
                        match_nuc = match_A[0]
                    else:
                        # always take the last entry -- this way if
                        # match_Z has multiple nuclei, we are taking
                        # the one with the highest Z (since we
                        # initially sorted by A and Z)
                        match_nuc = match_Z[-1]
                else:
                    match_nuc = match_A[0]

            if verbose:
                print(f"storing {old_n} as {match_nuc}")
            new_comp[match_nuc] += v

        return new_comp

    def plot(self, trace_threshold=0.1, hard_limit=None, size=(9, 5)):
        """Make a pie chart of Composition. group trace nuclei
        together and explode into bar chart

        Parameters
        ----------
        trace_threshold : float
            the threshold to consider a component to be trace.
        hard_limit : float
            limit below which an abundance will not be included
            in the trace nuclei wedget of the plot.
        size: tuple
            width, height of the plot in inches

        Returns
        -------
        matplotlib.figure.Figure
        """

        # find trace nuclei
        trace_keys = []
        trace_tot = 0.
        main_keys = []
        for k in self:
            # if below threshold, count as trace element
            if self[k] < trace_threshold:
                trace_keys.append(k)
                trace_tot += self[k]
            else:
                main_keys.append(k)

        # check if any trace nuclei
        if not trace_keys:
            # just do pie chart without including trace

            fig, ax = plt.subplots(1, 1, figsize=size)

            ax.pie(self.values(), labels=self.keys(), autopct=lambda p: f"{p/100:0.3f}")

        else:
            # find trace nuclei which contribute little to trace proportion
            if hard_limit is None:
                # make hardlimit proportional to trace abundance
                hard_limit = 0.05*trace_tot

            limited_trace_keys = []
            other_trace_tot = 0.
            for k in trace_keys:
                if self[k] < hard_limit:
                    other_trace_tot += self[k]
                else:
                    limited_trace_keys.append(k)

            # make figure and assign axis objects
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=size)
            fig.subplots_adjust(wspace=0)

            # pie chart parameters
            main_values = [trace_tot] + [self[k] for k in main_keys]
            main_labels = ['trace'] + main_keys
            explode = [0.2] + [0. for i in range(len(main_keys))]

            # rotate so that first wedge is split by the x-axis
            angle = -180 * main_values[0]
            wedges, *_ = ax1.pie(main_values, autopct=lambda p: f"{p/100:0.3f}", startangle=angle,
                                labels=main_labels, explode=explode)

            # bar chart parameters
            trace_values = [self[k] for k in limited_trace_keys] + [other_trace_tot]
            trace_labels = [f"{k}" for k in limited_trace_keys] + ['other']
            bottom = 1
            width = 0.1

            # Adding from the top matches the legend.
            alpha_list = np.linspace(0.1, 1, len(trace_values))
            trace_wedge_color = wedges[0].get_facecolor()

            for j, (height, label) in enumerate([*zip(trace_values, trace_labels)]):
                bottom -= height
                bc = ax2.bar(0, height, width, bottom=bottom, color=trace_wedge_color, label=label,
                            alpha=alpha_list[j])

                ax2.bar_label(bc, labels=[f"{height:.2e}"], label_type='center')
                ax2.bar_label(bc, labels=[f"{label:>30}"], label_type='center')

            ax2.set_title('Composition of Trace Nuclei')
            ax2.axis('off')
            ax2.set_xlim(- 2.5 * width, 2.5 * width)

            # use ConnectionPatch to draw lines between the two plots
            theta1, theta2 = wedges[0].theta1, wedges[0].theta2
            center, r = wedges[0].center, wedges[0].r
            bar_height = sum(trace_values)

            # draw top connecting line
            x = r * np.cos(np.pi / 180 * theta2) + center[0]
            y = r * np.sin(np.pi / 180 * theta2) + center[1]
            con = ConnectionPatch(xyA=(-width / 2, bar_height+bottom), coordsA=ax2.transData,
                                xyB=(x, y), coordsB=ax1.transData)
            con.set_color(trace_wedge_color)
            con.set_linewidth(4)
            ax2.add_artist(con)

            # draw bottom connecting line
            x = r * np.cos(np.pi / 180 * theta1) + center[0]
            y = r * np.sin(np.pi / 180 * theta1) + center[1]
            con = ConnectionPatch(xyA=(-width / 2, bottom), coordsA=ax2.transData,
                                xyB=(x, y), coordsB=ax1.transData)
            con.set_color(trace_wedge_color)
            ax2.add_artist(con)
            con.set_linewidth(4)

        plt.show()
        return fig


class RateCollection:
    """A collection of rates that together define a network.
    There are several arguments to the constructor -- any combination
    may be supplied.

    Parameters
    ----------
    rate_files : str, list, tuple
        a string or iterable of strings of file names that define valid
        rates. This can include Reaclib library files storing multiple
        rates.
    libraries : Library, list, tuple
        a Library or iterable of Library objects
    rates : Rate, list, tuple
        a Rate or iterable of Rate objects
    inert_nuclei : list, tuple
        an iterable of Nuclei that should be part of the collection but
        are not linked via reactions to the other Nuclei in the network.
    symmetric_screening : bool
        symmetric screening means that we screen the reverse rates
        using the same factor as the forward rates, for rates computed
        via detailed balance.
    do_screening : bool
        should we consider screening at all -- this mainly affects
        whether we build the screening map
    verbose : bool
        do we show informational messages?

    """

    pynucastro_dir = Path(__file__).parents[1]

    def __init__(self, rate_files=None, libraries=None, rates=None,
                 inert_nuclei=None,
                 symmetric_screening=False, do_screening=True,
                 verbose=False):

        self.rates = []
        combined_library = Library()

        self.inert_nuclei = Nucleus.cast_list(inert_nuclei, allow_None=True)

        self.symmetric_screening = symmetric_screening
        self.do_screening = do_screening

        self.verbose = verbose

        if rate_files:
            if isinstance(rate_files, str):
                rate_files = [rate_files]
            combined_library += self._read_rate_files(rate_files)

        if rates:
            if isinstance(rates, Rate):
                rates = [rates]
            for r in rates:
                if not isinstance(r, Rate):
                    raise ValueError('Expected Rate object or list of Rate objects passed as the rates argument.')
            rlib = Library(rates=rates)
            combined_library += rlib

        if libraries:
            if isinstance(libraries, Library):
                libraries = [libraries]
            for lib in libraries:
                if not isinstance(lib, Library):
                    raise ValueError('Expected Library object or list of Library objects passed as the libraries argument.')
            for lib in libraries:
                combined_library += lib

        self.rates = self.rates + combined_library.get_rates()

        self._build_collection()

    def _build_collection(self):

        # get the unique nuclei
        u = []
        for r in self.rates:
            t = set(r.reactants + r.products)
            u = set(list(u) + list(t))

        self.unique_nuclei = sorted(u)

        # approx nuclei are used in approximate rates
        self.approx_nuclei = []
        for r in self.rates:
            if isinstance(r, ApproximateRate):
                if r.intermediate_nucleus not in self.unique_nuclei + self.approx_nuclei:
                    self.approx_nuclei.append(r.intermediate_nucleus)

        if self.inert_nuclei is not None:
            for nuc in self.inert_nuclei:
                if nuc not in self.unique_nuclei:
                    self.unique_nuclei.append(nuc)

        # now make a list of each rate that touches each nucleus
        # we'll store this in a dictionary keyed on the nucleus
        self.nuclei_consumed = {}
        self.nuclei_produced = {}

        for n in self.unique_nuclei:
            self.nuclei_consumed[n] = [r for r in self.rates if n in r.reactants]
            self.nuclei_produced[n] = [r for r in self.rates if n in r.products]

        self.nuclei_rate_pairs = {}
        _rp = self.get_rate_pairs()

        for n in self.unique_nuclei:
            self.nuclei_rate_pairs[n] = \
                [rp for rp in _rp if rp.forward is not None and n in rp.forward.reactants + rp.forward.products or
                                     rp.reverse is not None and n in rp.reverse.reactants + rp.reverse.products]

        # Re-order self.rates so Reaclib rates come first,
        # followed by Tabular rates. This is needed if
        # reaclib coefficients are targets of a pointer array.
        # It is desired to avoid wasting array size
        # storing meaningless Tabular coefficient pointers.
        self.rates = sorted(self.rates,
                            key=lambda r: r.chapter == 't')

        self.tabular_rates = []
        self.reaclib_rates = []
        self.custom_rates = []
        self.approx_rates = []
        self.derived_rates = []
        self.modified_rates = []

        for r in self.rates:
            if isinstance(r, ApproximateRate):
                self.approx_rates.append(r)
                for cr in r.get_child_rates():
                    assert cr.chapter != "t"

                    # Check whether this child rate is removed or not.
                    # "removed" means that this rate is never used on
                    # its own to connect two nuclei in the network it
                    # is only used in one or more ApproximateRate or
                    # ModifiedRate
                    if cr not in self.rates:
                        cr.removed = True
                    else:
                        cr.removed = False

                    cr.fname = None
                    # pylint: disable-next=protected-access
                    cr._set_print_representation()

                    # child rates may be ReacLibRates, ModifiedRates,
                    # or DerivedRates.  Make sure we don't double
                    # count
                    if isinstance(cr, DerivedRate):
                        if cr not in self.derived_rates:
                            self.derived_rates.append(cr)
                    elif isinstance(cr, ModifiedRate):
                        if cr not in self.modified_rates:
                            self.modified_rates.append(cr)
                    else:
                        if cr not in self.reaclib_rates:
                            self.reaclib_rates.append(cr)

            elif isinstance(r, ModifiedRate):
                if r not in self.modified_rates:
                    self.modified_rates.append(r)

                cr = r.original_rate

                # Check whether this child rate is removed or not.
                # "removed" means that this rate is never used on
                # its own to connect two nuclei in the network it
                # is only used in one or more ApproximateRate or
                # ModifiedRate
                if cr not in self.rates:
                    cr.removed = True
                else:
                    cr.removed = False

                cr.fname = None
                # pylint: disable-next=protected-access
                cr._set_print_representation()

                # child rates may be ReacLibRates, ModifiedRates,
                # or DerivedRates.  Make sure we don't double
                # count
                if isinstance(cr, DerivedRate):
                    if cr not in self.derived_rates:
                        self.derived_rates.append(cr)
                elif isinstance(cr, ModifiedRate):
                    if cr not in self.modified_rates:
                        self.modified_rates.append(cr)
                else:
                    if cr not in self.reaclib_rates:
                        self.reaclib_rates.append(cr)

            elif r.chapter == 't':
                self.tabular_rates.append(r)
            elif r.chapter == "custom":
                self.custom_rates.append(r)
            elif isinstance(r, DerivedRate):
                if r not in self.derived_rates:
                    self.derived_rates.append(r)
            elif isinstance(r.chapter, int):
                if r not in self.reaclib_rates:
                    self.reaclib_rates.append(r)
                    if r.id == "n --> p <wc12_reaclib_weak_>":
                        msg = "ReacLib neutron decay rate (<n_to_p_weak_wc12>) does not account for degeneracy at high densities. Consider using tabular rate from Langanke."
                        warnings.warn(msg)
            else:
                raise NotImplementedError(f"Chapter type unknown for rate chapter {r.chapter}")

        self.all_rates = (self.reaclib_rates + self.custom_rates +
                          self.tabular_rates + self.approx_rates +
                          self.modified_rates + self.derived_rates)

        # finally check for duplicate rates -- these are not
        # allowed
        if dupes := self.find_duplicate_links():
            print(dupes)
            raise RateDuplicationError("Duplicate rates found")

    def _read_rate_files(self, rate_files):
        # get the rates
        combined_library = Library()
        for rf in rate_files:
            # create the appropriate rate object first
            try:
                rate = load_rate(rf)
            except RateFileError as ex:
                raise RateFileError(f"Error reading rate from file: {rf}") from ex

            # now create a library:
            rflib = Library(rates=[rate])
            combined_library += rflib
        return combined_library

    def get_forward_rates(self):
        """Return a list of the forward (exothermic) rates in the
        network

        Returns
        -------
        list(Rate)

        """

        # first handle the ones that have Q defined
        forward_rates = [r for r in self.rates if r.Q >= 0.0]

        return forward_rates

    def get_reverse_rates(self):
        """Return a list of the reverse (endothermic) rates).  Note
        these may not be the same as the reverse rates identified by
        ReacLib.

        Returns
        -------
        list(Rate)

        """

        # first handle the ones that have Q defined
        reverse_rates = [r for r in self.rates if r.Q < 0.0]

        return reverse_rates

    def find_reverse(self, forward_rate, reverse_rates=None):
        """Given a forward rate, locate the rate that is its reverse.

        Returns
        -------
        Rate
        """

        if reverse_rates is None:
            reverse_rates = self.get_reverse_rates()

        reverse = None

        for rr in reverse_rates:
            if sorted(forward_rate.reactants, key=lambda x: x.A) == sorted(rr.products, key=lambda x: x.A) and \
               sorted(forward_rate.products, key=lambda x: x.A) == sorted(rr.reactants, key=lambda x: x.A):
                reverse = rr
                break

        return reverse

    def get_rate_pairs(self):
        """Find pairs of forward (Q > 0) and reverse (Q < 0) rates for the
        same link between nuclei.

        Return
        ------
        list(RatePair)

        """

        rate_pairs = []

        reverse_rates = self.get_reverse_rates()

        # loop over all the forward rates and find the matching reverse rate
        # if it exists
        for fr in self.get_forward_rates():
            rp = RatePair(forward=fr)

            rr = self.find_reverse(fr, reverse_rates=reverse_rates)

            # since we found a match, remove the reverse rate we paired
            # from out list so no other forward rate can match with it
            if rr is not None:
                rp.reverse = rr
                reverse_rates.remove(rp.reverse)

            rate_pairs.append(rp)

        # we might have some reverse rates remaining for which there
        # were no forward rates -- add those now
        for rr in reverse_rates:
            rp = RatePair(reverse=rr)
            rate_pairs.append(rp)

        return rate_pairs

    def get_nuclei(self):
        """Get all the nuclei that are part of the network.

        Returns
        -------
        list(Nucleus)

        """
        return self.unique_nuclei

    def linking_nuclei(self, nuclei, return_type=None, **kwargs):
        """Return a new network containing only rates linking the
        given nuclei.

        Parameters
        ----------
        nuclei : list, tuple
            An iterable of Nucleus objects or string names of nuclei.
        return_type : Callable
            A different constructor (e.g., a superclass constructor)
            to use if the current class does not take a `libraries`
            keyword.
        kwargs : dict
            Additional arguments to pass onto the library linking_nuclei
            method.  See :py:mod:`pynucastro.rates.library.Library.linking_nuclei`

        Returns
        -------
        RateCollection

        """

        if return_type is None:
            return_type = self.__class__
        lib = Library(rates=self.rates)
        return return_type(libraries=lib.linking_nuclei(nuclei, **kwargs))

    def get_rates(self):
        """Get a list of the reaction rates in this network.

        Returns
        -------
        list(Rate)

        """
        return self.rates

    def get_hidden_rates(self):
        """Get a list of all of the rates approximated out of the
        network

        Returns
        -------
        list(Rate)

        """
        hidden_rates = []
        for r in self.get_rates():
            if isinstance(r, ApproximateRate):
                for c in r.get_child_rates():
                    if c.removed:
                        hidden_rates.append(c)
            elif isinstance(r, ModifiedRate):
                if r.original_rate.removed:
                    hidden_rates.append(r.original_rate)
        return set(hidden_rates)

    def get_rate(self, rid):
        """Return a rate matching the id provided.

        Parameters
        ----------
        rid : str
            The id of the rate, as returned by Rate.fname

        Returns
        -------
        Rate

        """
        try:
            rid_mod = capitalize_rid(rid, "_")
            return [r for r in self.rates if r.fname == rid_mod][0]
        except IndexError:
            raise LookupError(f"rate identifier {rid!r} does not match a rate in this network.") from None

    def get_rate_by_nuclei(self, reactants, products):
        """Given a list of reactants and products, return any matching rates

        Parameters
        ----------
        reactants : list(Nucleus), list(str)
            the reactants for the reaction.  These can either be string
            names or :py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`
            objects.
        products : list(Nucleus), list(str)
            the products for the reaction.  These can either be string
            names or :py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`
            objects.

        Returns
        -------
        rates : Rate, list(Rate)
            any matching rates

        """
        reactants = sorted(Nucleus.cast_list(reactants))
        products = sorted(Nucleus.cast_list(products))
        _tmp = [r for r in self.rates if
                sorted(r.reactants) == reactants and
                sorted(r.products) == products]

        if not _tmp:
            return None
        if len(_tmp) == 1:
            return _tmp[0]
        return _tmp

    def get_rate_by_name(self, name):
        """Given a rate in the form 'A(x,y)B' return the rate

        Parameters
        ----------
        name : str
            the name of the rate, in the form "A(x,y)B"

        Returns
        -------
        Rate

        """

        reactants, products = _rate_name_to_nuc(name)
        _r = self.get_rate_by_nuclei(reactants, products)
        if _r is None:
            return None
        return _r

    def get_nuclei_needing_partition_functions(self):
        """Return a list of nuclei that require partition functions
        for one or more :py:class:`DerivedRate
        <pynucastro.rates.derived_rate.DerivedRate>` in the collection

        Returns
        -------
        list(Nucleus)

        """

        nuclei_pfs = set()
        for r in self.all_rates:
            if isinstance(r, DerivedRate) and r.use_pf:
                for nuc in r.reactants + r.products:
                    if nuc.partition_function is not None:
                        nuclei_pfs.add(nuc)
        return sorted(nuclei_pfs)

    def dedupe_partition_function_temperatures(self):
        """Return a list of unique temperature arrays needed by
        partition function tables, along with a dictionary mapping
        each Nucleus to the corresponding index into that list

        Returns
        -------
        temp_arrays : list
            a list of NumPy ndarray specifying the temperature values
            for a particular partition function tabulation.
        temp_indices : dict
            a dictionary that keyed on Nucleus that maps a nucleus to
            the index in temp_arrays containing the temperature array
            for its partition function data.
        """

        nuclei = self.get_nuclei_needing_partition_functions()
        temp_arrays = []
        temp_indices = {}
        # nuclei must be sorted, so the output is deterministic
        for nuc in nuclei:
            nuc_temp = nuc.partition_function.temperature
            # do a sequential search on temp_arrays, since it should be short
            for i, temp in enumerate(temp_arrays):
                # np.array_equal handles comparing arrays of different shapes
                if np.array_equal(nuc_temp, temp):
                    temp_indices[nuc] = i
                    break
            else:
                # no match found, add a new entry
                temp_indices[nuc] = len(temp_arrays)
                temp_arrays.append(nuc_temp)

        return temp_arrays, temp_indices

    def remove_nuclei(self, nuc_list):
        """Remove the nuclei in nuc_list from the network along with
        any rates that directly involve them (this doesn't affect
        approximate rates that may have these nuclei as hidden
        intermediate links)

        """

        nuc_list = Nucleus.cast_list(nuc_list)
        rates_to_delete = []
        for nuc in nuc_list:
            for rate in self.rates:
                if nuc in rate.reactants + rate.products:
                    if self.verbose:
                        print(f"looking to remove {rate}")
                    rates_to_delete.append(rate)

        for rate in set(rates_to_delete):
            self.rates.remove(rate)

        self._build_collection()

    def remove_rates(self, rates):
        """Remove the Rate objects in rates from the network.  Note,
        if rate list is a dict, then the keys are assumed to be the
        rates to remove

        """

        if isinstance(rates, Rate):
            self.rates.remove(rates)
        else:
            for r in rates:
                self.rates.remove(r)

        self._build_collection()

    def add_rates(self, rates):
        """Add new rates to the network.  If the rate already exists,
        it will not be added.  The network is then regenerated using
        the updated rates

        Parameters
        ----------
        rates : Rate, list(Rate)
             a single Rate object or a list of Rate objects specifying the
             rates to be added to the network.

        """

        if isinstance(rates, Rate):
            if rates not in self.rates:
                self.rates.append(rates)

        else:
            for r in rates:
                if r not in self.rates:
                    self.rates.append(r)

        self._build_collection()

    def make_ap_pg_approx(self, intermediate_nuclei=None):
        """Combine the rates A(a,g)B and A(a,p)X(p,g)B (and the
        reverse) into a single effective approximate rate.  The new
        approximate rates will be added to the network and the original
        rates will be removed (although they are still carried by the
        ApproximateRate object.

        Parameters
        ----------
        intermediate_nuclei : list, tuple
            an iterable of Nucleus objects or string names representing
            the intermediate nucleus we wish to approximate out.

        """

        # make sure that the intermediate_nuclei list are Nuclei objects
        intermediate_nuclei = Nucleus.cast_list(intermediate_nuclei, allow_None=True)

        # find all of the (a,g) rates
        ag_rates = []
        for r in self.rates:
            if (len(r.reactants) == 2 and Nucleus("he4") in r.reactants and
                len(r.products) == 1):
                ag_rates.append(r)

        # for each (a,g), check to see if the remaining rates are present
        approx_rates = []

        for r_ag in ag_rates:
            prim_nuc = sorted(r_ag.reactants)[-1]
            prim_prod = sorted(r_ag.products)[-1]

            inter_nuc = Nucleus.from_Z_A(prim_nuc.Z+1, prim_nuc.A+3)

            if intermediate_nuclei and inter_nuc not in intermediate_nuclei:
                continue

            # look for A(a,p)X
            if not (r_ap := self.get_rate_by_nuclei([prim_nuc, Nucleus("he4")],
                                                    [inter_nuc, Nucleus("p")])):
                continue

            # look for X(p,g)B
            if not (r_pg := self.get_rate_by_nuclei([inter_nuc, Nucleus("p")],
                                                    [prim_prod])):
                continue

            # look for reverse B(g,a)A
            if not (r_ga := self.get_rate_by_nuclei([prim_prod],
                                                    [prim_nuc, Nucleus("he4")])):
                continue

            # look for reverse B(g,p)X
            if not (r_gp := self.get_rate_by_nuclei([prim_prod],
                                                    [inter_nuc, Nucleus("p")])):
                continue

            # look for reverse X(p,a)A
            if not (r_pa := self.get_rate_by_nuclei([inter_nuc, Nucleus("p")],
                                                    [Nucleus("he4"), prim_nuc])):
                continue

            # build the approximate rates
            ar = ApproximateRate(r_ag, [r_ap, r_pg], r_ga, [r_gp, r_pa], approx_type="ap_pg")
            ar_reverse = ApproximateRate(r_ag, [r_ap, r_pg], r_ga, [r_gp, r_pa], is_reverse=True, approx_type="ap_pg")

            if self.verbose:
                print(f"using approximate rate {ar}")
                print(f"using approximate rate {ar_reverse}")

            # approximate rates
            approx_rates += [ar, ar_reverse]

        # remove the old rates from the rate list and add the approximate rate
        for ar in approx_rates:
            for r in ar.get_child_rates():
                try:
                    self.rates.remove(r)
                    if self.verbose:
                        print(f"removing rate {r}")
                except ValueError:
                    pass

            # add the approximate rates
            self.rates.append(ar)

        # regenerate the links
        self._build_collection()

    def make_nn_g_approx(self, intermediate_nuclei=None):
        """Combine the rates A(n,g)X(n,g)B into a single effective
        rate. The new approximate rates will be added to the network
        and the original rates will be removed (although they are
        still carried by the ApproximateRate object.

        Parameters
        ----------
        intermediate_nuclei : list, tuple
            an iterable of `Nucleus <pynucastro.nucdata.nucleus.Nucleus>`
            or string names representing the intermediate nucleus we
            wish to approximate out.

        """

        # make sure that the intermediate_nuclei list are Nuclei objects
        intermediate_nuclei = Nucleus.cast_list(intermediate_nuclei, allow_None=True)

        # if we didn't pass in a list of nuclei, consider all as targets
        # for approximation
        if not intermediate_nuclei:
            intermediate_nuclei = self.unique_nuclei

        # for each intermediate nuclei X, look to see if we have A(n,g)X and X(n,g)B
        approx_rates = []
        nuclei_approximated_out = []

        for inter_nuc in intermediate_nuclei:

            if inter_nuc.A < 2:
                # can't approximate out protons or neutrons
                continue

            nuc_A = inter_nuc - Nucleus("n")
            nuc_B = inter_nuc + Nucleus("n")

            if nuc_A in nuclei_approximated_out:
                # don't try to approximate a rate sequence starting with
                # a nucleus that we already approximated out
                continue

            # look for A(n,g)X
            if not (rf1 := self.get_rate_by_nuclei([nuc_A, Nucleus("n")],
                                                   [inter_nuc])):
                continue

            # look for X(n,g)B
            if not (rf2 := self.get_rate_by_nuclei([inter_nuc, Nucleus("n")],
                                                   [nuc_B])):
                continue

            # look for reverse B(g,n)X
            if not (rr1 := self.get_rate_by_nuclei([nuc_B],
                                                   [inter_nuc, Nucleus("n")])):
                continue

            # look for reverse X(g,n)A
            if not (rr2 := self.get_rate_by_nuclei([inter_nuc],
                                                   [nuc_A, Nucleus("n")])):
                continue

            # build the approximate rates
            ar = ApproximateRate(None, [rf1, rf2],
                                 None, [rr1, rr2],
                                 approx_type="nn_g",
                                 use_identical_particle_factor=False)

            ar_reverse = ApproximateRate(None, [rf1, rf2],
                                         None, [rr1, rr2],
                                         is_reverse=True, approx_type="nn_g",
                                         use_identical_particle_factor=False)

            nuclei_approximated_out.append(inter_nuc)
            if self.verbose:
                print(f"approximating out {inter_nuc}")
                print(f"using approximate rate {ar}")
                print(f"using approximate rate {ar_reverse}")

            # approximate rates
            approx_rates += [ar, ar_reverse]

        # remove the old rates from the rate list and add the approximate rate
        for ar in approx_rates:
            for r in ar.get_child_rates():
                try:
                    self.rates.remove(r)
                    if self.verbose:
                        print(f"removing rate {r}")
                except ValueError:
                    pass

            # add the approximate rates
            self.rates.append(ar)

        # regenerate the links
        self._build_collection()

    def make_nse_protons(self, A):
        """Replace protons in rates involving nuclei with mass number
        >= A with NSE protons.  This will decouple these rates from
        the proton captures at lower mass number, simplifying the
        linear algebra.

        Parameters
        ----------
        A : int
            mass number above which to swap regular protons for
            NSE protons.

        """

        # we want to update both the forward and reverse rates,
        # so we are consistent
        for rp in self.get_rate_pairs():

            update = False
            if rp.forward is not None:
                heavy = [n for n in rp.forward.reactants + rp.forward.products
                         if n not in [Nucleus("p"), Nucleus("n"), Nucleus("he4")]]
                if heavy:
                    if (min(heavy, key=lambda x: x.A).A >= A and
                        Nucleus("p") in rp.forward.reactants + rp.forward.products):
                        update = True
            elif rp.reverse is not None:
                heavy = [n for n in rp.reverse.reactants + rp.reverse.products
                         if n not in [Nucleus("p"), Nucleus("n"), Nucleus("he4")]]
                if heavy:
                    if (min(heavy, key=lambda x: x.A).A >= A and
                        Nucleus("p") in rp.reverse.reactants + rp.reverse.products):
                        update = True

            if update:
                if rp.forward is not None:
                    if self.verbose:
                        print(f"modifying {rp.forward.fname} to use NSE protons")
                    rp.forward.swap_protons()
                if rp.reverse is not None:
                    if self.verbose:
                        print(f"modifying {rp.reverse.fname} to use NSE protons")
                    rp.reverse.swap_protons()

        self._build_collection()

    def summary(self):
        """Print a summary of the nuclei and rates for this network"""

        print("Network summary")
        print("---------------")
        print(f"  explicitly carried nuclei: {len(self.unique_nuclei)}")
        print(f"  approximated-out nuclei: {len(self.approx_nuclei)}")
        if self.inert_nuclei:
            print(f"  inert nuclei (included in carried): {len(self.inert_nuclei)}")
        else:
            print("  inert nuclei (included in carried): 0")

        print("")

        print(f"  total number of rates: {len(self.all_rates)}")
        print("")

        print(f"  rates explicitly connecting nuclei: {len(self.rates)}")
        print(f"  hidden rates: {len(self.get_hidden_rates())}")
        print("")

        print(f"  reaclib rates: {len(self.reaclib_rates)}")
        print(f"  tabular rates: {len(self.tabular_rates)}")
        print(f"  approximate rates: {len(self.approx_rates)}")
        print(f"  derived rates: {len(self.derived_rates)}")
        print(f"  modified rates: {len(self.modified_rates)}")
        print(f"  custom rates: {len(self.custom_rates)}")

    def evaluate_rates(self, rho, T, composition, screen_func=None):
        """Evaluate the rates for a specific density, temperature, and
        composition, with optional screening.  Note: this returns that
        rate as dY/dt, where Y is the molar fraction.  For a 2 body
        reaction, a + b, this will be of the form:

        ρ Y_a Y_b N_A <σv> / (1 + δ_{ab})

        where δ is the Kronecker delta that accounts for a = b.

        If you want dn/dt, where n is the number density (so you get
        n_a n_b <σv>), then you need to multiply the results here
        by ρ N_A (where N_A is Avogadro's number).

        Parameters
        ----------
        rho : float
            density used to evaluate rates
        T : float
            temperature used to evaluate rates
        composition : Composition
            composition used to evaluate rates
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.

        Returns
        -------
        dict(Rate)

        """

        rvals = {}
        ys = composition.get_molar()
        y_e = composition.ye

        if screen_func is not None:
            screen_factors = self.evaluate_screening(rho, T, composition, screen_func)
        else:
            screen_factors = {}

        for r in self.rates:
            val = r.prefactor * rho**r.dens_exp * r.eval(T, rho=rho, comp=composition)
            if (r.weak_type == 'electron_capture' and not isinstance(r, TabularRate)):
                val = val * y_e
            yfac = functools.reduce(mul, [ys[q] for q in r.reactants])
            rvals[r] = yfac * val * screen_factors.get(r, 1.0)

        return rvals

    def jacobian_mask(self):
        """Create an array of 0 and 1 indicating the sparsity pattern
        of the Jacobian.

        Returns
        -------
        numpy.ndarray

        """

        nnuc = len(self.unique_nuclei)
        jac = np.zeros((nnuc, nnuc), dtype=np.int8)

        for i, n_i in enumerate(self.unique_nuclei):
            for j, n_j in enumerate(self.unique_nuclei):

                jac[i, j] = 0

                # we are considering dYdot(n_i) / dY(n_j)
                for r in self.nuclei_consumed[n_i] + self.nuclei_produced[n_i]:
                    if n_j in r.reactants + r.products:
                        jac[i, j] = 1
                        break

        return jac

    def evaluate_jacobian(self, rho, T, comp, *,
                          screen_func=None, exclude_rates=None):
        """Return an array of the form J_ij = dYdot_i/dY_j for the
        network

        Parameters
        ----------
        rho : float
            density used to evaluate Jacobian terms
        T : float
            temperature used to evaluate Jacobian terms
        comp : Composition
            composition used to evaluate Jacobian terms
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.
        exclude_rates : Iterable(Rate)
            a list of rates to omit from the construction of the Jacobian.

        Returns
        -------
        numpy.ndarray

        """

        # the rate.eval_jacobian_term does not compute the screening,
        # so we multiply by the factors afterwards
        if screen_func is not None:
            screen_factors = self.evaluate_screening(rho, T, comp, screen_func)
        else:
            screen_factors = {}

        if exclude_rates is None:
            exclude_rates = []

        nnuc = len(self.unique_nuclei)
        jac = np.zeros((nnuc, nnuc), dtype=np.float64)

        for i, n_i in enumerate(self.unique_nuclei):
            for j, n_j in enumerate(self.unique_nuclei):

                # we are considering dYdot(n_i) / dY(n_j)

                jac[i, j] = 0.0

                for r in self.nuclei_consumed[n_i]:
                    if r in exclude_rates:
                        continue

                    # how many of n_i are destroyed by this reaction
                    c = r.reactant_count(n_i)
                    jac[i, j] -= c * screen_factors.get(r, 1.0) *\
                        r.eval_jacobian_term(T, rho, comp, n_j)

                for r in self.nuclei_produced[n_i]:
                    if r in exclude_rates:
                        continue

                    # how many of n_i are produced by this reaction
                    c = r.product_count(n_i)
                    jac[i, j] += c * screen_factors.get(r, 1.0) *\
                        r.eval_jacobian_term(T, rho, comp, n_j)

        return jac

    def spectral_radius(self, rho, T, comp, *,
                        screen_func=None, exclude_rates=None):
        """Compute the spectral radius of the Jacobian---this is the
        max{abs(e_i)}, where e_i are the eigenvalues of the Jacobian.

        Parameters
        ----------
        rho : float
            density used to evaluate Jacobian terms
        T : float
            temperature used to evaluate Jacobian terms
        comp : Composition
            composition used to evaluate Jacobian terms
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.
        exclude_rates : Iterable(Rate)
            a list of rates to omit from the calculation.  This is useful
            for testing how the spectral radius / stiffness is affected by
            the different rates.

        Returns
        -------
        float

        """

        J = self.evaluate_jacobian(rho, T, comp,
                                   screen_func=screen_func,
                                   exclude_rates=exclude_rates)
        e = eigvals(J)
        return np.max(np.abs(e))

    def validate(self, other_library, *, forward_only=True):
        """Perform various checks on the library, comparing to
        ``other_library``, to ensure that we are not missing important
        rates.  The idea is that the current library should be a
        reduced library (perhaps the result of filtering) and then we
        want to compare to the larger ``other_library`` to see if we
        missed something important.

        Parameters
        ----------
        other_library : Library
            the library to compare to
        forward_only : bool
            do we only check the forward rates?

        Returns
        -------
        bool

        """

        current_rates = sorted(self.get_rates())

        # check the forward rates to see if any of the products are
        # not consumed by other forward rates

        passed_validation = True

        for rate in current_rates:
            if rate.derived_from_inverse:
                continue
            for p in rate.products:
                found = False
                for orate in current_rates:
                    if orate == rate:
                        continue
                    if orate.derived_from_inverse:
                        continue
                    if p in orate.reactants:
                        found = True
                        break
                if not found:
                    passed_validation = False
                    msg = f"validation: {p} produced in {rate} never consumed."
                    print(msg)

        # now check if we are missing any rates from other_library with the exact same reactants

        other_by_reactants = collections.defaultdict(list)
        for rate in sorted(other_library.get_rates()):
            other_by_reactants[tuple(sorted(rate.reactants))].append(rate)

        for rate in current_rates:
            if forward_only and rate.derived_from_inverse:
                continue

            key = tuple(sorted(rate.reactants))
            for other_rate in other_by_reactants[key]:
                # check to see if other_rate is already in current_rates
                found = True
                if other_rate not in current_rates:
                    found = False

                if not found:
                    msg = f"validation: missing {other_rate} as alternative to {rate} (Q = {other_rate.Q} MeV)."
                    print(msg)

        return passed_validation

    def find_duplicate_links(self):
        """Check the network to see if there are multiple rates that
        share the same reactants and products.  These may not be the
        same Rate object (e.g., one could be tabular the other a
        simple decay), but they will present themselves in the network
        as the same link.

        We return a list, where each entry is a list of all the rates
        that share the same link.

        Returns
        -------
        list

        """

        duplicates = find_duplicate_rates(self.get_rates())

        # there are some allowed duplicates for special cases.  We
        # will now check for those
        dupe_to_remove = []
        for dupe in duplicates:
            if is_allowed_dupe(dupe):
                dupe_to_remove.append(dupe)

        for dupe in dupe_to_remove:
            duplicates.remove(dupe)

        return duplicates

    def find_unimportant_rates(self, states, cutoff_ratio, screen_func=None):
        """Evaluate the rates at multiple thermodynamic states, and
        find the rates that are always less than `cutoff_ratio` times
        the fastest rate for each state.  This returns a dict keyed by
        Rate giving the ratio of the rate to the largest rate.

        Parameters
        ----------
        states : list, tuple
             A tuple of the form (density, temperature, composition),
             where composition is a Composition object
        cutoff_ratio : float
             The ratio of a rate to the fastest rate, below which we
             consider this rate to be unimportant.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.

        Return
        ------
        dict(Rate)

        """
        largest_ratio = {r: 0 for r in self.rates}
        for rho, T, comp in states:
            rvals = self.evaluate_rates(rho, T, comp, screen_func)
            fastest = max(rvals.values())
            for r, value in rvals.items():
                largest_ratio[r] = max(largest_ratio[r], value / fastest)
        return {r: ratio for r, ratio in largest_ratio.items() if ratio < cutoff_ratio}

    def evaluate_screening(self, rho, T, composition, screen_func):
        """Evaluate the screening factors for each rate.

        Parameters
        ----------
        rho : float
            density used to evaluate screening
        T : float
            temperature used to evaluate screening
        composition : Composition
            composition used to evaluate screening
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`

        Returns
        -------
        dict(Rate)

        """
        # this follows the same logic as BaseCxxNetwork._compute_screening_factors()
        factors = {}
        ys = composition.get_molar()
        plasma_state = make_plasma_state(T, rho, ys)
        if not self.do_screening:
            screening_map = []
        else:
            screening_map = get_screening_map(self.get_rates(),
                                              symmetric_screening=self.symmetric_screening)

        for i, scr in enumerate(screening_map):
            if not (scr.n1.dummy or scr.n2.dummy):
                scn_fac = make_screen_factors(scr.n1, scr.n2)
                scor = screen_func(plasma_state, scn_fac)
            if scr.name == "He4_He4_He4":
                # we don't need to do anything here, but we want to avoid
                # immediately applying the screening
                pass
            elif scr.name == "He4_He4_He4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "He4_He4_He4"
                # handle the second part of the screening for 3-alpha
                scn_fac2 = make_screen_factors(scr.n1, scr.n2)
                scor2 = screen_func(plasma_state, scn_fac2)

                # there might be both the forward and reverse 3-alpha
                # if we are doing symmetric screening
                for r in scr.rates:
                    # use scor from the previous loop iteration
                    # pylint: disable-next=possibly-used-before-assignment
                    factors[r] = scor * scor2
            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now
                for r in scr.rates:
                    factors[r] = scor

        return factors

    def evaluate_ydots(self, rho, T, composition,
                       screen_func=None, rate_filter=None):
        """Evaluate net rate of change of molar abundance for each
        nucleus for a specific density, temperature, and composition

        Parameters
        ----------
        rho : float
            density used to evaluate rates
        T : float
            temperature used to evaluate rates
        composition : Composition
            composition used to evaluate rates
        screen_func : Callable
            a function from :py:mod:`pynucastro.screening` used to compute the
            screening enhancement for the rates.
        rate_filter : Callable
            a function that takes a :py:class:`Rate <pynucastro.rates.rate.Rate>`
            and returns `True` or `False` if it is to be evaluated.

        Returns
        -------
        dict(Nucleus)

        """

        rvals = self.evaluate_rates(rho, T, composition, screen_func)
        ydots = {}

        for nuc in self.unique_nuclei:

            # Rates that consume / produce nuc
            if rate_filter is None:
                consuming_rates = self.nuclei_consumed[nuc]
                producing_rates = self.nuclei_produced[nuc]
            else:
                consuming_rates = [r for r in self.nuclei_consumed[nuc] if rate_filter(r)]
                producing_rates = [r for r in self.nuclei_produced[nuc] if rate_filter(r)]

            # Number of nuclei consumed / produced
            nconsumed = (r.reactant_count(nuc) for r in consuming_rates)
            nproduced = (r.product_count(nuc) for r in producing_rates)

            # Multiply each rate by the count
            consumed = (c * rvals[r] for c, r in zip(nconsumed, consuming_rates))
            produced = (c * rvals[r] for c, r in zip(nproduced, producing_rates))

            # Net change is difference between produced and consumed
            ydots[nuc] = sum(produced) - sum(consumed)

        return ydots

    def evaluate_energy_generation(self, rho, T, composition,
                                   screen_func=None, return_enu=False):
        """Evaluate the specific energy generation rate of the network for a specific
        density, temperature and composition

        Parameters
        ----------
        rho : float
            density to evaluate the rates with
        T : float
            temperature to evaluate the rates with
        composition : Composition
            composition to evaluate the rates with
        screen_func : Callable
            a function from :py:mod:`pynucastro.screening` to
            call to compute the screening factor
        return_enu : bool
            return both enuc and enu -- the energy loss
            from neutrinos from weak reactions

        Returns
        -------
        enuc : float
            the energy generation rate
        enu : float
            the neutrino loss rate from weak reactions

        """

        ydots = self.evaluate_ydots(rho, T, composition, screen_func)
        enuc = 0.

        # compute constants and units

        # ion binding energy contributions. basically e=mc^2
        for nuc in self.unique_nuclei:
            enuc += ydots[nuc] * nuc.mass * constants.MeV2erg

        # convert from molar value to erg/g/s
        enuc *= -1*constants.N_A

        # subtract neutrino losses for tabular weak reactions
        enu = 0.0
        for r in self.rates:
            if isinstance(r, TabularRate):
                # get composition
                ys = composition.get_molar()

                # need to get reactant nucleus
                nuc = r.reactants[0]
                enu += constants.N_A * ys[nuc] * r.get_nu_loss(T, rho=rho, comp=composition)

        enuc -= enu
        if return_enu:
            return enuc, enu
        return enuc

    def evaluate_activity(self, rho, T, composition, screen_func=None):
        """Compute the activity for each nucleus--the sum of
        abs(creation rate) + abs(destruction rate), i.e., this neglects the
        sign of the terms.

        Parameters
        ----------
        rho : float
            density used to evaluate rates
        T : float
            temperature used to evaluate rates
        composition : Composition
            composition used to evaluate rates
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.

        Returns
        -------
        dict(Nucleus)

        """

        rvals = self.evaluate_rates(rho, T, composition, screen_func)
        act = {}

        for nuc in self.unique_nuclei:

            # Rates that consume / produce nuc
            consuming_rates = self.nuclei_consumed[nuc]
            producing_rates = self.nuclei_produced[nuc]
            # Number of nuclei consumed / produced
            nconsumed = (r.reactant_count(nuc) for r in consuming_rates)
            nproduced = (r.product_count(nuc) for r in producing_rates)
            # Multiply each rate by the count
            consumed = (c * rvals[r] for c, r in zip(nconsumed, consuming_rates))
            produced = (c * rvals[r] for c, r in zip(nproduced, producing_rates))
            # Net activity is sum of produced and consumed
            act[nuc] = sum(produced) + sum(consumed)

        return act

    def _get_network_chart(self, rho, T, composition):
        """Create a dict, keyed by rate that holds a list of tuples
        (Nucleus, ydot)

        """

        rvals = self.evaluate_rates(rho, T, composition)

        nc = {}

        for rate, rval in rvals.items():
            nucs = []
            for n in set(rate.reactants):
                nucs.append((n, -rate.reactant_count(n) * rval))
            for n in set(rate.products):
                nucs.append((n, rate.product_count(n) * rval))
            nc[rate] = nucs

        return nc

    def network_overview(self):
        """Return a verbose network overview showing for each nucleus
        which rates consume it and which produce it.

        Returns
        -------
        str

        """

        ostr = ""
        for n in self.unique_nuclei:
            ostr += f"{n}\n"
            ostr += "  consumed by:\n"
            for r in self.nuclei_consumed[n]:
                ostr += f"     {r.string}\n"

            ostr += "  produced by:\n"
            for r in self.nuclei_produced[n]:
                ostr += f"     {r.string}\n"

            ostr += "\n"
        return ostr

    def rate_pair_overview(self):
        """Return a verbose network overview in terms of
        forward-reverse pairs

        Returns
        -------
        str

        """

        ostr = ""
        for n in self.unique_nuclei:
            ostr += f"{n}\n"
            for rp in sorted(self.nuclei_rate_pairs[n]):
                ostr += f"     {rp}\n"
        return ostr

    def get_nuclei_latex_string(self):
        """Return a string listing the nuclei in latex format

        Returns
        -------
        str

        """

        ostr = ""
        for i, n in enumerate(self.unique_nuclei):
            ostr += f"${n.pretty}$"
            if i != len(self.unique_nuclei)-1:
                ostr += ", "
        return ostr

    def get_rates_latex_table_string(self):
        """Return a string giving the rows of a LaTeX table with
        forward rates in the first column and reverse rates in the
        second column.

        Returns
        -------
        str

        """

        ostr = ""
        for rp in sorted(self.get_rate_pairs()):
            if rp.forward:
                ostr += f"{rp.forward.pretty_string:38} & \n"
            else:
                ostr += f"{' ':38} \n &"

            if rp.reverse:
                ostr += rf"  {rp.reverse.pretty_string:38} \\"
            else:
                ostr += rf"  {' ':38} \\"

            ostr += "\n"

        return ostr

    def write_network(self, *args, **kwargs):
        """Write out the network.  For :py:class:`RateCollection`
        this is a no-op.  But derived classes will use this to
        create the network in a file (or files).

        We do a find check here that the rates in the network
        are distinguishable.

        """
        assert self._distinguishable_rates(), "ERROR: Rates not uniquely identified by Rate.fname"
        self._write_network(*args, **kwargs)

    def _distinguishable_rates(self):
        """Every Rate in this RateCollection should have a unique
        Rate.fname, as the network writers distinguish the rates on
        this basis.

        """
        names = [r.fname for r in self.rates]
        for n, r in zip(names, self.rates):
            k = names.count(n)
            if k > 1:
                print(f'Found rate {r} named {n} with {k} entries in the RateCollection.')
                print(f'Rate {r} has the original source:\n{r.original_source}')
                print(f'Rate {r} is in chapter {r.chapter}')
        return len(set(names)) == len(self.rates)

    def _write_network(self, *args, **kwargs):
        """Output the network.  This version is a stub that will be
        replaced by derived classes, as this is implementation
        dependent.

        """
        # pylint: disable=unused-argument
        print('To create network integration source code, use a class that implements a specific network type.')

    def create_network_graph(self, node_nuclei, *,
                             nuclei_custom_labels=None,
                             rotated=False,
                             rate_ydots=None, ydot_cutoff_value=None,
                             consuming_rate_threshold=None,
                             show_small_ydot=False,
                             hide_xalpha=False, hide_xp=False,
                             rate_filter_function=None,
                             highlight_filter_function=None):
        """Create a graph representation of the network using
        NetworkX.  This arranges the nuclei as nodes on a grid
        determined by their N and Z, and creates the edges that
        connect the nuclei.  The various parameters control how
        which edges are present and their weights.

        Parameters
        ----------
        node_nuclei : Iterable(Nucleus)
            the nuclei to represent as nodes in the graph
        nuclei_custom_labels : dict(Nucleus, str)
            a dictionary giving alternate labels for the nodes.  If not present
            the nuclei's isotope symbol is used.
        rotated : bool
            arrange the nodes as A - 2Z vs. Z or the default Z vs. N?
        rate_ydots : dict(Rate)
            the contribution of each rate to a nuclei's dY/dt evolution.
            This can be obtained from :py:meth:`.evaluate_rates`
        ydot_cutoff_value : float
            rate threshold below which we do not add an edge connecting
            nuclei.
        consuming_rate_threshold : float
            for a nucleus that has multiple rates that consume it, remove
            any rates that are ``consuming_rate_threshold`` smaller than
            the fastest rate consuming the nucleus.
        show_small_ydot : bool
            create edges for rates below ``ydot_cutoff_value``.  They will have
            the property "real" set to -1.
        hide_xalpha : bool
            don't create edges connecting alpha particles and heavy
            nuclei in reactions of the form A(alpha,X)B or A(X,alpha)B,
            except if alpha is the heaviest product.
        hide_xp : bool
            don't create edges connecting protons and heavy
            nuclei in reactions of the form A(p,X)B or A(X,p)B.
        rate_filter_function : Callable
            a function that takes a ``Rate`` object and returns True
            or False if an edge should be created for the nuclei
            it links.
        highlight_filter_function : Callable
            a function that takes a ``Rate`` object and returns True or
            False if we want to highlight the edge in the network.  This
            sets the "highlight" property of the edge.

        Returns
        -------
        networkx.classes.multidigraph.MultiDiGraph

        """

        G = nx.MultiDiGraph()
        G.position = {}
        G.labels = {}

        if nuclei_custom_labels is None:
            nuclei_custom_labels = {}

        for n in node_nuclei:
            G.add_node(n)
            if rotated:
                G.position[n] = (n.Z, n.A - 2*n.Z)
            else:
                G.position[n] = (n.N, n.Z)
            if n in nuclei_custom_labels:
                G.labels[n] = nuclei_custom_labels[n]
            else:
                G.labels[n] = fr"${n.pretty}$"

        # Do not show rates on the graph if their corresponding ydot
        # is less than ydot_cutoff_value
        invisible_rates = set()
        if ydot_cutoff_value is not None:
            for r in self.rates:
                if rate_ydots[r] < ydot_cutoff_value:
                    invisible_rates.add(r)

        # Consider each nucleus heavier than He and all the rates that
        # consume it.  If desired, only show rates that are within a
        # threshold of the fastest rate consuming that nucleus
        if consuming_rate_threshold is not None:
            assert consuming_rate_threshold > 0.0
            for n in node_nuclei:
                if n.Z <= 2.0:
                    continue
                consump_rates = [r for r in self.rates if n in r.reactants]
                if len(consump_rates) == 0:
                    continue
                max_rate = max(rate_ydots[r] for r in consump_rates)
                for r in consump_rates:
                    if rate_ydots[r] < consuming_rate_threshold * max_rate:
                        invisible_rates.add(r)

        # edges for the rates that are explicitly in the network
        for n in node_nuclei:
            if n not in self.nuclei_consumed:
                continue
            for r in self.nuclei_consumed[n]:
                if rate_filter_function is not None:
                    if not rate_filter_function(r):
                        continue

                highlight = False
                if highlight_filter_function is not None:
                    highlight = highlight_filter_function(r)

                for p in r.products:
                    if p not in node_nuclei:
                        continue

                    if hide_xalpha and _skip_xalpha(n, p, r):
                        continue

                    if hide_xp and _skip_xp(n, p, r):
                        continue

                    # networkx doesn't seem to keep the edges in
                    # any particular order, so we associate data
                    # to the edges here directly, in this case,
                    # the reaction rate, which will be used to
                    # color it
                    # here real means that it is not an approximate rate

                    if rate_ydots is None:
                        G.add_edges_from([(n, p)], weight=0.5,
                                         real=1, highlight=highlight)
                        continue

                    try:
                        rate_weight = math.log10(rate_ydots[r])
                    except ValueError:
                        # if rate_ydots[r] is zero, then set the
                        # weight to roughly the minimum exponent
                        # possible for python floats
                        rate_weight = -308

                    if r in invisible_rates:
                        if show_small_ydot:
                            # use real -1 for displaying rates that are below ydot_cutoff
                            G.add_edges_from([(n, p)], weight=rate_weight,
                                             real=-1, highlight=highlight)

                        continue

                    G.add_edges_from([(n, p)], weight=rate_weight,
                                     real=1, highlight=highlight)

        # now consider the rates that are approximated out of the network
        rate_seen = []
        for r in self.rates:
            if not isinstance(r, ApproximateRate):
                continue
            for sr in r.hidden_rates:
                if sr in rate_seen:
                    continue
                rate_seen.append(sr)

                highlight = False
                if highlight_filter_function is not None:
                    highlight = highlight_filter_function(sr)

                for n in sr.reactants:
                    if n not in node_nuclei:
                        continue
                    for p in sr.products:
                        if p not in node_nuclei:
                            continue

                        if hide_xalpha and _skip_xalpha(n, p, sr):
                            continue

                        if hide_xp and _skip_xp(n, p, sr):
                            continue

                        G.add_edges_from([(n, p)], weight=0, real=0, highlight=highlight)

        return G

    def plot(self, rho=None, T=None, comp=None, *,
             outfile=None,
             size=(800, 600), dpi=100, title=None,
             screen_func=None,
             ydot_cutoff_value=None, show_small_ydot=False,
             consuming_rate_threshold=None,
             node_size=1000, node_font_size=12,
             node_color="#444444", node_shape="o",
             color_nodes_by_abundance=False, node_abundance_cutoff=1.e-10,
             nuclei_custom_labels=None,
             curved_edges=False,
             N_range=None, Z_range=None, rotated=False,
             always_show_p=False, always_show_alpha=False,
             hide_xp=False, hide_xalpha=False,
             edge_labels=None,
             highlight_filter_function=None,
             nucleus_filter_function=None, rate_filter_function=None,
             legend_coord=None, plot_to_cbar_ratio=20,
             grid_spec=None):
        """Make a plot of the network structure showing the links between
        nuclei.  If a full set of thermodymamic conditions are
        provided (rho, T, comp), then the links are colored by rate
        strength.

        Parameters
        ----------
        rho : float
           density to evaluate rates with
        T : float
            temperature to evaluate rates with
        comp : Composition
            composition to evaluate rates with
        outfile : str
            output name of the plot (extension determines the type)
        size : (tuple, list)
            (width, height) of the plot in pixels
        dpi : int
            dots per inch used with size to set output image size
        title : str
            title to display on the plot
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.
        ydot_cutoff_value : float
            rate threshold below which we do not show a
            line corresponding to a rate
        show_small_ydot : bool
            show visible dashed lines for rates below ``ydot_cutoff_value``
        consuming_rate_threshold : float
            for a nucleus that has multiple rates that consume it, remove
            any rates that are ``consuming_rate_threshold`` smaller than
            the fastest rate consuming the nucleus.
        node_size : float
            size of a node (in networkx units)
        node_font_size : float
            size of the font used to write the isotope in the node
        node_color : str, Callable
            color to make the nodes. May be a callable that takes a Nucleus
            object and returns a color.
        node_shape : str
            shape of the node (using matplotlib marker names)
        nuclei_custom_labels : dict(Nucleus, str)
            a dict of the form {Nucleus: str} that provides alternate
            labels for nodes (instead of using the `pretty` attribute
            of the Nucleus.
        color_nodes_by_abundance : bool
            if true, the color of the nodes is set via
            log(X) for each nucleus.  Note: this cannot be
            used with a callable for ``node_color``.
        node_abundance_cutoff : float
            the lower cutoff value used for coloring nodes by abundance.
        curved_edges : bool
            do we use arcs to connect the nodes?
        N_range : Iterable
            range of neutron number to zoom in on
        Z_range : Iterable
            range of proton number to zoom in on
        rotate : bool
            plot A - 2Z vs. Z instead of the default Z vs. N
        always_show_p : bool
            include p as a node on the plot even if we
            don't have p+p reactions
        always_show_alpha : bool
            include He4 as a node on the plot even if
            we don't have 3-alpha
        hide_xalpha : bool
            don't connect the links to alpha for heavy
            nuclei reactions of the form A(alpha,X)B or A(X,alpha)B,
            except if alpha is the heaviest product.
        hide_xp : bool
            don't connect the links to p for heavy
            nuclei reactions of the form A(p,X)B or A(X,p)B.
        edge_labels : dict
            a dictionary of the form {(n1, n2): "label"}
            that gives labels for the edges in the network connecting
            nucleus n1 to n2.
        highlight_filter_function : Callable
            a function that takes a ``Rate`` object and returns True or
            False if we want to highlight the rate edge.
        nucleus_filter_function : Callable
            a function that takes a ``Nucleus`` object and returns
            True or False if it is to be shown as a node.
        rate_filter_function : Callable
            a function that takes a ``Rate`` object
            and returns True or False if it is to be shown as an edge.
        plot_to_cbar_ratio : float
            ratio of main axes to colorbar size
        grid_spec : matplotlib.gridspec.GridSpec
            a 2x2 matplotlib GridSpec to use in arranging the plot.
            If the colorbar is on the right, only the columns will be
            used.  If the colorbar is on the left, only the rows will
            be used.  This is only needed if you want to override the
            GridSpec created internally.

        Returns
        -------
        matplotlib.figure.Figure

        """

        if grid_spec is not None:
            fig = grid_spec.figure
        else:
            fig = plt.figure(constrained_layout=True,
                             figsize=(size[0]/dpi, size[1]/dpi))

        # we'll use a grid spec of 2 x 2.  We can merge columns / rows
        # as needed to give us the flexibility to have colorbars
        if rotated:
            if grid_spec is None:
                gs = mpl.gridspec.GridSpec(nrows=2, ncols=2,
                                           height_ratios=[plot_to_cbar_ratio, 1], figure=fig)
            else:
                gs = grid_spec
            # plot is the top row, colorbar(s) will be the bottom
            ax = fig.add_subplot(gs[0, :])
        else:
            if grid_spec is None:
                gs = mpl.gridspec.GridSpec(nrows=2, ncols=2,
                                           width_ratios=[plot_to_cbar_ratio, 1], figure=fig)
            else:
                gs = grid_spec
            # plot is the left column, colorbar(s) will be on the right
            ax = fig.add_subplot(gs[:, 0])

        # in general, we do not show p, n, alpha,
        # unless we have p + p, 3-a, etc.
        hidden_nuclei = ["n"]
        if not always_show_p:
            hidden_nuclei.append("p")
            hidden_nuclei.append("p_nse")
        if not always_show_alpha:
            hidden_nuclei.append("he4")

        # nodes -- the node nuclei will be all of the heavies
        # add all the nuclei into G.node
        node_nuclei = []
        colors = []

        if callable(node_color) and color_nodes_by_abundance:
            raise NotImplementedError("setting node_color to a callable and using color_nodes_by_abundance together is not supported.")

        if callable(node_color):
            get_node_color = node_color
        elif color_nodes_by_abundance:
            nuc_norm = mpl.colors.LogNorm(vmin=node_abundance_cutoff, vmax=1.0)
            nuc_sm = mpl.cm.ScalarMappable(norm=nuc_norm, cmap=mpl.colormaps.get_cmap("magma"))

            def get_node_color(_nuc):
                return nuc_sm.to_rgba(comp[_nuc])
        else:
            def get_node_color(_nuc):
                return node_color

        for n in self.unique_nuclei:
            if n.raw not in hidden_nuclei:
                node_nuclei.append(n)
                colors.append(get_node_color(n))
            else:
                # show hidden nuclei only if they react with themselves
                for r in self.rates:
                    if not isinstance(r, (ApproximateRate, ModifiedRate)) and r.reactant_count(n) > 1:
                        node_nuclei.append(n)
                        colors.append(get_node_color(n))
                        break

        # approx nuclei are given a different color
        for n in self.approx_nuclei:
            node_nuclei.append(n)
            colors.append("#888888")

        if nucleus_filter_function is not None:
            node_nuclei = list(filter(nucleus_filter_function, node_nuclei))
            # redo the colors:
            colors = []
            for n in node_nuclei:
                if n in self.approx_nuclei:
                    colors.append("#888888")
                else:
                    colors.append(get_node_color(n))

        # get the rates for each reaction
        if rho is not None and T is not None and comp is not None:
            rate_ydots = self.evaluate_rates(rho, T, comp,
                                             screen_func=screen_func)
        else:
            rate_ydots = None

        G = self.create_network_graph(node_nuclei,
                                      rate_ydots=rate_ydots,
                                      ydot_cutoff_value=ydot_cutoff_value,
                                      hide_xalpha=hide_xalpha, hide_xp=hide_xp,
                                      show_small_ydot=show_small_ydot,
                                      consuming_rate_threshold=consuming_rate_threshold,
                                      rate_filter_function=rate_filter_function,
                                      highlight_filter_function=highlight_filter_function,
                                      rotated=rotated,
                                      nuclei_custom_labels=nuclei_custom_labels)

        # It seems that networkx broke backwards compatibility, and 'zorder' is no longer a valid
        # keyword argument. The 'linewidth' argument has also changed to 'linewidths'.

        nx.draw_networkx_nodes(G, G.position,      # plot the element at the correct position
                               node_color=colors, alpha=1.0,
                               node_shape=node_shape, node_size=node_size, linewidths=2.0, ax=ax)

        if color_nodes_by_abundance:
            node_font_color = {}
            for n in node_nuclei:
                try:
                    r, g, b = get_node_color(n)[:3]
                    # simple rgb -> luminance conversion
                    # see: https://en.wikipedia.org/wiki/Luma_(video)
                    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
                    node_font_color[n] = "black" if luminance > 0.5 else "white"
                except KeyError:
                    # hidden nucleus
                    node_font_color[n] = "white"
        else:
            node_font_color = "w"

        nx.draw_networkx_labels(G, G.position, G.labels,   # label the name of element at the correct position
                                font_size=node_font_size, font_color=node_font_color, ax=ax)

        # now we'll draw edges in several groups

        if curved_edges:
            connectionstyle = "arc3, rad = 0.2"
            sort_reverse = False
        else:
            connectionstyle = "arc3"
            sort_reverse = True

        sorted_edges = sorted(G.edges(data=True), key=lambda edge: edge[-1].get("weight", 0),
                              reverse=sort_reverse)
        real_edges = [(u, v) for u, v, e in sorted_edges if e["real"] == 1]
        real_weights = [e["weight"] for u, v, e in sorted_edges if e["real"] == 1]

        if rate_ydots is None:
            edge_color = "C0"
        else:
            edge_color = real_weights
        ww = np.array(real_weights)
        min_weight = ww.min()
        max_weight = ww.max()
        dw = (max_weight - min_weight)/4
        widths = np.ones_like(ww)
        if dw > 0:
            widths[ww > min_weight + dw] = 1.5
            widths[ww > min_weight + 2*dw] = 2.5
            widths[ww > min_weight + 3*dw] = 4
        else:
            widths *= 2

        # draw the approximate rate edges
        approx_edges = [(u, v) for u, v, e in G.edges(data=True) if e["real"] == 0]

        _ = nx.draw_networkx_edges(G, G.position, width=1,
                                   edgelist=approx_edges, edge_color="0.5",
                                   connectionstyle=connectionstyle,
                                   style="dashed", node_size=node_size, ax=ax)

        # draw the edges for the rates that are below ydot_cutoff_value
        invis_edges = [(u, v) for u, v, e in G.edges(data=True) if e["real"] == -1]

        _ = nx.draw_networkx_edges(G, G.position, width=1,
                                   edgelist=invis_edges, edge_color="gray",
                                   connectionstyle=connectionstyle,
                                   style="dotted", node_size=node_size, ax=ax)

        # draw the edges that are real rates, and above any cutoffs
        real_edges_lc = nx.draw_networkx_edges(G, G.position, width=list(widths),
                                               edgelist=real_edges, edge_color=edge_color,
                                               connectionstyle=connectionstyle,
                                               node_size=node_size,
                                               edge_cmap=plt.cm.viridis, ax=ax)

        # highlight edges
        highlight_edges = [(u, v) for u, v, e in G.edges(data=True) if e["highlight"]]

        if rho is None:
            # we are not coloring edges by reaction rate, so highlight in yellow
            highlight_color = "yellow"
        else:
            # use C0 since it doesn't blend in with viridis
            highlight_color = "C0"

        _ = nx.draw_networkx_edges(G, G.position, width=5,
                                   edgelist=highlight_edges, edge_color=highlight_color, alpha=0.5,
                                   connectionstyle=connectionstyle,
                                   node_size=node_size, ax=ax)

        if edge_labels:
            nx.draw_networkx_edge_labels(G, G.position,
                                         connectionstyle=connectionstyle,
                                         font_size=node_font_size,
                                         edge_labels=edge_labels)

        # figure out the colorbar axes -- we have a single colorbar if we are doing the
        # rate_ydots.  We have 2 colorbars if we are also doing the color_nodes_by_abundance
        rate_cb_ax = None
        node_cb_ax = None
        if rate_ydots and color_nodes_by_abundance:
            if rotated:
                rate_cb_ax = fig.add_subplot(gs[1, 0])
                node_cb_ax = fig.add_subplot(gs[1, 1])
            else:
                rate_cb_ax = fig.add_subplot(gs[0, 1])
                node_cb_ax = fig.add_subplot(gs[1, 1])
        elif rate_ydots:
            if rotated:
                rate_cb_ax = fig.add_subplot(gs[1, :])
            else:
                rate_cb_ax = fig.add_subplot(gs[:, 1])

        orientation = "vertical"
        if rotated:
            orientation = "horizontal"

        if rate_ydots is not None:
            pc = mpl.collections.PatchCollection(real_edges_lc, cmap=plt.cm.viridis)
            pc.set_array(real_weights)
            fig.colorbar(pc, cax=rate_cb_ax, label="log10(rate)", orientation=orientation)

        if color_nodes_by_abundance:
            fig.colorbar(nuc_sm, cax=node_cb_ax, label="log10(X)",
                         orientation=orientation)

        if not rotated:
            ax.set_xlabel(r"$N$", fontsize="large")
            ax.set_ylabel(r"$Z$", fontsize="large")
        else:
            ax.set_xlabel(r"$Z$", fontsize="large")
            ax.set_ylabel(r"$A - 2Z$", fontsize="large")

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        if not rotated:
            if Z_range is not None and N_range is not None:
                ax.set_xlim(N_range[0], N_range[1])
                ax.set_ylim(Z_range[0], Z_range[1])
        else:
            if Z_range is not None:
                ax.set_xlim(Z_range[0], Z_range[1])

        # if we are rotated and all nuclei have Z = A, then make
        # the vertical axis symmetric
        if rotated:
            ZA = np.array([n.A - 2 * n.Z for n in node_nuclei])
            if ZA.min() == ZA.max():
                ax.set_ylim(ZA.min() - 0.5, ZA.min() + 0.5)

        if not rotated:
            ax.set_aspect("equal", "datalim")

        if legend_coord is not None:
            assert len(legend_coord) == 2
            eps = 0.1
            for label, dd in RATE_LINES.items():
                dZ = dd[0]
                dN = dd[1]
                if rotated:
                    ax.arrow(legend_coord[0], legend_coord[1],
                             dZ, dN-dZ, width=0.04,
                             length_includes_head=True)
                    ax.text(legend_coord[0]+dZ+eps, legend_coord[1]+dN-dZ+np.sign(dN-dZ)*eps,
                            label, fontsize="small")
                else:
                    ax.arrow(legend_coord[1], legend_coord[0],
                             dN, dZ, width=0.04,
                             length_includes_head=True)
                    ax.text(legend_coord[1]+dN+eps, legend_coord[0]+dZ+eps,
                            label, fontsize="small")

        if title is not None:
            fig.suptitle(title)

        if outfile is not None:
            fig.savefig(outfile, dpi=dpi)

        return fig

    def plot_jacobian(self, rho, T, comp, *,
                      outfile=None, screen_func=None,
                      rate_scaling=1.e10,
                      size=(800, 800), dpi=100):
        """Plot the Jacobian matrix of the system.

        Parameters
        ----------
        rho : float
            density used to evaluate terms
        T : float
            temperature used to evaluate terms
        comp : Composition
            composition used to evaluate terms
        outfile : str
            output file for plot (extension is used to specify file type)
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the evaluated rates will include the screening
            correction.
        rate_scaling : float
            the cutoff of values that we show, relative to the peak.  Any
            Jacobian element smaller than this will not be shown.
        size : (tuple, list)
            size in pixels for the output plot
        dpi : float
            dots per inch for the output plot

        Returns
        -------
        matplotlib.figure.Figure

        """

        jac = self.evaluate_jacobian(rho, T, comp, screen_func=screen_func)

        valid_max = np.abs(jac).max()

        # pylint: disable-next=redundant-keyword-arg
        norm = SymLogNorm(valid_max/rate_scaling, vmin=-valid_max, vmax=valid_max)

        fig, ax = plt.subplots()
        fig.set_size_inches(size[0]/dpi, size[1]/dpi)

        ax.set_xticks(np.arange(len(self.unique_nuclei)),
                      labels=[f"${n.pretty}$" for n in self.unique_nuclei], rotation=90)

        ax.set_yticks(np.arange(len(self.unique_nuclei)),
                      labels=[f"${n.pretty}$" for n in self.unique_nuclei])

        im = ax.imshow(jac, norm=norm, cmap=plt.cm.bwr)

        ax.set_aspect("equal")

        # Turn spines off and create white grid.
        #ax.spines[:].set_visible(False)

        ax.set_xticks(np.arange(jac.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(jac.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="w", linestyle='-', linewidth=2)
        ax.tick_params(which="minor", bottom=False, left=False)

        fig.colorbar(im, ax=ax, shrink=0.75)

        if outfile is not None:
            fig.savefig(outfile, bbox_inches="tight")

        return fig

    def plot_network_chart(self, rho=None, T=None, comp=None, *,
                           outfile=None,
                           size=(800, 800), dpi=100,
                           force_one_column=False,
                           max_ydot_ratio=1.e15,
                           plot_to_cbar_ratio=20):
        """Plot a heatmap showing which rates are affected by which
        nuclei.

        Parameters
        ----------
        rho : float
            density used to evaluate rates
        T : float
            temperature used to evaluate rates
        comp : Composition
            composition used to evaluate rates
        outfile : str
            filename to output image
        size : Iterable(int)
            image dimensions in pixels
        dpi : int
            dots per inch for physical size of image
        force_one_column : bool
            do we insist on a single column for the plot?
        max_ydot_ratio : float
            ratio between maximum ydot and minimum shown in the plot
        plot_to_cbar_ratio : float
            ratio of main axes to colorbar size

        Returns
        -------
        matplotlib.figure.Figure

        """

        nc = self._get_network_chart(rho, T, comp)

        # find the limits
        _ydot = []
        for r in self.rates:
            for _, y in nc[r]:
                _ydot.append(y)

        _ydot = np.asarray(_ydot)
        valid_max = np.abs(_ydot[_ydot != 0]).max()

        # pylint: disable-next=redundant-keyword-arg
        norm = SymLogNorm(valid_max/max_ydot_ratio, vmin=-valid_max, vmax=valid_max)

        # if there are a lot of rates, we split the network chart into
        # two side-by-side panes, with the first half of the rates on
        # the left and the second half of the rates on the right

        # how many panes?

        if len(self.rates) > 3 * len(self.unique_nuclei):
            npanes = 2
        else:
            npanes = 1

        if force_one_column:
            npanes = 1

        fig = plt.figure(constrained_layout=True,
                         figsize=(size[0]/dpi, size[1]/dpi))

        if npanes == 2:
            # we'll use a grid spec of 3x1 and make the main plot
            # areas and colorbar from it (colorbar on the right)
            gs = mpl.gridspec.GridSpec(figure=fig,
                                       nrows=1, ncols=3,
                                       width_ratios=[plot_to_cbar_ratio, plot_to_cbar_ratio, 1])
            drate = (len(self.rates) + 1) // 2
        else:
            # colorbar on the bottom
            drate = len(self.rates)
            gs = mpl.gridspec.GridSpec(figure=fig,
                                       nrows=2, ncols=1,
                                       height_ratios=[plot_to_cbar_ratio, 1])

        _rates = sorted(self.rates)

        for ipane in range(npanes):

            if npanes == 2:
                ax = fig.add_subplot(gs[0, ipane])
            else:
                ax = fig.add_subplot(gs[0, 0])

            istart = ipane * drate
            iend = min((ipane + 1) * drate - 1, len(self.rates)-1)

            nrates = iend - istart + 1

            data = np.zeros((nrates, len(self.unique_nuclei)), dtype=np.float64)

            # loop over rates -- each rate is a line in a grid of nuclei vs rate

            #ax = grid[ipane]

            for irate, r in enumerate(_rates):
                if istart <= irate <= iend:
                    irow = irate - istart
                    for n, ydot in nc[r]:
                        icol = self.unique_nuclei.index(n)
                        assert data[irow, icol] == 0.0
                        data[irow, icol] = ydot

            # each pane has all the nuclei
            ax.set_xticks(np.arange(len(self.unique_nuclei)),
                          labels=[f"${n.pretty}$"
                                  for n in self.unique_nuclei], rotation=90)

            # each pane only has its subset of rates
            ax.set_yticks(np.arange(nrates),
                          labels=[f"{r.pretty_string}"
                                  for irate, r in enumerate(_rates)
                                  if istart <= irate <= iend])

            im = ax.imshow(data, norm=norm, cmap=plt.cm.bwr)

            ax.set_aspect("equal")

            # Turn spines off and create white grid.
            ax.spines[:].set_visible(False)

            ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
            ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
            ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
            ax.tick_params(which="minor", bottom=False, left=False,
                           labelsize=8)

        if npanes == 1:
            cax = fig.add_subplot(gs[1, 0])
            fig.colorbar(im, cax=cax, orientation="horizontal")
        else:
            cax = fig.add_subplot(gs[0, 2])
            fig.colorbar(im, cax=cax, orientation="vertical")

        if outfile is not None:
            fig.savefig(outfile, bbox_inches="tight")

        return fig

    @staticmethod
    def _safelog(arr, small):

        arr = np.copy(arr)
        if np.any(arr < 0.0):
            raise ValueError("Negative values not allowed for logscale - try symlog instead.")
        zeros = arr == 0.0
        arr[zeros] = min(small, arr[~zeros].min() / 10)
        return np.log10(arr)

    @staticmethod
    def _symlog(arr, linthresh=1.0, linscale=1.0):

        symlog_transform = SymmetricalLogTransform(10, linthresh, linscale)
        arr = symlog_transform.transform_non_affine(arr)

        return arr

    @staticmethod
    def _scale(arr, minval=None, maxval=None):

        if minval is None:
            minval = arr.min()
        if maxval is None:
            maxval = arr.max()
        if minval != maxval:
            scaled = (arr - minval) / (maxval - minval)
        else:
            scaled = np.zeros_like(arr)
        scaled[scaled < 0.0] = 0.0
        scaled[scaled > 1.0] = 1.0
        return scaled

    def gridplot(self, rho=None, T=None, comp=None, color_field="X", **kwargs):
        """Plot nuclides as cells on a grid of Z vs. N, colored by `color_field`.
        If called without a composition, the function will just plot the grid
        with no color field.

        Parameters
        ----------
        rho : float
            density used to evaluate color_field
        T : float
            temperature used to evaluate color_field
        comp : Composition
            composition used to evaluate color_field
        color_field : str
            field to color by. Must be one of 'X' (mass fraction),
            'Y' (molar abundance), 'Xdot' (time derivative of X), 'Ydot' (time
            derivative of Y), or 'activity' (sum of contributions to Ydot of
            all rates, ignoring sign).
        scale : str
            One of 'linear', 'log', and 'symlog'. Linear by default.
        small : float
            If using logarithmic scaling, zeros will be replaced with
            this value. 1e-30 by default.
        linthresh : float
            Linearity threshold for symlog scaling.
        linscale : float
            The number of decades to use for each half of the linear
            range. Stretches linear range relative to the logarithmic range.
        filter_function : Callable
            A callable to filter :py:class:`Nucleus <pynucastro.nucdata.nucleus.Nucleus>`
            objects with. Should return `True` if the nuclide should be plotted.
        outfile : str
            Output file to save the plot to. The plot will be shown if
            not specified.
        dpi : float
            DPI to save the image file at.
        cmap : str
            Name of the matplotlib colormap to use. Default is 'magma'.
        edgecolor : str
            Color of grid cell edges.
        area : float
            Area of the figure without the colorbar, in square inches. 64
            by default.
        no_axes : bool
            Set to True to omit axis spines.
        no_ticks : bool
            Set to True to omit tickmarks.
        no_cbar : bool
            Set to True to omit colorbar.
        cbar_label : str
            Colorbar label.
        cbar_bounds : list, tuple
            Explicit colorbar bounds.
        cbar_format : str, matplotlib.ticker.Formatter
            Format string or formatter object for the colorbar ticks.
        cbar_ticks : int
            Number of ticks to use on the colorbar

        Returns
        -------
        matplotlib.figure.Figure
        """

        # Process kwargs
        outfile = kwargs.pop("outfile", None)
        scale = kwargs.pop("scale", "linear")
        cmap = kwargs.pop("cmap", "viridis")
        edgecolor = kwargs.pop("edgecolor", "grey")
        small = kwargs.pop("small", 1e-30)
        area = kwargs.pop("area", 64)
        no_axes = kwargs.pop("no_axes", False)
        no_ticks = kwargs.pop("no_ticks", False)
        no_cbar = kwargs.pop("no_cbar", False)
        cbar_label = kwargs.pop("cbar_label", None)
        cbar_format = kwargs.pop("cbar_format", None)
        cbar_bounds = kwargs.pop("cbar_bounds", None)
        filter_function = kwargs.pop("filter_function", None)
        dpi = kwargs.pop("dpi", 100)
        linthresh = kwargs.pop("linthresh", 1.0)
        linscale = kwargs.pop("linscale", 1.0)
        cbar_ticks = kwargs.pop("cbar_ticks", None)

        if kwargs:
            warnings.warn(f"Unrecognized keyword arguments: {kwargs.keys()}")

        # Get figure, colormap
        fig, ax = plt.subplots()
        cmap = mpl.colormaps.get_cmap(cmap)

        # Get nuclei and all 3 numbers
        nuclei = self.unique_nuclei
        if filter_function is not None:
            nuclei = list(filter(filter_function, nuclei))
        Ns = np.array([n.N for n in nuclei])
        Zs = np.array([n.Z for n in nuclei])
        As = Ns + Zs

        # Compute weights
        color_field = color_field.lower()
        if color_field not in {"x", "y", "ydot", "xdot", "activity"}:
            raise ValueError(f"Invalid color field: '{color_field}'")

        if comp is None:

            values = np.zeros(len(nuclei))

        elif color_field == "x":

            values = np.array([comp[nuc] for nuc in nuclei])

        elif color_field == "y":

            ys = comp.get_molar()
            values = np.array([ys[nuc] for nuc in nuclei])

        elif color_field in {"ydot", "xdot"}:

            if rho is None or T is None:
                raise ValueError("Need both rho and T to evaluate rates!")
            ydots = self.evaluate_ydots(rho, T, comp)
            values = np.array([ydots[nuc] for nuc in nuclei])
            if color_field == "xdot":
                values *= As

        elif color_field == "activity":

            if rho is None or T is None:
                raise ValueError("Need both rho and T to evaluate rates!")
            act = self.evaluate_activity(rho, T, comp)
            values = np.array([act[nuc] for nuc in nuclei])

        if scale == "log":
            values = self._safelog(values, small)
        elif scale == "symlog":
            values = self._symlog(values, linthresh, linscale)

        if cbar_bounds is None:
            cbar_bounds = values.min(), values.max()

        weights = self._scale(values, *cbar_bounds)

        # Plot a square for each nucleus
        for nuc, weight in zip(nuclei, weights):

            square = plt.Rectangle((nuc.N - 0.5, nuc.Z - 0.5), width=1, height=1,
                                   facecolor=cmap(weight), edgecolor=edgecolor)
            ax.add_patch(square)

        # Set limits
        maxN, minN = max(Ns), min(Ns)
        maxZ, minZ = max(Zs), min(Zs)

        plt.xlim(minN - 0.5, maxN + 0.6)
        plt.ylim(minZ - 0.5, maxZ + 0.6)

        # Set plot appearance
        rat = (maxN - minN) / (maxZ - minZ)
        width = np.sqrt(area * rat)
        height = area / width
        fig.set_size_inches(width, height)

        plt.xlabel(r"N $\rightarrow$")
        plt.ylabel(r"Z $\rightarrow$")

        if no_axes or no_ticks:

            plt.tick_params(
                axis='both',
                which='both',
                bottom=False,
                left=False,
                labelbottom=False,
                labelleft=False
            )

        else:

            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if no_axes:
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

        # Colorbar stuff
        if not no_cbar and comp is not None:

            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='3.5%', pad=0.1)

            if scale == "symlog":
                cbar_norm = mpl.colors.SymLogNorm(linthresh, linscale, *cbar_bounds)
            else:
                cbar_norm = mpl.colors.Normalize(*cbar_bounds)

            smap = mpl.cm.ScalarMappable(norm=cbar_norm, cmap=cmap)

            if not cbar_label:

                capfield = color_field.capitalize()
                if scale == "log":
                    cbar_label = f"log[{capfield}]"
                elif scale == "symlog":
                    cbar_label = f"symlog[{capfield}]"
                else:
                    cbar_label = capfield

            # set number of ticks
            if cbar_ticks is not None:
                tick_locator = mpl.ticker.MaxNLocator(nbins=cbar_ticks)
                tick_labels = tick_locator.tick_values(values.min(), values.max())

                # for some reason tick_locator doesn't give the label of the first tick
                # add them manually
                if scale == "symlog":
                    tick_labels = np.append(tick_labels, [linthresh, -linthresh])

                if cbar_format is None:
                    cbar_format = mpl.ticker.FormatStrFormatter("%.3g")

            else:
                tick_labels = None

            fig.colorbar(smap, cax=cax, orientation="vertical", ticks=tick_labels,
                         label=cbar_label, format=cbar_format)

        if outfile is not None:
            plt.tight_layout()
            plt.savefig(outfile, dpi=dpi)

        return fig

    def __repr__(self):
        string = ""
        for r in self.rates:
            string += f"{r.string}\n"
        return string


class Explorer:
    """A simple class that enables interactive exploration a
    RateCollection, presenting density and temperature sliders to
    update the reaction rate values.

    Parameters
    ----------
    rc : RateCollection
        The RateCollection we will visualize.
    comp : Composition
        A composition that will be used for evaluating the rates
    kwargs : dict
        Additional parameters that will be passed through to the
        RateCollection plot() function.  Note that "T" and "rho"
        will be ignored.

    """

    def __init__(self, rc, comp, **kwargs):
        self.rc = rc
        self.comp = comp
        self.kwargs = kwargs

        # we will override any T and rho passed in
        kwargs.pop("T", None)
        kwargs.pop("rho", None)

    def _make_plot(self, logrho, logT):
        self.rc.plot(rho=10.0**logrho, T=10.0**logT,
                     comp=self.comp, **self.kwargs)

    def explore(self, logrho=(2, 6, 0.1), logT=(7, 9, 0.1)):
        """Create the interactive visualization.  This uses ipywidgets.interact
        to create an interactive visualization.

        Parameters
        ----------
        logrho : list, tuple
            a tuple of (starting log(rho), ending log(rho), dlogrho) that
            defines the range of densities to explore with an interactive
            slider.
        logT : list, tuple
            a tuple of (starting log(T), ending log(T), dlogT) that
            defines the range of temperatures to explore with an interactive
            slider.
        """

        interact(self._make_plot, logrho=logrho, logT=logT)
