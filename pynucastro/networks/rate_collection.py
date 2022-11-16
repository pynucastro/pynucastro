"""A collection of classes and methods to deal with collections of
rates that together make up a network."""

import collections
import copy
import functools
import math
import os
import warnings
from operator import mul

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
from scipy import constants
from scipy.optimize import fsolve

# Import Rate
from pynucastro.nucdata import Nucleus, PeriodicTable
from pynucastro.rates import (ApproximateRate, DerivedRate, Library, Rate,
                              RatePair, TabularRate, load_rate)
from pynucastro.rates.library import _rate_name_to_nuc
from pynucastro.screening import make_plasma_state, make_screen_factors
from pynucastro.screening.screen import NseState

mpl.rcParams['figure.dpi'] = 100


def _skip_xalpha(n, p, r):
    """utility function to consider if we show an (a, x) or (x, a) rate.  Here, p is the
    product we want to link to"""

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
    """utility function to consider if we show an (p, x) or (x, p) rate.  Here, p is the
    product we want to link to"""

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


class Composition:
    """a composition holds the mass fractions of the nuclei in a network
    -- useful for evaluating the rates

    """
    def __init__(self, nuclei, small=1.e-16):
        """nuclei is an iterable of the nuclei (Nucleus objects) in the network"""
        if not isinstance(nuclei[0], Nucleus):
            raise ValueError("must supply an iterable of Nucleus objects")
        self.X = {k: small for k in nuclei}

    def set_solar_like(self, Z=0.02):
        """ approximate a solar abundance, setting p to 0.7, He4 to 0.3 - Z and
        the remainder evenly distributed with Z """
        num = len(self.X)
        rem = Z/(num-2)
        for k in self.X:
            if k == Nucleus("p"):
                self.X[k] = 0.7
            elif k.raw == "he4":
                self.X[k] = 0.3 - Z
            else:
                self.X[k] = rem

        self.normalize()

    def set_all(self, xval):
        """ set all species to a particular value """
        for k in self.X:
            self.X[k] = xval

    def set_nuc(self, name, xval):
        """ set nuclei name to the mass fraction xval """
        for k in self.X:
            if k.raw == name:
                self.X[k] = xval
                break

    def normalize(self):
        """ normalize the mass fractions to sum to 1 """
        X_sum = sum(self.X[k] for k in self.X)

        for k in self.X:
            self.X[k] /= X_sum

    def get_molar(self):
        """ return a dictionary of molar fractions"""
        molar_frac = {k: v/k.A for k, v in self.X.items()}
        return molar_frac

    def eval_ye(self):
        """ return the electron fraction """
        zvec = []
        avec = []
        xvec = []
        for n in self.X:
            zvec.append(n.Z)
            avec.append(n.A)
            xvec.append(self.X[n])
        zvec = np.array(zvec)
        avec = np.array(avec)
        xvec = np.array(xvec)
        electron_frac = np.sum(zvec*xvec/avec)/np.sum(xvec)
        return electron_frac

    def __str__(self):
        ostr = ""
        for k in self.X:
            ostr += f"  X({k}) : {self.X[k]}\n"
        return ostr

    def plot(self, trace_threshold=0.1, hard_limit=None, size=(9, 5)):
        """ Make a pie chart of Composition. group trace nuceli together and explode into bar chart

        parameters
        ----------

        trace_threshold : the threshold to consider a component to be trace.

        hard_limit : hard limit for nuclei to be labeled in the plot. Default is None which will set the limit to 5% of total trace abundance

        size: tuple giving width x height of the plot in inches

        """

        # find trace nuclei
        trace_keys = []
        trace_tot = 0.
        main_keys = []
        for k in self.X:
            # if below threshold, count as trace element
            if self.X[k] < trace_threshold:
                trace_keys.append(k)
                trace_tot += self.X[k]
            else:
                main_keys.append(k)

        # check if any trace nuclei
        if not trace_keys:
            # just do pie chart without including trace

            fig, ax = plt.subplots(1, 1, figsize=size)

            ax.pie(self.X.values(), labels=self.X.keys(), autopct=lambda p: f"{p/100:0.3f}")

        else:
            # find trace nuclei which contribute little to trace proportion
            if hard_limit is None:
                # make hardlimit proportional to trace abundance
                hard_limit = 0.05*trace_tot

            limited_trace_keys = []
            other_trace_tot = 0.
            for k in trace_keys:
                if self.X[k] < hard_limit:
                    other_trace_tot += self.X[k]
                else:
                    limited_trace_keys.append(k)

            # make figure and assign axis objects
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=size)
            fig.subplots_adjust(wspace=0)

            # pie chart parameters
            main_values = [trace_tot] + [self.X[k] for k in main_keys]
            main_labels = ['trace'] + main_keys
            explode = [0.2] + [0. for i in range(len(main_keys))]

            # rotate so that first wedge is split by the x-axis
            angle = -180 * main_values[0]
            wedges, *_ = ax1.pie(main_values, autopct=lambda p: f"{p/100:0.3f}", startangle=angle,
                                labels=main_labels, explode=explode)

            # bar chart parameters
            trace_values = [self.X[k] for k in limited_trace_keys] + [other_trace_tot]
            trace_labels = [k.short_spec_name for k in limited_trace_keys] + ['other']
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


class ScreeningPair:
    """a pair of nuclei that will have rate screening applied.  We store a
    list of all rates that match this pair of nuclei"""

    def __init__(self, name, nuc1, nuc2, rate=None):
        self.name = name
        self.n1 = nuc1
        self.n2 = nuc2

        if rate is None:
            self.rates = []
        else:
            self.rates = [rate]

    def add_rate(self, rate):
        if rate not in self.rates:
            self.rates.append(rate)

    def __str__(self):
        ostr = f"screening for {self.n1} + {self.n2}\n"
        ostr += "rates:\n"
        for r in self.rates:
            ostr += f"  {r}\n"
        return ostr

    def __eq__(self, other):
        """all we care about is whether the names are the same -- that conveys
        what the reaction is"""

        return self.name == other.name


class RateCollection:
    """ a collection of rates that together define a network """

    pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    def __init__(self, rate_files=None, libraries=None, rates=None, precedence=(),
                 inert_nuclei=None,
                 symmetric_screening=False, do_screening=True):
        """rate_files are the files that together define the network.  This
        can be any iterable or single string.

        This can include Reaclib library files storing multiple rates.

        If libraries is supplied, initialize a RateCollection using the rates
        in the Library object(s) in list 'libraries'.

        If rates is supplied, initialize a RateCollection using the
        Rate objects in the list 'rates'.

        Precedence should be sequence of rate labels (e.g. wc17) to be used to
        resolve name conflicts. If a nonempty sequence is provided, the rate
        collection will automatically be scanned for multiple rates with the
        same name. If all of their labels were given a ranking, the rate with
        the label that comes first in the sequence will be retained and the
        rest discarded.

        inert_nuclei is a list of nuclei (as Nucleus objects) that
        should be part of the collection but are not linked via reactions
        to the other nuclei in the network.

        symmetric_screening means that we screen the reverse rates
        using the same factor as the forward rates, for rates computed
        via detailed balance.

        Any combination of these options may be supplied.

        """

        self.files = []
        self.rates = []
        self.library = None

        self.inert_nuclei = inert_nuclei

        self.symmetric_screening = symmetric_screening
        self.do_screening = do_screening

        self.inert_nuclei = inert_nuclei

        if rate_files:
            if isinstance(rate_files, str):
                rate_files = [rate_files]
            self._read_rate_files(rate_files)

        if rates:
            if isinstance(rates, Rate):
                rates = [rates]
            for r in rates:
                if not isinstance(r, Rate):
                    raise ValueError('Expected Rate object or list of Rate objects passed as the rates argument.')
            rlib = Library(rates=rates)
            if not self.library:
                self.library = rlib
            else:
                self.library = self.library + rlib

        if libraries:
            if isinstance(libraries, Library):
                libraries = [libraries]
            for lib in libraries:
                if not isinstance(lib, Library):
                    raise ValueError('Expected Library object or list of Library objects passed as the libraries argument.')
            if not self.library:
                self.library = libraries.pop(0)
            for lib in libraries:
                self.library = self.library + lib

        if self.library:
            self.rates = self.rates + self.library.get_rates()

        if precedence:
            self._make_distinguishable(precedence)

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

        if self.inert_nuclei:
            for n in self.inert_nuclei:
                if isinstance(n, Nucleus):
                    nuc = n
                else:
                    nuc = Nucleus(n)
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
        self.approx_rates = []
        self.derived_rates = []

        for r in self.rates:
            if isinstance(r, ApproximateRate):
                self.approx_rates.append(r)
                for cr in r.get_child_rates():
                    assert cr.chapter != "t"
                    # child rates may be ReacLibRates or DerivedRates
                    # make sure we don't double count
                    if isinstance(cr, DerivedRate):

                        # Here we check whether this child rate is removed or not.
                        # removed means that this rate is never used on its own to connect two nuclei in the network
                        # it is only used in one or more ApproximateRate.
                        if cr not in self.rates:
                            cr.removed = True
                        else:
                            cr.removed = False

                        cr.fname = None
                        cr._set_print_representation()

                        if cr not in self.derived_rates:
                            self.derived_rates.append(cr)

                    else:
                        if cr not in self.rates:
                            cr.removed = True
                        else:
                            cr.removed = False

                        cr.fname = None
                        cr._set_print_representation()

                        if cr not in self.reaclib_rates:
                            self.reaclib_rates.append(cr)

            elif r.chapter == 't':
                self.tabular_rates.append(r)
            elif isinstance(r, DerivedRate):
                if r not in self.derived_rates:
                    self.derived_rates.append(r)
            elif isinstance(r.chapter, int):
                if r not in self.reaclib_rates:
                    self.reaclib_rates.append(r)
            else:
                raise NotImplementedError(f"Chapter type unknown for rate chapter {r.chapter}")

        self.all_rates = self.reaclib_rates + self.tabular_rates + self.approx_rates + self.derived_rates

    def _read_rate_files(self, rate_files):
        # get the rates
        self.files = rate_files
        for rf in self.files:
            # create the appropriate rate object first
            try:
                rate = load_rate(rf)
            except Exception as ex:
                raise Exception(f"Error reading rate from file: {rf}") from ex

            # now create a library:
            rflib = Library(rates=[rate])

            if not self.library:
                self.library = rflib
            else:
                self.library = self.library + rflib

    def get_forward_rates(self):
        """return a list of the forward (exothermic) rates"""

        # first handle the ones that have Q defined
        forward_rates = [r for r in self.rates if r.Q is not None and r.Q >= 0.0]

        # e-capture tabular rates don't have a Q defined, so just go off of the binding energy
        forward_rates += [r for r in self.rates if r.Q is None and r.reactants[0].nucbind <= r.products[0].nucbind]

        return forward_rates

    def get_reverse_rates(self):
        """return a list of the reverse (endothermic) rates)"""

        # first handle the ones that have Q defined
        reverse_rates = [r for r in self.rates if r.Q is not None and r.Q < 0.0]

        # e-capture tabular rates don't have a Q defined, so just go off of the binding energy
        reverse_rates += [r for r in self.rates if r.Q is None and r.reactants[0].nucbind > r.products[0].nucbind]

        return reverse_rates

    def find_reverse(self, forward_rate, reverse_rates=None):
        """given a forward rate, locate the rate that is its reverse"""

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
        """ return a list of RatePair objects, grouping the rates together
            by forward and reverse"""

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
        if reverse_rates:
            for rr in reverse_rates:
                rp = RatePair(reverse=rr)
                rate_pairs.append(rp)

        return rate_pairs

    def get_nuclei(self):
        """ get all the nuclei that are part of the network """
        return self.unique_nuclei

    def get_rates(self):
        """ get a list of the reaction rates in this network"""
        return self.rates

    def get_rate(self, rid):
        """ Return a rate matching the id provided.  Here rid should be
        the string return by Rate.fname"""
        try:
            return [r for r in self.rates if r.fname == rid][0]
        except IndexError:
            raise LookupError(f"rate identifier {rid!r} does not match a rate in this network.") from None

    def get_rate_by_nuclei(self, reactants, products):
        """given a list of reactants and products, return any matching rates"""
        _tmp = [r for r in self.rates if
                sorted(r.reactants) == sorted(reactants) and
                sorted(r.products) == sorted(products)]

        if not _tmp:
            return None
        return _tmp

    def get_rate_by_name(self, name):
        """given a rate in the form 'A(x,y)B' return the Rate"""

        reactants, products = _rate_name_to_nuc(name)
        _r = self.get_rate_by_nuclei(reactants, products)
        if _r is None:
            return None
        if len(_r) == 1:
            return _r[0]
        return _r

    def get_nuclei_needing_partition_functions(self):
        """return a list of the nuclei that require partition
        functions for one or more DerivedRates in the collection"""

        rates_with_pfs = [q for q in self.all_rates if isinstance(q, DerivedRate) and q.use_pf]

        if rates_with_pfs:
            nuclei_pfs = []
            for r in rates_with_pfs:
                nuclei_pfs += r.reactants + r.products
            return set(nuclei_pfs)
        return None

    def remove_nuclei(self, nuc_list):
        """remove the nuclei in nuc_list from the network along with any rates
        that directly involve them (this doesn't affect approximate rates that
        may have these nuclei as hidden intermediate links)"""

        rates_to_delete = []
        for nuc in nuc_list:
            nn = nuc
            if not isinstance(nuc, Nucleus):
                nn = Nucleus(nuc)
            for rate in self.rates:
                if nn in rate.reactants + rate.products:
                    print(f"looking to remove {rate}")
                    rates_to_delete.append(rate)

        for rate in set(rates_to_delete):
            self.rates.remove(rate)

        self._build_collection()

    def remove_rates(self, rates):
        """remove the Rate objects in rates from the network.  Note, if
        rate list is a dict, then the keys are assumed to be the rates
        to remove"""

        if isinstance(rates, Rate):
            self.rates.remove(rates)
        else:
            for r in rates:
                self.rates.remove(r)

        self._build_collection()

    def add_rates(self, rates):
        """add the Rate objects in rates from the network."""

        if isinstance(rates, Rate):
            self.rates.append(rates)

        else:
            for r in rates:
                self.rates.append(r)

        self._build_collection()

    def make_ap_pg_approx(self, intermediate_nuclei=None):
        """combine the rates A(a,g)B and A(a,p)X(p,g)B (and the reverse) into a single
        effective approximate rate."""

        # make sure that the intermediate_nuclei list are Nuclei objects
        _inter_nuclei_remove = []
        if intermediate_nuclei:
            for nn in intermediate_nuclei:
                if isinstance(nn, Nucleus):
                    _inter_nuclei_remove.append(nn)
                else:
                    _inter_nuclei_remove.append(Nucleus(nn))

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

            inter_nuc_Z = prim_nuc.Z + 1
            inter_nuc_A = prim_nuc.A + 3

            element = PeriodicTable.lookup_Z(inter_nuc_Z)

            inter_nuc = Nucleus(f"{element.abbreviation}{inter_nuc_A}")

            if intermediate_nuclei and inter_nuc not in _inter_nuclei_remove:
                continue

            # look for A(a,p)X
            _r = self.get_rate_by_nuclei([prim_nuc, Nucleus("he4")], [inter_nuc, Nucleus("p")])

            if _r:
                r_ap = _r[-1]
            else:
                continue

            # look for X(p,g)B
            _r = self.get_rate_by_nuclei([inter_nuc, Nucleus("p")], [prim_prod])

            if _r:
                r_pg = _r[-1]
            else:
                continue

            # look for reverse B(g,a)A
            _r = self.get_rate_by_nuclei([prim_prod], [prim_nuc, Nucleus("he4")])

            if _r:
                r_ga = _r[-1]
            else:
                continue

            # look for reverse B(g,p)X
            _r = self.get_rate_by_nuclei([prim_prod], [inter_nuc, Nucleus("p")])

            if _r:
                r_gp = _r[-1]
            else:
                continue

            # look for reverse X(p,a)A
            _r = self.get_rate_by_nuclei([inter_nuc, Nucleus("p")], [Nucleus("he4"), prim_nuc])

            if _r:
                r_pa = _r[-1]
            else:
                continue

            # build the approximate rates
            ar = ApproximateRate(r_ag, [r_ap, r_pg], r_ga, [r_gp, r_pa], approx_type="ap_pg")
            ar_reverse = ApproximateRate(r_ag, [r_ap, r_pg], r_ga, [r_gp, r_pa], is_reverse=True, approx_type="ap_pg")

            print(f"using approximate rate {ar}")
            print(f"using approximate rate {ar_reverse}")

            # approximate rates
            approx_rates += [ar, ar_reverse]

        # remove the old rates from the rate list and add the approximate rate
        for ar in approx_rates:
            for r in ar.get_child_rates():
                try:
                    self.rates.remove(r)

                    print(f"removing rate {r}")
                except ValueError:
                    pass

            # add the approximate rates
            self.rates.append(ar)

        # regenerate the links
        self._build_collection()

    def evaluate_rates(self, rho, T, composition, screen_func=None):
        """evaluate the rates for a specific density, temperature, and
        composition, with optional screening"""
        rvals = {}
        ys = composition.get_molar()
        y_e = composition.eval_ye()

        if screen_func is not None:
            screen_factors = self.evaluate_screening(rho, T, composition, screen_func)
        else:
            screen_factors = {}

        for r in self.rates:
            val = r.prefactor * rho**r.dens_exp * r.eval(T, rho * y_e)
            if (r.weak_type == 'electron_capture' and not isinstance(r, TabularRate)):
                val = val * y_e
            yfac = functools.reduce(mul, [ys[q] for q in r.reactants])
            rvals[r] = yfac * val * screen_factors.get(r, 1.0)

        return rvals

    def evaluate_jacobian(self, rho, T, comp, screen_func=None):
        """return an array of the form J_ij = dYdot_i/dY_j for the network"""

        # the rate.eval_jacobian_term does not compute the screening,
        # so we multiply by the factors afterwards
        if screen_func is not None:
            screen_factors = self.evaluate_screening(rho, T, comp, screen_func)
        else:
            screen_factors = {}

        nnuc = len(self.unique_nuclei)
        jac = np.zeros((nnuc, nnuc), dtype=np.float64)

        for i, n_i in enumerate(self.unique_nuclei):
            for j, n_j in enumerate(self.unique_nuclei):

                # we are considering dYdot(n_i) / dY(n_j)

                jac[i, j] = 0.0

                for r in self.nuclei_consumed[n_i]:
                    # how many of n_i are destroyed by this reaction
                    c = r.reactants.count(n_i)
                    jac[i, j] -= c * screen_factors.get(r, 1.0) *\
                        r.eval_jacobian_term(T, rho, comp, n_j)

                for r in self.nuclei_produced[n_i]:
                    # how many of n_i are produced by this reaction
                    c = r.products.count(n_i)
                    jac[i, j] += c * screen_factors.get(r, 1.0) *\
                        r.eval_jacobian_term(T, rho, comp, n_j)

        return jac

    def validate(self, other_library, forward_only=True, ostream=None):
        """perform various checks on the library, comparing to other_library,
        to ensure that we are not missing important rates.  The idea
        is that self should be a reduced library where we filtered out
        a few rates and then we want to compare to the larger
        other_library to see if we missed something important.

        ostream is the I/O stream to send output to (for instance, a
        file object or StringIO object).  If it is None, then output
        is to stdout.

        """

        current_rates = sorted(self.get_rates())

        # check the forward rates to see if any of the products are
        # not consumed by other forward rates

        passed_validation = True

        for rate in current_rates:
            if rate.reverse:
                continue
            for p in rate.products:
                found = False
                for orate in current_rates:
                    if orate == rate:
                        continue
                    if orate.reverse:
                        continue
                    if p in orate.reactants:
                        found = True
                        break
                if not found:
                    passed_validation = False
                    msg = f"validation: {p} produced in {rate} never consumed."
                    if ostream is None:
                        print(msg)
                    else:
                        ostream.write(msg + "\n")

        # now check if we are missing any rates from other_library with the exact same reactants

        other_by_reactants = collections.defaultdict(list)
        for rate in sorted(other_library.get_rates()):
            other_by_reactants[tuple(sorted(rate.reactants))].append(rate)

        for rate in current_rates:
            if forward_only and rate.reverse:
                continue

            key = tuple(sorted(rate.reactants))
            for other_rate in other_by_reactants[key]:
                # check to see if other_rate is already in current_rates
                found = True
                if other_rate not in current_rates:
                    found = False

                if not found:
                    msg = f"validation: missing {other_rate} as alternative to {rate} (Q = {other_rate.Q} MeV)."
                    if ostream is None:
                        print(msg)
                    else:
                        ostream.write(msg + "\n")

        return passed_validation

    def find_unimportant_rates(self, states, cutoff_ratio, screen_func=None):
        """evaluate the rates at multiple thermodynamic states, and find the
        rates that are always less than `cutoff_ratio` times the fastest rate
        for each state

        Here, states is a list of tuple of the form (density, temperature, composition),
        where composition is of type `Composition`.
        """
        largest_ratio = {r: 0 for r in self.rates}
        for rho, T, comp in states:
            rvals = self.evaluate_rates(rho, T, comp, screen_func)
            fastest = max(rvals.values())
            for r, value in rvals.items():
                largest_ratio[r] = max(largest_ratio[r], value / fastest)
        return {r: ratio for r, ratio in largest_ratio.items() if ratio < cutoff_ratio}

    def evaluate_screening(self, rho, T, composition, screen_func):
        """Evaluate the screening factors for each rate, using one of the
        methods in :py:mod:`pynucastro.screening`"""
        # this follows the same logic as BaseCxxNetwork._compute_screening_factors()
        factors = {}
        ys = composition.get_molar()
        plasma_state = make_plasma_state(T, rho, ys)
        screening_map = self.get_screening_map()

        for i, scr in enumerate(screening_map):
            if not (scr.n1.dummy or scr.n2.dummy):
                scn_fac = make_screen_factors(scr.n1, scr.n2)
                scor = screen_func(plasma_state, scn_fac)
            if scr.name == "he4_he4_he4":
                # we don't need to do anything here, but we want to avoid
                # immediately applying the screening
                pass
            elif scr.name == "he4_he4_he4_dummy":
                # make sure the previous iteration was the first part of 3-alpha
                assert screening_map[i - 1].name == "he4_he4_he4"
                # handle the second part of the screening for 3-alpha
                scn_fac2 = make_screen_factors(scr.n1, scr.n2)
                scor2 = screen_func(plasma_state, scn_fac2)

                # there might be both the forward and reverse 3-alpha
                # if we are doing symmetric screening
                for r in scr.rates:
                    # use scor from the previous loop iteration
                    factors[r] = scor * scor2
            else:
                # there might be several rates that have the same
                # reactants and therefore the same screening applies
                # -- handle them all now
                for r in scr.rates:
                    factors[r] = scor

        return factors

    def _evaluate_n_e(self, state, Xs):

        m_u = constants.value("unified atomic mass unit") * 1.0e3  # atomic unit mass in g

        n_e = 0.0
        for nuc in self.unique_nuclei:
            n_e += nuc.Z * state.dens * Xs[nuc] / (nuc.A * m_u)

        return n_e

    def _evaluate_mu_c(self, n_e, state, use_coulomb_corr=True):
        """ A helper equation that finds the mass fraction of each nuclide in NSE state,
        u[0] is chemical potential of proton  while u[1] is chemical potential of neutron"""

        # If there is proton included in network, upper limit of ye is 1
        # And if neutron is included in network, lower limit of ye is 0.
        # However, there are networks where either of them are included
        # So here I add a general check to find the upper and lower limit of ye
        # so the input doesn't go outside of the scope and the solver won't be able to converge if it did
        ye_low = min(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        ye_max = max(nuc.Z/nuc.A for nuc in self.unique_nuclei)
        assert state.ye >= ye_low and state.ye <= ye_max, "input electron fraction goes outside of scope for current network"

        # Seting up the constants needed to compute mu_c
        k = constants.value("Boltzmann constant") * 1.0e7          # boltzmann in erg/K
        Erg2MeV = 624151.0

        # These are three constants for calculating coulomb corrections of chemical energy, see Calders paper: iopscience 510709, appendix
        if use_coulomb_corr:
            A_1 = -0.9052
            A_2 = 0.6322
            A_3 = -0.5 * np.sqrt(3.0) - A_1 / np.sqrt(A_2)

        # u_c is the coulomb correction term for NSE
        # Calculate the composition at NSE, equations found in appendix of Calder paper
        u_c = {}
        for nuc in set(list(self.unique_nuclei) + [Nucleus("p")]):
            if use_coulomb_corr:
                Gamma = state.gamma_e_fac * n_e ** (1.0/3.0) * nuc.Z ** (5.0 / 3.0) / state.temp
                u_c[nuc] = Erg2MeV * k * state.temp * (A_1 * (np.sqrt(Gamma * (A_2 + Gamma)) - A_2 * np.log(np.sqrt(Gamma / A_2) +
                                      np.sqrt(1.0 + Gamma / A_2))) + 2.0 * A_3 * (np.sqrt(Gamma) - np.arctan(np.sqrt(Gamma))))
            else:
                u_c[nuc] = 0.0

        return u_c

    def _nucleon_fraction_nse(self, u, u_c, state):

        # Define constants: amu, boltzmann, planck, and electron charge
        m_u = constants.value("unified atomic mass unit") * 1.0e3  # atomic unit mass in g
        k = constants.value("Boltzmann constant") * 1.0e7          # boltzmann in erg/K
        h = constants.value("Planck constant") * 1.0e7             # in cgs
        Erg2MeV = 624151.0

        Xs = {}
        up_c = u_c[Nucleus("p")]

        for nuc in self.unique_nuclei:
            if nuc.partition_function:
                pf = nuc.partition_function(state.temp)
            else:
                pf = 1.0

            if not nuc.spin_states:
                raise ValueError(f"The spin of {nuc} is not implemented for now.")

            nse_exponent = (nuc.Z * u[0] + nuc.N * u[1] - u_c[nuc] + nuc.Z * up_c + nuc.nucbind * nuc.A) / k / state.temp / Erg2MeV
            nse_exponent = min(500.0, nse_exponent)

            Xs[nuc] = m_u * nuc.A_nuc * pf * nuc.spin_states / state.dens * (2.0 * np.pi * m_u * nuc.A_nuc * k * state.temp / h**2) ** (3. / 2.) \
                    * np.exp(nse_exponent)

        return Xs

    def _constraint_eq(self, u, u_c, state):
        """ Constraint Equations used to evaluate chemical potential for proton and neutron,
        which is used when evaluating composition at NSE"""

        # Don't use eval_ye() since it does automatic mass fraction normalization.
        # However, we should force normalization through constraint eq1.

        # m_u = constants.value("unified atomic mass unit") * 1.0e3  # atomic unit mass in g

        Xs = self._nucleon_fraction_nse(u, u_c, state)

        eq1 = 0.0
        for nuc in self.unique_nuclei:
            eq1 += Xs[nuc]

        eq1 += - 1.0

        eq2 = 0.0
        for nuc in self.unique_nuclei:
            eq2 += Xs[nuc] * nuc.Z / nuc.A
        # for nuc in self.unique_nuclei:
        #     eq2 += ((state.ye - 1.0) * nuc.Z + state.ye * nuc.N ) * Xs[nuc] / (nuc.A_nuc * m_u)

        eq2 += - state.ye
        return [eq1, eq2]

    def get_comp_nse(self, rho, T, ye, init_guess=(-3.5, -15), tol=1.0e-11, use_coulomb_corr=False, return_sol=False):
        """
        Returns the NSE composition given density, temperature and prescribed electron fraction
        using scipy.fsolve, `tol` is an optional parameter for the tolerance of scipy.fsolve.
        init_guess is optional, however one should change init_guess accordingly if unable or
        taking long time to find solution.
        One can enable printing the actual guess that found the solution after fine-tuning, which is useful
        when calling this method multiple times such as making a plot.
        """

        #From here we convert the init_guess list into a np.array object:

        init_guess = np.array(init_guess)
        state = NseState(T, rho, ye)

        u_c = {}
        for nuc in set(list(self.unique_nuclei) + [Nucleus("p")]):
            u_c[nuc] = 0.0

        Xs = {}

        j = 0
        init_guess = np.array(init_guess)
        is_pos_old = False
        found_sol = False

        # Filter out runtimewarnings from fsolve, here we check convergence by np.isclose
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        # This nested loops should fine-tune the initial guess if fsolve is unable to find a solution
        while (j < 20):
            i = 0
            guess = copy.deepcopy(init_guess)
            init_dx = 0.5

            while (i < 20):
                u = fsolve(self._constraint_eq, guess, args=(u_c, state), xtol=tol, maxfev=800)
                Xs = self._nucleon_fraction_nse(u, u_c, state)
                n_e = self._evaluate_n_e(state, Xs)
                u_c = self._evaluate_mu_c(n_e, state, use_coulomb_corr)

                res = self._constraint_eq(u, u_c, state)
                is_pos_new = all(k > 0 for k in res)
                found_sol = np.all(np.isclose(res, [0.0, 0.0], rtol=1.0e-7, atol=1.0e-7))

                if found_sol:
                    Xs = self._nucleon_fraction_nse(u, u_c, state)
                    comp = Composition(self.unique_nuclei)
                    comp.X = Xs

                    if return_sol:
                        return comp, u

                    return comp

                if is_pos_old != is_pos_new:
                    init_dx *= 0.8

                if is_pos_new:
                    guess -= init_dx
                else:
                    guess += init_dx

                is_pos_old = is_pos_new
                i += 1

            j += 1
            init_guess[0] -= 0.5

        raise ValueError("Unable to find a solution, try to adjust initial guess manually")

    def evaluate_ydots(self, rho, T, composition, screen_func=None):
        """evaluate net rate of change of molar abundance for each nucleus
        for a specific density, temperature, and composition"""

        rvals = self.evaluate_rates(rho, T, composition, screen_func)
        ydots = {}

        for nuc in self.unique_nuclei:

            # Rates that consume / produce nuc
            consuming_rates = self.nuclei_consumed[nuc]
            producing_rates = self.nuclei_produced[nuc]
            # Number of nuclei consumed / produced
            nconsumed = (r.reactants.count(nuc) for r in consuming_rates)
            nproduced = (r.products.count(nuc) for r in producing_rates)
            # Multiply each rate by the count
            consumed = (c * rvals[r] for c, r in zip(nconsumed, consuming_rates))
            produced = (c * rvals[r] for c, r in zip(nproduced, producing_rates))
            # Net change is difference between produced and consumed
            ydots[nuc] = sum(produced) - sum(consumed)

        return ydots

    def evaluate_energy_generation(self, rho, T, composition, screen_func=None):
        """evaluate the specific energy generation rate of the network for a specific
        density, temperature and composition"""

        ydots = self.evaluate_ydots(rho, T, composition, screen_func)
        enuc = 0.

        # compute constants and units
        m_n_MeV = constants.value('neutron mass energy equivalent in MeV')
        m_p_MeV = constants.value('proton mass energy equivalent in MeV')
        m_e_MeV = constants.value('electron mass energy equivalent in MeV')
        MeV2erg = (constants.eV * constants.mega) / constants.erg

        # ion binding energy contributions. basically e=mc^2
        for nuc in self.unique_nuclei:
            # add up mass in MeV then convert to erg
            mass = ((nuc.A - nuc.Z) * m_n_MeV + nuc.Z * (m_p_MeV + m_e_MeV) - nuc.A * nuc.nucbind) * MeV2erg
            enuc += ydots[nuc] * mass

        #convert from molar value to erg/g/s
        enuc *= -1*constants.Avogadro

        #subtract neutrino losses for tabular weak reactions
        for r in self.rates:
            if isinstance(r, TabularRate):
                # get composition
                ys = composition.get_molar()
                y_e = composition.eval_ye()

                # need to get reactant nucleus
                nuc = r.reactants[0]
                enuc -= constants.Avogadro * ys[nuc] * r.get_nu_loss(T, rho * y_e)

        return enuc

    def evaluate_activity(self, rho, T, composition, screen_func=None):
        """sum over all of the terms contributing to ydot,
        neglecting sign"""

        rvals = self.evaluate_rates(rho, T, composition, screen_func)
        act = {}

        for nuc in self.unique_nuclei:

            # Rates that consume / produce nuc
            consuming_rates = self.nuclei_consumed[nuc]
            producing_rates = self.nuclei_produced[nuc]
            # Number of nuclei consumed / produced
            nconsumed = (r.reactants.count(nuc) for r in consuming_rates)
            nproduced = (r.products.count(nuc) for r in producing_rates)
            # Multiply each rate by the count
            consumed = (c * rvals[r] for c, r in zip(nconsumed, consuming_rates))
            produced = (c * rvals[r] for c, r in zip(nproduced, producing_rates))
            # Net activity is sum of produced and consumed
            act[nuc] = sum(produced) + sum(consumed)

        return act

    def _get_network_chart(self, rho, T, composition):
        """a network chart is a dict, keyed by rate that holds a list of tuples (Nucleus, ydot)"""

        rvals = self.evaluate_rates(rho, T, composition)

        nc = {}

        for rate, rval in rvals.items():
            nucs = []
            for n in set(rate.reactants):
                nucs.append((n, -rate.reactants.count(n) * rval))
            for n in set(rate.products):
                nucs.append((n, rate.products.count(n) * rval))
            nc[rate] = nucs

        return nc

    def network_overview(self):
        """ return a verbose network overview """
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
        """ return a verbose network overview in terms of forward-reverse pairs"""
        ostr = ""
        for n in self.unique_nuclei:
            ostr += f"{n}\n"
            for rp in sorted(self.nuclei_rate_pairs[n]):
                ostr += f"     {rp}\n"
        return ostr

    def get_nuclei_latex_string(self):
        """return a string listing the nuclei in latex format"""

        ostr = ""
        for i, n in enumerate(self.unique_nuclei):
            ostr += f"${n.pretty}$"
            if i != len(self.unique_nuclei)-1:
                ostr += ", "
        return ostr

    def get_rates_latex_table_string(self):
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

    def get_screening_map(self):
        """a screening map is just a list of tuples containing the information
        about nuclei pairs for screening: (descriptive name of nuclei,
        nucleus 1, nucleus 2, rate, 1-based index of rate).  If symmetric_screening=True,
        then for reverse rates, we screen using the forward rate nuclei (assuming that we
        got here via detailed balance).

        """
        screening_map = []
        if not self.do_screening:
            return screening_map

        # we need to consider the child rates that come with ApproximateRate
        all_rates = []
        for r in self.rates:
            if isinstance(r, ApproximateRate):
                all_rates += r.get_child_rates()
            else:
                all_rates.append(r)

        for r in all_rates:
            screen_nuclei = r.ion_screen
            if self.symmetric_screening:
                screen_nuclei = r.symmetric_screen

            # screen_nuclei may be [] if it is a decay, gamma-capture, or neutron-capture
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

                # if we got here because nuc == "he4_he4_he4",
                # then we also have to add to "he4_he4_he4_dummy"

                if nucs == "he4_he4_he4":
                    scr2 = [q for q in screening_map if q.name == nucs + "_dummy"]
                    assert len(scr2) == 1

                    scr2[0].add_rate(r)

            else:

                # we handle 3-alpha specially -- we actually need
                # 2 screening factors for it

                if nucs == "he4_he4_he4":
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

    def write_network(self, *args, **kwargs):
        """Before writing the network, check to make sure the rates
        are distinguishable by name."""
        assert self._distinguishable_rates(), "ERROR: Rates not uniquely identified by Rate.fname"
        self._write_network(*args, **kwargs)

    def _distinguishable_rates(self):
        """Every Rate in this RateCollection should have a unique Rate.fname,
        as the network writers distinguish the rates on this basis."""
        names = [r.fname for r in self.rates]
        for n, r in zip(names, self.rates):
            k = names.count(n)
            if k > 1:
                print(f'Found rate {r} named {n} with {k} entries in the RateCollection.')
                print(f'Rate {r} has the original source:\n{r.original_source}')
                print(f'Rate {r} is in chapter {r.chapter}')
        return len(set(names)) == len(self.rates)

    def _make_distinguishable(self, precedence):
        """If multiple rates have the same name, eliminate the extraneous ones according to their
        labels' positions in the precedence list. Only do this if all of the labels have
        rankings in the list."""

        nameset = {r.fname for r in self.rates}
        precedence = {lab: i for i, lab in enumerate(precedence)}

        def sorting_key(i):
            return precedence[self.rates[i].label]

        for n in nameset:

            # Count instances of name, and cycle if there is only one
            ind = [i for i, r in enumerate(self.rates) if r.fname == n]
            k = len(ind)
            if k <= 1:
                continue

            # If there were multiple instances, use the precedence settings to delete extraneous
            # rates
            labels = [self.rates[i].label for i in ind]

            if all(lab in precedence for lab in labels):

                sorted_ind = sorted(ind, key=sorting_key)
                r = self.rates[sorted_ind[0]]
                for i in sorted(sorted_ind[1:], reverse=True):
                    del self.rates[i]
                print(f'Found rate {r} named {n} with {k} entries in the RateCollection.')
                print(f'Kept only entry with label {r.label} out of {labels}.')

    def _write_network(self, *args, **kwargs):
        """A stub for function to output the network -- this is implementation
        dependent."""
        # pylint: disable=unused-argument
        print('To create network integration source code, use a class that implements a specific network type.')

    def plot(self, outfile=None, rho=None, T=None, comp=None,
             size=(800, 600), dpi=100, title=None,
             ydot_cutoff_value=None, show_small_ydot=False,
             node_size=1000, node_font_size=13, node_color="#A0CBE2", node_shape="o",
             curved_edges=False,
             N_range=None, Z_range=None, rotated=False,
             always_show_p=False, always_show_alpha=False,
             hide_xp=False, hide_xalpha=False,
             highlight_filter_function=None,
             nucleus_filter_function=None, rate_filter_function=None):
        """Make a plot of the network structure showing the links between
        nuclei.  If a full set of thermodymamic conditions are
        provided (rho, T, comp), then the links are colored by rate
        strength.


        parameters
        ----------

        outfile: output name of the plot -- extension determines the type

        rho: density to evaluate rates with

        T: temperature to evaluate rates with

        comp: composition to evaluate rates with

        size: tuple giving width x height of the plot in pixels

        dpi: pixels per inch used by matplotlib in rendering bitmap

        title: title to display on the plot

        ydot_cutoff_value: rate threshold below which we do not show a
        line corresponding to a rate

        show_small_ydot: if true, then show visibile lines for rates below
        ydot_cutoff_value

        node_size: size of a node

        node_font_size: size of the font used to write the isotope in the node

        node_color: color to make the nodes

        node_shape: shape of the node (using matplotlib marker names)

        curved_edges: do we use arcs to connect the nodes?

        N_range: range of neutron number to zoom in on

        Z_range: range of proton number to zoom in on

        rotate: if True, we plot A - 2Z vs. Z instead of the default Z vs. N

        always_show_p: include p as a node on the plot even if we
        don't have p+p reactions

        always_show_alpha: include He4 as a node on the plot even if
        we don't have 3-alpha

        hide_xalpha=False: dont connect the links to alpha for heavy
        nuclei reactions of the form A(alpha,X)B or A(X,alpha)B,
        except if alpha is the heaviest product.

        hide_xp=False: dont connect the links to p for heavy
        nuclei reactions of the form A(p,X)B or A(X,p)B.

        highlight_filter_function: name of a custom function that
        takes a Rate object and returns true or false if we want
        to highlight the rate edge.

        nucleus_filter_funcion: name of a custom function that takes a
        Nucleus object and returns true or false if it is to be shown
        as a node.

        rate_filter_funcion: name of a custom function that takes a Rate
        object and returns true or false if it is to be shown as an edge.

        """

        G = nx.MultiDiGraph()
        G.position = {}
        G.labels = {}

        fig, ax = plt.subplots()
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes('right', size='15%', pad=0.05)

        #ax.plot([0, 0], [8, 8], 'b-')

        # in general, we do not show p, n, alpha,
        # unless we have p + p, 3-a, etc.
        hidden_nuclei = ["n"]
        if not always_show_p:
            hidden_nuclei.append("p")
        if not always_show_alpha:
            hidden_nuclei.append("he4")

        # nodes -- the node nuclei will be all of the heavies
        # add all the nuclei into G.node
        node_nuclei = []
        colors = []
        for n in self.unique_nuclei:
            if n.raw not in hidden_nuclei:
                node_nuclei.append(n)
                colors.append(node_color)
            else:
                for r in self.rates:
                    if r.reactants.count(n) > 1:
                        node_nuclei.append(n)
                        colors.append(node_color)
                        break

        # approx nuclei are given a different color
        for n in self.approx_nuclei:
            node_nuclei.append(n)
            colors.append("#555555")

        if nucleus_filter_function is not None:
            node_nuclei = list(filter(nucleus_filter_function, node_nuclei))
            # redo the colors:
            colors = []
            for n in node_nuclei:
                if n in self.approx_nuclei:
                    colors.append("#555555")
                else:
                    colors.append(node_color)

        for n in node_nuclei:
            G.add_node(n)
            if rotated:
                G.position[n] = (n.Z, n.A - 2*n.Z)
            else:
                G.position[n] = (n.N, n.Z)
            G.labels[n] = fr"${n.pretty}$"

        # get the rates for each reaction
        if rho is not None and T is not None and comp is not None:
            ydots = self.evaluate_rates(rho, T, comp)
        else:
            ydots = None

        # Do not show rates on the graph if their corresponding ydot is less than ydot_cutoff_value
        invisible_rates = set()
        if ydot_cutoff_value is not None:
            for r in self.rates:
                if ydots[r] < ydot_cutoff_value:
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

                    if ydots is None:
                        G.add_edges_from([(n, p)], weight=0.5,
                                         real=1, highlight=highlight)
                        continue

                    try:
                        rate_weight = math.log10(ydots[r])
                    except ValueError:
                        # if ydots[r] is zero, then set the weight
                        # to roughly the minimum exponent possible
                        # for python floats
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
            for sr in r.secondary_rates + r.secondary_reverse:
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

        # It seems that networkx broke backwards compatability, and 'zorder' is no longer a valid
        # keyword argument. The 'linewidth' argument has also changed to 'linewidths'.

        nx.draw_networkx_nodes(G, G.position,      # plot the element at the correct position
                               node_color=colors, alpha=1.0,
                               node_shape=node_shape, node_size=node_size, linewidths=2.0, ax=ax)

        nx.draw_networkx_labels(G, G.position, G.labels,   # label the name of element at the correct position
                                font_size=node_font_size, font_color="w", ax=ax)

        # now we'll draw edges in two groups -- real links and approximate links

        if curved_edges:
            connectionstyle = "arc3, rad = 0.2"
        else:
            connectionstyle = "arc3"

        real_edges = [(u, v) for u, v, e in G.edges(data=True) if e["real"] == 1]
        real_weights = [e["weight"] for u, v, e in G.edges(data=True) if e["real"] == 1]

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

        # plot the arrow of reaction
        real_edges_lc = nx.draw_networkx_edges(G, G.position, width=list(widths),
                                               edgelist=real_edges, edge_color=edge_color,
                                               connectionstyle=connectionstyle,
                                               node_size=node_size,
                                               edge_cmap=plt.cm.viridis, ax=ax)

        approx_edges = [(u, v) for u, v, e in G.edges(data=True) if e["real"] == 0]

        _ = nx.draw_networkx_edges(G, G.position, width=1,
                                   edgelist=approx_edges, edge_color="0.5",
                                   connectionstyle=connectionstyle,
                                   style="dotted", node_size=node_size, ax=ax)

        # plot invisible rates, rates that are below ydot_cutoff_value
        invis_edges = [(u, v) for u, v, e in G.edges(data=True) if e["real"] == -1]

        _ = nx.draw_networkx_edges(G, G.position, width=1,
                                   edgelist=invis_edges, edge_color="gray",
                                   connectionstyle=connectionstyle,
                                   style="dashed", node_size=node_size, ax=ax)

        # highlight edges
        highlight_edges = [(u, v) for u, v, e in G.edges(data=True) if e["highlight"]]

        _ = nx.draw_networkx_edges(G, G.position, width=5,
                                   edgelist=highlight_edges, edge_color="C0", alpha="0.25",
                                   connectionstyle=connectionstyle,
                                   node_size=node_size, ax=ax)

        if ydots is not None:
            pc = mpl.collections.PatchCollection(real_edges_lc, cmap=plt.cm.viridis)
            pc.set_array(real_weights)
            if not rotated:
                plt.colorbar(pc, ax=ax, label="log10(rate)")
            else:
                plt.colorbar(pc, ax=ax, label="log10(rate)", orientation="horizontal", fraction=0.05)

        Ns = [n.N for n in node_nuclei]
        Zs = [n.Z for n in node_nuclei]

        if not rotated:
            ax.set_xlim(min(Ns)-1, max(Ns)+1)
        else:
            ax.set_xlim(min(Zs)-1, max(Zs)+1)

        #plt.ylim(min(Zs)-1, max(Zs)+1)

        if not rotated:
            plt.xlabel(r"$N$", fontsize="large")
            plt.ylabel(r"$Z$", fontsize="large")
        else:
            plt.xlabel(r"$Z$", fontsize="large")
            plt.ylabel(r"$A - 2Z$", fontsize="large")

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

        if Z_range is not None and N_range is not None:
            if not rotated:
                ax.set_xlim(N_range[0], N_range[1])
                ax.set_ylim(Z_range[0], Z_range[1])
            else:
                ax.set_xlim(Z_range[0], Z_range[1])

        if not rotated:
            ax.set_aspect("equal", "datalim")

        fig.set_size_inches(size[0]/dpi, size[1]/dpi)

        if title is not None:
            fig.suptitle(title)

        if outfile is not None:
            plt.tight_layout()
            plt.savefig(outfile, dpi=dpi)

        return fig

    def plot_jacobian(self, outfile=None, rho=None, T=None, comp=None,
                      screen_func=None,
                      size=(800, 800), dpi=100):

        jac = self.evaluate_jacobian(rho, T, comp, screen_func=screen_func)

        valid_max = np.abs(jac).max()

        # pylint: disable-next=redundant-keyword-arg
        norm = SymLogNorm(valid_max/1.e10, vmin=-valid_max, vmax=valid_max)

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

    def plot_network_chart(self, outfile=None, rho=None, T=None, comp=None,
                           size=(800, 800), dpi=100, force_one_column=False):

        nc = self._get_network_chart(rho, T, comp)

        # find the limits
        _ydot = []
        for r in self.rates:
            for _, y in nc[r]:
                _ydot.append(y)

        _ydot = np.asarray(_ydot)
        valid_max = np.abs(_ydot[_ydot != 0]).max()

        # pylint: disable-next=redundant-keyword-arg
        norm = SymLogNorm(valid_max/1.e15, vmin=-valid_max, vmax=valid_max)

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

        fig, _ax = plt.subplots(1, npanes, constrained_layout=True)

        fig.set_size_inches(size[0]/dpi, size[1]/dpi)

        if npanes == 1:
            drate = len(self.rates)
        else:
            drate = (len(self.rates) + 1) // 2

        _rates = sorted(self.rates)

        for ipane in range(npanes):

            if npanes == 2:
                ax = _ax[ipane]
            else:
                ax = _ax

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
            ax.set_xticks(np.arange(len(self.unique_nuclei)), labels=[f"${n.pretty}$" for n in self.unique_nuclei], rotation=90)

            # each pane only has its subset of rates
            ax.set_yticks(np.arange(nrates), labels=[f"{r.pretty_string}" for irate, r in enumerate(_rates) if istart <= irate <= iend])

            im = ax.imshow(data, norm=norm, cmap=plt.cm.bwr)

            ax.set_aspect("equal")

            # Turn spines off and create white grid.
            ax.spines[:].set_visible(False)

            ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
            ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
            ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
            ax.tick_params(which="minor", bottom=False, left=False)

        if npanes == 1:
            fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.75)
        else:
            fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.25)

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

    def gridplot(self, comp=None, color_field="X", rho=None, T=None, **kwargs):
        """
        Plot nuclides as cells on a grid of Z vs. N, colored by *color_field*. If called
        without a composition, the function will just plot the grid with no color field.

        :param comp: Composition of the environment.
        :param color_field: Field to color by. Must be one of 'X' (mass fraction),
            'Y' (molar abundance), 'Xdot' (time derivative of X), 'Ydot' (time
            derivative of Y), or 'activity' (sum of contributions to Ydot of
            all rates, ignoring sign).
        :param rho: Density to evaluate rates at. Needed for fields involving time
            derivatives.
        :param T: Temperature to evaluate rates at. Needed for fields involving time
            derivatives.

        :Keyword Arguments:

            - *scale* -- One of 'linear', 'log', and 'symlog'. Linear by default.
            - *small* -- If using logarithmic scaling, zeros will be replaced with
              this value. 1e-30 by default.
            - *linthresh* -- Linearity threshold for symlog scaling.
            - *linscale* --  The number of decades to use for each half of the linear
              range. Stretches linear range relative to the logarithmic range.
            - *filter_function* -- A callable to filter Nucleus objects with. Should
              return *True* if the nuclide should be plotted.
            - *outfile* -- Output file to save the plot to. The plot will be shown if
              not specified.
            - *dpi* -- DPI to save the image file at.
            - *cmap* -- Name of the matplotlib colormap to use. Default is 'magma'.
            - *edgecolor* -- Color of grid cell edges.
            - *area* -- Area of the figure without the colorbar, in square inches. 64
              by default.
            - *no_axes* -- Set to *True* to omit axis spines.
            - *no_ticks* -- Set to *True* to omit tickmarks.
            - *no_cbar* -- Set to *True* to omit colorbar.
            - *cbar_label* -- Colorbar label.
            - *cbar_bounds* -- Explicit colorbar bounds.
            - *cbar_format* -- Format string or Formatter object for the colorbar ticks.
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
        cmap = mpl.cm.get_cmap(cmap)

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

            values = np.array([comp.X[nuc] for nuc in nuclei])

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
    """ interactively explore a rate collection """
    def __init__(self, rc, comp, **kwargs):
        """ take a RateCollection and a composition """
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
        """Perform interactive exploration of the network structure."""
        interact(self._make_plot, logrho=logrho, logT=logT)
