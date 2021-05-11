"""A collection of classes and methods to deal with collections of
rates that together make up a network."""

# Common Imports
from __future__ import print_function
import warnings
import functools
import math
import os

from operator import mul
from collections import OrderedDict

from ipywidgets import interact

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
import networkx as nx

# Import Rate
from pynucastro.rates import Rate, Nucleus, Library

mpl.rcParams['figure.dpi'] = 100

class Composition:
    """a composition holds the mass fractions of the nuclei in a network
    -- useful for evaluating the rates

    """
    def __init__(self, nuclei, small=1.e-16):
        """nuclei is an iterable of the nuclei (Nucleus objects) in the network"""
        if not isinstance(nuclei[0], Nucleus):
            raise ValueError("must supply an iterable of Nucleus objects")
        else:
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
        X_sum = sum([self.X[k] for k in self.X])

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
        
    def eval_abar(self):
        """ return the mean molecular weight """

        avec = np.zeros(len(self.X), dtype=np.int32)
        xvec = np.zeros(len(self.X), dtype=np.float64)
        for i, n in enumerate(self.X):
            avec[i] = n.A
            xvec[i] = self.X[n]
        return 1. / np.sum(xvec / avec)

    def __str__(self):
        ostr = ""
        for k in self.X:
            ostr += "  X({}) : {}\n".format(k, self.X[k])
        return ostr

class RateCollection:
    """ a collection of rates that together define a network """

    pynucastro_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    def __init__(self, rate_files=None, libraries=None, rates=None, precedence=()):
        """
        rate_files are the files that together define the network.  This
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

        Any combination of these options may be supplied.
        """

        self.files = []
        self.rates = []
        self.library = None

        if rate_files:
            if isinstance(rate_files, str):
                rate_files = [rate_files]
            self._read_rate_files(rate_files)

        if rates:
            if isinstance(rates, Rate):
                rates = [rates]
            try:
                for r in rates:
                    assert isinstance(r, Rate)
            except:
                print('Expected Rate object or list of Rate objects passed as the rates argument.')
                raise
            else:
                rlib = Library(rates=rates)
                if not self.library:
                    self.library = rlib
                else:
                    self.library = self.library + rlib

        if libraries:
            if isinstance(libraries, Library):
                libraries = [libraries]
            try:
                for lib in libraries:
                    assert isinstance(lib, Library)
            except:
                print('Expected Library object or list of Library objects passed as the libraries argument.')
                raise
            else:
                if not self.library:
                    self.library = libraries.pop(0)
                for lib in libraries:
                    self.library = self.library + lib

        if self.library:
            self.rates = self.rates + self.library.get_rates()
            
        if precedence:
            self._make_distinguishable(precedence)

        # get the unique nuclei
        u = []
        for r in self.rates:
            t = set(r.reactants + r.products)
            u = set(list(u) + list(t))

        self.unique_nuclei = sorted(u)

        # now make a list of each rate that touches each nucleus
        # we'll store this in a dictionary keyed on the nucleus
        self.nuclei_consumed = OrderedDict()
        self.nuclei_produced = OrderedDict()

        for n in self.unique_nuclei:
            self.nuclei_consumed[n] = [r for r in self.rates if n in r.reactants]
            self.nuclei_produced[n] = [r for r in self.rates if n in r.products]

        # Re-order self.rates so Reaclib rates come first,
        # followed by Tabular rates. This is needed if
        # reaclib coefficients are targets of a pointer array
        # in the Fortran network.
        # It is desired to avoid wasting array size
        # storing meaningless Tabular coefficient pointers.
        self.rates = sorted(self.rates,
                            key=lambda r: r.chapter == 't')

        self.tabular_rates = []
        self.reaclib_rates = []
        for n, r in enumerate(self.rates):
            if r.chapter == 't':
                self.tabular_rates.append(n)
            elif isinstance(r.chapter, int):
                self.reaclib_rates.append(n)
            else:
                print('ERROR: Chapter type unknown for rate chapter {}'.format(
                    str(r.chapter)))
                exit()

    def _read_rate_files(self, rate_files):
        # get the rates
        self.files = rate_files
        for rf in self.files:
            try:
                rflib = Library(rf)
            except:
                print("Error reading library from file: {}".format(rf))
                raise
            else:
                if not self.library:
                    self.library = rflib
                else:
                    self.library = self.library + rflib

    def add_nucleus(self, nuc):
        """Add a nucleus to this network without adding any rates that utilize it."""
        
        self.unique_nuclei.append(nuc)
        self.nuclei_consumed[nuc] = []
        self.nuclei_produced[nuc] = []
        
    def get_nuclei(self):
        """ get all the nuclei that are part of the network """
        return self.unique_nuclei
        
    def linking_nuclei(self, nuclei):
        """
        Return a new RateCollection/Network object containing only rates linking the
        given nuclei (parameter *nuclei*). Nuclei can be provided as an iterable of Nucleus
        objects or a list of abbreviations.
        """
        
        return self.__class__(libraries=self.library.linking_nuclei(nuclei))      

    def evaluate_rates(self, rho, T, composition):
        """evaluate the rates for a specific density, temperature, and
        composition"""
        rvals = OrderedDict()
        ys = composition.get_molar()
        y_e = composition.eval_ye()

        for r in self.rates:
            val = r.prefactor * rho**r.dens_exp * r.eval(T, rho * y_e)
            if (r.weak_type == 'electron_capture' and not r.tabular):
                val = val * y_e
            yfac = functools.reduce(mul, [ys[q] for q in r.reactants])
            rvals[r] = yfac * val

        return rvals
        
    def evaluate_ydots(self, rho, T, composition):
        """evaluate net rate of change of molar abundance for each nucleus
        for a specific density, temperature, and composition"""
        
        rvals = self.evaluate_rates(rho, T, composition)
        ydots = dict()
        
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
        
    def evaluate_activity(self, rho, T, composition):
        """sum over all of the terms contributing to ydot,
        neglecting sign"""
        
        rvals = self.evaluate_rates(rho, T, composition)
        act = dict()
        
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

    def network_overview(self):
        """ return a verbose network overview """
        ostr = ""
        for n in self.unique_nuclei:
            ostr += "{}\n".format(n)
            ostr += "  consumed by:\n"
            for r in self.nuclei_consumed[n]:
                ostr += "     {}\n".format(r.string)

            ostr += "  produced by:\n"
            for r in self.nuclei_produced[n]:
                ostr += "     {}\n".format(r.string)

            ostr += "\n"
        return ostr

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
                print('Found rate {} named {} with {} entries in the RateCollection.'.format(r, n, k))
                print('Rate {} has the original source:\n{}'.format(r, r.original_source))
                print('Rate {} is in chapter {}'.format(r, r.chapter))
        return len(set(names)) == len(self.rates)
        
    def _make_distinguishable(self, precedence):
        """If multiple rates have the same name, eliminate the extraneous ones according to their
        labels' positions in the precedence list. Only do this if all of the labels have
        rankings in the list."""
        
        nameset = set([r.fname for r in self.rates])
        precedence = {lab: i for i, lab in enumerate(precedence)}
        def sorting_key(i): return precedence[self.rates[i].label]
        
        for n in nameset:
            
            # Count instances of name, and cycle if there is only one
            ind = [i for i, r in enumerate(self.rates) if r.fname == n]
            k = len(ind)
            if k <= 1: continue
            
            # If there were multiple instances, use the precedence settings to delete extraneous
            # rates
            labels = [self.rates[i].label for i in ind]
            
            if all((lab in precedence for lab in labels)):
                
                sorted_ind = sorted(ind, key=sorting_key)
                r = self.rates[sorted_ind[0]]
                for i in sorted(sorted_ind[1:], reverse=True): del self.rates[i]
                print('Found rate {} named {} with {} entries in the RateCollection.'.format(r, n, k))
                print('Kept only entry with label {} out of {}.'.format(r.label, labels))

    def _write_network(self, *args, **kwargs):
        """A stub for function to output the network -- this is implementation
        dependent."""
        print('To create network integration source code, use a class that implements a specific network type.')
        return

    def plot(self, outfile=None, rho=None, T=None, comp=None, size=(800, 600), dpi=100, title=None,
             ydot_cutoff_value=None, always_show_p=False, always_show_alpha=False, filter_function=None):
        """Make a plot of the network structure showing the links between nuclei"""

        G = nx.MultiDiGraph()
        G.position = {}
        G.labels = {}

        fig, ax = plt.subplots()
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes('right', size='15%', pad=0.05)

        ax.plot([0, 0], [8, 8], 'b-')

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
        for n in self.unique_nuclei:
            if n.raw not in hidden_nuclei:
                node_nuclei.append(n)
            else:
                for r in self.rates:
                    if r.reactants.count(n) > 1:
                        node_nuclei.append(n)
                        break

        if filter_function is not None:
            node_nuclei = list(filter(filter_function, node_nuclei))
        
        for n in node_nuclei:
            G.add_node(n)
            G.position[n] = (n.N, n.Z)
            G.labels[n] = r"${}$".format(n.pretty)

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

        # edges
        for n in node_nuclei:
            for r in self.nuclei_consumed[n]:
                for p in r.products:
                    if p in node_nuclei:
                        # networkx doesn't seem to keep the edges in
                        # any particular order, so we associate data
                        # to the edges here directly, in this case,
                        # the reaction rate, which will be used to
                        # color it
                        if ydots is None:
                            G.add_edges_from([(n, p)], weight=0.5)
                        else:
                            if r in invisible_rates:
                                continue
                            try:
                                rate_weight = math.log10(ydots[r])
                            except ValueError:
                                # if ydots[r] is zero, then set the weight
                                # to roughly the minimum exponent possible
                                # for python floats
                                rate_weight = -308
                            except:
                                raise
                            G.add_edges_from([(n, p)], weight=rate_weight)
                            
        # It seems that networkx broke backwards compatability, and 'zorder' is no longer a valid
        # keyword argument. The 'linewidth' argument has also changed to 'linewidths'.

        try:
            nx.draw_networkx_nodes(G, G.position,      # plot the element at the correct position
                                   node_color="#A0CBE2", alpha=1.0,
                                   node_shape="o", node_size=1000, linewidth=2.0, zorder=10, ax=ax)
        except TypeError:
            nx.draw_networkx_nodes(G, G.position,      # plot the element at the correct position
                                   node_color="#A0CBE2", alpha=1.0,
                                   node_shape="o", node_size=1000, linewidths=2.0, ax=ax)

        try:
            nx.draw_networkx_labels(G, G.position, G.labels,   # label the name of element at the correct position
                                    font_size=13, font_color="w", zorder=100, ax=ax)
        except TypeError:
            nx.draw_networkx_labels(G, G.position, G.labels,   # label the name of element at the correct position
                                    font_size=13, font_color="w", ax=ax)

        # get the edges and weights coupled in the same order
        edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())

        # plot the arrow of reaction
        try:
            edges_lc = nx.draw_networkx_edges(G, G.position, width=3,    # plot the arrow of reaction
                                              edgelist=edges, edge_color=weights,
                                              node_size=1000,
                                              edge_cmap=plt.cm.viridis, zorder=1, ax=ax)
        except TypeError:
            edges_lc = nx.draw_networkx_edges(G, G.position, width=3,    # plot the arrow of reaction
                                              edgelist=edges, edge_color=weights,
                                              node_size=1000,
                                              edge_cmap=plt.cm.viridis, ax=ax)

        # for networkx <= 2.0 draw_networkx_edges returns a
        # LineCollection matplotlib type which we can use for the
        # colorbar directly.  For networkx >= 2.1, it is a collection
        # of FancyArrowPatch-s, which we need to run through a
        # PatchCollection.  See:
        # https://stackoverflow.com/questions/18658047/adding-a-matplotlib-colorbar-from-a-patchcollection

        if ydots is not None:
            pc = mpl.collections.PatchCollection(edges_lc, cmap=plt.cm.viridis)
            pc.set_array(weights)
            plt.colorbar(pc, label="log10(rate)")

        Ns = [n.N for n in node_nuclei]
        Zs = [n.Z for n in node_nuclei]

        plt.xlim(min(Ns)-1, max(Ns)+1)
        #plt.ylim(min(Zs)-1, max(Zs)+1)
        plt.xlabel(r"$N$", fontsize="large")
        plt.ylabel(r"$Z$", fontsize="large")

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

        ax.set_aspect("equal", "datalim")

        fig.set_size_inches(size[0]/dpi, size[1]/dpi)

        if title is not None:
            fig.suptitle(title)

        if outfile is None:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(outfile, dpi=dpi)

    def get_reaction_network_graph(self):
        G = nx.DiGraph()
        G.position = {}
        G.labels = {}

        # nodes -- the node nuclei will be all of the heavies
        # add all the nuclei into G.node
        node_nuclei = []
        for n in self.unique_nuclei:
            node_nuclei.append(n)

        for n in node_nuclei:
            G.add_node(n)
            G.position[n] = (n.N, n.Z)
            G.labels[n] = r"{}".format(n.raw)

        # edges
        for n in node_nuclei:
            for r in self.nuclei_consumed[n]:
                for p in r.products:
                    if p in node_nuclei:
                        G.add_edges_from([(n, p)], weight=1)

        return G
    
    @staticmethod        
    def _safelog(arr, small):
        
        arr = np.copy(arr)
        if np.any(arr < 0.0):
            raise ValueError("Negative values not allowed for logscale - try symlog instead.")
        zeros = arr == 0.0
        arr[zeros] = min(small, arr[~zeros].min() / 10)
        return np.log10(arr)
        
    @staticmethod
    def _symlog(arr, linthresh=1.0):
        
        assert linthresh >= 1.0
        neg = arr < 0.0
        arr = np.abs(arr)
        needslog = arr > linthresh
        
        arr[needslog] = np.log10(arr[needslog]) + linthresh
        arr[neg] *= -1
        return arr
    
    @staticmethod    
    def _scale(arr, minval=None, maxval=None):
        
        if minval is None: minval = arr.min()
        if maxval is None: maxval = arr.max()
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
        
        if kwargs: warnings.warn("Unrecognized keyword arguments: {}".format(kwargs.keys()))
        
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
            raise ValueError("Invalid color field: '{}'".format(color_field))
        
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
            if color_field == "xdot": values *= As
            
        elif color_field == "activity":
            
            if rho is None or T is None:
                raise ValueError("Need both rho and T to evaluate rates!")
            act = self.evaluate_activity(rho, T, comp)
            values = np.array([act[nuc] for nuc in nuclei])
            
        if scale == "log": values = self._safelog(values, small)
        elif scale == "symlog": values = self._symlog(values, linthresh)
        
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
            
            plt.tick_params \
            (
                axis = 'both',  
                which = 'both',      
                bottom = False,      
                left = False,
                labelbottom = False,
                labelleft = False
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
            cbar_norm = mpl.colors.Normalize(*cbar_bounds)
            smap = mpl.cm.ScalarMappable(norm=cbar_norm, cmap=cmap)
            
            if not cbar_label:
                
                capfield = color_field.capitalize()
                if scale == "log":
                    cbar_label = "log[{}]".format(capfield)
                elif scale == "symlog":
                    cbar_label = "symlog[{}]".format(capfield)
                else:
                    cbar_label = capfield
                
            fig.colorbar(smap, cax=cax, orientation="vertical",
                    label=cbar_label, format=cbar_format)
        
        # Show or save
        if outfile is None:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(outfile, dpi=dpi)

    def __repr__(self):
        string = ""
        for r in self.rates:
            string += "{}\n".format(r.string)
        return string


class Explorer(object):
    """ interactively explore a rate collection """
    def __init__(self, rc, comp, size=(800, 600),
                 ydot_cutoff_value=None,
                 always_show_p=False, always_show_alpha=False):
        """ take a RateCollection and a composition """
        self.rc = rc
        self.comp = comp
        self.size = size
        self.ydot_cutoff_value = ydot_cutoff_value
        self.always_show_p = always_show_p
        self.always_show_alpha = always_show_alpha

    def _make_plot(self, logrho, logT):
        self.rc.plot(rho=10.0**logrho, T=10.0**logT,
                     comp=self.comp, size=self.size,
                     ydot_cutoff_value=self.ydot_cutoff_value,
                     always_show_p=self.always_show_p,
                     always_show_alpha=self.always_show_alpha)

    def explore(self, logrho=(2, 6, 0.1), logT=(7, 9, 0.1)):
        """Perform interactive exploration of the network structure."""
        interact(self._make_plot, logrho=logrho, logT=logT)
