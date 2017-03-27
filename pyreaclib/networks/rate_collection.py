# Common Imports
from __future__ import print_function

import glob
import os

from ipywidgets import interact

import functools
from operator import mul

import math

import networkx as nx

import matplotlib
matplotlib.rcParams['figure.dpi'] = 100

import matplotlib.pyplot as plt

# Import Rate
from pyreaclib.rates import Rate, Nucleus

class Composition(object):
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
            if k.raw == "p":
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

    def __str__(self):
        ostr = ""
        for k in self.X:
            ostr += "  X({}) : {}\n".format(k, self.X[k])
        return ostr

class RateCollection(object):
    """ a collection of rates that together define a network """

    def __init__(self, rate_files, use_cse=False):
        """
        rate_files are the files that together define the network.  This
        can be any iterable or single string, and can include
        wildcards
        """

        self.pyreaclib_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.files = []
        self.rates = []
        self.use_cse = use_cse

        if type(rate_files) is str:
            rate_files = [rate_files]

        # get the rates
        self.pyreaclib_rates_dir = os.path.join(self.pyreaclib_dir,
                                                'rates')
        exit_program = False
        for p in rate_files:
            # check to see if the rate file is in the working dir
            fp = glob.glob(p)
            if fp:
                self.files += fp
            else:
                # check to see if the rate file is in pyreaclib/reaclib-rates
                fp = glob.glob(os.path.join(self.pyreaclib_rates_dir, p))
                if fp:
                    self.files += fp
                else: # Notify of all missing files before exiting
                    print('ERROR: File {} not found in {} or the working directory!'.format(
                        p,self.pyreaclib_rates_dir))
                    exit_program = True
        if exit_program:
            exit()

        for rf in self.files:
            try:
                self.rates.append(Rate(rf))
            except:
                print("Error with file: {}".format(rf))
                raise

        # get the unique nuclei
        u = []
        for r in self.rates:
            t = set(r.reactants + r.products)
            u = set(list(u) + list(t))

        self.unique_nuclei = sorted(u)

        # now make a list of each rate that touches each nucleus
        # we'll store this in a dictionary keyed on the nucleus
        self.nuclei_consumed = {}
        self.nuclei_produced = {}

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
                            key = lambda r: r.chapter=='t')

        self.tabular_rates = []
        self.reaclib_rates = []
        for n, r in enumerate(self.rates):
            if r.chapter == 't':
                self.tabular_rates.append(n)
            elif type(r.chapter)==int:
                self.reaclib_rates.append(n)
            else:
                print('ERROR: Chapter type unknown for rate chapter {}'.format(
                    str(r.chapter)))
                exit()

    def get_nuclei(self):
        """ get all the nuclei that are part of the network """
        return self.unique_nuclei


    def evaluate_rates(self, rho, T, composition):
        """evaluate the rates for a specific density, temperature, and
        composition"""
        rvals = {}
        ys = composition.get_molar()

        for r in self.rates:
            val = r.prefactor * rho**r.dens_exp * r.eval(T)
            yfac = functools.reduce(mul, [ys[q] for q in r.reactants])
            rvals[r] = yfac * val

        return rvals

    def print_network_overview(self):
        for n in self.unique_nuclei:
            print(n)
            print("  consumed by: ")
            for r in self.nuclei_consumed[n]:
                print("     {}".format(r.string))

            print("  produced by: ")
            for r in self.nuclei_produced[n]:
                print("     {}".format(r.string))

            print(" ")

    def write_network(self):
        print('To create network integration source code, use a class that implements a specific network type.')
        return

    def plot(self, outfile=None, rho=None, T=None, comp=None):
        G = nx.MultiDiGraph()
        G.position={}
        G.labels = {}

        plt.plot([0,0], [8,8], 'b-')

        # nodes -- the node nuclei will be all of the heavies, but not
        # p, n, alpha, unless we have p + p, 3-a, etc.
        node_nuclei = []
        for n in self.unique_nuclei:
            if n.raw not in ["p", "n", "he4"]:
                node_nuclei.append(n)
            else:
                for r in self.rates:
                    if r.reactants.count(n) > 1:
                        node_nuclei.append(n)
                        break

        for n in node_nuclei:
            G.add_node(n)
            G.position[n] = (n.N, n.Z)
            G.labels[n] = r"${}$".format(n.pretty)

        if rho is not None and T is not None and comp is not None:
            ydots = self.evaluate_rates(rho, T, comp)
        else:
            ydots = None

        #for rr in ydots:
        #    print("{}: {}".format(rr, ydots[rr]))

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
                            G.add_edges_from([(n, p)], weight=math.log10(ydots[r]))

        nx.draw_networkx_nodes(G, G.position,
                               node_color="#A0CBE2", alpha=1.0,
                               node_shape="o", node_size=1000, linewidth=2.0, zorder=10)

        nx.draw_networkx_labels(G, G.position, G.labels,
                                font_size=13, font_color="w", zorder=100)

        # get the edges and weights coupled in the same order
        edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())

        edges_lc = nx.draw_networkx_edges(G, G.position, width=3,
                                          edgelist=edges, edge_color=weights,
                                          edge_cmap=plt.cm.viridis, zorder=1)

        # draw_networkx_edges returns a LineCollection matplotlib type
        # which we can use for the colorbar
        if ydots is not None:
            plt.colorbar(edges_lc)

        Zs = [n.Z for n in node_nuclei]
        Ns = [n.N for n in node_nuclei]

        plt.xlim(min(Zs)-1, max(Zs)+1)
        plt.xlabel(r"$N$", fontsize="large")
        plt.ylabel(r"$Z$", fontsize="large")

        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        ax = plt.gca()
        ax.set_aspect("equal", "datalim")

        if outfile is None:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(outfile, dpi=100)

    def __repr__(self):
        string = ""
        for r in self.rates:
            string += "{}\n".format(r.string)
        return string


class Explorer(object):
    """ interactively explore a rate collection """
    def __init__(self, rc, comp):
        """ take a RateCollection and a composition """
        self.rc = rc
        self.comp = comp

    def _make_plot(self, logrho, logT):
        self.rc.plot(rho=10.0**logrho, T=10.0**logT, comp=self.comp)

    def explore(self):
        interact(self._make_plot, logrho=(2, 6, 0.1), logT=(7, 9, 0.1))

