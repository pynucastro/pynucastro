# parse the reaclib stuff

from __future__ import print_function

import glob
import os
import re

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from periodictable import elements


class Nucleus(object):
    """
    a nucleus that participates in a reaction -- we store it in a
    class to hold its properties, define a sorting, and give it a
    pretty printing string

    """
    def __init__(self, name):
        self.raw = name

        # element symbol and atomic weight
        if name == "p":
            self.el = "H"
            self.A = 1
        else:
            e = re.match("([a-zA-Z]*)(\d*)", name)
            self.el = e.group(1).title()  # chemical symbol
            self.A = int(e.group(2))

        # atomic number comes from periodtable
        i = elements.isotope("{}-{}".format(self.A, self.el))
        self.Z = i.number
        self.N = self.A - self.Z

        # latex formatted style
        self.pretty = r"{{}}^{{{}}}\mathrm{{{}}}".format(self.A, self.el)
 
    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        return self.raw == other.raw

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        else:
            return self.A < other.A


class Tfactors(object):
    """ precompute temperature factors for speed """

    def __init__(self, T):
        """ return the Tfactors object.  Here, T is temperature in Kelvin """
        self.T9 = T/1.e9
        self.T9i = 1.0/self.T9
        self.T913i = self.T9i**(1./3.)
        self.T913 = self.T9**(1./3.)
        self.T953 = self.T9**(5./3.)
        self.lnT9 = np.log(self.T9)


class SingleSet(object):
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

        A single rate in Reaclib can be composed of multiple sets
    """

    def __init__(self, a, label=None):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.label = label


    def f(self):
        """
        return a function for this set -- note: Tf here is a Tfactors
        object
        """
        return lambda tf: np.exp(self.a[0] +
                                 self.a[1]*tf.T9i +
                                 self.a[2]*tf.T913i +
                                 self.a[3]*tf.T913 +
                                 self.a[4]*tf.T9 +
                                 self.a[5]*tf.T953 +
                                 self.a[6]*tf.lnT9)


    def set_string(self, prefix="set", plus_equal=False):
        """
        return a string containing the python code for this set
        """
        if plus_equal:
            string =  "{} += np.exp( ".format(prefix)
        else:
            string =  "{} = np.exp( ".format(prefix)
        string += " {}".format(self.a[0])
        if not self.a[1] == 0.0: string += " + {}*tf.T9i".format(self.a[1])
        if not self.a[2] == 0.0: string += " + {}*tf.T913i".format(self.a[2])
        if not self.a[3] == 0.0: string += " + {}*tf.T913".format(self.a[3])
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            string += "\n{}         ".format(len(prefix)*" ")
        if not self.a[4] == 0.0: string += " + {}*tf.T9".format(self.a[4])
        if not self.a[5] == 0.0: string += " + {}*tf.T953".format(self.a[5])
        if not self.a[6] == 0.0: string += " + {}*tf.lnT9".format(self.a[6])
        string += ")"
        return string


class Rate(object):
    """ a single Reaclib rate, which can be composed of multiple sets """

    def __init__(self, file):
        self.file = os.path.basename(file)
        self.chapter = None    # the Reaclib chapter for this reaction
        self.original_source = None   # the contents of the original rate file
        self.reactants = []
        self.products = []
        self.sets = []

        self.dens_exp = 1
        self.prefactor = 1.0    # this is 1/2 for rates like a + a (double counting)

        idx = self.file.rfind("-")
        self.fname = self.file[:idx].replace("--","-").replace("-","_")

        self.Q = 0.0

        # read in the file, parse the different sets and store them as
        # SingleSet objects in sets[]
        f = open(file, "r")
        lines = f.readlines()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = int(lines[0])

        # remove any black lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        # the rest is the sets
        first = 1
        while len(set_lines) > 0:
            # sets are 3 lines long
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)

            # first line of a set has up to 6 nuclei, then the label,
            # and finally the Q value
            f = s1.split()
            Q = f.pop()
            label = f.pop()

            if first:
                self.Q = Q

                # what's left are the nuclei -- their interpretation
                # depends on the chapter
                if self.chapter == 1:
                    # e1 -> e2
                    self.reactants.append(Nucleus(f[0]))
                    self.products.append(Nucleus(f[1]))

                    self.string = "{} -> {}".format(*(self.reactants + self.products))
                    self.dens_exp = 0

                elif self.chapter == 2:
                    # e1 -> e2 + e3
                    self.reactants.append(Nucleus(f[0]))
                    self.products += [Nucleus(f[1]), Nucleus(f[2])]

                    self.string = "{} -> {} + {}".format(*(self.reactants + self.products))
                    self.dens_exp = 0

                elif self.chapter == 3:
                    # e1 -> e2 + e3 + e4
                    self.reactants.append(Nucleus(f[0]))
                    self.products += [Nucleus(f[1]), Nucleus(f[2]), Nucleus(f[3])]

                    self.string = "{} -> {} + {} + {}".format(*(self.reactants + self.products))
                    self.dens_exp = 0

                elif self.chapter == 4:
                    # e1 + e2 -> e3
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                    self.products.append(Nucleus(f[2]))

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 1

                elif self.chapter == 5:
                    # e1 + e2 -> e3 + e4
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                    self.products += [Nucleus(f[2]), Nucleus(f[3])]

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 1

                elif self.chapter == 6:
                    # e1 + e2 -> e3 + e4 + e5
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                    self.products += [Nucleus(f[2]), Nucleus(f[3]), Nucleus(f[4])]

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 1

                elif self.chapter == 7:
                    # e1 + e2 -> e3 + e4 + e5 + e6
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                    self.products += [Nucleus(f[2]), Nucleus(f[3]), Nucleus(f[4]), Nucleus(f[5])]

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./2.

                    self.dens_exp = 1

                elif self.chapter == 8:
                    # e1 + e2 + e3 -> e4
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                    self.products.append(Nucleus(f[3]))

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./6.  # 1/3!
                    elif len(set(self.reactants)) == 2:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                elif self.chapter == 9:
                    # e1 + e2 + e3 -> e4 + e5
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                    self.products += [Nucleus(f[3]), Nucleus(f[4])]

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./6.  # 1/3!
                    elif len(set(self.reactants)) == 2:
                        self.prefactor = 1./2.

                    self.dens_exp = 2

                elif self.chapter == 10:
                    # e1 + e2 + e3 + e4 -> e5 + e6
                    self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2]), Nucleus(f[3])]
                    self.products += [Nucleus(f[4]), Nucleus(f[5])]

                    if len(set(self.reactants)) == 1:
                        self.prefactor = 1./24.  # 1/4!
                    elif len(set(self.reactants)) == 2:
                        # there may be some instances where we have a + a + b + b,
                        # so prefactor = 1/4?
                        self.prefactor = 1./6. # 1/3!
                    elif len(set(self.reactants)) == 3:
                        self.prefactor = 1./2.

                    self.dens_exp = 3

                elif self.chapter == 11:
                    # e1 -> e2 + e3 + e4 + e5
                    self.reactants.append(Nucleus(f[0]))
                    self.products += [Nucleus(f[1]), Nucleus(f[2]), Nucleus(f[3]), Nucleus(f[4])]

                    self.dens_exp = 0


                self.string = ""
                self.pretty_string = r"$"
                for n, r in enumerate(self.reactants):
                    self.string += "{}".format(r)
                    self.pretty_string += r"{}".format(r.pretty)
                    if not n == len(self.reactants)-1:
                        self.string += " + "
                        self.pretty_string += r" + "

                self.string += " --> "
                self.pretty_string += r" \rightarrow "

                for n, p in enumerate(self.products):
                    self.string += "{}".format(p)
                    self.pretty_string += r"{}".format(p.pretty)
                    if not n == len(self.products)-1:
                        self.string += " + "
                        self.pretty_string += r" + "

                self.pretty_string += r"$"

                first = 0

            # the second line contains the first 4 coefficients
            # the third lines contains the final 3
            # we can't just use split() here, since the fields run into one another
            n = 13  # length of the field
            a = [s2[i:i+n] for i in range(0, len(s2), n)]
            a += [s3[i:i+n] for i in range(0, len(s3), n)]

            a = [float(e) for e in a if not e.strip() == ""]
            self.sets.append(SingleSet(a, label=label))


    def __repr__(self):
        return self.string

    def eval(self, T):
        """ evauate the reaction rate for temperature T """
        tf = Tfactors(T)
        r = 0.0
        for s in self.sets:
            f = s.f()
            r += f(tf)

        return r

    def get_rate_exponent(self, T0):
        """ 
        for a rate written as a power law, r = r_0 (T/T0)**nu, return
        nu corresponding to T0
        """
        
        # nu = dln r /dln T, so we need dr/dT
        r1 = self.eval(T0)
        dT = 1.e-8*T0
        r2 = self.eval(T0 + dT)

        drdT = (r2 - r1)/dT
        return (T0/r1)*drdT
        
    def plot(self, Tmin=1.e7, Tmax=1.e10):

        T = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(T)

        for n in range(len(T)):
            r[n] = self.eval(T[n])

        plt.loglog(T, r)

        plt.xlabel(r"$T$")

        if self.dens_exp == 0:
            plt.ylabel(r"\tau")
        elif self.dens_exp == 1:
            plt.ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            plt.ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        plt.title(r"{}".format(self.pretty_string))

        plt.show()


    def rate_string(self, indent=0, prefix="rate"):
        """
        return the functional form of rate as a function of
        the temperature (as Tfactors)
        """

        tstring = "# {}\n".format(self.string)
        tstring += "{} = 0.0\n\n".format(prefix)

        for s in self.sets:
            tstring += "# {}\n".format(s.label)
            tstring += "{}\n".format(s.set_string(prefix=prefix, plus_equal=True))


        string = ""
        for t in tstring.split("\n"):
            string += indent*" " + t + "\n"

        return string


    def function_string(self):
        """
        return a string containing python function that computes the
        rate
        """

        string = ""
        string += "def {}(tf):\n".format(self.fname)
        string += "{}".format(self.rate_string(indent=4))
        string += "    return rate\n\n"

        return string


    def ydot_string(self):
        """
        return a string containing the term in a dY/dt equation
        in a reaction network corresponding to this rate
        """

        # composition dependence
        Y_string = ""
        for n, r in enumerate(set(self.reactants)):
            c = self.reactants.count(r)
            if c > 1:
                Y_string += "Y[i{}]**{}".format(r, c)
            else:
                Y_string += "Y[i{}]".format(r, c)

            if n < len(set(self.reactants))-1:
                Y_string += "*"

        # density dependence
        if self.dens_exp == 0:
            dens_string = ""
        elif self.dens_exp == 1:
            dens_string = "rho*"
        else:
            dens_string = "rho**{}*".format(self.dens_exp)

        # prefactor
        if not self.prefactor == 1.0:
            prefactor_string = "{}*".format(self.prefactor)
        else:
            prefactor_string = ""

        return "{}{}{}*lambda_{}".format(prefactor_string, dens_string, Y_string, self.fname)


class RateCollection(object):
    """ a collection of rates that together define a network """

    def __init__(self, rate_files):
        """
        rate_files are the files that together define the network.  This
        can be any iterable or single string, and can include
        wildcards

        """

        self.files = []
        self.rates = []

        if type(rate_files) is str:
            rate_files = [rate_files]

        # get the rates
        for p in rate_files:
            self.files += glob.glob(p)

        for rf in self.files:
            self.rates.append(Rate(rf))


        # get the unique nuclei
        u = []
        for r in self.rates:
            t = list(set(r.reactants + r.products))
            u = list(set(u + t))

        self.unique_nuclei = sorted(u)

        # now make a list of each rate that touches each nucleus
        # we'll store this in a dictionary keyed on the nucleus
        self.nuclei_consumed = {}
        self.nuclei_produced = {}

        for n in self.unique_nuclei:
            self.nuclei_consumed[n] = []
            for r in self.rates:
                if n in r.reactants:
                    self.nuclei_consumed[n].append(r)

            self.nuclei_produced[n] = []
            for r in self.rates:
                if n in r.products:
                    self.nuclei_produced[n].append(r)


    def print_network_overview(self):
        for n in self.unique_nuclei:
            print(n)
            print("  consumed by: ")
            for r in self.nuclei_consumed[n]:
                print("     {} : {}".format(r.string, r.ydot_string()))

            print("  produced by: ")
            for r in self.nuclei_produced[n]:
                print("     {} : {}".format(r.string, r.ydot_string()))

            print(" ")

    
    def make_network(self, outfile="net.py"):
        """
        this is the actual RHS for the system of ODEs that
        this network describes
        """

        try: of = open(outfile, "w")
        except: raise

        of.write("import numpy as np\n")
        of.write("import reaclib\n\n")

        # integer keys
        for i, n in enumerate(self.unique_nuclei):
            of.write("i{} = {}\n".format(n, i))

        of.write("nnuc = {}\n\n".format(len(self.unique_nuclei)))

        of.write("A = np.zeros((nnuc), dtype=np.int32)\n\n")
        for n in self.unique_nuclei:
            of.write("A[i{}] = {}\n".format(n, n.A))

        of.write("\n")

        for r in self.rates:
            of.write(r.function_string())

        of.write("def rhs(t, Y, rho, T):\n\n")

        indent = 4*" "

        # get the rates
        of.write("{}tf = reaclib.Tfactors(T)\n\n".format(indent))
        for r in self.rates:
            of.write("{}lambda_{} = {}(tf)\n".format(indent, r.fname, r.fname))

        of.write("\n")

        of.write("{}dYdt = np.zeros((nnuc), dtype=np.float64)\n\n".format(indent))

        # now make the RHSs
        for n in self.unique_nuclei:
            of.write("{}dYdt[i{}] = (\n".format(indent, n))
            for r in self.nuclei_consumed[n]:
                c = r.reactants.count(n)
                if c == 1:
                    of.write("{}   -{}\n".format(indent, r.ydot_string()))
                else:
                    of.write("{}   -{}*{}\n".format(indent, c, r.ydot_string()))
            for r in self.nuclei_produced[n]:
                of.write("{}   +{}\n".format(indent, r.ydot_string()))
            of.write("{}   )\n\n".format(indent))

        of.write("{}return dYdt\n".format(indent))


    def plot(self):
        G = nx.Graph()
        G.position={}
        G.labels = {}

        plt.plot([0,0], [8,8], 'b-')

        # nodes
        for n in self.unique_nuclei:
            G.add_node(n)
            G.position[n] = (n.N, n.Z)
            G.labels[n] = r"${}$".format(n.pretty)

        # edges
        for n in self.unique_nuclei:
            for r in self.nuclei_consumed[n]:
                for p in r.products:
                    G.add_edge(n, p)

        nx.draw_networkx_nodes(G, G.position, 
                               node_color="blue", alpha=0.8, 
                               node_shape="s", node_size=1000)
        nx.draw_networkx_edges(G, G.position, edge_color="0.5")
        nx.draw_networkx_labels(G, G.position, G.labels, font_size=14, font_color="1.0")
        plt.show()

    def __repr__(self):
        string = ""
        for r in self.rates:
            string += "{}\n".format(r.string)
        return string



if __name__ == "__main__":
    r = Rate("examples/CNO/c13-pg-n14-nacr")
    print(r.rate_string(indent=3))
    print(r.eval(1.0e9))
    print(r.get_rate_exponent(2.e7))


    rc = RateCollection("examples/CNO/*-*")
    print(rc)
    rc.print_network_overview()
    rc.make_network()

    
