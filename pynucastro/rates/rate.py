""" Classes and methods to deal with a single nuclear reaction rate."""

import os
import re
import numpy as np
import matplotlib.pyplot as plt

from periodictable import elements

def list_known_rates():
    """ list the rates found in the library """

    lib_path = "{}/../library/".format(os.path.dirname(__file__))

    for _, _, filenames in os.walk(lib_path):
        for f in filenames:
            try:
                r = Rate(f)
            except:
                continue
            else:
                print("{:32} : {}".format(f, r))


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
            string = "{} += np.exp( ".format(prefix)
        else:
            string = "{} = np.exp( ".format(prefix)
        string += " {}".format(self.a[0])
        if not self.a[1] == 0.0:
            string += " + {}*tf.T9i".format(self.a[1])
        if not self.a[2] == 0.0:
            string += " + {}*tf.T913i".format(self.a[2])
        if not self.a[3] == 0.0:
            string += " + {}*tf.T913".format(self.a[3])
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            string += "\n{}         ".format(len(prefix)*" ")
        if not self.a[4] == 0.0:
            string += " + {}*tf.T9".format(self.a[4])
        if not self.a[5] == 0.0:
            string += " + {}*tf.T953".format(self.a[5])
        if not self.a[6] == 0.0:
            string += " + {}*tf.lnT9".format(self.a[6])
        string += ")"
        return string

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
            self.short_spec_name = "h1"
        elif name == "d":
            self.el = "H"
            self.A = 2
            self.short_spec_name = "h2"
        elif name == "t":
            self.el = "H"
            self.A = 3
            self.short_spec_name = "h3"
        elif name == "n":
            self.el = "n"
            self.A = 1
            self.short_spec_name = "n"
        else:
            e = re.match(r"([a-zA-Z]*)(\d*)", name)
            self.el = e.group(1).title()  # chemical symbol
            assert(self.el)
            self.A = int(e.group(2))
            assert(self.A >= 0)
            self.short_spec_name = name

        # atomic number comes from periodtable
        i = elements.isotope("{}-{}".format(self.A, self.el))
        self.Z = i.number
        assert(type(self.Z)==int)
        assert(self.Z >= 0)
        self.N = self.A - self.Z
        assert(type(self.N)==int)
        assert(self.N >= 0)

        # long name
        if i.name == 'neutron':
            self.spec_name = i.name
        else:
            self.spec_name = '{}-{}'.format(i.name, self.A)

        # latex formatted style
        self.pretty = r"{{}}^{{{}}}\mathrm{{{}}}".format(self.A, self.el)

    def __repr__(self):
        return self.raw

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        return self.el.lower() == other.el.lower() and \
               self.Z == other.Z and self.A == other.A

    def __lt__(self, other):
        if not self.Z == other.Z:
            return self.Z < other.Z
        else:
            return self.A < other.A

class LibraryFile(object):
    """ a single file containing one or many Reaclib rates,
    possibly containing multiple sets per rate. """

    def __init__(self, libfile):
        self._library_file = libfile
        self._rate_list = []
        self._library_source_lines = []

        # loop through library file, read lines
        try:
            flib = open(self._library_file, 'r')
        except:
            print('Could not open file {}'.format(self._library_file))
            raise
        for line in flib:
            ls = line.strip()
            if ls:
                self._library_source_lines.append(ls)
        flib.close()

        # identify distinct rates from library lines
        chapter = 0
        for i, line in enumerate(self._library_source_lines):
            # detect chapter if it's supplied
            try:
                chapter = int(line)
                continue
            except:
                if line == 't' or line == 'T':
                    chapter = 't'
                    continue
            # line was not a chapter, so see if it contains a nuclide
            ls = line.split()
            try:
                scratch = Nucleus(ls[0])
            except:
                # ls[0] is not a Nucleus, so line must be a Set line
                # do nothing, continue to the next line
                continue
            # ls[0] is a Nucleus, so create a Rate from this and
            # the following two lines.
            try:
                #NEED TO CONVERT STRINGS INTO A STRINGIO OBJECT?
                self._rate_list.append(Rate(self._library_source_lines[i:i+3],
                                           chapter=chapter))
            except:
                print('Could not extract Set from {}, line {} and following.'.format(
                    self._library_file, i))
                raise

    def rates(self):
        """ return the rates in this library file. """
        return self._rate_list[:]

    def text(self):
        """ return the source text of this library file. """
        return '\n'.join(self._library_source_lines)

class Rate(object):
    """ a single Reaclib rate, which can be composed of multiple sets """
    #NEED TO REVISE API FOR LIBRARY
    def __init__(self, rfile):
        self.rfile_path = rfile
        self.rfile = os.path.basename(rfile)
        self.chapter = None    # the Reaclib chapter for this reaction
        self.original_source = None   # the contents of the original rate file
        self.reactants = []
        self.products = []
        self.sets = []

        # Tells if this rate is eligible for screening
        # using screenz.f90 provided by StarKiller Microphysics.
        # If not eligible for screening, set to None
        # If eligible for screening, then
        # Rate.ion_screen is a 2-element list of Nucleus objects for screening
        self.ion_screen = None

        idx = self.rfile.rfind("-")
        self.fname = self.rfile[:idx].replace("--", "-").replace("-", "_")

        self.Q = 0.0

        # read in the file, parse the different sets and store them as
        # SingleSet objects in sets[]
        f = open(self.rfile_path, "r")

        lines = f.readlines()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        # catch table prescription
        if self.chapter != "t":
            self.chapter = int(self.chapter)

        # remove any black lines
        set_lines = [l for l in lines[1:] if not l.strip() == ""]

        if self.chapter == "t":
            # e1 -> e2, Tabulated
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)
            s4 = set_lines.pop(0)
            s5 = set_lines.pop(0)
            f = s1.split()
            self.reactants.append(Nucleus(f[0]))
            self.products.append(Nucleus(f[1]))

            self.table_file = s2.strip()
            self.table_header_lines = int(s3.strip())
            self.table_rhoy_lines = int(s4.strip())
            self.table_temp_lines = int(s5.strip())
            self.table_num_vars = 6 # Hard-coded number of variables in tables for now.
            self.table_index_name = 'j_{}_{}'.format(self.reactants[0], self.products[0])

        else:
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

                    elif self.chapter == 2:
                        # e1 -> e2 + e3
                        self.reactants.append(Nucleus(f[0]))
                        self.products += [Nucleus(f[1]), Nucleus(f[2])]

                    elif self.chapter == 3:
                        # e1 -> e2 + e3 + e4
                        self.reactants.append(Nucleus(f[0]))
                        self.products += [Nucleus(f[1]), Nucleus(f[2]), Nucleus(f[3])]

                    elif self.chapter == 4:
                        # e1 + e2 -> e3
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                        self.products.append(Nucleus(f[2]))

                    elif self.chapter == 5:
                        # e1 + e2 -> e3 + e4
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                        self.products += [Nucleus(f[2]), Nucleus(f[3])]

                    elif self.chapter == 6:
                        # e1 + e2 -> e3 + e4 + e5
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                        self.products += [Nucleus(f[2]), Nucleus(f[3]), Nucleus(f[4])]

                    elif self.chapter == 7:
                        # e1 + e2 -> e3 + e4 + e5 + e6
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1])]
                        self.products += [Nucleus(f[2]), Nucleus(f[3]),
                                          Nucleus(f[4]), Nucleus(f[5])]

                    elif self.chapter == 8:
                        # e1 + e2 + e3 -> e4
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                        self.products.append(Nucleus(f[3]))

                    elif self.chapter == 9:
                        # e1 + e2 + e3 -> e4 + e5
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1]), Nucleus(f[2])]
                        self.products += [Nucleus(f[3]), Nucleus(f[4])]

                    elif self.chapter == 10:
                        # e1 + e2 + e3 + e4 -> e5 + e6
                        self.reactants += [Nucleus(f[0]), Nucleus(f[1]),
                                           Nucleus(f[2]), Nucleus(f[3])]
                        self.products += [Nucleus(f[4]), Nucleus(f[5])]

                    elif self.chapter == 11:
                        # e1 -> e2 + e3 + e4 + e5
                        self.reactants.append(Nucleus(f[0]))
                        self.products += [Nucleus(f[1]), Nucleus(f[2]),
                                          Nucleus(f[3]), Nucleus(f[4])]

                    first = 0

                # the second line contains the first 4 coefficients
                # the third lines contains the final 3
                # we can't just use split() here, since the fields run into one another
                n = 13  # length of the field
                a = [s2[i:i+n] for i in range(0, len(s2), n)]
                a += [s3[i:i+n] for i in range(0, len(s3), n)]

                a = [float(e) for e in a if not e.strip() == ""]
                self.sets.append(SingleSet(a, label=label))

        # compute self.prefactor and self.dens_exp from the reactants
        self.prefactor = 1.0  # this is 1/2 for rates like a + a (double counting)
        self.inv_prefactor = 1
        for r in set(self.reactants):
            self.inv_prefactor = self.inv_prefactor * np.math.factorial(self.reactants.count(r))
        self.prefactor = self.prefactor/float(self.inv_prefactor)
        self.dens_exp = len(self.reactants)-1

        # determine if this rate is eligible for screening
        nucz = []
        for parent in self.reactants:
            if parent.Z != 0:
                nucz.append(parent)
        if len(nucz) > 1:
            nucz.sort(key=lambda x: x.Z)
            self.ion_screen = []
            self.ion_screen.append(nucz[0])
            self.ion_screen.append(nucz[1])

        self.string = ""
        self.pretty_string = r"$"

        # put p, n, and alpha second
        treactants = []
        for n in self.reactants:
            if n.raw not in ["p", "he4", "n"]:
                treactants.insert(0, n)
            else:
                treactants.append(n)

        for n, r in enumerate(treactants):
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
        """plot the rate's temperature sensitivity vs temperature"""

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T)

        plt.loglog(temps, r)

        plt.xlabel(r"$T$")

        if self.dens_exp == 0:
            plt.ylabel(r"\tau")
        elif self.dens_exp == 1:
            plt.ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            plt.ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        plt.title(r"{}".format(self.pretty_string))

        plt.show()
