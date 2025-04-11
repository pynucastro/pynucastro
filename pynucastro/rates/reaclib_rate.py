import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pynucastro.nucdata import Nucleus
from pynucastro.rates.files import RateFileError, _find_rate_file
from pynucastro.rates.rate import Rate, RateSource, Tfactors


class SingleSet:
    """ a set in Reaclib is one piece of a rate, in the form

        lambda = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

    A single rate in Reaclib can be composed of multiple sets

    :param a: the coefficients of the exponential fit
    :param labelprops: a collection of flags that classify a ReacLib rate

    """

    def __init__(self, a, labelprops):
        """here a is iterable (e.g., list or numpy array), storing the
           coefficients, a0, ..., a6

        """
        self.a = a
        self.labelprops = labelprops
        self.label = None
        self.resonant = None
        self.weak = None
        self.derived_from_inverse = None

        self._update_label_properties()

    def _update_label_properties(self):
        """ Set label and flags indicating Set is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        assert len(self.labelprops) == 6

        self.label = self.labelprops[0:4]
        self.resonant = self.labelprops[4] == 'r'
        self.weak = self.labelprops[4] == 'w'
        self.derived_from_inverse = self.labelprops[5] == 'v'

    def __eq__(self, other):
        """ Determine whether two SingleSet objects are equal to each other. """
        x = True

        for ai, aj in zip(self.a, other.a):
            x = x and (ai == aj)

        x = x and (self.label == other.label)
        x = x and (self.resonant == other.resonant)
        x = x and (self.weak == other.weak)
        x = x and (self.derived_from_inverse == other.derived_from_inverse)
        return x

    def f(self):
        """ return a function for rate(tf) where tf is a Tfactors
        object """
        return lambda tf: float(np.exp(self.a[0] +
                                       self.a[1]*tf.T9i +
                                       self.a[2]*tf.T913i +
                                       self.a[3]*tf.T913 +
                                       self.a[4]*tf.T9 +
                                       self.a[5]*tf.T953 +
                                       self.a[6]*tf.lnT9))

    def dfdT(self):
        """ return a function for this dratedT(tf), where tf is a
        Tfactors object """

        # we have lambda = exp(f(T_9))
        # so dlambda/dT9 = lambda * df/dT9
        # and dlambda/dT = dlambda/dT9 / 1.e9

        return lambda tf: self.f()(tf) * (-self.a[1] * tf.T9i * tf.T9i +
                                          -(1./3.) * self.a[2] * tf.T913i * tf.T9i +
                                          (1./3.) * self.a[3] * tf.T913i * tf.T913i +
                                          self.a[4] +
                                          (5./3.) * self.a[5] * tf.T913 * tf.T913 +
                                          self.a[6] * tf.T9i) / 1.e9

    def set_string_py(self, prefix="set", plus_equal=False):
        """
        return a string containing the python code for this set
        """
        if plus_equal:
            string = f"{prefix} += np.exp( "
        else:
            string = f"{prefix} = np.exp( "
        string += f" {self.a[0]}"
        if not self.a[1] == 0.0:
            string += f" + {self.a[1]}*tf.T9i"
        if not self.a[2] == 0.0:
            string += f" + {self.a[2]}*tf.T913i"
        if not self.a[3] == 0.0:
            string += f" + {self.a[3]}*tf.T913"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]}*tf.T9"
        if not self.a[5] == 0.0:
            string += f" + {self.a[5]}*tf.T953"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]}*tf.lnT9"
        string += ")"
        return string

    def set_string_cxx(self, prefix="set", plus_equal=False, with_exp=True):
        """
        return a string containing the C++ code for this set
        """
        if plus_equal:
            string = f"{prefix} += "
        else:
            string = f"{prefix} = "
        if with_exp:
            string += "std::exp( "
        string += f" {self.a[0]}"
        if not self.a[1] == 0.0:
            string += f" + {self.a[1]} * tfactors.T9i"
        if not self.a[2] == 0.0:
            string += f" + {self.a[2]} * tfactors.T913i"
        if not self.a[3] == 0.0:
            string += f" + {self.a[3]} * tfactors.T913"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]} * tfactors.T9"
        if not self.a[5] == 0.0:
            string += f" + {self.a[5]} * tfactors.T953"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]} * tfactors.lnT9"
        if with_exp:
            string += ");"
        else:
            string += ";"
        if all(q == 0.0 for q in self.a[1:]):
            string += "\namrex::ignore_unused(tfactors);"
        return string

    def dln_set_string_dT9_cxx(self, prefix="dset_dT", plus_equal=False):
        """
        return a string containing the C++ code for d/dT9 ln(set)
        """
        if plus_equal:
            string = f"{prefix} += "
        else:
            string = f"{prefix} = "

        if all(q == 0.0 for q in self.a[1:]):
            string += "0.0;"
            return string

        if not self.a[1] == 0.0:
            string += f" {-self.a[1]} * tfactors.T9i * tfactors.T9i"
        if not self.a[2] == 0.0:
            string += f" + -(1.0/3.0) * {self.a[2]} * tfactors.T943i"
        if not self.a[3] == 0.0:
            string += f" + (1.0/3.0) * {self.a[3]} * tfactors.T923i"
        if not (self.a[4] == 0.0 and self.a[5] == 0.0 and self.a[6] == 0.0):
            indent = len(prefix)*" "
            string += f"\n{indent}         "
        if not self.a[4] == 0.0:
            string += f" + {self.a[4]}"
        if not self.a[5] == 0.0:
            string += f" + (5.0/3.0) * {self.a[5]} * tfactors.T923"
        if not self.a[6] == 0.0:
            string += f" + {self.a[6]} * tfactors.T9i"
        string += ";"
        return string


class ReacLibRate(Rate):
    """A single reaction rate.  Currently, this is a ReacLib rate, which
    can be composed of multiple sets, or a tabulated electron capture
    rate.

    :raises: :class:`.RateFileError`, :class:`.UnsupportedNucleus`
    """
    def __init__(self, rfile=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
        """ rfile can be either a string specifying the path to a rate file or
        an io.StringIO object from which to read rate information. """
        # pylint: disable=super-init-not-called

        self.rfile_path = None
        self.rfile = None
        self.source = None

        if isinstance(rfile, (str, Path)):
            rfile = Path(rfile)
            self.rfile_path = _find_rate_file(rfile)
            self.rfile = rfile.name

        self.chapter = chapter    # the Reaclib chapter for this reaction
        self.original_source = original_source   # the contents of the original rate file
        self.fname = None

        if reactants:
            self.reactants = Nucleus.cast_list(reactants)
        else:
            self.reactants = []

        if products:
            self.products = Nucleus.cast_list(products)
        else:
            self.products = []

        if sets:
            self.sets = sets
        else:
            self.sets = []

        # a modified rate is one where we manually changed some of its
        # properties

        self.modified = False

        self.labelprops = labelprops

        self.approx = self.labelprops == "approx"

        self.derived = self.labelprops == "derived"

        self.label = None
        self.resonant = None
        self.weak = None
        self.weak_type = None
        self.derived_from_inverse = None

        self.removed = None

        self.Q = Q

        self.tabular = False

        self.use_identical_particle_factor = True

        self.rate_eval_needs_rho = False
        self.rate_eval_needs_comp = False

        if isinstance(rfile, Path):
            # read in the file, parse the different sets and store them as
            # SingleSet objects in sets[]
            f = self.rfile_path.open()
        elif isinstance(rfile, io.StringIO):
            # Set f to the io.StringIO object
            f = rfile
        else:
            f = None

        if f:
            self._read_from_file(f)
            f.close()
        else:
            self._set_label_properties()

        self._set_rhs_properties()
        self._set_screening()
        self._set_print_representation()

    def _set_print_representation(self):
        """ compose the string representations of this Rate. """

        super()._set_print_representation()

        # This is used to determine which rates to detect as the same reaction
        # from multiple sources in a Library file, so it should not be unique
        # to a given source, e.g. wc12, but only unique to the reaction.
        reactants_str = '_'.join([repr(nuc) for nuc in self.reactants])
        products_str = '_'.join([repr(nuc) for nuc in self.products])
        self.fname = f'{reactants_str}__{products_str}'

        if self.weak:
            self.fname += f'__weak__{self.weak_type}'
        if self.modified:
            self.fname += "__modified"
        if self.approx:
            self.fname += "__approx"
        if self.derived:
            self.fname += "__derived"
        if self.removed:
            self.fname += "__removed"

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """ Determine whether two Rate objects are equal.
        They are equal if they contain identical reactants and products and
        if they contain the same SingleSet sets and if their chapters are equal."""

        if not isinstance(other, ReacLibRate):
            return False

        x = (self.chapter == other.chapter) and (self.products == other.products) and \
                (self.reactants == other.reactants)
        if not x:
            return x
        x = len(self.sets) == len(other.sets)
        if not x:
            return x

        for si in self.sets:
            scomp = False
            for sj in other.sets:
                if si == sj:
                    scomp = True
                    break
            x = x and scomp

        return x

    def __add__(self, other):
        """Combine the sets of two Rate objects if they describe the same
           reaction. Must be Reaclib rates."""
        assert self.reactants == other.reactants
        assert self.products == other.products
        assert self.chapter == other.chapter
        assert isinstance(self.chapter, int)
        assert self.label == other.label
        assert self.weak == other.weak
        assert self.weak_type == other.weak_type
        assert self.derived_from_inverse == other.derived_from_inverse

        if self.resonant != other.resonant:
            self._labelprops_combine_resonance()
        new_rate = ReacLibRate(chapter=self.chapter,
                               original_source='\n'.join([self.original_source,
                                                          other.original_source]),
                               reactants=self.reactants,
                               products=self.products,
                               sets=self.sets + other.sets,
                               labelprops=self.labelprops,
                               Q=self.Q)
        return new_rate

    def _set_label_properties(self, labelprops=None):
        """ Calls _update_resonance_combined and then
            _update_label_properties. """
        if labelprops:
            self.labelprops = labelprops

        # Update labelprops based on the Sets in this Rate
        # to set the resonance_combined flag properly
        self._update_resonance_combined()
        self._update_label_properties()

    def _update_resonance_combined(self):
        """ Checks the Sets in this Rate and updates the
            resonance_combined flag as well as
            self.labelprops[4] """
        sres = [s.resonant for s in self.sets]
        if True in sres and False in sres:
            self._labelprops_combine_resonance()

    def _labelprops_combine_resonance(self):
        """ Update self.labelprops[4] = 'c'"""
        llp = list(self.labelprops)
        llp[4] = 'c'
        self.labelprops = ''.join(llp)

    def _update_label_properties(self):
        """ Set label and flags indicating Rate is resonant,
            weak, or reverse. """
        assert isinstance(self.labelprops, str)
        if self.labelprops == "approx":
            self.label = "approx"
            self.resonant = False
            self.weak = False
            self.weak_type = None
            self.derived_from_inverse = False
        elif self.labelprops == "derived":
            self.label = "derived"
            self.resonant = False  # Derived may be resonant in some cases
            self.weak = False
            self.weak_type = None
            self.derived_from_inverse = False
        else:
            assert len(self.labelprops) == 6
            self.label = self.labelprops[0:4]
            self.resonant = self.labelprops[4] == 'r'
            self.weak = self.labelprops[4] == 'w'
            if self.weak:
                if self.label.strip() == 'ec' or self.label.strip() == 'bec':
                    self.weak_type = 'electron_capture'
                else:
                    self.weak_type = self.label.strip().replace('+', '_pos_').replace('-', '_neg_')
            else:
                self.weak_type = None
            self.derived_from_inverse = self.labelprops[5] == 'v'
            self.source = RateSource.source(self.label)

    def _read_from_file(self, f):
        """ given a file object, read rate data from the file. """
        lines = f.readlines()
        f.close()

        self.original_source = "".join(lines)

        # first line is the chapter
        self.chapter = lines[0].strip()
        self.chapter = int(self.chapter)

        # remove any blank lines
        set_lines = [line for line in lines[1:] if not line.strip() == ""]

        # the rest is the sets
        first = 1
        while len(set_lines) > 0:
            # check for a new chapter id in case of Reaclib v2 format
            check_chapter = set_lines[0].strip()
            try:
                # see if there is a chapter number preceding the set
                check_chapter = int(check_chapter)
                # check that the chapter number is the same as the first
                # set in this rate file
                if check_chapter != self.chapter:
                    raise RateFileError(f'read chapter {check_chapter}, expected chapter {self.chapter} for this rate set.')
                # get rid of chapter number so we can read a rate set
                set_lines.pop(0)
            except (TypeError, ValueError):
                # there was no chapter number, proceed reading a set
                pass

            # sets are 3 lines long
            s1 = set_lines.pop(0)
            s2 = set_lines.pop(0)
            s3 = set_lines.pop(0)

            # first line of a set has up to 6 nuclei, then the label,
            # and finally the Q value

            # get rid of first 5 spaces
            s1 = s1[5:]

            # next follows 6 fields of 5 characters containing nuclei
            # the 6 fields are padded with spaces
            f = []
            for i in range(6):
                ni = s1[:5]
                s1 = s1[5:]
                ni = ni.strip()
                if ni:
                    f.append(ni)

            # next come 8 spaces, so get rid of them
            s1 = s1[8:]

            # next is a 4-character set label and 2 character flags
            labelprops = s1[:6]
            s1 = s1[6:]

            # next come 3 spaces
            s1 = s1[3:]

            # next comes a 12 character Q value followed by 10 spaces
            Q = float(s1.strip())

            if first:
                self.Q = Q

                # what's left are the nuclei -- their interpretation
                # depends on the chapter

                chapter_dict = {
                    1: ((1,), (2,)),  # e1 -> e2
                    2: ((1,), (2, 3)),  # e1 -> e2 + e3
                    3: ((1,), (2, 3, 4)),  # e1 -> e2 + e3 + e4
                    4: ((1, 2), (3,)),  # e1 + e2 -> e3
                    5: ((1, 2), (3, 4)),  # e1 + e2 -> e3 + e4
                    6: ((1, 2), (3, 4, 5)),  # e1 + e2 -> e3 + e4 + e5
                    7: ((1, 2), (3, 4, 5, 6)),  # e1 + e2 -> e3 + e4 + e5 + e6
                    8: ((1, 2, 3), (4,)),  # e1 + e2 + e3 -> e4
                    9: ((1, 2, 3), (4, 5)),  # e1 + e2 + e3 -> e4 + e5
                    10: ((1, 2, 3, 4), (5, 6)),  # e1 + e2 + e3 + e4 -> e5 + e6
                    11: ((1,), (2, 3, 4, 5))  # e1 -> e2 + e3 + e4 + e5
                }

                try:
                    r, p = chapter_dict[self.chapter]
                    self.reactants += [Nucleus.from_cache(f[i-1]) for i in r]
                    self.products += [Nucleus.from_cache(f[j-1]) for j in p]

                    # support historical format, where chapter 8 also handles what are
                    # now chapter 9 rates
                    if self.chapter == 8 and len(f) == 5:
                        self.products.append(Nucleus.from_cache(f[4]))

                except KeyError as exc:
                    raise RateFileError(f'Chapter {self.chapter} could not be identified in {self.original_source}') from exc

                first = 0

            # the second line contains the first 4 coefficients
            # the third lines contains the final 3
            # we can't just use split() here, since the fields run into one another
            n = 13  # length of the field
            a = [s2[i:i+n] for i in range(0, len(s2), n)]
            a += [s3[i:i+n] for i in range(0, len(s3), n)]

            a = [float(e) for e in a if not e.strip() == ""]
            self.sets.append(SingleSet(a, labelprops=labelprops))
            self._set_label_properties(labelprops)

    def write_to_file(self, f):
        """ Given a file object, write rate data to the file. """

        if self.original_source is None:
            raise NotImplementedError(
                f"Original source is not stored for this rate ({self})."
                " At present, we cannot reconstruct the rate representation without"
                " storing the original source."
            )

        print(self.original_source, file=f)

    def get_rate_id(self):
        """ Get an identifying string for this rate.
        Don't include resonance state since we combine resonant and
        non-resonant versions of reactions. """

        srev = ''
        if self.derived_from_inverse:
            srev = 'derived_from_inverse'

        sweak = ''
        if self.weak:
            sweak = 'weak'

        ssrc = 'reaclib'

        return f'{self.rid} <{self.label.strip()}_{ssrc}_{sweak}_{srev}>'

    def function_string_py(self):
        """
        Return a string containing python function that computes the
        rate
        """

        fstring = ""
        fstring += "@numba.njit()\n"
        fstring += f"def {self.fname}(rate_eval, tf):\n"
        fstring += f"    # {self.rid}\n"
        fstring += "    rate = 0.0\n\n"

        for s in self.sets:
            fstring += f"    # {s.labelprops[0:5]}\n"
            set_string = s.set_string_py(prefix="rate", plus_equal=True)
            for t in set_string.split("\n"):
                fstring += "    " + t + "\n"

        fstring += "\n"
        fstring += f"    rate_eval.{self.fname} = rate\n\n"
        return fstring

    def function_string_cxx(self, dtype="double", specifiers="inline", leave_open=False, extra_args=()):
        """
        Return a string containing C++ function that computes the
        rate
        """

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.cname()}({', '.join(args)}) {{\n\n"
        fstring += f"    // {self.rid}\n\n"
        fstring += "    rate = 0.0;\n"
        fstring += "    drate_dT = 0.0;\n\n"
        fstring += f"    {dtype} ln_set_rate{{0.0}};\n"
        fstring += f"    {dtype} dln_set_rate_dT9{{0.0}};\n"
        fstring += f"    {dtype} set_rate{{0.0}};\n\n"

        for s in self.sets:
            fstring += f"    // {s.labelprops[0:5]}\n"
            set_string = s.set_string_cxx(prefix="ln_set_rate", plus_equal=False, with_exp=False)
            for t in set_string.split("\n"):
                fstring += "    " + t + "\n"
            fstring += "\n"

            fstring += "    if constexpr (do_T_derivatives) {\n"
            dln_set_string_dT9 = s.dln_set_string_dT9_cxx(prefix="dln_set_rate_dT9", plus_equal=False)
            for t in dln_set_string_dT9.split("\n"):
                fstring += "        " + t + "\n"
            fstring += "    }\n"
            fstring += "\n"

            fstring += "    // avoid underflows by zeroing rates in [0.0, 1.e-100]\n"
            fstring += "    ln_set_rate = std::max(ln_set_rate, -230.0);\n"
            fstring += "    set_rate = std::exp(ln_set_rate);\n"

            fstring += "    rate += set_rate;\n"

            fstring += "    if constexpr (do_T_derivatives) {\n"
            fstring += "        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;\n"
            fstring += "    }\n\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def eval(self, T, *, rho=None, comp=None):
        """ evauate the reaction rate for temperature T """

        tf = Tfactors(T)
        r = 0.0
        for s in self.sets:
            f = s.f()
            r += f(tf)

        return r

    def eval_deriv(self, T, *, rho=None, comp=None):
        """ evauate the derivative of reaction rate with respect to T """
        _ = rho  # unused by this subclass
        _ = comp  # unused by this subclass

        tf = Tfactors(T)
        drdT = 0.0
        for s in self.sets:
            dfdT = s.dfdT()
            drdT += dfdT(tf)

        return drdT

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

    def plot(self, Tmin=1.e8, Tmax=1.6e9, rhoYmin=3.9e8, rhoYmax=2.e9,
             figsize=(10, 10)):
        """plot the rate's temperature sensitivity vs temperature

        :param float Tmin:    minimum temperature for plot
        :param float Tmax:    maximum temperature for plot
        :param float rhoYmin: minimum electron density to plot (e-capture rates only)
        :param float rhoYmax: maximum electron density to plot (e-capture rates only)
        :param tuple figsize: figure size specification for matplotlib

        :return: a matplotlib figure object
        :rtype: matplotlib.figure.Figure

        """
        _ = (rhoYmin, rhoYmax)  # unused by this subclass

        fig, ax = plt.subplots(figsize=figsize)

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T)

        ax.loglog(temps, r)
        ax.set_xlabel(r"$T$")

        if self.dens_exp == 0:
            ax.set_ylabel(r"$\tau$")
        elif self.dens_exp == 1:
            ax.set_ylabel(r"$N_A <\sigma v>$")
        elif self.dens_exp == 2:
            ax.set_ylabel(r"$N_A^2 <n_a n_b n_c v>$")

        ax.set_title(fr"{self.pretty_string}")

        return fig
