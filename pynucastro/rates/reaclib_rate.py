"""Classes and methods for working with rates from the ReacLib
library.

"""

import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pynucastro.nucdata import Nucleus
from pynucastro.rates.files import RateFileError, _find_rate_file
from pynucastro.rates.rate import Rate, RateSource, Tfactors


class SingleSet:
    """A single ReacLib set for a reaction in the form:

    λ = exp[ a_0 + sum_{i=1}^5  a_i T_9**(2i-5)/3  + a_6 log T_9]

    A single rate in Reaclib can be composed of multiple sets.

    Parameters
    ----------
    a : list, numpy.ndarray
        the coefficients of the exponential fit
    labelprops : str
        a collection of flags that classify a ReacLib rate

    """

    def __init__(self, a, labelprops):
        self.a = a
        self.labelprops = labelprops
        self.label = None
        self.resonant = None
        self.weak = None
        self.derived_from_inverse = None

        self._update_label_properties()

    def _update_label_properties(self):
        """Set label and flags indicating Set is resonant, weak, or
        reverse.

        """
        assert isinstance(self.labelprops, str)
        assert len(self.labelprops) == 6

        self.label = self.labelprops[0:4]
        self.resonant = self.labelprops[4] == 'r'
        self.weak = self.labelprops[4] == 'w'
        self.derived_from_inverse = self.labelprops[5] == 'v'

    def __eq__(self, other):
        x = True

        for ai, aj in zip(self.a, other.a):
            x = x and (ai == aj)

        x = x and (self.label == other.label)
        x = x and (self.resonant == other.resonant)
        x = x and (self.weak == other.weak)
        x = x and (self.derived_from_inverse == other.derived_from_inverse)
        return x

    def f(self):
        """Return a function for ``rate(tf)`` where ``tf`` is a
        :py:class:`Tfactors <pynucastro.rates.rate.Tfactors>` object

        Returns
        -------
        Callable

        """
        return lambda tf: float(np.exp(self.a[0] +
                                       self.a[1]*tf.T9i +
                                       self.a[2]*tf.T913i +
                                       self.a[3]*tf.T913 +
                                       self.a[4]*tf.T9 +
                                       self.a[5]*tf.T953 +
                                       self.a[6]*tf.lnT9))

    def dfdT(self):
        """Return a function for the temperature derivative of the
        set, ``dratedT(tf)``, where ``tf`` is a :py:class:`Tfactors
        <pynucastro.rates.rate.Tfactors>` object

        """

        # we have lambda = exp(f(T_9))
        # so dlambda/dT9 = lambda * df/dT9
        # and dlambda/dT = dlambda/dT9 / 1.e9

        return lambda tf: self.f()(tf) * (-self.a[1] * tf.T9i * tf.T9i +
                                          -(1./3.) * self.a[2] * tf.T913i * tf.T9i +
                                          (1./3.) * self.a[3] * tf.T913i * tf.T913i +
                                          self.a[4] +
                                          (5./3.) * self.a[5] * tf.T913 * tf.T913 +
                                          self.a[6] * tf.T9i) / 1.e9

    def set_string_py(self, *, prefix="set", plus_equal=False):
        """Generate the python code needed to evaluate the set.

        Parameters
        ----------
        prefix : str
            variable name used to store the set
        plus_equal : bool
            do we add to the existing set? or create a new
            variable and initialize it to this set?

        Returns
        -------
        str

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

    def set_string_cxx(self, *, prefix="set", plus_equal=False,
                       with_exp=True):
        """
        Generate the C++ code needed to evaluate the set.

        Parameters
        ----------
        prefix : str
            variable name used to store the set
        plus_equal : bool
            do we add to the existing set? or create a new
            variable and initialize it to this set?
        with_exp : bool
            do we compute the set (``True``) or the log of the
            set (``False``)?  The later is useful if we also
            are computing the derivative.

        Returns
        -------
        str

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

    def dln_set_string_dT9_cxx(self, *, prefix="dset_dT",
                               plus_equal=False):
        """Generate the C++ code to evaluate d/dT9 ln(set).

        Parameters
        ----------
        prefix : str
            variable name used to store the set
        plus_equal : bool
            do we add to the existing set? or create a new
            variable and initialize it to this set?

        Returns
        -------
        str

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
    """A single reaction rate from the ReacLib library, which
    can be composed of multiple sets.

    Parameters
    ----------
    rfile : str, pathlib.Path, io.StringIO
        the data file or string containing the rate in ReacLib format.
    chapter : int
        the ReacLib chapter describing the number of reactants and products
    original_source : str
        the original source.  This is usually set automatically when
        reading ``rfile``, but can be manually provided when adding
        rates together.
    reactants : list(str), list(Nucleus)
        the reactants for the reaction
    products : list(str), list(Nucleus)
        the products for the reaction
    sets : list(SingleSet)
        the sets that make up the rate
    labelprops : str
        a collection of flags that classify a ReacLib rate
    Q : float
        the energy release (in MeV)

    Raises
    ------
    RateFileError
        If the rate file is not correctly formatted.
    UnsupportedNucleus
        If the nucleus is unknown to pynucastro

    """

    def __init__(self, rfile=None, chapter=None, original_source=None,
                 reactants=None, products=None, sets=None, labelprops=None, Q=None):
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

        self.tabular = False

        self.use_identical_particle_factor = True

        self.rate_eval_needs_rho = False
        self.rate_eval_needs_comp = False

        # some subclasses might define a stoichmetry as a dict{Nucleus}
        # that gives the numbers for the dY/dt equations
        self.stoichiometry = None

        if Q is None:
            self._set_q()
        else:
            self.Q = Q

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

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other):
        """Determine whether two Rate objects are equal.  They are
        equal if they contain identical reactants and products and if
        they contain the same SingleSet sets and if their chapters are
        equal.

        """

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
        """Combine the sets of two Rate objects if they describe the
        same reaction. Must be Reaclib rates.

        """
        assert self.reactants == other.reactants
        assert self.products == other.products
        assert self.chapter == other.chapter
        assert isinstance(self.chapter, int)
        assert self.label == other.label
        assert self.src == other.src
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
        """Call _update_resonance_combined and then
        _update_label_properties.

        """
        if labelprops:
            self.labelprops = labelprops

        # Update labelprops based on the Sets in this Rate
        # to set the resonance_combined flag properly
        self._update_resonance_combined()
        self._update_label_properties()

    def _update_resonance_combined(self):
        """Check the Sets in this Rate and updates the
        resonance_combined flag as well as self.labelprops[4]

        """
        sres = [s.resonant for s in self.sets]
        if True in sres and False in sres:
            self._labelprops_combine_resonance()

    def _labelprops_combine_resonance(self):
        """Update self.labelprops[4] = 'c'"""
        llp = list(self.labelprops)
        llp[4] = 'c'
        self.labelprops = ''.join(llp)

    def _update_label_properties(self):
        """Set label and flags indicating Rate is resonant, weak, or
        reverse.

        """
        assert isinstance(self.labelprops, str)
        if self.labelprops == "approx":
            self.label = "approx"
            self.src = ""
            self.resonant = False
            self.weak = False
            self.weak_type = None
            self.derived_from_inverse = False
        elif self.labelprops == "derived":
            self.label = "derived"
            self.src = ""
            self.resonant = False  # Derived may be resonant in some cases
            self.weak = False
            self.weak_type = None
            self.derived_from_inverse = False
        else:
            assert len(self.labelprops) == 6
            self.label = "reaclib"
            self.src = self.labelprops[0:4]
            self.source = RateSource.source(self.src)
            self.resonant = self.labelprops[4] == 'r'
            self.weak = self.labelprops[4] == 'w'
            if self.weak:
                if self.src.strip() == 'ec' or self.src.strip() == 'bec':
                    self.weak_type = 'electron_capture'
                else:
                    self.weak_type = self.src.strip().replace('+', '_pos_').replace('-', '_neg_')
            else:
                self.weak_type = None
            self.derived_from_inverse = self.labelprops[5] == 'v'

    def _read_from_file(self, f):
        """Given a file object, read rate data from the file.

        Parameters
        ----------
        f : io.TextIOWrapper, io.StringIO

        """
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
        """Given a file object, write rate data to the file.

        Parameters
        ----------
        f : io.TextIOWrapper, io.StringIO

        """

        if self.original_source is None:
            raise NotImplementedError(
                f"Original source is not stored for this rate ({self})."
                " At present, we cannot reconstruct the rate representation without"
                " storing the original source."
            )

        print(self.original_source, file=f)

    def function_string_py(self):
        """Return a string containing the python function that
        computes the rate.

        Returns
        -------
        str

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

    def function_string_cxx(self, dtype="double", specifiers="inline",
                            leave_open=False, extra_args=None):
        """Return a string containing the C++ function that computes
        the rate

        Parameters
        ----------
        dtype : str
            The C++ datatype to use for all declarations
        specifiers : str
            C++ specifiers to add before each function declaration
            (i.e. "inline")
        leave_open : bool
            If ``true``, then we leave the function unclosed (no "}"
            at the end).  This can allow additional functions to add
            to this output.
        extra_args : list, tuple
            A list of strings representing additional arguments that
            should be appended to the argument list when defining the
            function interface.

        Returns
        -------
        str

        """

        if extra_args is None:
            extra_args = ()

        args = ["const tf_t& tfactors", f"{dtype}& rate", f"{dtype}& drate_dT", *extra_args]
        fstring = ""
        fstring += "template <int do_T_derivatives>\n"
        fstring += f"{specifiers}\n"
        fstring += f"void rate_{self.fname}({', '.join(args)}) {{\n\n"
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
            fstring += "        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;\n"
            fstring += "    }\n\n"

        if not leave_open:
            fstring += "}\n\n"

        return fstring

    def eval(self, T, *, rho=None, comp=None,
             screen_func=None):
        """Evaluate the reaction rate for temperature T

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at (not needed for ReacLib
            rates), but needed for evaluating screening effects.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with (not needed for ReacLib rates),
            but needed for evaluating screening effects.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include screening correction.

        Returns
        -------
        float

        """

        tf = Tfactors(T)
        r = 0.0
        for s in self.sets:
            f = s.f()
            r += f(tf)

        scor = 1.0
        if screen_func is not None:
            if rho is None or comp is None:
                raise ValueError("rho (density) and comp (Composition) needs to be defined when applying electron screening.")
            scor = self.evaluate_screening(rho, T, comp, screen_func)

        r *= scor

        return r

    def eval_deriv(self, T, *, rho=None, comp=None):
        """Evaluate the derivative of reaction rate with respect to T.
        This currently does NOT consider electron screening effects.

        Parameters
        ----------
        T : float
            the temperature to evaluate the rate at
        rho : float
            the density to evaluate the rate at (not needed for ReacLib
            rates).
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with (not needed for ReacLib rates).

        Returns
        -------
        float

        """

        _ = rho  # unused by this subclass
        _ = comp  # unused by this subclass

        tf = Tfactors(T)
        drdT = 0.0
        for s in self.sets:
            dfdT = s.dfdT()
            drdT += dfdT(tf)

        return drdT

    def get_rate_exponent(self, T0, *, rho=None, comp=None,
                          screen_func=None):
        """For a rate written as a power law, r = r_0 (T/T0)**nu,
        return nu corresponding to T0. This also considers electron
        screening effect if screen_func is passed in.

        Parameters
        ----------
        T0 : float
            the temperature to base the power law from
        rho : float
            the density to evaluate the rate at (not needed for ReacLib
            rates), but needed for evaluating screening effects.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the rate with (not needed for ReacLib rates),
            but needed for evaluating screening effects.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate exponent will include screening correction.

        Returns
        -------
        float

        """

        # nu = dln r /dln T, so we need dr/dT
        r1 = self.eval(T0, rho=rho, comp=comp, screen_func=screen_func)
        dT = 1.e-8*T0
        r2 = self.eval(T0 + dT, rho=rho, comp=comp, screen_func=screen_func)

        drdT = (r2 - r1)/dT
        return (T0/r1)*drdT

    def plot(self, Tmin=1.e8, Tmax=1.6e9, rhoYmin=3.9e8, rhoYmax=2.e9,
             figsize=(10, 10), *, rho=None, comp=None, screen_func=None):
        """Plot the rate's temperature sensitivity vs temperature

        Parameters
        ----------
        Tmin : float
            minimum temperature for the plot
        Tmax : float
            maximum temperature for the plot
        rhoYmin : float
            unused for ReacLib rates
        rhoYmax : float
            unused for ReacLib rates
        figsize : tuple
            the horizontal, vertical size (in inches) for the plot
        rho : float
            the density to evaluate the screening effect.
        comp : float
            the composition (of type
            :py:class:`Composition <pynucastro.networks.rate_collection.Composition>`)
            to evaluate the screening effect.
        screen_func : Callable
            one of the screening functions from :py:mod:`pynucastro.screening`
            -- if provided, then the rate will include the screening correction.

        Returns
        -------
        matplotlib.figure.Figure

        """
        _ = (rhoYmin, rhoYmax)  # unused by this subclass

        fig, ax = plt.subplots(figsize=figsize)

        temps = np.logspace(np.log10(Tmin), np.log10(Tmax), 100)
        r = np.zeros_like(temps)

        for n, T in enumerate(temps):
            r[n] = self.eval(T, rho=rho, comp=comp, screen_func=screen_func)

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
