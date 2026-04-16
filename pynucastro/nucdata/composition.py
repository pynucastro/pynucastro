"""Classes and methods for describing a collection of nuclei and their
abundances.

"""

import collections
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch

from pynucastro.nucdata.nucleus import Nucleus

# Current present day isotope mass fractions as obatained by Lodders
# et al. 2020 & 2021 The dataset is avialble at
# https://sites.wustl.edu/planetarychemistrylaboratory/data-tables/

LODDERS_DATA = {'h1': 0.7462, 'h2': 0.0, 'he3': 7.677e-05, 'he4': 0.2388, 'li6': 6.197e-10, 'li7': 8.842e-09,
            'be9': 1.375e-10, 'b10': 8.637e-10, 'b11': 3.798e-09, 'c12': 0.002596, 'c13': 2.941e-05, 'n14': 0.0007314,
            'n15': 1.775e-06, 'o16': 0.006362, 'o17': 2.362e-06, 'o18': 1.35e-05, 'f19': 5.767e-07, 'ne20': 0.00195,
            'ne21': 4.909e-06, 'ne22': 0.0001528, 'na23': 3.184e-05, 'mg24': 0.0004655, 'mg25': 6.166e-05, 'mg26': 7.035e-05,
            'al27': 5.288e-05, 'si28': 0.000618, 'si29': 3.249e-05, 'si30': 2.219e-05, 'p31': 6.13e-06, 's32': 0.0003184,
            's33': 2.583e-06, 's34': 1.489e-05, 's36': 5.515e-08, 'cl35': 3.36e-06, 'cl37': 1.134e-06, 'ar36': 7.124e-05,
            'ar38': 1.367e-05, 'ar40': 2.295e-08, 'k39': 3.14e-06, 'k40': 4.022e-10, 'k41': 2.382e-07, 'ca40': 5.313e-05,
            'ca42': 3.72e-07, 'ca43': 7.926e-08, 'ca44': 1.257e-06, 'ca46': 2.202e-09, 'ca48': 1.229e-07, 'sc45': 3.63e-08,
            'ti46': 2.235e-07, 'ti47': 2.059e-07, 'ti48': 2.079e-06, 'ti49': 1.56e-07, 'ti50': 1.532e-07, 'v50': 8.377e-10,
            'v51': 3.351e-07, 'cr50': 6.833e-07, 'cr52': 1.369e-05, 'cr53': 1.586e-06, 'cr54': 4.006e-07, 'mn55': 1.197e-05,
            'fe54': 6.591e-05, 'fe56': 0.001074, 'fe57': 2.524e-05, 'fe58': 3.415e-06, 'co59': 3.191e-06, 'ni58': 4.595e-05,
            'ni60': 1.838e-05, 'ni61': 8.103e-07, 'ni62': 2.626e-06, 'ni64': 6.893e-07, 'cu63': 5.579e-07, 'cu65': 2.567e-07,
            'zn64': 9.497e-07, 'zn66': 5.513e-07, 'zn67': 8.178e-08, 'zn68': 3.792e-07, 'zn70': 1.34e-08, 'ga69': 3.584e-08,
            'ga71': 2.447e-08, 'ge70': 4.138e-08, 'ge72': 5.687e-08, 'ge73': 1.625e-08, 'ge74': 7.775e-08, 'ge76': 1.692e-08,
            'as75': 1.09e-08, 'se74': 1.027e-09, 'se76': 1.133e-08, 'se77': 9.455e-09, 'se78': 2.987e-08, 'se80': 6.445e-08,
            'se82': 1.17e-08, 'br79': 1.184e-08, 'br81': 1.181e-08, 'kr78': 3.547e-10, 'kr80': 2.298e-09, 'kr82': 1.176e-08,
            'kr83': 1.18e-08, 'kr84': 5.867e-08, 'kr86': 1.809e-08, 'rb85': 1.034e-08, 'rb87': 4.081e-09, 'sr84': 2.614e-10,
            'sr86': 4.734e-09, 'sr87': 3.353e-09, 'sr88': 4.052e-08, 'y89': 9.266e-09, 'zr90': 1.211e-08, 'zr91': 2.67e-09,
            'zr92': 4.124e-09, 'zr94': 4.273e-09, 'zr96': 7.032e-10, 'nb93': 1.736e-09, 'mo92': 8.368e-10, 'mo94': 5.355e-10,
            'mo95': 9.369e-10, 'mo96': 9.95e-10, 'mo97': 5.781e-10, 'mo98': 1.478e-09, 'mo100': 6.056e-10, 'ru96': 2.298e-10,
            'ru98': 7.976e-11, 'ru99': 5.451e-10, 'ru100': 5.458e-10, 'ru101': 7.447e-10, 'ru102': 1.392e-09, 'ru104': 8.365e-10,
            'rh103': 8.334e-10, 'pd102': 3.443e-11, 'pd104': 3.824e-10, 'pd105': 7.739e-10, 'pd106': 9.566e-10, 'pd108': 9.437e-10,
            'pd110': 4.266e-10, 'ag107': 6.609e-10, 'ag109': 6.236e-10, 'cd106': 5.075e-11, 'cd108': 3.62e-11, 'cd110': 5.188e-10,
            'cd111': 5.368e-10, 'cd112': 1.022e-09, 'cd113': 5.221e-10, 'cd114': 1.239e-09, 'cd116': 3.305e-10, 'in113': 2.164e-11,
            'in115': 4.708e-10, 'sn112': 9.384e-11, 'sn114': 6.55e-11, 'sn115': 3.304e-11, 'sn116': 1.45e-09, 'sn117': 7.731e-10,
            'sn118': 2.458e-09, 'sn119': 8.775e-10, 'sn120': 3.364e-09, 'sn122': 4.849e-10, 'sn124': 6.175e-10, 'sb121': 5.939e-10,
            'sb123': 4.535e-10, 'te120': 1.436e-11, 'te122': 3.593e-10, 'te123': 1.266e-10, 'te124': 6.739e-10, 'te125': 1.009e-09,
            'te126': 2.697e-09, 'te128': 4.579e-09, 'te130': 4.961e-09, 'i127': 4.835e-09, 'xe124': 2.078e-11, 'xe126': 1.81e-11,
            'xe128': 3.739e-10, 'xe129': 4.654e-09, 'xe130': 7.439e-10, 'xe131': 3.751e-09, 'xe132': 4.579e-09, 'xe134': 1.716e-09,
            'xe136': 1.417e-09, 'cs133': 1.172e-09, 'ba130': 1.556e-11, 'ba132': 1.58e-11, 'ba134': 3.529e-10, 'ba135': 9.697e-10,
            'ba136': 1.162e-09, 'ba137': 1.676e-09, 'ba138': 1.078e-08, 'la138': 1.322e-12, 'la139': 1.528e-09, 'ce136': 6.513e-12,
            'ce138': 9.913e-12, 'ce140': 3.453e-09, 'ce142': 4.386e-10, 'pr141': 5.908e-10, 'nd142': 7.956e-10, 'nd143': 3.595e-09,
            'nd144': 7.069e-10, 'nd145': 2.639e-10, 'nd146': 5.174e-10, 'nd148': 1.737e-10, 'nd150': 1.724e-10, 'sm144': 2.896e-11,
            'sm147': 1.433e-10, 'sm148': 1.081e-10, 'sm149': 1.338e-10, 'sm150': 7.184e-11, 'sm152': 2.643e-10, 'sm154': 2.275e-10,
            'eu151': 1.732e-10, 'eu153': 1.916e-10, 'gd152': 4.619e-13, 'gd154': 2.783e-11, 'gd155': 1.901e-10, 'gd156': 2.645e-10,
            'gd157': 2.036e-10, 'gd158': 3.251e-10, 'gd160': 2.899e-10, 'tb159': 2.38e-10, 'dy156': 1.868e-13, 'dy158': 1.514e-12,
            'dy160': 3.64e-11, 'dy161': 2.965e-10, 'dy162': 4.027e-10, 'dy163': 3.958e-10, 'dy164': 4.521e-10, 'ho165': 3.521e-10,
            'er162': 1.383e-12, 'er164': 1.613e-11, 'er166': 3.416e-10, 'er167': 2.346e-10, 'er168': 2.784e-10, 'er170': 1.557e-10,
            'tm169': 1.631e-10, 'yb168': 1.207e-12, 'yb170': 3.054e-11, 'yb171': 1.458e-10, 'yb172': 2.253e-10, 'yb173': 1.682e-10,
            'yb174': 3.367e-10, 'yb176': 1.383e-10, 'lu175': 1.551e-10, 'lu176': 4.216e-12, 'hf174': 1.042e-12, 'hf176': 3.457e-11,
            'hf177': 1.225e-10, 'hf178': 1.808e-10, 'hf179': 9.09e-11, 'hf180': 2.35e-10, 'ta180': 1.121e-14, 'ta181': 9.321e-11,
            'w180': 8.623e-13, 'w182': 1.661e-10, 'w183': 9.03e-11, 'w184': 1.939e-10, 'w186': 1.818e-10, 're185': 8.641e-11,
            're187': 1.46e-10, 'os184': 5.73e-13, 'os186': 4.638e-11, 'os187': 4.927e-11, 'os188': 3.912e-09, 'os189': 4.797e-10,
            'os190': 7.843e-10, 'os192': 1.232e-09, 'ir191': 1.08e-09, 'ir193': 1.835e-09, 'pt190': 7.282e-13, 'pt192': 4.507e-11,
            'pt194': 1.891e-09, 'pt195': 1.958e-09, 'pt196': 1.473e-09, 'pt198': 4.297e-10, 'au197': 9.203e-10, 'hg196': 4.695e-12,
            'hg198': 1.802e-10, 'hg199': 3.051e-10, 'hg200': 4.168e-10, 'hg201': 2.359e-10, 'hg202': 5.42e-10, 'hg204': 1.271e-10,
            'tl203': 2.578e-10, 'tl205': 6.188e-10, 'pb204': 3.226e-10, 'pb206': 3.079e-09, 'pb207': 3.401e-09, 'pb208': 9.743e-09,
            'bi209': 7.08e-10, 'th232': 2.341e-10, 'u234': 2.725e-15, 'u235': 3.638e-13, 'u238': 5.076e-11}

class Composition(collections.UserDict):
    """A composition holds the mass fractions of the nuclei in a
    network.

    Parameters
    ----------
    nuclei : list, tuple
        an iterable of Nucleus objects
    small : float
        a floor for nuclei mass fractions, used as the default value.
    init : str
        Different modes to set up the initial composition. Valid choices
        are [`uniform`, `random`, `solar`]. If "solar" is selected, assume
        we work with metallicity, Z=0.02. If init=None, then no initial composition
        will be set.
    """

    def __init__(self, nuclei, small=1.e-16, init=None):
        try:
            super().__init__({Nucleus.cast(k): small for k in nuclei})
        except TypeError:
            raise ValueError("must supply an iterable of Nucleus objects or strings") from None

        if init == "uniform":
            self.set_equal()
        elif init == "random":
            self.set_random()
        elif init == "solar":
            self.set_solar_like(Z=0.02)
        elif init is None:
            pass
        else:
            raise ValueError(f"{init} is not valid. Choose from ['uniform','random','solar']")

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


class LoddersComposition(Composition):
    """A class to use present day solar abundances from Lodders et al. 2020 & 2021.

    Parameters
    ----------
    Z : float
        Target metallicity :math:`Z` to scale the Lodders solar
        mixture to. If ``None`` (the default), the unscaled Lodders
        abundances are used.
    half_life_thresh : float
        The half life value below which to zero the mass fraction
        of a nucleus.  This prevents us from making a composition
        that is not really stable.

    Notes
    -----
    The Lodders abundances are stored in ``LODDERS_DATA`` as mass
    fractions for individual isotopes.
    """

    def __init__(self, Z=None, half_life_thresh=None):

        nuclei = [Nucleus(name) for name in LODDERS_DATA]

        # base composition initialize in Composition
        super().__init__(nuclei)

        # now give Lodders abundances with scaling if needed (raw data)
        self._get_from_lodders(Z, half_life_thresh=half_life_thresh)

    def _get_from_lodders(self, Z=None, half_life_thresh=None):

        for name, X_raw in LODDERS_DATA.items():
            self[Nucleus(name)] = X_raw

        # normalize Lodders raw data to sum upto 1
        self.normalize(half_life_thresh=None)

        X_solar = 0.0
        Y_solar = 0.0

        for nuc, X_i in self.items():
            if nuc.Z == 1:
                X_solar += X_i
            elif nuc.Z == 2:
                Y_solar += X_i

        # default metallicity Z is what Lodders has 1 - X -Y
        Z_solar = 1.0 - X_solar - Y_solar

        # if no Z is given, we don't need to do anything else
        # if desired Z is given then we scale as following

        if Z is not None:
            Z_target = float(Z)
            Z_scale = Z_target / Z_solar   # scaling factor for metals

            sum_XY = X_solar + Y_solar
            XY_scale = (1.0 - Z_target) / sum_XY   # scaling factor for H and He

            for nuc in self:
                if nuc.Z in (1, 2):
                    self[nuc] *= XY_scale
                elif nuc.Z >= 3:
                    self[nuc] *= Z_scale

        self.normalize(half_life_thresh=half_life_thresh)
