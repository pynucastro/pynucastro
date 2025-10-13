# Changelog

## 2.8.0

  * `network_heloper` now takes a `verbose` flag (#1134)

  * numerous redundant examples (from the `examples/` directory)
    were removed (#1131)

  * He burning example has been updated (#1129)

  * `RateCollection.plot` can now use screening and take a
    matplotlib `GridSpec` that defines the axes (instead of
    creating them internally. (#1126)

  * `RateCollection` now takes a `verbose` option, and most
    informational prints are disabled by default (#1130)

  * Fix LaTeX rate rendering (#1128)

  * fix a mutability issue in `Library` addition (#1112)

  * add a method to `Library` to eliminate duplicate rates (#1111)

  * load all nuclear spins by default, even those that are
    "unreliable" (#1115)

  * fix a CI issue on `DerivedRate` (#1124)

  * for `ReacLibRate`, use `r.derived_from_inverse` instead of
    `r.reverse` (#938)

  * an electron-positron EOS was added.  This directly integrates the
    Fermi-Dirac integrals and constructs n, p, and e. (#1076)

  * `NumpyNetwork` has been removed -- it did not give much of a
    performance boost (#1109)

  * MPI utils have been moved from `reduction/` to a higher-level
    for reuse (#1104)

  * Clean-up of the reduction code and new example documentation
    (#1098, #1103)

  * update links docs CI block list (#1113)

  * `Composition` can now zero out short-lived species (#1107)

  * Weak rates from Pruet & Fuller for A = 65 to 80 have been added
    (#1091, #1093)

  * `RateFilter` can now filter on an endpoint (#1106)

  * The `write_to_file` method has been moved from `Library` to
    `ReacLibLibrary` since it only works with those rates. (#1105)

  * A bounds issue has been fixed in the partition function C++ code
    (#1100)

  * CI fixes due to python or Microphysics library changes (#1086,
    #1092, #1094, #1099, #1102)

  * added an adaptived-difference method for differencing (#1075)

  * `AmrexAstroCxxNetwork` templates now use the new neutrino cooling
    method (#1069)

  * Binding energy example has been updated with a better example (#1087)

  * Doc and API docstring improvements (#1079, #1082, #1083, #1084,
    #1116, #1117, #1127, #1132)

  * `pyproject.toml` updated with license info (#1081)

  * `plot_network_chart` now uses a `GridSpec` for a better layout
    (#1060)


## 2.7.1

  * fix a numba import for conda-forge (#1078)

## 2.7.0

  * add a `color_nodes_by_abundance` option to `RateCollection.plot()`
    (#1058)

  * Sort the edges before drawing network plots (#1045)

  * Include the pynucastro version information in
    `AmrexAstroCxxNetwork` output (#1028)

  * Updated `README.md` (#1065)

  * Improved error messages for `DerivedRate` (#1064)

  * Removed `setup.py` and `setup.cfg` (#1056)

  * Added pydocstyle CI (#1052)

  * Make library warnings more informative (#1049)

  * add a class that computes the Fermi-Dirac integrals (#1039, #1040,
    #1041, #1046, #1047, #1057, #1068, #1070)

  * added a new `network_helper` function to make simple networks
    more easily (#1025)

  * fix some Sphinx issues (#1017, #1024)

  * `AmrexAstroCxxNetwork` C++ cleaning and syncing with Microphysics
    (#993, #1003, #1009, #1010, #1020, #1022, #1036, #1038)

  * update CI due to Microphysics changes (#1021, #1023, #1029, #1050)

  * `DerivedRate` now works on a `ModifiedRate` (#1011, #1012).

  * `DerivedRate` is now consistent with NSE to machine precision
    (#1007)

  * `get_nuclei_in_range` is now more flexible in inputs (#1004)

  * The weak rates from Oda et al. have been added (#1005)

  * Creating a `Nucleus` now works with strings of the form "56ni",
    and also fixes some regex issues (#998, #1016)

  * A new method for evaluating the stiffness in a network was added
    (#1006, #1008)

  * docstring updates (#978, #989, #994, #1042, #1067)

  * a new method for exporting a NetworkX graph from a network was
    added (#990, #991)

  * general documentation fixes and API links (#980, #982, #984, #985,
    #987, #988, #1026, #1030, #1031, #1032, #1033, #1034, #1044,
    #1061, #1062, #1063, #1071, #1072, #1073, #1074)

  * add CI testing on MacOS (#1035)

## 2.6.0

  * Added a `summary()` method to `RateCollection` that gives details
    about the number of rates and nuclei (#954)

  * Added a new `ModifiedRate` class that allows for altering
    reactants, products, and stoichiometry while using an existing
    underlying rate (#950)

  * Drop requirement for `sphinx_rtd_theme` (#969)

  * `SimpleCxxNetwork` now stores and computes Y_e so it can work with
    weak rates (#965), and it also can take commandline arguments for
    rho and T (#970)

  * Fix `README.md` logo for PyPI (#973)

  * Add a new `BaryonConservationError` exception to `Rate` that
    ensures a rate is balanced with respect to baryons. (#975)

  * Added a new CI test that directly compares a `PythonNetwork` and
    `SimpleCxxNetwork`  evaluation of dY/dt (#964)

  * `SimpleCxxNetwork` now allows for screening to be disabled
    (compile with `DISABLE_SCREENING=TRUE`).  It also prints more
    precision in the test (#963)

  * Node colors in the `RateCollection` `plot()` can now be customized
    (#962)

  * Some helper functions for working with yt datasets have been added
    (#961)

  * update to the latest ReacLib database (#917)

  * don't use `std::pow` in constexpr in a `SimpleCxxNetwork` as that is
    C++23 and not yet widely supported.  (#953)

  * allow the nodes in a network plot to have custom names (#946)

  * create a `Nucleus` `.summary()` method to display the information
    we know about a nucleus. (#949)

  * allow rates to support their own stoichiometry, instead of simply
    relying on the number of nuclei in the reactants or products
    list. (#947, #952)

  * Added a new network method, `get_hidden_rates()` that returns the
    rates that are not explicitly linked in the network but instead
    are used internally by `ApproximateRate` rates. (#943)

  * Clean up `Library` (#932)

  * Fix the binding energy notebook example -- we don't want to show
    negative binding energies (#939)

  * `TabularLibrary()` now takes a keyword argument `ordering` that
    can specify which sources to use and the order to read
    them. (#922)

  * Docstring updates (#892, #899, #926, #927, #928, #931, #933, #934,
    #935, #937, #956, #957, #960, #966, #967, #972, #977) and more API
    coverage (#929) and reorganization of the API docs (#955)

  * New codespell github action (#930)

  * mailmap updates (#924)

  * Add ruff settings to `pyproject.toml` (#923)

  * Added the Fuller, Fowler, & Newman weak rates (#910)

  * Updated the colors used for `RateCollection` plots (#921)

  * Added the ability to show a legend on a `RateCollection` network
    plot (#817, #942)

  * Documentation updates (#901, #902, #907, #913, #919, #920) and
    fixes (#941, #976)

  * `PythonNetwork` now pairs the forward and reverse rates when
    summing the ydot term, mirroring what we do in C++ (#784)

  * Highlighting of edges now uses better colors in `RateCollection`
    network plots (#904)

  * `AmrexAstroCxxNetwork` code cleaning (#909)

  * `DerivedRate` now sets `reverse=True` (#915)

  * Weak rate tables have been renamed to include the source in the
    filename (#912)

  * CI fixes (#905, #906)

  * `RateCoillection` class cleaning (#900)

## 2.5.0

  * fix the zenodo authors (#898) and automate the zenodo bibtex in
    the docs (#880)

  * codespell fixes (#897)

  * `AmrexAstroCxxNetwork` now ifdefs out neutrinos (#884)

  * officially support python 3.13 (#860)

  * update package requirements (#896)

  * Implement Debye-Huckel screening (#870)

  * Fix a CODATA constants inconsistency between Nubase and SciPy
    (#877) and compatibility with SciPy < 1.15.0 (#889)

  * Address a new warning in matplotlib 3.10.0 (#879, #885)

  * Bump up jinja2 version (#871)

  * Added the ability to approximate A(n,gamma)X(n,gamma)B into an
    effective rate A(nn,gamma)B (#818)

  * `SimpleCxxNetwork` now includes electron screening using the
    prescription in Chugunov 2007 (#844)

  * added "NSE protons" that evaluate the same as a proton but can be
    used to decouple light and heavy reactions (#853, #854, #862, #863)

  * `modify_products` was moved to the base `Rate` class all rate classes
    can use it (#851)

  * Fix an issue with AMD HIP/ROCm memory access on the GPU in the
    partition functions (#859)

  * A `FortranNetwork` is now available that provides a wrapper around
    `SimpleCxxNetwork` with Fortran interfaces. (#841, #847)

  * `Nucleus` objects can now be added and subtracted (#835)

  * Fix a compiler issue with `SimpleCxxNetwork` and ReacLib weak rates
    (#846)

  * Documentation improvements (#782, #831, #832, #843, #848, #852,
    #855, #858, #861, #864, #866, #867, #874, #881, #887, #888),
    including new examples (#783, #800, #826) and build enhancements
    (#893)

  * Clean-up unused class data in `SimpleCxxNetwork` (#845)

  * Clean up the `rate` module by separating it into submodules (#840, #842)

  * Add the ability to label edges in network plots (#837)

  * Clean-up the `make_ap_pg_approx` code (#830, #850)

  * `PartitionFunction.eval` now returns a python float instead of a
    numpy float (#839)

  * Rates now store their source from publications (#822, #834)

  * Change some methods to properties in `Composition`, `Library`, and
    `Rate` (#829)

  * `PythonNetwork` rates now can take rho and Y. (#828)

  * CI fixes and updates (#827, #890)

  * Child rates in `ApproximateRate` are now stored as a dictionary
    (#824) + additional cleaning (#812)

  * Allow the identical particle factor in `Rate` to be disabled (#821)

  * `Rate.eval()` now take rho and comp (#820)

  * Update Zenodo authors (#815)

  * Use python's `pathlib.Path` throughout (#811)

  * `Nucleus` now stores the halflife (#801)

## 2.4.0

  * expanded documentation on exporting networks (#805) and physical
    constants (#793)

  * some doc test fixes (#806, #807, #809)

  * expanded unit testing (#802)

  * update pylint disabled features (#794)

  * dropped support for python 3.9 (#787)

  * we now compute the Q value for rates from masses instead of
    binding energies if it is not available (#797)

  * updated the nuclear properties to all be based on Nubase 2020,
    with the mass excesses being the primary quantity used and binding
    energies now derived from the mass excesses. (#788, #792, #796)

  * added a `get_all_nuclei()` method that returns every nucleus for
    which we have a mass (#798)

  * removed unused scripts in the nucdata module (#789)

  * fixed the labels for trace nuclei in a `Composition` pie chart
    (#781)

## 2.3.0

  * added a warning if we use the ReacLib n decay rate (#778)

  * doc updates, including an example of dealing with duplicate rates
    (#777), new theme (#775), axis label fixes (#774), some examples
    of nuclear astrophysics concepts (#762, #766)

  * fix the limits for `TabularRate` plots (#779)

  * created a new `alternate_rates` module and added the de Boer et
    al. 2017 C12(a,g)O16 rate (#769)

  * the minimum version of python supported is now 3.8 (#768)

  * some test fixes (#755, #757, #760, #763)

  * C++ autodiff changes / cleanups (#752, #759)

  * `AmrexAstroCxxNetwork` now create a ydot_weak function that
    evaluates the weak rate contributions to compute dYe/dt (#739)

## 2.2.1

  * numpy 2.0 support added (#748)

  * new logo! (#744, #746)

  * pynucastro can now be installed by conda (#741, #743)

  * minor test fixes (#729, #735, #740, #745)

  * doc fixes (#737)

  * pylint fixes (#733)

  * `Composition` now allows for direct key-value indexing (#731),
    and has a new `set_random()` method (#728)

  * The `screen5` screening routine has been added (#730)

  * C++ partition function interpolation is now cached (#726)

## 2.2.0

  * Switch to using `math.fsum()` instead of `sum()` for better python
    3.12 compatibility (#720)

  * Some modernization of the AMReX-Astro C++ code output (#709, #714,
    #717)

  * Documentation improvements, including references and contributing
    instructions, and NSE table (#669, #679, #695, #697, #706, #707)

  * Clean-up up some plot interfaces to always take rho, T, comp in
    that order (#693)

  * Performance improvements (#680, #685, #692)

  * A new `constants` module was added to provide fundamental constants
    (#688)

  * A new `NumpyNetwork` was added for use with the rate reduction
    algorithm (#684)

  * Neutrino cooling routines based on Itoh et al. 1996 were added
    (#628)

  * A new method for generating an NSE table was added (#676)

  * A new `NSENetwork` class was added for dealing with NSE (#675)

  * Moved some screening utilities out of `RateCollection` (#666)

  * A `PythonNetwork` now write out a function to compute the energy
    release (#640)

## 2.1.1

  * `RateCollection.evaluate_energy_generation()` can now return the
    neutrino losses from weak rates separately. (#665)

  * `AmrexAstroCxxNetwork` now capitalizes the species names to be
    consistent with other AMReX-Astro Microphysics networks (#652)

  * The pynucastro project configuration was updated and reorganized
    to `pyproject.toml` (#657, #660, #661)

  * `Library.get_rate_by_name` no longer returns None if one of
    the input rate names is not found (#659)

  * A C++ indexing error at the boundary of a weak rate table was
    fixed. (#654)

  * The C++ `AmrexAstroCxxNetwork` for evaluating weak rate neutrino
    losses was cleaned up to fix a bounds error (#655)

  * We now explicitly check in the C++ code if a weak rate table
    if found and contains data (#654)

  * a `RateCollection` now enforces that there are no duplicate rates.  A
    `find_duplicate_links()` method has been added to `Library` to allow
    one to construct a `Library` without duplicates and then create a
    `RateCollection` or other network from that. (#651)

  * Some documentation updates (#643, #650)

  * The positron rates in the Langanke weak rate tables have been
    merged with the corresponding electron rates and the n <-> p rate
    was added (#648)

  * The `Composition.bin_as()` method now takes a list of nuclei
    to exclude from being binned into (#642)

  * `RateCollection.evaluate_ydots()` now can take a filter to
    only evaluate a subset of the rates. (#641)

## 2.1.0

  * add `eval_zbar()` to Composition (#632)

  * fix `get_rate_by_name` to work with "pp" reactions (#632)

  * created a method to reduce a `Composition` from one set
    of nuclei to another based on the nuclei masses and
    charge number (#625)

  * switch `AmrexAstroCxxNetwork` to do bilinear interpolation
    in terms of log(T) and log(rhoY) (#592, #611)

  * tabular rates in python now do linear interpolation (#602)

  * added an example of creating a custom rate (#615)

  * Rate now calls `_set_q()` to set the Q value if it is not
    passed in (#617)

  * added a `TableInterpolator` that works both for interactive
    python and `PythonNetwork` (#612, 610, 609)

  * added a `RateCollection` method to find duplicate links (#601)

  * python networks with tabular rates now copy the table files (#605)

  * added a `get_nuclei_in_range` method to return a range of
    nuclei (#593)

  * we now do a binary search in the C++ partition function
    interpolation (#581)

  * added a simple C++ network (`SimpleCxxNetwork`) (#591, #585)

  * added the weak rates from Langanke 2001 (#536)

  * cleaned up partition functions in C++ (#578, 573, 569, 565)

  * converted the Suzuki tabular rates to be in terms of log (#550)

  * fixed a bounds issue in C++ table interpolation (#566)

  * eliminated a variable length array in the C++ table interpolation
    (#567)

  * added rate indices to the C++ networks (#563)

  * added a network reduction algorithm (#529, #528, #527, #526, #525,
    #523)

  * added a molar fraction method to `Composition` (#546)

  * added examples of interfacing with Julia (#539)

  * added a code of conduct (#504)

  * added gamma heating to tabular weak rates (#502)

## 2.0.3

  * add functions to numerically evaluate the Jacobian (#467)

  * Fix C++ compilation with partition functions and approximate
    rates and derived rates (#407, #494)

  * The Explorer class now takes all of the options as
    `RateCollection.plot()` (#496)

  * Rates edges can be highlighted via a user-defined function in the
    `RateCollection` plot (#486)

  * All of the plot routines in `RateCollection` now return a
    matplotlib `Figure` object (#490)

  * add `get_rate_by_nuclei` and `add_rate` to `Library()` (#484)

  * `RateCollection.remove_rates()` can now act on a single rate and a
    new method, `add_rates()`, will add a rate directly to the
    collection.

  * Fixed printing of `TabularLibrary` (#479)

  * Documentation updates (#473)

  * added `get_rate_by_name` to `RateCollection` (#468)

  * added tabular rates to `PythonNetwork` (#451)

  * fixed NSE state for network without protons (#463)

## 2.0.2

  * Fix the installation of the C++ templates via `setup.py` (#460)

## 2.0.1

  * Fix some `setup.py` metadata for PyPI (#453)

  * Change how tabular rate tables are stored and read (#450)

  * Support historical ReacLib format (#449)

## 2.0.0

  * Added a method to get a rate by short name: A(x,y)B (#438)

  * Added the option to show rate links hidden by ydot <
    `ydot_cutoff_value` (#436)

  * Added a `TabularLibrary` (#429)

  * Moved `validate` into `RateCollection` (#433)

  * renamed `StarKillerCxxNetwork` to `AmrexAstroCxxNetwork` (#426)

  * Fixed formatting of tabular rate strings (#421)

  * Fixed normalization of colorbar in `grid_plot` (#424)

  * Added the ability to plot neutrino loss for tabular rates (#422)

  * Added support for partition functions in C++ networks (#404)

  * NSE solver now can return chemical potentials (#415)

  * Add spins to NSE calculations (#405)

  * Separate the Coulomb calculation from the NSE screening, make it
    switchable (#398, #380)

  * Fixed `StarKillerCxxNetwork` output directory (#400)

  * Added Potekhin screening (#385)

  * Added screening to python networks (#383)

  * Numba accelerated screening (#373)

  * Fixed C++ approximate rate screening (#384)

  * Added `RateCollection` `remove_rates()` (#368), and allow
    `remove_rates` to operate on dictionary keys (#375)

  * Added NSE solver (#339, #364, #377)

  * added `find_unimportant_rates` method (#367, #369, #374)

  * added spins to partition functions (#371)

  * Split rate up into several classes, including `ReacLibRate` (#362)
    and `TabularRate` (#360)

  * added partition function support to python networks (#358)

  * fixed definition of inverse rate for symmetric screening (#363)

  * Moved `Nucleus` into `nucdata` (#346)

  * Added screening to `RateCollection` (#344)

  * Added approximate rate support to C++ (#333)

  * C++ networks now hardcode the coefficients in a function to
    compute the rates instead of storing a metadata table that is
    read at runtime (#329)


## 1.7.0

  * the `Rate` class now knows how to make its function string in
    python and C++ (#328)

  * C++ networks now have a `std::vector<std::string>` of rate names
    (#326)

  * support for Fortran networks was removed (#324)

  * numerous optimizations (#263, #264, #321, #331)

  * a `DerivedRate` class was added (#272)

  * approximate rate support was added to python networks (#316, #315,
    #313, #271)

  * energy generation was added to python networks (#301)

  * support for inert nuclei was added (#307)

  * the ability to disable rates in C++ networks was added (#304, #290)

  * methods for finding reverse rates were added (#303)

  * a method to find neutrino loss energy from tabulated rates was
    added (#302)

  * the ability to run without Numba was added (#300)

  * weak rate table units in the header were fixed (#297)

  * python nuclei variable indices now begin with 'j' (#296)

  * python implementations of screening were added (#281)

  * the network chart plot was added (#202)

  * a rate filter for the network plot was added (#187)

  * the `Explorer` class was expanded to have more options (#251)

  * the rate plot now returns a matplotlib `Figure` object (#231)

  * added the ability to modify rate products (#252)

  * allow for the edges between nodes to be curved (#257)

  * add a `RatePair` object that groups forward and reverse rates (#212)

  * updated the ReacLib rate library (#248)

  * added nuclear spin tables (#238)

  * added partition function tables (#241, #204)

  * a `Nucleus` now knows its binding energy (#220)

  * many improvements to C++ output (#246, #214, #185, #183)

  * added a diff method to a `Library` (#194)

  * fixed rate sorting so it is more deterministic (#216)

  * added `forward()` and `backward()` methods to `Library` (#207)

  * added a default `ReacLibLibrary` function (#206)

  * added a `validate()` method for a library to find potentially
    important missing rates (#188, #172)

  * added a method to get the number of rates in a library (#173)

  * add a method to remove a rate from a library (#199)

  * libraries are now sorted when printed, with more information
    shown (#195, #168)

  * added a symmetric screening option (#178)

  * a "rotated" plot type for the network structure was added (#161)

  * versioning is now managed by `setuptools_scm` (#158)

## 1.6.0

  * added support for C++ StarKiller / AMReX Microphysics networks
    (#152, #151, #149)

  * centralized sympy code generation and removed common
    subexpression support from network generation (#145)

  * added an example on integrating a python network

## 1.5.0

  * Added `gridplot` function for plotting nuclides on a grid and
    coloring the cells by X, Y, Xdot, Ydot or activity

  * Created a notebook and a script for generating rp-process
    networks. The script allows an arbitrary endpoint to be
    specified.

  * Added a `filter_function` option to `RateFilter` and to the plotting
    functions that allows for filtering rates or nuclei using a
    user-defined Boolean map.

  * Fixed a `write_network` crash if the `RateCollection` contains
    multiple rates with the same name.

  * Deleted unused BLAS and VODE files previously used for the
    standalone Fortran outputs. These have been deprecated in favor
    of the StarKiller Microphysics network output as of v1.3.0.

  * fixed screening for the 3-alpha rate

## 1.4.1

  * Improvements for locating rate files in directories also
    containing python scripts or Jupyter notebooks

  * Fixed a warning when using Numba

  * Updates to the StarKiller Microphysics format for compile time
    variables

  * Enhancements for the `RateCollection` `Explorer` to use NetworkX 2.5

  * Updated the requirements

## 1.4.0

  * Added general support for tabulated weak rates from Suzuki, et
    al. 2016

  * Added a simple pp network example

  * Updated StarKiller Microphysics output to track latest changes in
    Microphysics network API

  * Added a core developers policy in the Readme

## 1.3.0

  * Replace double precision reals in StarKiller networks with custom
    real type `_rt`

  * Incorporate modifications to thermal neutrino Jacobian
    contributions from Microphysics (#210) in StarKiller networks

  * Simplify rate evaluation code in StarKiller networks

  * Deprecated standalone Fortran output with VODE integration.

  * `BaseFortranNetwork` is now an abstract class to require users to
    either use StarKillerNetwork or provide templates

## 1.2.0

  * Fix tabular rate bug where electron chemical potential
    contributions were not included

  * Update documentation and code comments

  * Add Numba support for Python networks

  * Enable sparse Jacobian (CSR) for StarKiller networks

  * Incorporate CUDA Fortran port for StarKiller networks including
    tabulated rates

  * Optimize rate screening evaluation for StarKiller networks

  * Fix bug to include the electron fraction and correct density
    dependence for Reaclib electron capture rates

  * Allow a nuclide to be both reactant and product

  * Updates an error message for multiple degenerate rates

  * Add example code generating a 160-isotope network

  * Fix table compilation issue for StarKiller networks

  * Moved energy generation subroutine into RHS module for StarKiller
    networks

## 1.1.1

  * pynucastro published on JOSS is archived to Zenodo

## 1.1.0

  * JOSS reviewer changes merged and pynucastro accepted to JOSS
