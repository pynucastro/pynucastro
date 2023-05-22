# 2.0.3

  * add functions to numerically evaluate the Jacobian (#467)

  * Fix C++ compilation with partition functions and approximate
    rates and derived rates (#407, #494)

  * The Explorer class now takes all of the options as
    RateCollection.plot() (#496)

  * Rates edges can be highlighted via a user-defined function in the
    RateCollection plot (#486)

  * All of the plot routines in RateCollection now return a
    matplotlib Figure object (#490)

  * add get_rate_by_nuclei and add_rate to Library() (#484)

  * RateCollection.remove_rates() can now act on a single rate and a
    new method, add_rates(), will add a rate directly to the
    collection.

  * Fixed printing of TabularLibrary (#479)

  * Documentation updates (#473)

  * added get_rate_by_name to RateCollection (#468)

  * added tabular rates to PythonNetwork (#451)

  * fixed NSE state for network without protons (#463)

# 2.0.2

  * Fix the installation of the C++ templates via setup.py (#460)

# 2.0.1

  * Fix some setup.py metadata for PyPI (#453)

  * Change how tabular rate tables are stored and read (#450)

  * Support historial ReacLib format (#449)

# 2.0.0

  * Added a method to get a rate by short name: A(x,y)B (#438)

  * Added the option to show rate links hidden by ydot <
    ydot_cutoff_value (#436)

  * Added a TabularLibrary (#429)

  * Moved validate into RateCollection (#433)

  * renamed StarKillerCxxNetwork to AmrexAstroCxxNetwork (#426)

  * Fixed formatting of tabular rate strings (#421)

  * Fixed normalization of colorbar in grid_plot (#424)

  * Added the ability to plot neutrino loss for tabular rates (#422)

  * Added support for partition functions in C++ networks (#404)

  * NSE solver now can return chemical potentials (#415)

  * Add spins to NSE calculatons (#405)

  * Separate the Coulomb calculation from the NSE screening, make it
    switchable (#398, #380)

  * Fixed StarKillerCxxNetwork output directory (#400)

  * Added Potekhin screening (#385)

  * Added screening to python networks (#383)

  * Numba accelerated screening (#373)

  * Fixed C++ approximate rate screening (#384)

  * Added RateCollection remove_rates() (#368), and allow
    remove_rates to operate on dictionary keys (#375)

  * Added NSE solver (#339, #364, #377)

  * added find_unimportant_rates method (#367, #369, #374)

  * added spins to partition functions (#371)

  * Split rate up into several classes, including ReacLibRate (#362)
    and TabularRate (#360)

  * added partition function support to python networks (#358)

  * fixed definition of inverse rate for symmetric screening (#363)

  * Moved Nucleus into nucdata (#346)

  * Added screening to RateCollection (#344)

  * Added approximate rate support to C++ (#333)

  * C++ networks now hardcode the coefficients in a function to
    compute the rates instead of storing a metadata table that is
    read at runtime (#329)


# 1.7.0

  * the Rate class now knows how to make its function string in
    python and C++ (#328)

  * C++ networks now have a std::vector<std::string> of rate names
    (#326)

  * support for Fortran networks was removed (#324)

  * numerous optimizations (#263, #264, #321, #331)

  * a DerivedRate class was added (#272)

  * approximate rate support was added to python networks (#316, #315, #313, #271)

  * energy generation was added to python networks (#301)

  * support for inert nuclei was added (#307)

  * the ability to disable rates in C++ networks was added (#304, #290)

  * methods for finding reverse rates were added (#303)

  * a method to find neutrino loss energy from tabulated rates was added (#302)

  * the ability to run without Numba was added (#300)

  * weak rate table units in the header were fixed (#297)

  * python nuclei variable indicies now begin with 'j' (#296)

  * python implementations of screening were added (#281)

  * the network chart plot was added (#202)

  * a rate filter for the network plot was added (#187)

  * the Explorer class was expanded to have more options (#251)

  * the rate plot now returns a matplotlib Figure object (#231)

  * added the ability to modify rate products (#252)

  * allow for the edges between nodes to be curved (#257)

  * add a RatePair object that groups forward and reverse rates (#212)

  * updated the ReacLib rate library (#248)

  * added nuclear spin tables (#238)

  * added partition function tables (#241, #204)

  * a Nucleus no knows its binding energy (#220)

  * many improvements to C++ output (#246, #214, #185, #183)

  * added a diff method to a Library (#194)

  * fixed rate sorting so it is more deterministic (#216)

  * added forward() and backward() methods to Library (#207)

  * added a default ReacLibLibrary function (#206)

  * added a validate() method for a library to find potentially
    important missing rates (#188, #172)

  * added a method to get the number of rates in a library (#173)

  * add a method to remove a rate from a library (#199)

  * libraries are now sorted when printed, with more information
    shown (#195, #168)

  * added a symmetric screening option (#178)

  * a "rotated" plot type for the network structure was added (#161)

  * versioning is now managed by setuptools_scm (#158)

# 1.6.0

  * added support for C++ StarKiller / AMReX Microphysics networks
    (#152, #151, #149)

  * centralized sympy code generation and removed common
    subexpression support from network generation (#145)

  * added an example on integrating a python network

# 1.5.0

  * Added gridplot function for plotting nuclides on a grid and
    coloring the cells by X, Y, Xdot, Ydot or activity

  * Created a notebook and a script for generating rp-process
    networks. The script allows an arbitrary endpoint to be
    specified.

  * Added a filter_function option to RateFilter and to the plotting
    functions that allows for filtering rates or nuclei using a
    user-defined Boolean map.

  * Fixed a write_network crash if the RateCollection contains
    multiple rates with the same name.

  * Deleted unused BLAS and VODE files previously used for the
    standalone Fortran outputs. These have been deprecated in favor
    of the StarKiller Microphysics network output as of v1.3.0.

  * fixed screening for the 3-alpha rate

# 1.4.1

  * Improvements for locating rate files in directories also
    containing python scripts or Jupyter notebooks

  * Fixed a warning when using Numba

  * Updates to the StarKiller Microphysics format for compile time
    variables

  * Enhancements for the RateCollection Explorer to use NetworkX 2.5

  * Updated the requirements

# 1.4.0

  * Added general support for tabulated weak rates from Suzuki, et
    al. 2016

  * Added a simple pp network example

  * Updated StarKiller Microphysics output to track latest changes in
    Microphysics network API

  * Added a core developers policy in the Readme

# 1.3.0

  * Replace double precision reals in StarKiller networks with custom
    real type `_rt`

  * Incorporate modifications to thermal neutrino Jacobian
    contributions from Microphysics (#210) in StarKiller networks

  * Simplify rate evaluation code in StarKiller networks

  * Deprecated standalone Fortran output with VODE integration.

  * BaseFortranNetwork is now an abstract class to require users to
    either use StarKillerNetwork or provide templates

# 1.2.0

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

# 1.1.1

  * pynucastro published on JOSS is archived to Zenodo

# 1.1.0

  * JOSS reviewer changes merged and pynucastro accepted to JOSS
