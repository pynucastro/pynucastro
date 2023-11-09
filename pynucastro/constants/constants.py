import scipy.constants as physical_constants

# atomic unit mass in g
m_u = physical_constants.value("unified atomic mass unit") * physical_constants.kilo
k = physical_constants.value("Boltzmann constant") / physical_constants.erg
h = physical_constants.value("Planck constant") / physical_constants.erg

erg2MeV = physical_constants.erg / (physical_constants.eV * physical_constants.mega)
