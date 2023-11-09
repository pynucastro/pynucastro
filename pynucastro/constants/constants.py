import scipy.constants as physical_constants

# atomic unit mass in g
m_u = physical_constants.value("unified atomic mass unit") * physical_constants.kilo
m_u_MeV = physical_constants.value('atomic mass constant energy equivalent in MeV')
k = physical_constants.value("Boltzmann constant") / physical_constants.erg
k_MeV = physical_constants.value("Boltzmann constant in eV/K") / physical_constants.mega
h = physical_constants.value("Planck constant") / physical_constants.erg
hbar = physical_constants.value("reduced Planck constant") / physical_constants.erg
N_A = physical_constants.value("Avogadro constant")
q_e = physical_constants.value("elementary charge") * (physical_constants.c * 100) / 10  # C to statC (esu)

erg2MeV = physical_constants.erg / (physical_constants.eV * physical_constants.mega)
MeV2erg = (physical_constants.eV * physical_constants.mega) / physical_constants.erg
