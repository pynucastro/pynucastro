import scipy.constants as physical_constants

# atomic unit mass in g and MeV
m_u = physical_constants.value("unified atomic mass unit") * physical_constants.kilo
m_u_MeV = physical_constants.value("atomic mass constant energy equivalent in MeV")

# proton mass in g and MeV
m_p = physical_constants.value("proton mass") * physical_constants.kilo
m_p_MeV = physical_constants.value("proton mass energy equivalent in MeV")

# neutron mass in g and MeV
m_n = physical_constants.value("neutron mass") * physical_constants.kilo
m_n_MeV = physical_constants.value("neutron mass energy equivalent in MeV")

# electron mass in g and MeV
m_e = physical_constants.value("electron mass") * physical_constants.kilo
m_e_MeV = physical_constants.value("electron mass energy equivalent in MeV")

N_A = physical_constants.value("Avogadro constant")

# Boltzmann constant in erg/K and MeV/K
k = physical_constants.value("Boltzmann constant") / physical_constants.erg
k_MeV = physical_constants.value("Boltzmann constant in eV/K") / physical_constants.mega

h = physical_constants.value("Planck constant") / physical_constants.erg
hbar = physical_constants.value("reduced Planck constant") / physical_constants.erg

q_e = physical_constants.value("elementary charge") * (physical_constants.c * 100) / 10  # C to statC (esu)
c_light = physical_constants.value("speed of light in vacuum") / physical_constants.centi

# conversions
erg2MeV = physical_constants.erg / (physical_constants.eV * physical_constants.mega)
MeV2erg = (physical_constants.eV * physical_constants.mega) / physical_constants.erg

# CODATA 2018 values of masses needed for Nubase 2020 consistency
# these should no longer be needed when the next Nubase is released

m_u_MeV_C18 = 931.49410242
m_n_MeV_C18 = 939.56542052

m_u_C18 = 931.49410242 * MeV2erg / c_light**2
m_n_C18 = 939.56542052 * MeV2erg / c_light**2
