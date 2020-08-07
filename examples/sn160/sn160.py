# This builds a subch network with all the
# Reaclib rates linking the specified nuclei.

import pynucastro as pyna
from pynucastro.networks import StarKillerNetwork

library_file = "20180319default2"
mylibrary = pyna.rates.Library(library_file)

all_nuclei = ['n',
              'h1', 'h2',
              'he3', 'he4',
              'li6', 'li7',
              'be7', 'be9',
              'b8', 'b10', 'b11',
              'c12', 'c13', 'c14',
              'n13', 'n14', 'n15',
              'o14', 'o15', 'o16', 'o17', 'o18',
              'f17', 'f18', 'f19',
              'ne18', 'ne19', 'ne20', 'ne21', 'ne22',
              'na21', 'na22', 'na23',
              'mg23', 'mg24', 'mg25', 'mg26',
              'al25', 'al26', 'al27',
              'si28', 'si29', 'si30', 'si31', 'si32',
              'p29', 'p30', 'p31', 'p32', 'p33',
              's32', 's33', 's34', 's35', 's36',
              'cl33', 'cl34', 'cl35', 'cl36', 'cl37',
              'ar36', 'ar37', 'ar38', 'ar39', 'ar40',
              'k37', 'k38', 'k39', 'k40', 'k41',
              'ca40', 'ca41', 'ca42', 'ca43', 'ca44', 'ca45', 'ca46', 'ca47', 'ca48',
              'sc43', 'sc44', 'sc45', 'sc46', 'sc47', 'sc48', 'sc49',
              'ti44', 'ti45', 'ti46', 'ti47', 'ti48', 'ti49', 'ti50', 'ti51',
              'v46', 'v47', 'v48', 'v49', 'v50', 'v51', 'v52',
              'cr48', 'cr49', 'cr50', 'cr51', 'cr52', 'cr53', 'cr54',
              'mn50', 'mn51', 'mn52', 'mn53', 'mn54', 'mn55',
              'fe52', 'fe53', 'fe54', 'fe55', 'fe56', 'fe57', 'fe58',
              'co53', 'co54', 'co55', 'co56', 'co57', 'co58', 'co59',
              'ni56', 'ni57', 'ni58', 'ni59', 'ni60', 'ni61', 'ni62', 'ni63', 'ni64',
              'cu57', 'cu58', 'cu59', 'cu60', 'cu61', 'cu62', 'cu63', 'cu64', 'cu65',
              'zn59', 'zn60', 'zn61', 'zn62', 'zn63', 'zn64', 'zn65', 'zn66',
              'ga62', 'ga63', 'ga64',
              'ge63', 'ge64']

sn160 = mylibrary.linking_nuclei(all_nuclei, with_reverse=True)

net = StarKillerNetwork(libraries=[sn160])
net.write_network()
