class Element:
    """An element---this holds the name, abbreviation, and proton
    number (Z)

    Parameters
    ----------
    abbreviation : str
        The element abbreviation.
    name : str
        The full name of the element.
    Z : int
        The proton number.

    """
    def __init__(self, abbreviation: str, name: str, Z: int) -> None:
        self.abbreviation = abbreviation
        self.name = name
        self.Z = Z


class UnidentifiedElement(Exception):
    pass


class PeriodicTable:
    """The periodic table of elements.

    """

    _table = {'h':  Element('h',  'hydrogen', 1),
              'he': Element('he', 'helium', 2),
              'li': Element('li', 'lithium', 3),
              'be': Element('be', 'beryllium', 4),
              'b':  Element('b',  'boron', 5),
              'c':  Element('c',  'carbon', 6),
              'n':  Element('n',  'nitrogen', 7),
              'o':  Element('o',  'oxygen', 8),
              'f':  Element('f',  'fluorine', 9),
              'ne': Element('ne', 'neon', 10),
              'na': Element('na', 'sodium', 11),
              'mg': Element('mg', 'magnesium', 12),
              'al': Element('al', 'aluminum', 13),
              'si': Element('si', 'silicon', 14),
              'p':  Element('p',  'phosphorus', 15),
              's':  Element('s',  'sulfur', 16),
              'cl': Element('cl', 'chlorine', 17),
              'ar': Element('ar', 'argon', 18),
              'k':  Element('k',  'potassium', 19),
              'ca': Element('ca', 'calcium', 20),
              'sc': Element('sc', 'scandium', 21),
              'ti': Element('ti', 'titanium', 22),
              'v':  Element('v',  'vanadium', 23),
              'cr': Element('cr', 'chromium', 24),
              'mn': Element('mn', 'manganese', 25),
              'fe': Element('fe', 'iron', 26),
              'co': Element('co', 'cobalt', 27),
              'ni': Element('ni', 'nickel', 28),
              'cu': Element('cu', 'copper', 29),
              'zn': Element('zn', 'zinc', 30),
              'ga': Element('ga', 'gallium', 31),
              'ge': Element('ge', 'germanium', 32),
              'as': Element('as', 'arsenic', 33),
              'se': Element('se', 'selenium', 34),
              'br': Element('br', 'bromine', 35),
              'kr': Element('kr', 'krypton', 36),
              'rb': Element('rb', 'rubidium', 37),
              'sr': Element('sr', 'strontium', 38),
              'y':  Element('y',  'yttrium', 39),
              'zr': Element('zr', 'zirconium', 40),
              'nb': Element('nb', 'niobium', 41),
              'mo': Element('mo', 'molybdenum', 42),
              'tc': Element('tc', 'technetium', 43),
              'ru': Element('ru', 'ruthenium', 44),
              'rh': Element('rh', 'rhodium', 45),
              'pd': Element('pd', 'palladium', 46),
              'ag': Element('ag', 'silver', 47),
              'cd': Element('cd', 'cadmium', 48),
              'in': Element('in', 'indium', 49),
              'sn': Element('sn', 'tin', 50),
              'sb': Element('sb', 'antimony', 51),
              'te': Element('te', 'tellurium', 52),
              'i':  Element('i',  'iodine', 53),
              'xe': Element('xe', 'xenon', 54),
              'cs': Element('cs', 'cesium', 55),
              'ba': Element('ba', 'barium', 56),
              'la': Element('la', 'lanthanum', 57),
              'ce': Element('ce', 'cerium', 58),
              'pr': Element('pr', 'praseodymium', 59),
              'nd': Element('nd', 'neodymium', 60),
              'pm': Element('pm', 'promethium', 61),
              'sm': Element('sm', 'samarium', 62),
              'eu': Element('eu', 'europium', 63),
              'gd': Element('gd', 'gadolinium', 64),
              'tb': Element('tb', 'terbium', 65),
              'dy': Element('dy', 'dysprosium', 66),
              'ho': Element('ho', 'holmium', 67),
              'er': Element('er', 'erbium', 68),
              'tm': Element('tm', 'thulium', 69),
              'yb': Element('yb', 'ytterbium', 70),
              'lu': Element('lu', 'lutetium', 71),
              'hf': Element('hf', 'hafnium', 72),
              'ta': Element('ta', 'tantalum', 73),
              'w':  Element('w',  'tungsten', 74),
              're': Element('re', 'rhenium', 75),
              'os': Element('os', 'osmium', 76),
              'ir': Element('ir', 'iridium', 77),
              'pt': Element('pt', 'platinum', 78),
              'au': Element('au', 'gold', 79),
              'hg': Element('hg', 'mercury', 80),
              'tl': Element('tl', 'thallium', 81),
              'pb': Element('pb', 'lead', 82),
              'bi': Element('bi', 'bismuth', 83),
              'po': Element('po', 'polonium', 84),
              'at': Element('at', 'astatine', 85),
              'rn': Element('rn', 'radon', 86),
              'fr': Element('fr', 'francium', 87),
              'ra': Element('ra', 'radium', 88),
              'ac': Element('ac', 'actinium', 89),
              'th': Element('th', 'thorium', 90),
              'pa': Element('pa', 'protactinium', 91),
              'u':  Element('u',  'uranium', 92),
              'np': Element('np', 'neptunium', 93),
              'pu': Element('pu', 'plutonium', 94),
              'am': Element('am', 'americium', 95),
              'cm': Element('cm', 'curium', 96),
              'bk': Element('bk', 'berkelium', 97),
              'cf': Element('cf', 'californium', 98),
              'es': Element('es', 'einsteinium', 99),
              'fm': Element('fm', 'fermium', 100),
              'md': Element('md', 'mendelevium', 101),
              'no': Element('no', 'nobelium', 102),
              'lr': Element('lr', 'lawrencium', 103),
              'rf': Element('rf', 'rutherfordium', 104),
              'db': Element('db', 'dubnium', 105),
              'sg': Element('sg', 'seaborgium', 106),
              'bh': Element('bh', 'bohrium', 107),
              'hs': Element('hs', 'hassium', 108),
              'mt': Element('mt', 'meitnerium', 109),
              'ds': Element('ds', 'darmstadtium', 110),
              'rg': Element('rg', 'roentgenium', 111),
              'cn': Element('cn', 'copernicium', 112),
              'nh': Element('nh', 'nihonium', 113),
              'fl': Element('fl', 'flerovium', 114),
              'mc': Element('mc', 'moscovium', 115),
              'lv': Element('lv', 'livermorium', 116),
              'ts': Element('ts', 'tennessine', 117),
              'og': Element('og', 'oganesson', 118)}

    @classmethod
    def lookup_abbreviation(cls, abbrev: str) -> Element:
        """Given an abbreviation, return the :py:obj:`Element`.

        Parameters
        ----------
        abbrev : str
            The element abbreviation.

        Returns
        -------
        Element

        Raises
        ------
        UnidentifiedElement
            If the element abbreviation is not found.

        """

        try:
            return cls._table[abbrev.lower()]
        except IndexError:
            raise UnidentifiedElement(f'Could not identify element: {abbrev}') from None

    @classmethod
    def lookup_Z(cls, Z: int) -> Element | None:
        """Given the proton number, return the :py:obj:`Element`.

        Parameters
        ----------
        Z : int
            The proton number of the element.

        Returns
        -------
        Element

        """

        for element in cls._table.values():
            if element.Z == Z:
                return element
        return None
