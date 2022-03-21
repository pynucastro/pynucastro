# unit tests for rates
import pynucastro.networks as networks
import pynucastro.rates as rates

import io

ANSWER = \
"""validation: ni56 produced in fe52 + he4 --> ni56 never consumed.
validation: ni56 produced in co55 + p --> ni56 never consumed.
validation: missing he4 + he4 + he4 --> n + c11 as alternative to he4 + he4 + he4 --> c12 (Q = -11.4466 MeV).
validation: missing he4 + he4 + he4 --> p + b11 as alternative to he4 + he4 + he4 --> c12 (Q = -8.682 MeV).
validation: missing c12 + he4 --> n + o15 as alternative to c12 + he4 --> o16 (Q = -8.502 MeV).
validation: missing c12 + he4 --> p + n15 as alternative to c12 + he4 --> o16 (Q = -4.966 MeV).
validation: missing n13 + n --> p + c13 as alternative to n13 + n --> n14 (Q = 3.003 MeV).
validation: missing n13 + n --> he4 + b10 as alternative to n13 + n --> n14 (Q = -1.059 MeV).
validation: missing n13 + he4 --> n + f16 as alternative to n13 + he4 --> p + o16 (Q = -10.981 MeV).
validation: missing n14 + n --> n15 as alternative to n14 + n --> p + c14 (Q = 10.8333 MeV).
validation: missing n14 + n --> d + c13 as alternative to n14 + n --> p + c14 (Q = -5.32663 MeV).
validation: missing n14 + n --> he4 + b11 as alternative to n14 + n --> p + c14 (Q = -0.158 MeV).
validation: missing c14 + he4 --> n + o17 as alternative to c14 + he4 --> o18 (Q = -1.81745 MeV).
validation: missing n14 + he4 --> n + f17 as alternative to n14 + he4 --> f18 (Q = -4.735 MeV).
validation: missing n14 + he4 --> p + o17 as alternative to n14 + he4 --> f18 (Q = -1.19182 MeV).
validation: missing c14 --> he4 + be10 as alternative to c14 --> n14 (Q = -12.0119 MeV).
validation: missing c14 --> n + c13 as alternative to c14 --> n14 (Q = -8.1763 MeV).
validation: missing o16 + he4 --> n + ne19 as alternative to o16 + he4 --> ne20 (Q = -12.1348 MeV).
validation: missing o16 + he4 --> p + f19 as alternative to o16 + he4 --> ne20 (Q = -8.114 MeV).
validation: missing f18 + n --> f19 as alternative to f18 + n --> p + o18 (Q = 10.432 MeV).
validation: missing f18 + n --> he4 + n15 as alternative to f18 + n --> p + o18 (Q = 6.41745 MeV).
validation: missing o18 + he4 --> ne22 as alternative to o18 + he4 --> n + ne21 (Q = 9.6681 MeV).
validation: missing o18 + he4 --> p + f21 as alternative to o18 + he4 --> n + ne21 (Q = -5.599 MeV).
validation: missing f18 + he4 --> na22 as alternative to f18 + he4 --> p + ne21 (Q = 8.48 MeV).
validation: missing f18 + he4 --> n + na21 as alternative to f18 + he4 --> p + ne21 (Q = -2.58854 MeV).
validation: missing f18 --> p + o17 as alternative to f18 --> o18 (Q = -5.6065 MeV).
validation: missing f18 --> n + f17 as alternative to f18 --> o18 (Q = -9.15 MeV).
validation: missing ne20 + n --> p + f20 as alternative to ne20 + n --> ne21 (Q = -6.242 MeV).
validation: missing ne20 + n --> he4 + o17 as alternative to ne20 + n --> ne21 (Q = -0.586 MeV).
validation: missing ne21 + he4 --> mg25 as alternative to ne21 + he4 --> n + mg24 (Q = 9.882 MeV).
validation: missing ne21 + he4 --> p + na24 as alternative to ne21 + he4 --> n + mg24 (Q = -2.178 MeV).
validation: missing na23 + he4 --> n + al26 as alternative to na23 + he4 --> al27 (Q = -2.96595 MeV).
validation: missing na23 + he4 --> p + mg26 as alternative to na23 + he4 --> al27 (Q = 1.82067 MeV).
validation: missing mg23 + he4 --> n + si26 as alternative to mg23 + he4 --> si27 (Q = -3.97554 MeV).
validation: missing mg23 + he4 --> p + al26 as alternative to mg23 + he4 --> si27 (Q = 1.87249 MeV).
validation: missing mg23 --> he4 + ne19 as alternative to mg23 --> na23 (Q = -9.65012 MeV).
validation: missing mg23 --> p + na22 as alternative to mg23 --> na23 (Q = -7.5803 MeV).
validation: missing mg23 --> n + mg22 as alternative to mg23 --> na23 (Q = -13.1481 MeV).
validation: missing al27 + he4 --> n + p30 as alternative to al27 + he4 --> p31 (Q = -2.643 MeV).
validation: missing al27 + he4 --> p + si30 as alternative to al27 + he4 --> p31 (Q = 2.37222 MeV).
validation: missing si27 + he4 --> n + s30 as alternative to si27 + he4 --> s31 (Q = -3.96817 MeV).
validation: missing si27 + he4 --> p + p30 as alternative to si27 + he4 --> s31 (Q = 2.95222 MeV).
validation: missing si27 --> p + al26 as alternative to si27 --> al27 (Q = -7.464 MeV).
validation: missing si27 --> n + si26 as alternative to si27 --> al27 (Q = -13.311 MeV).
validation: missing p31 + he4 --> n + cl34 as alternative to p31 + he4 --> cl35 (Q = -5.64751 MeV).
validation: missing p31 + he4 --> p + s34 as alternative to p31 + he4 --> cl35 (Q = 0.626848 MeV).
validation: missing s31 --> p + p30 as alternative to s31 --> p31 (Q = -6.13304 MeV).
validation: missing s31 --> n + s30 as alternative to s31 --> p31 (Q = -13.0534 MeV).
validation: missing s32 + he4 --> n + ar35 as alternative to s32 + he4 --> ar36 (Q = -8.61469 MeV).
validation: missing cl35 + p --> n + ar35 as alternative to cl35 + p --> ar36 (Q = -6.74847 MeV).
validation: missing cl35 + p --> n + ar35 as alternative to cl35 + p --> he4 + s32 (Q = -6.74847 MeV).
validation: missing cl35 + he4 --> n + k38 as alternative to cl35 + he4 --> k39 (Q = -5.85925 MeV).
validation: missing cl35 + he4 --> p + ar38 as alternative to cl35 + he4 --> k39 (Q = 0.836956 MeV).
validation: missing ar36 + he4 --> n + ca39 as alternative to ar36 + he4 --> ca40 (Q = -8.60354 MeV).
validation: missing k39 + p --> n + ca39 as alternative to k39 + p --> ca40 (Q = -7.31496 MeV).
validation: missing k39 + p --> n + ca39 as alternative to k39 + p --> he4 + ar36 (Q = -7.31496 MeV).
validation: missing k39 + he4 --> n + sc42 as alternative to k39 + he4 --> sc43 (Q = -7.33217 MeV).
validation: missing k39 + he4 --> p + ca42 as alternative to k39 + he4 --> sc43 (Q = -0.123994 MeV).
validation: missing ca40 + he4 --> n + ti43 as alternative to ca40 + he4 --> ti44 (Q = -11.1716 MeV).
validation: missing sc43 + p --> n + ti43 as alternative to sc43 + p --> ti44 (Q = -7.64917 MeV).
validation: missing sc43 + p --> n + ti43 as alternative to sc43 + p --> he4 + ca40 (Q = -7.64917 MeV).
validation: missing sc43 + he4 --> n + v46 as alternative to sc43 + he4 --> v47 (Q = -4.76132 MeV).
validation: missing sc43 + he4 --> p + ti46 as alternative to sc43 + he4 --> v47 (Q = 3.07144 MeV).
validation: missing ti44 + he4 --> n + cr47 as alternative to ti44 + he4 --> cr48 (Q = -8.63648 MeV).
validation: missing ti44 + he4 --> n + cr47 as alternative to ti44 + he4 --> p + v47 (Q = -8.63648 MeV).
validation: missing v47 + p --> n + cr47 as alternative to v47 + p --> cr48 (Q = -8.22601 MeV).
validation: missing v47 + he4 --> n + mn50 as alternative to v47 + he4 --> mn51 (Q = -5.02164 MeV).
validation: missing v47 + he4 --> p + cr50 as alternative to v47 + he4 --> mn51 (Q = 3.39339 MeV).
validation: missing cr48 + he4 --> n + fe51 as alternative to cr48 + he4 --> fe52 (Q = -8.24324 MeV).
validation: missing cr48 + he4 --> n + fe51 as alternative to cr48 + he4 --> p + mn51 (Q = -8.24324 MeV).
validation: missing mn51 + p --> n + fe51 as alternative to mn51 + p --> fe52 (Q = -8.80135 MeV).
validation: missing mn51 + he4 --> n + co54 as alternative to mn51 + he4 --> co55 (Q = -5.8782 MeV).
validation: missing mn51 + he4 --> p + fe54 as alternative to mn51 + he4 --> co55 (Q = 3.14706 MeV).
validation: missing fe52 + he4 --> n + ni55 as alternative to fe52 + he4 --> ni56 (Q = -8.64244 MeV).
validation: missing fe52 + he4 --> n + ni55 as alternative to fe52 + he4 --> p + co55 (Q = -8.64244 MeV).
validation: missing co55 + p --> n + ni55 as alternative to co55 + p --> ni56 (Q = -9.47432 MeV).
"""

class TestValidate:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """
        pass

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        pass

    def setup_method(self):
        """ this is run before each test """

        self.reaclib_library = rates.ReacLibLibrary()

        all_reactants = ["n", "p",
                         "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                         "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                         "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                         "c14", "n13", "n14", "o18", "f18", "ne21", "mg23", "na23", "si27", "s31"]

        self.reduced_library = self.reaclib_library.linking_nuclei(all_reactants)

    def teardown_method(self):
        """ this is run after each test """
        self.reaclib_library = None
        self.reduced_library = None

    def test_validate(self):
        output = io.StringIO()
        self.reduced_library.validate(self.reaclib_library, ostream=output)
        assert ANSWER == output.getvalue()
