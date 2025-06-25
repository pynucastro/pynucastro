# unit tests for rates
import pytest

import pynucastro as pyna

ANSWER = \
"""validation: Ni56 produced in Fe52 + He4 ⟶ Ni56 + 𝛾 never consumed.
validation: Ni56 produced in Co55 + p ⟶ Ni56 + 𝛾 never consumed.
validation: missing 3 He4 ⟶ p + B11 as alternative to 3 He4 ⟶ C12 + 𝛾 (Q = -8.682 MeV).
validation: missing 3 He4 ⟶ n + C11 as alternative to 3 He4 ⟶ C12 + 𝛾 (Q = -11.4466 MeV).
validation: missing C12 + He4 ⟶ p + N15 as alternative to C12 + He4 ⟶ O16 + 𝛾 (Q = -4.966 MeV).
validation: missing C12 + He4 ⟶ n + O15 as alternative to C12 + He4 ⟶ O16 + 𝛾 (Q = -8.502 MeV).
validation: missing C14 + He4 ⟶ n + O17 as alternative to C14 + He4 ⟶ O18 + 𝛾 (Q = -1.81745 MeV).
validation: missing C14 ⟶ He4 + Be10 as alternative to C14 ⟶ N14 + e⁻ + 𝜈 (Q = -12.0119 MeV).
validation: missing C14 ⟶ n + C13 as alternative to C14 ⟶ N14 + e⁻ + 𝜈 (Q = -8.1763 MeV).
validation: missing N13 + n ⟶ He4 + B10 as alternative to N13 + n ⟶ N14 + 𝛾 (Q = -1.059 MeV).
validation: missing N13 + n ⟶ p + C13 as alternative to N13 + n ⟶ N14 + 𝛾 (Q = 3.003 MeV).
validation: missing N13 + He4 ⟶ n + F16 as alternative to N13 + He4 ⟶ p + O16 (Q = -10.981 MeV).
validation: missing N14 + n ⟶ He4 + B11 as alternative to N14 + n ⟶ p + C14 (Q = -0.158 MeV).
validation: missing N14 + n ⟶ H2 + C13 as alternative to N14 + n ⟶ p + C14 (Q = -5.32663 MeV).
validation: missing N14 + n ⟶ N15 + 𝛾 as alternative to N14 + n ⟶ p + C14 (Q = 10.8333 MeV).
validation: missing N14 + He4 ⟶ p + O17 as alternative to N14 + He4 ⟶ F18 + 𝛾 (Q = -1.19182 MeV).
validation: missing N14 + He4 ⟶ n + F17 as alternative to N14 + He4 ⟶ F18 + 𝛾 (Q = -4.735 MeV).
validation: missing O16 + He4 ⟶ p + F19 as alternative to O16 + He4 ⟶ Ne20 + 𝛾 (Q = -8.114 MeV).
validation: missing O16 + He4 ⟶ n + Ne19 as alternative to O16 + He4 ⟶ Ne20 + 𝛾 (Q = -12.1348 MeV).
validation: missing O18 + He4 ⟶ p + F21 as alternative to O18 + He4 ⟶ n + Ne21 (Q = -5.599 MeV).
validation: missing O18 + He4 ⟶ Ne22 + 𝛾 as alternative to O18 + He4 ⟶ n + Ne21 (Q = 9.6681 MeV).
validation: missing F18 + n ⟶ He4 + N15 as alternative to F18 + n ⟶ p + O18 (Q = 6.41745 MeV).
validation: missing F18 + n ⟶ F19 + 𝛾 as alternative to F18 + n ⟶ p + O18 (Q = 10.432 MeV).
validation: missing F18 + He4 ⟶ n + Na21 as alternative to F18 + He4 ⟶ p + Ne21 (Q = -2.58854 MeV).
validation: missing F18 + He4 ⟶ Na22 + 𝛾 as alternative to F18 + He4 ⟶ p + Ne21 (Q = 8.48 MeV).
validation: missing F18 ⟶ p + O17 as alternative to F18 ⟶ O18 + e⁺ + 𝜈 (Q = -5.6065 MeV).
validation: missing F18 ⟶ n + F17 as alternative to F18 ⟶ O18 + e⁺ + 𝜈 (Q = -9.15 MeV).
validation: missing Ne20 + n ⟶ He4 + O17 as alternative to Ne20 + n ⟶ Ne21 + 𝛾 (Q = -0.586 MeV).
validation: missing Ne20 + n ⟶ p + F20 as alternative to Ne20 + n ⟶ Ne21 + 𝛾 (Q = -6.242 MeV).
validation: missing Ne21 + He4 ⟶ p + Na24 as alternative to Ne21 + He4 ⟶ n + Mg24 (Q = -2.178 MeV).
validation: missing Ne21 + He4 ⟶ Mg25 + 𝛾 as alternative to Ne21 + He4 ⟶ n + Mg24 (Q = 9.882 MeV).
validation: missing Na23 + He4 ⟶ p + Mg26 as alternative to Na23 + He4 ⟶ Al27 + 𝛾 (Q = 1.82067 MeV).
validation: missing Na23 + He4 ⟶ n + Al26 as alternative to Na23 + He4 ⟶ Al27 + 𝛾 (Q = -2.96595 MeV).
validation: missing Mg23 + He4 ⟶ p + Al26 as alternative to Mg23 + He4 ⟶ Si27 + 𝛾 (Q = 1.87249 MeV).
validation: missing Mg23 + He4 ⟶ n + Si26 as alternative to Mg23 + He4 ⟶ Si27 + 𝛾 (Q = -3.97554 MeV).
validation: missing Mg23 ⟶ He4 + Ne19 as alternative to Mg23 ⟶ Na23 + e⁺ + 𝜈 (Q = -9.65012 MeV).
validation: missing Mg23 ⟶ p + Na22 as alternative to Mg23 ⟶ Na23 + e⁺ + 𝜈 (Q = -7.5803 MeV).
validation: missing Mg23 ⟶ n + Mg22 as alternative to Mg23 ⟶ Na23 + e⁺ + 𝜈 (Q = -13.1481 MeV).
validation: missing Al27 + He4 ⟶ p + Si30 as alternative to Al27 + He4 ⟶ P31 + 𝛾 (Q = 2.37222 MeV).
validation: missing Al27 + He4 ⟶ n + P30 as alternative to Al27 + He4 ⟶ P31 + 𝛾 (Q = -2.643 MeV).
validation: missing Si27 + He4 ⟶ p + P30 as alternative to Si27 + He4 ⟶ S31 + 𝛾 (Q = 2.95222 MeV).
validation: missing Si27 + He4 ⟶ n + S30 as alternative to Si27 + He4 ⟶ S31 + 𝛾 (Q = -3.96817 MeV).
validation: missing Si27 ⟶ p + Al26 as alternative to Si27 ⟶ Al27 + e⁺ + 𝜈 (Q = -7.464 MeV).
validation: missing Si27 ⟶ n + Si26 as alternative to Si27 ⟶ Al27 + e⁺ + 𝜈 (Q = -13.311 MeV).
validation: missing P31 + He4 ⟶ p + S34 as alternative to P31 + He4 ⟶ Cl35 + 𝛾 (Q = 0.626848 MeV).
validation: missing P31 + He4 ⟶ n + Cl34 as alternative to P31 + He4 ⟶ Cl35 + 𝛾 (Q = -5.64751 MeV).
validation: missing S31 ⟶ p + P30 as alternative to S31 ⟶ P31 + e⁺ + 𝜈 (Q = -6.13304 MeV).
validation: missing S31 ⟶ n + S30 as alternative to S31 ⟶ P31 + e⁺ + 𝜈 (Q = -13.0534 MeV).
validation: missing S32 + He4 ⟶ n + Ar35 as alternative to S32 + He4 ⟶ Ar36 + 𝛾 (Q = -8.61469 MeV).
validation: missing Cl35 + p ⟶ n + Ar35 as alternative to Cl35 + p ⟶ He4 + S32 (Q = -6.74847 MeV).
validation: missing Cl35 + p ⟶ n + Ar35 as alternative to Cl35 + p ⟶ Ar36 + 𝛾 (Q = -6.74847 MeV).
validation: missing Cl35 + He4 ⟶ p + Ar38 as alternative to Cl35 + He4 ⟶ K39 + 𝛾 (Q = 0.836956 MeV).
validation: missing Cl35 + He4 ⟶ n + K38 as alternative to Cl35 + He4 ⟶ K39 + 𝛾 (Q = -5.85925 MeV).
validation: missing Ar36 + He4 ⟶ n + Ca39 as alternative to Ar36 + He4 ⟶ Ca40 + 𝛾 (Q = -8.60354 MeV).
validation: missing K39 + p ⟶ n + Ca39 as alternative to K39 + p ⟶ He4 + Ar36 (Q = -7.31496 MeV).
validation: missing K39 + p ⟶ n + Ca39 as alternative to K39 + p ⟶ Ca40 + 𝛾 (Q = -7.31496 MeV).
validation: missing K39 + He4 ⟶ p + Ca42 as alternative to K39 + He4 ⟶ Sc43 + 𝛾 (Q = -0.123994 MeV).
validation: missing K39 + He4 ⟶ n + Sc42 as alternative to K39 + He4 ⟶ Sc43 + 𝛾 (Q = -7.33217 MeV).
validation: missing Ca40 + He4 ⟶ n + Ti43 as alternative to Ca40 + He4 ⟶ Ti44 + 𝛾 (Q = -11.1716 MeV).
validation: missing Sc43 + p ⟶ n + Ti43 as alternative to Sc43 + p ⟶ He4 + Ca40 (Q = -7.64917 MeV).
validation: missing Sc43 + p ⟶ n + Ti43 as alternative to Sc43 + p ⟶ Ti44 + 𝛾 (Q = -7.64917 MeV).
validation: missing Sc43 + He4 ⟶ p + Ti46 as alternative to Sc43 + He4 ⟶ V47 + 𝛾 (Q = 3.07144 MeV).
validation: missing Sc43 + He4 ⟶ n + V46 as alternative to Sc43 + He4 ⟶ V47 + 𝛾 (Q = -4.76132 MeV).
validation: missing Ti44 + He4 ⟶ n + Cr47 as alternative to Ti44 + He4 ⟶ p + V47 (Q = -8.63648 MeV).
validation: missing Ti44 + He4 ⟶ n + Cr47 as alternative to Ti44 + He4 ⟶ Cr48 + 𝛾 (Q = -8.63648 MeV).
validation: missing V47 + p ⟶ n + Cr47 as alternative to V47 + p ⟶ Cr48 + 𝛾 (Q = -8.22601 MeV).
validation: missing V47 + He4 ⟶ p + Cr50 as alternative to V47 + He4 ⟶ Mn51 + 𝛾 (Q = 3.39339 MeV).
validation: missing V47 + He4 ⟶ n + Mn50 as alternative to V47 + He4 ⟶ Mn51 + 𝛾 (Q = -5.02164 MeV).
validation: missing Cr48 + He4 ⟶ n + Fe51 as alternative to Cr48 + He4 ⟶ p + Mn51 (Q = -8.24324 MeV).
validation: missing Cr48 + He4 ⟶ n + Fe51 as alternative to Cr48 + He4 ⟶ Fe52 + 𝛾 (Q = -8.24324 MeV).
validation: missing Mn51 + p ⟶ n + Fe51 as alternative to Mn51 + p ⟶ Fe52 + 𝛾 (Q = -8.80135 MeV).
validation: missing Mn51 + He4 ⟶ p + Fe54 as alternative to Mn51 + He4 ⟶ Co55 + 𝛾 (Q = 3.14706 MeV).
validation: missing Mn51 + He4 ⟶ n + Co54 as alternative to Mn51 + He4 ⟶ Co55 + 𝛾 (Q = -5.8782 MeV).
validation: missing Fe52 + He4 ⟶ n + Ni55 as alternative to Fe52 + He4 ⟶ p + Co55 (Q = -8.64244 MeV).
validation: missing Fe52 + He4 ⟶ n + Ni55 as alternative to Fe52 + He4 ⟶ Ni56 + 𝛾 (Q = -8.64244 MeV).
validation: missing Co55 + p ⟶ n + Ni55 as alternative to Co55 + p ⟶ Ni56 + 𝛾 (Q = -9.47432 MeV).
"""


class TestValidate:

    @pytest.fixture(scope="class")
    def reduced_library(self, reaclib_library):
        all_reactants = ["n", "p",
                         "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                         "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                         "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                         "c14", "n13", "n14", "o18", "f18", "ne21", "mg23", "na23", "si27", "s31"]

        return reaclib_library.linking_nuclei(all_reactants)

    def test_validate(self, reduced_library, reaclib_library, capsys):
        rc = pyna.RateCollection(libraries=[reduced_library])
        rc.validate(reaclib_library)
        captured = capsys.readouterr()
        assert ANSWER == captured.out
