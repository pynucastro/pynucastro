# unit tests for rates
import pynucastro as pyna


class TestPythonNetwork2:
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

        reaclib_library = pyna.ReacLibLibrary()
        mylib = reaclib_library.linking_nuclei(["he4", "c12", "o16"])
        self.pynet = pyna.PythonNetwork(libraries=[mylib], inert_nuclei=["ne20"])

    def teardown_method(self):
        """ this is run after each test """
        self.pynet = None

    def test_full_ydot_string(self):

        dyodt = \
"""dYdt[jo16] = (
   -Y[jo16]*lambda_o16__he4_c12
   +rho*Y[jhe4]*Y[jc12]*lambda_he4_c12__o16
   )

"""

        dynedt = \
"""dYdt[jne20] = 0.0

"""

        assert self.pynet.full_ydot_string(pyna.Nucleus("o16")) == dyodt
        assert self.pynet.full_ydot_string(pyna.Nucleus("ne20")) == dynedt
