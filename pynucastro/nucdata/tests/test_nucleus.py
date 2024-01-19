from pytest import approx, raises

from pynucastro.nucdata import Nucleus, get_nuclei_in_range


class TestNucleus:
    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """

    def setup_method(self):
        """ this is run before each test """

        self.p = Nucleus("p")
        self.n = Nucleus("n")
        self.h1 = Nucleus("H1")
        self.d = Nucleus("d")
        self.he4 = Nucleus("He4")
        self.c12 = Nucleus("C12")
        self.o16 = Nucleus("O16")
        self.ni56 = Nucleus("Ni56")
        self.u238 = Nucleus("U238")
        self.he4_also = Nucleus("he4")
        self.ne41 = Nucleus("ne41")
        self.ni61 = Nucleus("ni61")
        self.pb237 = Nucleus("pb237")

    def teardown_method(self):
        """ this is run after each test """

    def test_atomic_weights(self):
        assert self.d.A == 2
        assert self.he4.A == 4
        assert self.c12.A == 12
        assert self.ni56.A == 56
        assert self.u238.A == 238

    def test_atomic_numbers(self):
        assert self.d.Z == 1
        assert self.he4.Z == 2
        assert self.c12.Z == 6
        assert self.ni56.Z == 28
        assert self.u238.Z == 92

    def test_comparisons(self):
        assert self.p == self.h1
        assert self.d != self.he4
        assert self.d < self.he4
        assert self.ni56 > self.c12
        assert self.he4_also == self.he4

    def test_names(self):
        assert self.ni56.c() == "Ni56"
        assert self.ni56.pretty == "{}^{56}\\mathrm{Ni}"

    def test_binding(self):
        assert self.ni56.nucbind == approx(8.642779)

    def test_spin(self):

        assert int(self.p.spin_states) == 2
        assert int(self.he4.spin_states) == 1
        assert int(self.d.spin_states) == 3
        assert int(self.c12.spin_states) == 1
        assert int(self.ni56.spin_states) == 1

    def test_partition_low_temp(self):

        assert not self.p.partition_function
        assert not self.h1.partition_function
        assert self.ne41.partition_function.eval(0.35e9) == approx(1.0121446711436666)
        assert self.ni61.partition_function.eval(0.35e9) == approx(1.160524742683722)
        assert self.pb237.partition_function.eval(0.35e9) == approx(1.4410114805045504)

    def test_partition_high_temp(self):

        assert not self.p.partition_function
        assert not self.h1.partition_function
        assert self.ne41.partition_function.eval(32.0e9) == approx(4.901052000000001)
        assert self.ni61.partition_function.eval(32.0e9) == approx(1927800.437886083)
        assert self.pb237.partition_function.eval(32.0e9) == approx(5.05620611030359e+28)

    def test_mass(self):

        assert self.p.A_nuc == approx(1.0078250307554963)
        assert self.n.A_nuc == approx(1.0086649179839473)
        assert self.o16.A_nuc == approx(15.994914621587304)
        assert self.c12.A_nuc == 12.0

    def test_range(self):

        nuc_list = get_nuclei_in_range(6, 8, 12, 16)

        assert len(nuc_list) == 15
        assert nuc_list[0] == Nucleus("c12")
        assert nuc_list[1] == Nucleus("c13")
        assert nuc_list[2] == Nucleus("c14")
        assert nuc_list[3] == Nucleus("c15")
        assert nuc_list[4] == Nucleus("c16")
        assert nuc_list[5] == Nucleus("n12")
        assert nuc_list[6] == Nucleus("n13")
        assert nuc_list[7] == Nucleus("n14")
        assert nuc_list[8] == Nucleus("n15")
        assert nuc_list[9] == Nucleus("n16")
        assert nuc_list[10] == Nucleus("o12")
        assert nuc_list[11] == Nucleus("o13")
        assert nuc_list[12] == Nucleus("o14")
        assert nuc_list[13] == Nucleus("o15")
        assert nuc_list[14] == Nucleus("o16")

    def test_cast(self):
        assert Nucleus.cast("c12") == self.c12
        assert Nucleus.cast("C12") == self.c12
        assert Nucleus.cast(self.c12) == self.c12
        assert Nucleus.cast("n") == Nucleus("n")

    def test_cast_list(self):
        expected = [self.p, self.n, self.he4, self.c12, self.pb237]
        assert Nucleus.cast_list(["p", self.n, "a", "c12", self.pb237]) == expected
        assert Nucleus.cast_list(["p", "n", "he4", "c12", "pb237"]) == expected
        assert Nucleus.cast_list(expected) == expected
        assert Nucleus.cast_list([self.p, "n", "he4", self.c12, self.pb237], allow_single=True) == expected
        assert Nucleus.cast_list(expected, allow_None=True, allow_single=True) == expected
        assert Nucleus.cast_list(expected, allow_None=True, allow_single=True) == expected

        with raises(ValueError):
            Nucleus.cast_list("he4")
        with raises(ValueError):
            Nucleus.cast_list(self.he4)
        assert Nucleus.cast_list(self.he4, allow_single=True) == [self.he4]
        assert Nucleus.cast_list("he4", allow_single=True) == [self.he4]
        assert Nucleus.cast_list([self.he4], allow_single=True) == [self.he4]
        assert Nucleus.cast_list(["he4"]) == [self.he4]

        with raises(TypeError):
            Nucleus.cast_list(None)
        assert Nucleus.cast_list(None, allow_None=True) is None
        assert Nucleus.cast_list([]) == []
