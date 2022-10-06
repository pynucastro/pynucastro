from pytest import approx

from pynucastro.nucdata import Nucleus


class TestNucleus:
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
        pass

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
        assert self.ne41.partition_function(0.35e9) == approx(1.0121446711436666)
        assert self.ni61.partition_function(0.35e9) == approx(1.160524742683722)
        assert self.pb237.partition_function(0.35e9) == approx(1.4410114805045504)

    def test_partition_high_temp(self):

        assert not self.p.partition_function
        assert not self.h1.partition_function
        assert self.ne41.partition_function(32.0e9) == approx(4.901052000000001)
        assert self.ni61.partition_function(32.0e9) == approx(1927800.437886083)
        assert self.pb237.partition_function(32.0e9) == approx(5.05620611030359e+28)

    def test_mass(self):

        assert self.p.A_nuc == approx(1.0078250307554963)
        assert self.n.A_nuc == approx(1.0086649179839473)
        assert self.o16.A_nuc == approx(15.994914621587304)
        assert self.c12.A_nuc == 12.0
