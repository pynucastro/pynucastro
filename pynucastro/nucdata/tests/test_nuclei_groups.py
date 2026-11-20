from pynucastro.nucdata import Nucleus, parse_nuclei_group


class TestNucleiGroup:

    def test_range(self):

        nucstr = "ni56-58"
        nuc = parse_nuclei_group(nucstr)

        assert len(nuc) == 3
        assert Nucleus("ni56") in nuc
        assert Nucleus("ni57") in nuc
        assert Nucleus("ni58") in nuc

    def test_range2(self):

        nucstr = "ni56-58,60"
        nuc = parse_nuclei_group(nucstr)

        assert len(nuc) == 4
        assert Nucleus("ni56") in nuc
        assert Nucleus("ni57") in nuc
        assert Nucleus("ni58") in nuc
        assert Nucleus("ni60") in nuc

    def test_range3(self):

        nucstr = "ni56-58,60,62-63"
        nuc = parse_nuclei_group(nucstr)

        assert len(nuc) == 6
        assert Nucleus("ni56") in nuc
        assert Nucleus("ni57") in nuc
        assert Nucleus("ni58") in nuc
        assert Nucleus("ni60") in nuc
        assert Nucleus("ni62") in nuc
        assert Nucleus("ni63") in nuc
