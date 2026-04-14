import numpy as np
import pytest

import pynucastro as pyna


class TestStarLibLibrary:
    @pytest.fixture(scope="class")
    def sl_median(self, starlib_library):
        return starlib_library

    @pytest.fixture(scope="class")
    def median_eval_1GK(self, sl_median):
        T = 1.0e9
        return np.array([r.eval(T) for r in sl_median.get_rates()])

    @pytest.fixture(scope="class")
    def sl_sampled(self):
        # use arbitrary yet fixed seed
        return pyna.StarLibLibrary(seed=142)

    def test_num_rates(self, sl_median, sl_sampled):
        # number of rates will need to be adjusted once Al26 isomers
        # are included
        assert sl_median.num_rates == 73905
        assert sl_sampled.num_rates == 73905

    def test_sampling_against_median_rates(self, median_eval_1GK, sl_sampled):
        # Some rates have near zero uncertainties. Thus, instead of
        # forcing the temp eval of all rates to significantly differ
        # from their median, we enforce that the temp eval for most
        # rates differs from their median.

        T = 1.0e9
        sampled_evals = np.array([r.eval(T) for r in sl_sampled.get_rates()])

        # note: atol=0 forces numpy to use rtol
        nchanged = np.sum(~np.isclose(median_eval_1GK, sampled_evals,
                                      rtol=1.0e-8, atol=0))

        assert nchanged >= 0.9 * len(median_eval_1GK)

    def test_sampling_gives_same_rates_for_same_seeds(self, sl_sampled):
        T = 1.0e9
        s1_evals = np.array([r.eval(T) for r in sl_sampled.get_rates()])

        sl_sampled.resample(seed=142)  # same seed as initialization
        s2_evals = np.array([r.eval(T) for r in sl_sampled.get_rates()])

        assert np.array_equal(s1_evals, s2_evals)

    def test_sampling_gives_diff_rates_for_diff_seeds(self, sl_sampled):
        T = 1.0e9

        sl_sampled.resample(seed=487)  # fresh seed
        s1_evals = [r.eval(T) for r in sl_sampled.get_rates()]

        sl_sampled.resample(seed=932)  # fresh seed
        s2_evals = [r.eval(T) for r in sl_sampled.get_rates()]

        nchanged = np.sum(~np.isclose(s1_evals, s2_evals, rtol=1.0e-8, atol=0))
        # once again we check that a meaningful fraction of rates have different evals

        assert nchanged >= 0.9 * len(s1_evals)

    def test_unsample(self, median_eval_1GK, sl_sampled):
        T = 1.0e9

        sl_sampled.unsample()
        unsamp_evals = [r.eval(T) for r in sl_sampled.get_rates()]

        assert np.array_equal(unsamp_evals, median_eval_1GK)

    def test_pp_d(self, sl_median):
        # p + p -> d and p + e + p -> d should be distinguishable
        lib = sl_median.linking_nuclei(["p", "d"])
        rr = lib.get_rates()
        assert len(rr) == 2
        assert rr[0] != rr[1]
