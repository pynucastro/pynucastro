# unit tests for rates

from pytest import approx

from pynucastro.neutrino_cooling import sneut5


class TestNeutrinos:

    def test_sneut5(self):

        # these values were obtained from the C++ Microphysics
        # test_neutrino_cooling unit test

        # first state

        rho = 435964.0405
        T = 11659144.01
        abar = 1.230769231
        zbar = 1.076923077

        spair = 0
        splas = 1.200819134e-06
        sphot = 1.80613144e-16
        sbrem = 4.244860688e-07
        sreco = 0

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # second state

        rho = 91028210.15
        T = 39810717.06
        abar = 1.255100538
        zbar = 1.08784409

        spair = 4.550779002e-274
        splas = 2.941922305e-10
        sphot = 4.284345605e-35
        sbrem = 0.0007036893859
        sreco = 0

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # third state

        rho = 549.2802717
        T = 857695898.6
        abar = 1.255100538
        zbar = 1.08784409

        spair = 1.241049006e+10
        splas = 0.1260851805
        sphot = 4400976.955
        sbrem = 25.38526796
        sreco = 0.001846613691

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # fourth state

        rho = 1315341948
        T = 5411695265
        abar = 1.255100538
        zbar = 1.08784409

        spair = 2.190629927e+11
        splas = 2.766747688e+11
        sphot = 5.544103242e+10
        sbrem = -5486545486
        sreco = 0

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # fifth state

        rho = 23946644.66
        T = 5411695265
        abar = 1.526742695
        zbar = 1.209769625

        spair = 4.83775822e+14
        splas = 3827484288
        sphot = 1.827462615e+12
        sbrem = 941456652.3
        sreco = 2831.153079

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # sixth state

        rho = 91028210.15
        T = 39810717.06
        abar = 1.564362233
        zbar = 1.22665501

        spair = 5.623151799e-269
        splas = 7.550402286e-10
        sphot = 4.619933319e-34
        sbrem = 0.0006869215968
        sreco = 0

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)

        # seventh state

        rho = 6299605.249
        T = 2928644565
        abar = 1.735406634
        zbar = 1.30342763

        spair = 5.333457343e+12
        splas = 94616482.83
        sphot = 2.217659785e+10
        sbrem = 22237302.51
        sreco = 262.7413606

        snu, comps = sneut5(rho, T, abar=abar, zbar=zbar, full_output=True)

        assert comps.splas == approx(splas, rel=1.e-7, abs=1.e-100)
        assert comps.spair == approx(spair, rel=1.e-7, abs=1.e-100)
        assert comps.sphot == approx(sphot, rel=1.e-7, abs=1.e-100)
        assert comps.sbrem == approx(sbrem, rel=1.e-7, abs=1.e-100)
        assert comps.sreco == approx(sreco, rel=1.e-7, abs=1.e-100)
        assert snu == approx(spair + splas + sphot + sbrem + sreco)
