from pathlib import Path

from numpy import array

from pynucastro.nucdata import (PartitionFunctionCollection,
                                PartitionFunctionTable)

nucdata_dir = Path(__file__).parents[1]
pf_dir = nucdata_dir/'PartitionFunction'

dir_etfsiq_low = pf_dir/'etfsiq_low.txt'
dir_frdm_low = pf_dir/'frdm_low.txt'
dir_etfsiq_high = pf_dir/'etfsiq_high.txt'
dir_frdm_high = pf_dir/'frdm_high.txt'

ANSWER_ETFSIQ_LOW = array([1.000271, 1.002656, 1.009124, 1.035543, 1.076750, 1.128518, 1.187847, 1.252797,
                           1.322103, 1.394926, 1.470693, 1.883077, 2.339548, 2.835353, 3.371056, 3.949365,
                           4.574281, 5.250894, 5.985411, 7.659520, 9.675912, 12.147961, 15.237089, 19.172457])

ANSWER_FRDM_LOW = array([1.000157, 1.001534, 1.005265, 1.020486, 1.044185, 1.073899, 1.107886, 1.145013,
                         1.184544, 1.225988, 1.269010, 1.501456, 1.755494, 2.027655, 2.317420, 2.625339,
                         2.952505, 3.300375, 3.670713, 4.487378, 5.423159, 6.505528, 7.771334, 9.270601])

ANSWER_ETFSIQ_HIGH = array([5.79E+000, 1.07E+001, 2.13E+001, 4.38E+001, 9.23E+001, 1.97E+002, 4.23E+002,
                             9.12E+002, 1.97E+003, 4.25E+003, 2.92E+004, 2.00E+005, 1.36E+006, 9.31E+006,
                             6.34E+007, 4.31E+008, 2.92E+009, 1.97E+010, 1.33E+011, 8.93E+011, 5.98E+012,
                             3.99E+013, 2.65E+014, 1.76E+015, 1.16E+016, 7.66E+016, 5.03E+017, 3.30E+018,
                             2.16E+019, 1.41E+020, 9.21E+020, 6.00E+021, 3.91E+022, 2.54E+023, 1.65E+024,
                             1.07E+025, 6.97E+025, 4.52E+026, 2.94E+027, 1.91E+028, 8.07E+029, 3.42E+031,
                             1.46E+033, 6.23E+034, 2.68E+036, 1.16E+038, 5.03E+039, 6.45E+043])

ANSWER_FRDM_HIGH = array([9.40e+007, 2.81e+009, 4.93e010, 1.95e+012, 8.84e+013, 3.66e+015, 1.44e+017,
                          5.48e+018, 2.04e+020, 7.48e+021, 5.72e+025, 4.07e+029, 2.69e+033, 1.66e+037,
                          9.60e+040, 5.20e+044, 2.65e+048, 1.28e+052, 5.85e+055, 2.55e+059, 1.06e+063,
                          4.27e+066, 1.65e+070, 6.16e+073, 2.23e+077, 7.87e+080, 2.71e+084, 9.15e+087,
                          3.03e+091, 9.86e+094, 3.17e+098, 1.00e+102, 3.14e+105, 9.77e+108, 3.01e+112,
                          9.23e+115, 2.82e+119, 8.56e+122, 2.59e+126, 7.85e+129, 7.18e+136, 6.59e+143,
                          6.11e+150, 5.74e+157, 5.48e+164, 5.35e+171, 5.34e+178, 1.88e+196])

TEMPERATURES_LOW = array([0.01E+9, 0.15E+9, 0.2E+9, 0.3E+9, 0.4E+9, 0.5E+9, 0.6E+9,
                          0.7E+9, 0.8E+9, 0.9E+9, 1.0E+9, 1.5E+9, 2.0E+9, 2.5E+9,
                          3.0E+9, 3.5E+9, 4.0E+9, 4.5E+9, 5.0E+9, 6.0E+9, 7.0E+9,
                          8.0E+9, 9.0E+9, 10.0E+9])

TEMPERATURES_HIGH = array([12.0E+9, 14.0E+9, 16.0E+9, 18.0E+9, 20.0E+9, 22.0E+9, 24.0E+9,
                           26.0E+9, 28.0E+9, 30.0E+9, 35.0E+9, 40.0E+9, 45.0E+9, 50.0E+9,
                           55.0E+9, 60.0E+9, 65.0E+9, 70.0E+9, 75.0E+9, 80.0E+9, 85.0E+9,
                           90.0E+9, 95.0E+9, 100.0E+9, 105.0E+9, 110.0E+9, 115.0E+9, 120.0E+9,
                           125.0E+9, 130.0E+9, 135.0E+9, 140.0E+9, 145.0E+9, 150.0E+9, 155.0E+9,
                           160.0E+9, 165.0E+9, 170.0E+9, 175.0E+9, 180.0E+9, 190.0E+9, 200.0E+9,
                           210.0E+9, 220.0E+9, 230.0E+9, 240.0E+9, 250.0E+9, 275.0E+9])


class TestPartition:

    @classmethod
    def setup_class(cls):
        """ this is run once for each class before any tests """

        cls.pf_table_etfsiq_low = PartitionFunctionTable(dir_etfsiq_low)
        cls.pf_table_frdm_low = PartitionFunctionTable(dir_frdm_low)
        cls.pf_table_etfsiq_high = PartitionFunctionTable(dir_etfsiq_high)
        cls.pf_table_frdm_high = PartitionFunctionTable(dir_frdm_high)

        cls.pf_collection_frdm = PartitionFunctionCollection(use_high_temperatures=True, use_set='frdm')
        cls.pf_collection_etfsiq = PartitionFunctionCollection(use_high_temperatures=True, use_set='etfsiq')

    @classmethod
    def teardown_class(cls):
        """ this is run once for each class after all tests """
        del cls.pf_table_etfsiq_low
        del cls.pf_table_frdm_low
        del cls.pf_table_etfsiq_high
        del cls.pf_table_frdm_high

        del cls.pf_collection_frdm
        del cls.pf_collection_etfsiq

    def setup_method(self):
        """ this is run before each test """

    def teardown_method(self):
        """ this is run after each test """

    def test_pf_table(self):

        co46_pf_etfsiq_low = self.pf_table_etfsiq_low.get_partition_function('co46')
        ne37_pf_frdm_low = self.pf_table_frdm_low.get_partition_function('ne37')
        fe47_pf_etfsiq_high = self.pf_table_etfsiq_high.get_partition_function('fe47')
        po188_pf_frdm_high = self.pf_table_frdm_high.get_partition_function('po188')

        assert all(co46_pf_etfsiq_low.partition_function == ANSWER_ETFSIQ_LOW)
        assert all(co46_pf_etfsiq_low.temperature == TEMPERATURES_LOW)

        assert all(ne37_pf_frdm_low.partition_function == ANSWER_FRDM_LOW)
        assert all(ne37_pf_frdm_low.temperature == TEMPERATURES_LOW)

        assert all(fe47_pf_etfsiq_high.partition_function == ANSWER_ETFSIQ_HIGH)
        assert all(fe47_pf_etfsiq_high.temperature == TEMPERATURES_HIGH)

        assert all(po188_pf_frdm_high.partition_function == ANSWER_FRDM_HIGH)
        assert all(po188_pf_frdm_high.temperature == TEMPERATURES_HIGH)

    def test_pfsum(self):

        ne19_pf_frdm_low = self.pf_table_frdm_low.get_partition_function('ne19')
        ne19_pf_frdm_high = self.pf_table_frdm_high.get_partition_function('ne19')

        co60_pf_etfsiq_low = self.pf_table_etfsiq_low.get_partition_function('co60')
        co60_pf_etfsiq_high = self.pf_table_etfsiq_high.get_partition_function('co60')

        assert self.pf_collection_frdm.get_partition_function('ne19') == ne19_pf_frdm_high + ne19_pf_frdm_low
        assert self.pf_collection_etfsiq.get_partition_function('co60') == co60_pf_etfsiq_low + co60_pf_etfsiq_high
