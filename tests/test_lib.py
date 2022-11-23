import unittest
import pathlib
import os

import numpy as np
import pandas

import pauNy


DATA = "data"
FILE = "scerevisiae_05x_flye.ctg.fa"
REF = "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
PATH = f"{DATA}/{FILE}"
GZPATH = f"{DATA}/{FILE}.gz"
REFPATH = f"{DATA}/{REF}"
GSIZE = 12000000


class TestAssembly(unittest.TestCase):

    def setUp(self):
        self.asm = pauNy.Assembly(path=PATH)
        self.asm_gzip = pauNy.Assembly(path=GZPATH)
        self.asm_gsize = pauNy.Assembly(path=PATH, gsize=GSIZE)

    def test_init_path_conversion(self):
        self.assertTrue(isinstance(self.asm.path, pathlib.PosixPath))

    def test_init_name(self):
        self.assertEqual(self.asm.name, FILE)

    def test_init_nogsize(self):
        self.assertEqual(self.asm.gsize, 0)

    def test_init_gzip(self):
        self.assertFalse(self.asm.gzipped)

    def test_init_gzip_true(self):
        self.assertTrue(self.asm_gzip.gzipped)

    def test_sequence_lengths(self):
        asm_length = self.asm._get_sequence_lengths()
        self.assertEqual(np.sum(asm_length), 9156495)

    def test_sequence_lengths_gzip(self):
        asm_length = self.asm_gzip._get_sequence_lengths()
        self.assertEqual(np.sum(asm_length), 9156495)

    def test_sequence_length_type(self):
        asm_length = self.asm._get_sequence_lengths()
        self.assertTrue(isinstance(asm_length, np.ndarray))

    def test_calc_nx_nog(self):
        nx_gt = np.array([217632, 217632, 165609, 165609, 153671, 150046, 150046, 143359,
         143359, 139628, 133354, 133354, 129089, 127360, 123995, 123995,
         116616, 111672, 107487, 101126, 99284, 99284, 97680, 96648, 96390,
         94620, 93798, 92750, 92613, 89105, 88740, 88665, 86365, 85123,
         83774, 80865, 80829, 80500, 80076, 78740, 77226, 70549, 69691, 69227,
         68079, 67870, 66963, 66020, 65365, 64957, 61564, 57893, 57127, 56663,
         55813, 54854, 54100, 53660, 52587, 50276, 48681, 47929, 47429, 46490,
         46219, 45195, 45065, 43519, 43093, 42789, 41483, 39846, 39505, 38885,
         38122, 37981, 37556, 37240, 36882, 35991, 35453, 34401, 33231, 30391,
         29844, 28875, 28208, 27490, 27185, 25344, 23752, 22332, 19853, 18936,
         16843, 15489, 14098, 11770, 9709, 587])
        nx = self.asm.calculate_Nx()
        self.assertTrue(np.allclose(nx_gt, nx))


    def test_calc_nx_g(self):
        nx_gt = np.array([217632, 165609, 165609, 153671, 150046, 143359, 139628,
         139628, 133354, 129089, 127360, 123995, 116616, 111672, 107487, 101126,
         97680, 96648, 96390, 94620, 92750, 92613, 89105, 88665, 86365, 85123,
         82591, 80865, 80500, 80076, 77226, 74292, 69691, 68550, 67870, 66963,
         66020, 64957, 61564, 57822, 56833, 55813, 54854, 53660, 52587, 50276,
         47929, 47429, 46269, 46076, 45065, 43381, 42948, 41658, 39846, 39487,
         38225, 37981, 37306, 37090, 35991, 34767, 33732, 30391, 29608, 28295,
         27659, 27104, 25195, 23009, 19853, 18028, 16520, 14098, 10946, 8294,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        nx = self.asm_gsize.calculate_Nx()
        self.assertTrue(np.allclose(nx_gt, nx))


    def test_calc_aun_nog(self):
        aun = self.asm.calculate_auN()
        self.assertAlmostEqual(aun, 72722.85214517127)


    def test_calc_aun_g(self):
        aun = self.asm_gsize.calculate_auN()
        self.assertAlmostEqual(aun, 55490.536004416674)



class TestReference(unittest.TestCase):

    def setUp(self):
        self.ref = pauNy.Reference(path=PATH)


    def test_genome_size(self):
        gsize = self.ref.genome_size()
        self.assertEqual(gsize, 9156495)



class TestAssemblyCollection(unittest.TestCase):

    def setUp(self):
        self.ac = pauNy.AssemblyCollection(paths=[PATH])
        self.ac_ref = pauNy.AssemblyCollection(paths=[PATH], reference_path=REFPATH)
        self.ac_gsize = pauNy.AssemblyCollection(paths=[PATH], genome_size=GSIZE)


    def test_init_conversion(self):
        check01 = all([isinstance(p, pathlib.PosixPath) for p in self.ac._assembly_paths])
        check02 = all([isinstance(p, pathlib.PosixPath) for p in self.ac_ref._assembly_paths])
        self.assertTrue(check01)
        self.assertTrue(check02)


    def test_ref_inits(self):
        self.assertEqual(self.ac_ref.ref_path, pathlib.Path(REFPATH))
        self.assertEqual(self.ac_ref.ref_name, REF)

    def test_ref_path_addition(self):
        refpath_in_asm = self.ac_ref.ref_path in self.ac_ref._assembly_paths
        self.assertTrue(refpath_in_asm)

    def test_ref_size(self):
        self.assertEqual(self.ac_ref.genome_size, 12157105)

    def test_ref_size_input(self):
        self.assertEqual(self.ac_gsize.genome_size, GSIZE)


    def test_initialised_assemblies(self):
        check01 = all([isinstance(a, pauNy.Assembly) for a in self.ac.assemblies])
        check02 = all([isinstance(a, pauNy.Assembly) for a in self.ac_ref.assemblies])
        self.assertTrue(check01)
        self.assertTrue(check02)


    def test_calculate_metrics_ref(self):
        nx, aun = self.ac_ref.calculate_metrics()
        self.assertEqual(len(nx), 2)
        self.assertEqual(len(aun), 2)
        self.assertTrue(isinstance(nx, dict))
        self.assertTrue(isinstance(aun, dict))


    def test_calculate_metrics(self):
        nx, aun = self.ac.calculate_metrics()
        self.assertEqual(len(nx), 1)
        self.assertEqual(len(aun), 1)
        self.assertTrue(isinstance(nx, dict))
        self.assertTrue(isinstance(aun, dict))


    def test_generate_df(self):
        nx, aun = self.ac.calculate_metrics()
        nxf, aunf = self.ac.generate_dataframes(nx, aun)
        self.assertTrue(isinstance(nxf, pandas.core.frame.DataFrame))
        self.assertTrue(isinstance(aunf, pandas.core.frame.DataFrame))
        self.assertEqual(nxf.shape, (100, 4))
        self.assertEqual(aunf.shape, (1, 3))


    def test_generate_df_ref(self):
        nx, aun = self.ac_ref.calculate_metrics()
        nxf, aunf = self.ac_ref.generate_dataframes(nx, aun)
        self.assertTrue(isinstance(nxf, pandas.core.frame.DataFrame))
        self.assertTrue(isinstance(aunf, pandas.core.frame.DataFrame))
        self.assertEqual(nxf.shape, (200, 4))
        self.assertEqual(aunf.shape, (2, 3))


    def test_metric_frames(self):
        # only ref case
        self.ac_ref.metric_dataframes(out_name="testing_frame")
        check_exist01 = os.path.isfile('testing_frame.nx.csv')
        check_exist02 = os.path.isfile('testing_frame.aun.csv')
        check_notempty01 = os.path.getsize('testing_frame.nx.csv')
        check_notempty02 = os.path.getsize('testing_frame.aun.csv')
        self.assertTrue(check_exist01)
        self.assertTrue(check_exist02)
        self.assertTrue(check_notempty01 > 0)
        self.assertTrue(check_notempty02 > 0)


    def test_metric_plots(self):
        self.ac_ref.metric_plots(out_name="testing")
        check_exist = os.path.isfile('testing.nx.pdf')
        check_notempty = os.path.getsize('testing.nx.pdf')
        self.assertTrue(check_exist)
        self.assertTrue(check_notempty > 0)

        check_exist = os.path.isfile('testing.aun.pdf')
        check_notempty = os.path.getsize('testing.aun.pdf')
        self.assertTrue(check_exist)
        self.assertTrue(check_notempty > 0)




if __name__ == '__main__':
    unittest.main()
