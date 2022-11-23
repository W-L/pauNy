import unittest
import os

import plotnine

import pauNy
import pauny  # runscript import

PATH_ALL = "data"
REFPATH = "data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"



class TestPlots(unittest.TestCase):

    def setUp(self):
        paths = pauny.parse_input_arg([PATH_ALL])
        ac = pauNy.AssemblyCollection(paths=paths, reference_path=REFPATH)
        nx, aun = ac.calculate_metrics()
        nxf, aunf = ac.generate_dataframes(nx, aun)
        self.nx_frame = nxf
        self.aun_frame = aunf


    def test_nx_plot(self):
        p = pauNy.plot_nx(self.nx_frame)
        self.assertTrue(isinstance(p, plotnine.ggplot))

    def test_aun_plot(self):
        q = pauNy.plot_aun(self.aun_frame)
        self.assertTrue(isinstance(q, plotnine.ggplot))

    def test_save_plot(self):
        p = pauNy.plot_nx(self.nx_frame)
        pauNy.save_plot(p, type="nx", out_name="testing", format="pdf")
        check_exist = os.path.isfile('testing.nx.pdf')
        check_notempty = os.path.getsize('testing.nx.pdf')
        self.assertTrue(check_exist)
        self.assertTrue(check_notempty > 0)



if __name__ == '__main__':
    unittest.main()