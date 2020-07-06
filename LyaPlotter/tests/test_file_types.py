from LyaPlotter.file_types import *
from LyaPlotter.sims import CoLoReSim, LyaCoLoReSim
import numpy as np
import unittest
import os

class TestCoLoRe(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/colore_output'

    def setUp(self):
        self.sim = CoLoReSim(0,self.test_data)
        self.simfiles = self.sim.get_Sources(ifiles=[0],lr_max=1200.)

    def test_id(self):
        self.assertEqual(self.sim.id_,0)

    def test_qsos(self):
        self.assertAlmostEqual(np.mean(self.simfiles.z),2.0320902)
        self.assertEqual(self.simfiles.N_obj,1000)

    def test_file_path(self):
        self.assertEqual(self.simfiles.file_paths[0],self.test_data+'/out_srcs_s1_0.fits')

    def test_RA(self):
        self.assertAlmostEqual(np.mean(self.simfiles.RA),92.35906,5)
    
    def test_DEC(self):
        self.assertAlmostEqual(np.mean(self.simfiles.DEC),-33.71186,5)

    def test_z_skewer(self):
        self.assertAlmostEqual(np.mean(np.mean(self.simfiles.z_skewer)), 1.3284745)

    def test_wavelength(self):
        self.assertAlmostEqual(np.mean(np.mean(self.simfiles.wavelength)),2830.6567,places=4)

    def test_vrad(self):
        self.assertAlmostEqual(np.mean(np.mean(self.simfiles.vrad)),4.2440406e-05)

    def test_mask(self):
        self.assertAlmostEqual(np.mean(self.simfiles.mask), 0.7251377)

    def test_delta_skewers(self):
        self.assertAlmostEqual(np.mean(np.mean(self.simfiles.delta_skewers)),0.01087344)

class TestLyaCoLoRe(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/lyacolore_output_standard'

    def setUp(self):
        self.sim = LyaCoLoReSim(1, self.test_data)
        self.transmission = self.sim.get_Transmission(pixels=[1,2,3])

    def test_qsos(self):
        self.assertAlmostEqual(np.mean(self.transmission.z),2.4256012)
        self.assertEqual(self.transmission.N_obj,38)

    def test_z_skewer(self):
        self.assertAlmostEqual(np.mean(np.mean(self.transmission.z_skewer)), 3.100537,6)

    def test_wavelength(self):
        self.assertAlmostEqual(np.mean(np.mean(self.transmission.wavelength)),4984.9,places=1)

    def test_mask(self):
        self.assertAlmostEqual(np.mean(self.transmission.mask), 0.21329512)

    def test_lya_absorption(self):
        self.assertAlmostEqual(np.mean(np.mean(self.transmission.lya_absorption)),0.9593648)
    
    def test_lyb_absorption(self):
        self.assertAlmostEqual(np.mean(np.mean(self.transmission.lyb_absorption)),0.9910432)

    def test_delta_lya_absorption(self):
        self.assertAlmostEqual(np.mean(self.transmission.delta_lya_absorption),0.9999999)

    def test_delta_lyb_absorption(self):
        self.assertAlmostEqual(np.mean(self.transmission.delta_lyb_absorption),0.99999976)

    def test_id(self):
        self.assertEqual(self.transmission.id[0],208)
    
    def test_qsos_norsd(self):
        self.assertAlmostEqual(np.mean(self.transmission.z_noRSD),2.4263303)