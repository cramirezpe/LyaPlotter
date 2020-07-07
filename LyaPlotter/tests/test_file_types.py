from LyaPlotter.file_types import *
from LyaPlotter.sims import CoLoReSim, LyaCoLoReSim, QuickQuasarsSim
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

    def test_downsampling_all(self):
        self.simfiles.downsampling = 0.999999
        self.assertNotEqual(self.simfiles.downsampling, 1)
        self.assertAlmostEqual(np.mean(self.simfiles.DEC), -33.71186, 5)

    def test_downsampling_partial(self):
        self.simfiles.downsampling = 0.4
        self.assertAlmostEqual(np.mean(self.simfiles.DEC), -34.545387, 6)
        self.assertAlmostEqual(np.mean(self.simfiles.delta_skewers), 0.009134396)
    # def test_get_multiple_data(self):
    #     z_skewer, RA, delta_skewers = self.simfiles.get_multiple_data(['z_skewer','RA','delta_skewers'])

    #     self.assertAlmostEqual(np.mean(np.mean(z_skewer)), 1.3284745)
    #     self.assertAlmostEqual(np.mean(RA),92.35906,5)
    #     self.assertAlmostEqual(np.mean(np.mean(delta_skewers)),0.01087344)

class TestLyaCoLoReTransmission(unittest.TestCase):
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

    def test_vrad(self):
        self.assertAlmostEqual(np.mean(self.transmission.vrad), 0.9593648)

class TestGaussianCoLoReFiles(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/lyacolore_output_standard'

    def setUp(self):
        self.sim = LyaCoLoReSim(1, self.test_data)
        self.gaussian = self.sim.get_GaussianCoLoRe(pixels=[1])

    def test_z(self):   
        self.assertAlmostEqual(np.mean(self.gaussian.z), 2.3428938)

    def test_id(self):
        self.assertAlmostEqual(np.mean(self.gaussian.id),  52144.301656151416)
    
    def test_wavelength(self):
        self.assertAlmostEqual(np.mean(self.gaussian.wavelength),  2830.6567, 4)

    def test_z_skewer(self):
        self.assertAlmostEqual(np.mean(self.gaussian.z_skewer), 1.3284745)
    
    def test_mask(self):
        self.assertAlmostEqual(np.mean(self.gaussian.mask), 0.791356, 6)
        #self.assertAlmostEqual(np.mean(self.gaussian.mask), 0.46572945)

    def test_delta_skewers(self):
        self.assertAlmostEqual(np.mean(np.mean(self.gaussian.delta_skewers)), 0.001748952)
    
    def test_vrad(self):
        self.assertAlmostEqual(np.mean(np.mean(self.gaussian.vrad)), 9.936127e-05)

class TestPiccaStyleFiles(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/lyacolore_output_standard'

    def setUp(self):
        self.sim = LyaCoLoReSim(1, self.test_data)
        self.piccastyle = self.sim.get_PiccaStyleFiles('picca-gaussian-colorecell', pixels=[1])

    def test_z(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.z), 2.3484433)

    def test_id(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.id),  52129.51813471503)
    
    def test_wavelength(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.wavelength),  4435.072, 3)

    def test_z_skewer(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.z_skewer), 2.6482527)
    
    def test_mask(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.mask),  0.33045572)

    def test_values(self):
        self.assertAlmostEqual(np.mean(np.mean(self.piccastyle.values)),  0.013509239)
    
    def test_RA(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.RA), 47.832645, 6)
    
    def test_DEC(self):
        self.assertAlmostEqual(np.mean(self.piccastyle.DEC), 4.8258276)

class TestSpectra(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/quickquasars_output'

    def setUp(self):
        self.sim = QuickQuasarsSim(1, self.test_data, compression=False)
        self.spectra = self.sim.get_spectra('R',pixels=[1])

    def test_z(self):
        self.assertAlmostEqual(np.mean(self.spectra.z), 2.6764663863182068)
    
    def test_id(self):
        self.assertAlmostEqual(np.mean(self.spectra.id), 34836.9, 6)
    
    def test_RA(self):
        self.assertAlmostEqual(np.mean(self.spectra.RA), 47.757553, 6)
    
    def test_DEC(self):
        self.assertAlmostEqual(np.mean(self.spectra.DEC), 5.8909206)

    def test_wavelength(self):
        self.assertAlmostEqual(np.mean(self.spectra.wavelength),  5649.89990234375)

    def test_z_skewer(self):
        self.assertAlmostEqual(np.mean(self.spectra.z_skewer),  3.6475605241091325)
    
    def test_mask(self):
        self.assertAlmostEqual(np.mean(self.spectra.mask),  0.0, 2)

    def test_flux(self):
        self.assertAlmostEqual(np.mean(self.spectra.flux),  0.80128753)
    

class TestSpectraTruth(unittest.TestCase):
    test_data = os.path.dirname(os.path.realpath(__file__)) + '/data/quickquasars_output'

    def setUp(self):
        self.sim = QuickQuasarsSim(1, self.test_data, compression=False)
        self.spectra = self.sim.get_spectra('R',pixels=[1], redshift_to_use='truth')
    
    def test_z(self):
        self.assertAlmostEqual(np.mean(self.spectra.z), 2.6764665)