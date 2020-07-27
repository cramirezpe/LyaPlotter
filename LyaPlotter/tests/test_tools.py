from LyaPlotter.tools import *

import unittest
import os
from astropy.io import fits
import warnings
from unittest.mock import patch

class TestCatGenerator(unittest.TestCase):
    catalog_data = os.path.dirname(os.path.realpath(__file__)) +  '/data/catalogs'

    master_file = catalog_data + '/master.fits'
    zcat_file   = catalog_data + '/zcat_0.5.fits'
    test_zcat_file = catalog_data + '/zcat_test.fits'

    def setUp(self):
        warnings.simplefilter('ignore', category=RuntimeWarning)

    def tearDown(self):
        if os.path.exists(self.test_zcat_file):
            os.remove(self.test_zcat_file)

    @patch('builtins.print')
    def test_cat_gen(self, mocked_print):
        master_to_qso_cat(in_file= self.master_file, out_file=self.test_zcat_file, min_zcat = 0,nside=16, downsampling = 0.5, downsampling_seed = 2)


        diff = fits.FITSDiff(self.zcat_file, self.test_zcat_file)
        self.assertEqual(diff.identical, True)

        