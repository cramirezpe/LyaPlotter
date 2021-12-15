from LyaPlotter.computations import *
import numpy as np
import unittest
import os

class TestComputations(unittest.TestCase):
    
    def test_overall_mean(self):
        a = [1,2,3,4]
        mask = [True,True,True,False]

        self.assertAlmostEqual(Computations.overall_mean(a,mask), 2)

    def test_overall_sigma(self):
        a = np.asarray([1,2,3,4])
        mask = [False,True,False,True]
        
        self.assertAlmostEqual(Computations.overall_sigma(a,mask),1)

    def test_mean_per_pixel(self):
        values = np.asarray([[1,2,3],[4,5,6],[7,8,9]])

        np.testing.assert_equal(Computations.mean_per_pixel(values),np.asarray([4,5,6]))
    
    def test_mean_per_pixel_masked(self):
        values = np.asarray([[1,2,3],[4,5,6],[7,8,9]])
        mask = np.asarray([[1,0,1],[1,0,1],[1,0,0]])

        np.testing.assert_equal(Computations.mean_per_pixel(values,mask), np.asarray([4,4.5]))

    def test_std_per_pixel_masked(self):
        values = np.asarray([[1,2,3],[4,5,6],[7,8,9]])
        mask = np.asarray([[1,0,1],[1,0,1],[1,0,0]])

        np.testing.assert_equal(Computations.std_per_pixel(values,mask), np.asarray([ 2.4494897427831780982, 1.5]))


    def test_std_per_pixel(self):
        values = np.asarray([[1,2,8],[4,5,6],[7,8,9]])

        np.testing.assert_almost_equal(Computations.std_per_pixel(values), np.asarray([ 2.4494897427831780982, 2.4494897427831780982, 1.24721912892464712853]))