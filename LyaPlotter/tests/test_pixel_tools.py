from LyaPlotter.pixel_tools import *

import unittest

class TestGetNeighbours(unittest.TestCase):
    def test_neighbours(self):
        target = {870, 871, 1126, 872, 1061, 1062, 1319, 807, 808, 810, 1063, 1064, 1065, 809, 811, 1321, 552, 553, 1066, 1067, 1068, 1385, 1190, 869, 1125, 615, 616, 617, 1383, 875, 1128, 1129, 874, 1130, 1384, 1131, 873, 1448, 1194, 1449, 1195, 743, 998, 999, 933, 934, 1191, 935, 937, 682, 936, 938, 679, 680, 681, 1000, 1192, 1193, 939, 940, 1512, 1258, 1127, 806, 996, 997, 742, 1254, 488, 745, 1002, 1256, 1003, 1257, 1004, 746, 1255, 744, 1001, 1320, 1322}

        self.assertEqual(pixel_get_all_neighbours(16,1000,4, nest=False), target)

        target = {872, 936, 937, 999, 1001, 1064, 1065, 1128}
        self.assertEqual(pixel_get_all_neighbours(16,1000,3,True, nest=False)[1],target)