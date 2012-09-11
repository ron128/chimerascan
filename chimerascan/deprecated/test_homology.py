'''
Created on Jul 21, 2011

@author: mkiyer
'''
import unittest

from chimerascan.lib.seq import calc_homology

class TestLibraries(unittest.TestCase):

    def testHomology(self):
        a = "AAAAGGGGTTTTCCCC"
        b = "AAAAGGGGTTTTCCCC"
        self.assertEquals(calc_homology(a, b, 0), 16)
        b = "AAAAGGGGTTTTCCCG"
        self.assertEquals(calc_homology(a, b, 0), 15)
        b = "AAATTTGGTTTTCCCC"
        self.assertEquals(calc_homology(a, b, 0), 3)
        self.assertEquals(calc_homology(a, b, 1), 4)
        self.assertEquals(calc_homology(a, b, 2), 5)
        self.assertEquals(calc_homology(a, b, 3), 16)

    



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()