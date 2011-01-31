'''
Created on Jan 13, 2011

@author: mkiyer
'''
import unittest


from ..join_segmented_pe_sam import get_contiguous_indexes, gen_segment_intervals, find_missing_indexes

class Test(unittest.TestCase):

    def testJoinContiguousIndexes(self):
        
        #print get_contiguous_indexes([0,2,3,5,6])
        d = {(0,1): 'hi',
             (2,): 'missing1',
             (3,4,5): 'hi2',
             (4,5,6,7): 'missing2',
             (8,): 'missing3'}
        
        print find_missing_indexes([2,4,5,6,7,8], d)
        #print list(gen_segment_intervals((0,1,2,3)))        
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()