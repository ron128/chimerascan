'''
Created on Jan 13, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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