'''
Created on Jun 11, 2011

@author: mkiyer
'''


class Breakpoint(object):
    
    def __init__(self):
        self.name = None
        self.seq5p = None
        self.seq3p = None
        self.chimera_names = []
        self.chrom5p = None
        self.pos5p = 0
        self.strand5p = 0
        self.chrom3p = None
        self.pos3p = 0
        self.strand3p = 0

    @staticmethod
    def from_list(fields):
        b = Breakpoint()
        b.name = fields[0]
        b.seq5p = fields[1]
        b.seq3p = fields[2]
        b.chimera_names = fields[3].split(',')
        
    def to_list(self):
        fields = [self.name, self.seq5p, self.seq3p]
        fields.append(','.join(self.chimera_names))
        return fields