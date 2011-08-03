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

    @property
    def pos(self):
        """
        return position of break along sequence measured from 5' -> 3'
        """
        return len(self.seq5p)

    @staticmethod
    def from_list(fields):
        b = Breakpoint()
        b.name = fields[0]
        b.seq5p = fields[1]
        b.seq3p = fields[2]
        b.chimera_names = fields[3].split(',')
        return b
        
    def to_list(self):
        fields = [self.name, self.seq5p, self.seq3p]
        fields.append(','.join(self.chimera_names))
        return fields