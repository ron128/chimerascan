'''
Created on Jan 24, 2011

@author: mkiyer
'''
from collections import defaultdict

class EmpiricalCdf3D(object):
    
    def prob(self, x, y, z):
        if self.n == 0:
            return 0.0
        # find prob(X = x) by summing all y's and z'a
        nx = 0
        ydict = self.D[x]
        for zdict in ydict.itervalues():        
            nz_given_y = sum(zdict.itervalues())
            nx += nz_given_y
        if nx == 0:
            return 0.0
        px = nx / float(self.n)        
        # find prob(Y = y | X = x)
        ny_given_x = sum(self.D[x][y].itervalues())
        if ny_given_x == 0:
            return 0.0
        py_given_x = ny_given_x / float(nx)
        # find prob(Z = z | Y=y, X=x)
        nz_given_xy = self.D[x][y][z]
        if nz_given_xy == 0:
            return 0.0
        pz_given_xy = nz_given_xy / float(ny_given_x) 
        # multiply together
        return pz_given_xy * py_given_x * px

    def _count(self, x, y, z):
        total = 0
        xkeys = sorted(self.D.iterkeys())
        for xval in xkeys:
            if xval > x:
                break
            ykeys = sorted(self.D[xval].iterkeys())
            for yval in ykeys:
                if yval > y:
                    break
                zkeys = sorted(self.D[xval][yval].iterkeys())
                for zval in zkeys:
                    if zval > z:
                        break
                    total += self.D[xval][yval][zval]
        return total

    def __call__(self, x, y, z):
        return self.CDF[x][y][z] / float(self.n)

    def __init__(self, data_iter):
        # use dict as sparse matrix for now
        self.D = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
        self.n = 0
        for x,y,z in data_iter:
            self.n += 1
            self.D[x][y][z] += 1
        # turn into dicts
        for xval, ydict in self.D.iteritems():
            self.D[xval] = dict(ydict)
            for yval, zdict in ydict.iteritems():
                self.D[xval][yval] = dict(zdict)
        self.CDF = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))                
        # turn into cumulative counts
        xkeys = sorted(self.D.iterkeys())
        for xval in xkeys:        
            ykeys = sorted(self.D[xval].iterkeys())
            for yval in ykeys:
                zkeys = sorted(self.D[xval][yval].iterkeys())
                for zval in zkeys:
                    c = self._count(xval, yval, zval)
                    print xval, yval, zval, c                     
                    self.CDF[xval][yval][zval] = c 


if __name__ == '__main__':
    import random
    X = [random.randrange(0, 5) for x in xrange(100)]
    Y = [random.randrange(0, 5) for y in xrange(100)]
    Z = [random.randrange(0, 5) for z in xrange(100)]
    import itertools
    x = EmpiricalCdf3D(itertools.izip(X,Y,Z))
    print x.n    
    print x.cdf(4, 4, 4)

    