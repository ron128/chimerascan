'''
Created on Jan 30, 2011

@author: mkiyer
'''
import math
from math import log
from collections import defaultdict

def kl_divergence(arr):
    t = sum(arr)
    if t == 0:
        return 0
    expected = t / float(len(arr))
    kldiv = sum((x/float(t))*log(x/expected) for x in arr
                if x > 0)
    return kldiv

def poisson(m):
    '''
    courtesy (http://telliott99.blogspot.com/2010/02/replot-poisson-example-with-python.html)
    '''
    def f(k):
        e = math.e**(-m)
        f = math.factorial(k)
        g = m**k
        return g*e/f
    return f

def std(a):
    # find the mean
    n = len(a)
    mean = mean(a)
    # find the standard deviation
    std = sum((x - mean)**2 for x in a)
    std = (std / float(n-1))**0.5
    return std

def normmeanCI(p, xbar, sd, n):
    """
    Computes a p x 100 CI for the given arguments
    p    - confidence coefficient, common values are 0.99, 0.95, 0.90
    xbar - sample point estimate of unknown pop. mean.
    sd   - standard deviation
    n    - sample size
    """
    se    = sd / (n ** 0.5)
    alphadiv2 = (1.0- p)/2.0
    z2    = stat.norm. ppf(1-alphadiv2)
    a     = xbar - z2 * se
    b     = xbar + z2 * se
    return (a, b)

def median(a):
    b = sorted(a)
    ind,odd = divmod(len(b),2)
    median = (b[ind] + b[ind+odd]) / 2.0

def mean(a):
    return sum(a)/float(len(a))

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
                    self.CDF[xval][yval][zval] = c 

    def __call__(self, x, y, z):
        return self.CDF[x][y][z] / float(self.n)

if __name__ == '__main__':
    import random
    X = [random.randrange(0, 5) for x in xrange(100)]
    Y = [random.randrange(0, 5) for y in xrange(100)]
    Z = [random.randrange(0, 5) for z in xrange(100)]
    import itertools
    x = EmpiricalCdf3D(itertools.izip(X,Y,Z))
    print x.n    
    print x(4, 4, 4)

    

