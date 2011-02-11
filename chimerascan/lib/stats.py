'''
Created on Jan 30, 2011

@author: mkiyer
'''
import math
from math import log
from collections import defaultdict

def comb(N,k):
    """
    This function was taken from scipy 0.9.0rc1
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    POSSIBILITY OF SUCH DAMAGE.
    
    The number of combinations of N things taken k at a time.
    This is often expressed as "N choose k".

    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.

    Returns
    -------
    val : int, array
        The total number of combinations.

    Notes
    -----
    - Array arguments accepted only for exact=0 case.
    - If k > N, N < 0, or k < 0, then a 0 is returned.

    Examples
    --------
    >>> k = np.array([3, 4])
    >>> n = np.array([10, 10])
    >>> comb(n, k, exact=False)
    array([ 120.,  210.])
    >>> comb(10, 3, exact=True)
    120L
    """
    if (k > N) or (N < 0) or (k < 0):
        return 0L
    val = 1L
    for j in xrange(min(k, N-k)):
        val = (val*(N-j))//(j+1)
    return val

def normal_pdf(x, m, v):
    return 1.0/math.sqrt(2*math.pi*v) * math.exp(-(x-m)**2/(2*v))

def binomial_pdf(p, n, k):
    if n < 100:
        return comb(n, k) * p**k * p**(n-k)  # Fall back to your current method
    return normal_pdf(k, n*p, n*p*(1.0-p))

def binomial_cdf(p, n, k):
    return sum(binomial_pdf(p,n,x) for x in xrange(k+1))

def _interpolate(a, b, fraction):
    """
    This function was taken from scipy 0.9.0rc1
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    POSSIBILITY OF SUCH DAMAGE.

    Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a)*fraction;

def scoreatpercentile(values, p):
    """
    This function was taken from scipy 0.9.0rc1
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
    COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    POSSIBILITY OF SUCH DAMAGE.    
    
    Calculate the score at the given `per` percentile of the sequence `a`.

    For example, the score at per=50 is the median. If the desired quantile
    lies between two data points, we interpolate between them. If the parameter
    `limit` is provided, it should be a tuple (lower, upper) of two values.
    Values of `a` outside this (closed) interval will be ignored.

    Parameters
    ----------
    a : ndarray
        Values from which to extract score.
    per : int or float
        Percentile at which to extract score.
    limit : tuple, optional
        Tuple of two scalars, the lower and upper limits within which to
        compute the percentile.

    Returns
    -------
    score : float
        Score at percentile.

    See Also
    --------
    percentileofscore

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5

    """
    idx = p * (values.shape[0] - 1)
    if (idx % 1 == 0):
        return values[idx]
    else:
        return _interpolate(values[int(idx)], values[int(idx) + 1], idx % 1)

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

    

