'''
Created on Jan 30, 2011

@author: mkiyer
'''
import math

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



