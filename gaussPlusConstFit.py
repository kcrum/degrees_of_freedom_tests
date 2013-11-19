# Step 1: Pull nbkg events from Norm(nbkgT,nbkgerr) from flat dist. 
# Step 2: Pull nsig events from Pois(nsigT) in hist range from expo dist. 
#         ----> ensure bin contents large so Gaussian approx. holds true
# Step 3: Construct cov. mat.: M = M_stat
# Step 4: Construct chi^2 (1 free param for expo. scale, 1 for flat line):
#         chi^2 = [D-P(scale,b)]*M^-1*[D-P(scale,b)] + (b-nbkgT)^2/nbkgerr
# Step 4: Minimize w.r.t scale, b. Find delta_chi^2 of 1 and 4 for both.
#         (maybe also try in 2-D?)
# Step 5: Verify coverage of scaleT and nbkgT, make appropriate plots.
#

import sys 
import numpy as np
from matplotlib import pyplot as plt

def generateBackground(nbkg, sigma, xlow, xupper):
    simnbkg = np.random.normal(nbkg,sigma)
    bkgset = np.random.uniform(xlow, xupper, simbkg)
    return bkgset, simnbkg

def generateSignal(scale, nsig, xlow, xupper):
    simnsig = np.random.poisson(nsig)
    sigset = np.random.exponential(scale,simnsig)
    return sigset, simnsig

        
if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print 'You must pass a true time constant, number of signal events, true number of background events, and an uncertainty on the number of background events. Exiting.'
        sys.exit(1)

    import scipy.stats as stats

    scaleT = float(sys.argv[1]); nsigT = int(sys.argv[2])
    nbkgT = int(sys.argv[3]); nbkgerr = float(sys.argv[4])
    nbins = 10; xlow=0.; xupper = 10.

    print 'scaleT: %s   nsigT: %s   nbkgT: %s   bkgerr: %s' % \
        (scaleT,nsigT,nbkgT,nbkgerr)

    # Calculate fraction of signal events above xupper
    x = np.linspace(xlow,xupper,1000)
    expcdf = stats.expon.cdf(x,scale=scaleT)
    print 'Fraction of signal events above %s is %s.' % (xupper,1-expcdf[1000-1])
    print 'Number of expected signal events above %s in one experiment is %s.' \
        % (xupper, nsigT*(1-expcdf[1000-1]) )
