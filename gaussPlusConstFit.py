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

#############################################################################
def generateBackground(nbkgd, sigma, xlow, xupper):
    """
    This generates an integer number of flat background events beween the limits
    of xlow and xupper. The number of events is pulled from a gaussian centered
    at nbkgd with a std dev of sigma. The resultant pull is rounded, yielding the
    true number of background events in the simulation, simnbkgd. Note: I rounded
    the output of random.normal() since the 3rd arg of uniform() should be an
    integer. This will likely alter the confidence interval's coverage!
    """
    simnbkgd = np.round(np.random.normal(nbkgd,sigma))
    bkgset = np.random.uniform(xlow, xupper, simnbkgd)
    return bkgset, simnbkgd

#############################################################################
def generateSignal(scale, nsigtot, xupper):
    """
    This pulls nsigtot events from f(x) = (1/scale)*exp[-x/scale]. The total
    number of signal events falling below xupper is returned as nsigtot - nover.
    """
    sigset = np.random.exponential(scale,nsigtot)
    nover = 0
    for evt in sigset: 
        if evt >= xupper: nover += 1
    return sigset, nsigtot - nover

#############################################################################
def binEventsEvenly(eventset, xlower, xupper, nbins):
    """
    This bins eventset into a histogram of nbins of even width within the range
    [xlower, xupper). The histogram is returned as an array; the number of 
    overflow events is also returned. eventset will probably be filled by the 
    first value returned by generateSignal(...) or generateBackground(...).
    """
    # Create histogram filled with zeros, calculate bin width
    hist = np.zeros(nbins)
    xwidth = (float(xupper) - float(xlower))/float(nbins)    

    eventset.sort()
    currbin = 1
    noverflow = 0
    # Loop over sorted events and fill into histogram
    for evt in eventset:
        # If currbin is the overflow bin, automatically fill into overflow
        if currbin == nbins + 1: noverflow += 1
        else: 
            # If event falls below current bin's upper edge, fill current bin
            if evt < currbin*xwidth: hist[currbin-1] += 1
            else: 
                # Increment currbin until event falls into currbin, or once we've
                # reached the overflow bin.
                while evt >= currbin*xwidth and currbin != nbins + 1: 
                    currbin += 1
                if currbin != nbins + 1: hist[currbin-1] += 1
                else: noverflow += 1
    return hist, noverflow

#############################################################################
def signalFractionVector(scale, xlower, xupper, nbins):
    """
    For a given even binning--specified by xlower, xupper, and nbins--this 
    returns an array containing the fraction of an exponential exp[-x/scale] in 
    each bin. This is found by integrating (1/scale)*exp[-x/scale] from x0 to x1.
    This definite integral equals:  exp[-x0/scale] - exp[-x1/scale]
    """
    # Create array of bin fractions (filled with zeros to start)
    fracs = np.zeros(nbins)
    xwidth = (float(xupper) - float(xlower))/float(nbins)

    # Loop over bins
    for ibin in xrange(nbins):
        x0, x1 = xwidth*ibin, xwidth*(ibin+1)
        fracs[ibin] = np.exp(-x0/scale) - np.exp(-x1/scale)
    return fracs

#############################################################################
def backgroundFractionVector(xlower, xupper, nbins):
    """
    For a given even binning--specified by xlower, xupper, and nbins--this 
    returns an array containing the fraction of flat background in each bin.
    This is simply 1./nbins.
    """
    # Create empty array of bin fractions and fill 
    fracs = np.ndarray(nbins)
    fracs.fill(1./nbins)
    return fracs

#############################################################################
def chisquare(datavec, predvec, bkgdrate, bkgdsigma, bkgdnorm, covmat):
    """
    Calculates chi^2 for a fit to exponential + flat background. We're fitting 
    one free parameter for the signal normalization (assume exponential constant
    is known) and one free parameter for the background normalization. Our 
    chi^2 function is defined as:
        chi^2 = [D-P]*M^-1*[D-P] + (bkgdnorm-bkgdrate)^2/(bkgdsigma^2)
    Where: 
       'D' is "datavec"
       'P' is "predvec" (sum of signal and background predicitons)
       'M' is "covmat" (the covariance matrix, which only has Poisson errors)
       'bkgdrate' is the true value for the number of background events
       'bkgdnorm' is the fitted value for the number of background events
       'bkgdsigma' is the 1-sigma uncertainty on bkgdnorm.
    """
    
    # Pass covmat to numpy matrix and invert it.
    covmat = np.matrix(covmat).I
    
    chisq0 = ((datavec-predvec).dot(covmat)).dot((datavec-predvec))
    chisq1 = (bkgdnorm-bkgdrate)*(bkgdnorm-bkgdrate)/(bkgdsigma*bkgdsigma)
    return chisq0 + chisq1

#############################################################################
def chisquareFCN(params, hyperparams):
    """
    This is the function that gets passed to the optimizer. It basically handles 
    inputs to the chisquare(...) function, which actually builds and calculates
    the chi^2. "params" will be filled with the two values we want to minimize 
    over (signal normalization and background rate), and "hyperparams" will hold 
    all of the other necessary chi^2 inputs. All of these are described in 
    chisquare(...) except "Pearson".

    If "Pearson" is true, we will use the Poisson errors of the prediction vector
    to construct the covariance matrix. If False, we'll use the Neyman treatment,
    where the Poisson errors of the data vector will be used to construct the 
    covariance matrix.
    """
    if len(params) != 2 or len(hyperparams) != 6:
        print 'Did not pass the right number of arguments to chisquareFCN(...)'
        print 'Exiting!'
        sys.exit()
    signorm, bkgdnorm = params
    datavec, sigfracvec, bkgdfracvec, bkgdrate, bkgdsigma, Pearson = hyperparams

    predvec = sigfracvec*signorm + bkgdfracvec*bkgdnorm

    if len(datavec) != len(predvec):
        print 'Datavec and predvec have different dimensions. Exiting!'
        sys.exit()
    
    # Create covariance matrix (ndarray filled with zeros), then fill appropriate
    # vector into diagonal. 
    covmat = np.zeros((len(datavec),len(datavec)))
    if Pearson:
        for (i,ent) in enumerate(predvec): covmat[i,i] = ent
    else:
        for (i,ent) in enumerate(datavec): covmat[i,i] = ent

    return chisquare(datavec, predvec, bkgdrate, bkgdsigma, bkgdnorm, covmat)

#############################################################################
def fitExperiment(binning, signal, background)
    xlower, xupper, nbins = binning
    nsigtot, scale, sigfracvec = signal
    nbkgdCV, bkgdsigma, bkgdfracvec = background

    bkgddata, nbkgdtrue = generateBackground(nbkgdCV, bkgdsigma, xlower, xupper)
    sigdata, nsigtrue = generateSignal(scale, nsigtot, xupper)

#############################################################################
def dumbtest():
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

#############################################################################
if __name__ == "__main__":
    dumbtest()
