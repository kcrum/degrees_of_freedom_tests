import sys
import numpy as np
from matplotlib import pyplot as plt

##########################################################################
def four_bin_chisq_fit(tmean, tsigma, sample, usetruesig = False, 
                       plothist = False):
    """
    This returns the chi-square of a fitted mean and standard deviation from a
    the array of normal draws called 'sample' which is placed in a 4-bin 
    histogram. This histogram  itself is also returned. The bins are set to be of
    equal probability, so the bin cutpoints are at:
    {tmean-0.6745*tsigma, tmean, tmean+0.6745*tsigma}, where 'tmean' and 'tsigma'
    are the sample's true mean and true standard deviation, respectively. 

    If 'usetruesig' is set True, the denominator of the chisquare sum is set to 
    the true poisson variance (sample size/4). Otherwise the poisson sample 
    variance (= bin content) is used. The 'plothist' bool is mainly for
    debugging; if true it outputs the sample and plots it in a histogram.
    """
    sampsize = sample.shape[0] # Number of entries in sample

    hist = np.zeros(shape=4)    
    for x in sample:
#        if plothist: print 'x: %s' % x
        if x <= tmean - 0.6745*tsigma: hist[0] += 1
        elif x <= tmean: hist[1] += 1
        elif x <= tmean + 0.6745*tsigma: hist[2] += 1
        else: hist[3] += 1

    chisq = 0
    for b in hist:
        sigsqrd = 0
        if usetruesig: sigsqrd = sampsize/4.
        else: sigsqrd = b

        if sigsqrd == 0: # Don't divide by zero!
            print 'Warning: one of your bins had zero events. Gaussian approximation for your chi square fit has definitely been violated, and this bin will be ignored'
            continue 
        else: chisq += ((b - (sampsize/4.))**2) / sigsqrd

    if plothist:        
        for x in hist:
            print 'bin content: %s' % x
        print 'chisq: %s' % chisq

        binedges = [tmean - 2*0.6745*tsigma, tmean - 0.6745*tsigma, tmean, 
                    tmean + 0.6745*tsigma]
        plt.bar(binedges, hist, width = 0.6745*tsigma)
        plt.axhline(y=sampsize/4., linewidth=2, color='r')
        plt.show()

    return chisq, hist

##########################################################################
def dist_of_chi_vals(mean, sigma, nsamples, nevtsPerSamp):

    if nevtsPerSamp < 20:
        print 'Warning: your array length is shorter than 20 elements. The gaussian approximation needed for a chi-square fit may fail in this instance!'
        
    chisq = 0
    for i in range(0, nsamples): 
        sample = np.random.normal(mean, sigma, nevtsPerSamp) # This fills an 
        # array with 'nevtsPerSamp' draws from a gaussian.

        chisq = four_bin_chisq_fit(mean, sigma, sample, True)[0]
        sampchisq = four_bin_chisq_fit(mean, sigma, sample, False)[0]

        print 'chisq: %s \t sampchisq: %s' % (chisq, sampchisq)
        

##########################################################################
def truth_vs_sample_sigma(nbins, xmax, nsamples, nevtsPerSamp, mean=0, sigma=1):
    """
    Calculate the chi-squares of 'nsamples' collections of 'nevtsPerSamp' 
    gaussian draws from a distribution with 'mean' and 'sigma'. Plot the 
    histograms of calculated chi-square values, where in one case the chi-square
    is calculated using the true sigma, and in the other case it is calculated 
    with the sample sigma. 
    """

    truesigchis = np.zeros(shape=nsamples)
    sampsigchis = np.zeros(shape=nsamples)

    for i in range(nsamples):
        sample = np.random.normal(mean, sigma, nevtsPerSamp) # This fills an 
        # array with 'nevtsPerSamp' draws from a gaussian.
        truesigchis[i] = four_bin_chisq_fit(mean, sigma, sample, True)[0]
        sampsigchis[i] = four_bin_chisq_fit(mean, sigma, sample, False)[0]
#        if i%100 == 0: print 'truesigchis[%s] = %s' % (i,truesigchis[i])
#        print 'truesigchis[%s] = %s' % (i,truesigchis[i])

    bins = np.linspace(0, xmax, nbins)
    plt.hist(truesigchis, bins, alpha=0.7, hatch="\\") # alpha sets transparency
    plt.hist(sampsigchis, bins, alpha=0.7, hatch="/") # alpha sets transparency
    
    plt.show()

##########################################################################
if __name__ == "__main__":
    mean = 5
    sigma = 1
    
    truth_vs_sample_sigma(10,5,1000,100)

#    sample = np.random.normal(mean, sigma, 24) # This fills an array
#    # with 24 draws from a gaussian.

#    chisq, hist = four_bin_chisq_fit(mean, sigma, sample, False, True)
