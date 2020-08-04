# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:50:16 2020

@author: Ludovic.spaeth
"""

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


def noiseFit(gridOfNoise,bins, plot=True):
    '''
    Computes gaussian fit on noise distribution to extract sigma value for Zscore calculation
    
    INPUTS:
    gridOfNoise (array or nd-array) : the noise distribution
    bins (positive integer) : the number of bins
    plot (bool) : if True, will display an independent plot with the distribution and the fit
    
    OUTPUTS
    root_x (float) : the HWHM value of the fit (used as sigme for Zscore)
    '''
    
    import numpy as np 
    import matplotlib.pyplot as plt
    import seaborn as sn
    from scipy import stats

    #Curve fit for noise 
    noiseToFit  = np.abs(gridOfNoise).ravel()
    
    #Gaussian func
    def gaussian(x, mu, sigma):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))
   
    #Compute noise histogram
    n,bins = np.histogram(noiseToFit, bins,density=1)

    #Determine fit parameters to normal function
    (mu,sigma) = stats.norm.fit(noiseToFit)
    print ("mu={0}, sigma={1}".format(mu, sigma))
    
    #Do fit 
    x_dummy = np.linspace(np.min(noiseToFit),np.max(noiseToFit), 100)
    gaussianFit = gaussian(x_dummy,mu,sigma)
    
    #Determine scaling factor for gaussian fit 
    maxGauss = np.nanmean(gaussianFit)
    maxDist = np.nanmean(n)
    factor = maxGauss/maxDist
    
    root_x = mu+sigma
    
    #Plot if needed
    if plot == True:
        noisefig, noiseplot = plt.subplots(1,2)
        n,bins,patches, = noiseplot[0].hist(noiseToFit, bins,density=1,facecolor='gray',alpha=0.5, label='noise distrubution',zorder=1) #Histogram
        noiseplot[0].plot(x_dummy,gaussianFit/factor,label='Gaussian Fit',zorder=2)
        noiseplot[0].scatter(mu+sigma,gaussian(mu+sigma,mu,sigma)/factor,label='Sigma={}'.format(round(mu+sigma,2)),zorder=3)
        noiseplot[0].legend(loc='best') 
        sn.distplot(noiseToFit, ax=noiseplot[1], fit=stats.norm)
    
    return root_x


if __name__ == '__main__':
    
    import numpy as np
    
    gridOfNoise = np.genfromtxt('E:/000_PAPER/Inhibition/ENR/Amplitude/191219(1)_Amp_Noisemap_OK.csv',delimiter=',')

    number_bin_fit = 20
    
    a = noiseFit(gridOfNoise,bins=number_bin_fit)
    
    
