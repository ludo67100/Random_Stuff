def fit_noise(noise, bins, polarity='excitation', method='half_gauss',closefig=True,title=None,
                        savefig=False,url=None):
    
    '''
    Computes gaussian, half_gauss or gauss(multi-factor) fit on input distribution
    
    INPUTS
    noise (1D-array) : the distribution, canonically noise array from uncaging experiment
    bins (int) : number of bins for histogram
    polarity (str) : 'excitation' for EPSC based noise, 'inhibitory' for IPSC based noise
    method (str) : 'half_gauss', 'gaussian' or 'gauss' fit function 
    closefig (bool) : if True, will not display the figure
    title (str, optional) : title for the figure
    savefig (bool) : if True, saves the figure at url+name location 
    url (str, optional) : path for save plot 
    
    OUTPUTS 
    root_x (float) : value of half-width at half maximum of distribution fit 
    
    
    '''

    import numpy as np 
    from scipy.interpolate import UnivariateSpline
    from matplotlib import pyplot as plt 
    from scipy.optimize import curve_fit

    def gaussian(x, mu, sigma,A):
        return A*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))
    
    def gauss(x, a, b, c, d):
        return a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d
    
    def half_gauss(x,sigma):
        return 2*sigma/np.pi*np.exp(-x**2*sigma**2/np.pi)
    
    def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
        return gaussian(x,mu1,sigma1,A1)+gaussian(x,mu2,sigma2,A2)        
    
    
    
    if polarity == 'excitation': #Extract negative values only
        reverse_noise = [item*-1 for item in noise if item<=0]
    
    else: #Extract positive values only
        reverse_noise = [item for item in noise if item>=0]
    
    
    plt.figure()
    plt.title('Gaussian fit {}'.format(title))
    n,bins,patches, = plt.hist(reverse_noise, bins=bins,density=0,facecolor='gray',alpha=0.5, label='noise distrubution') #Histogram
    
    x = bins[:-1]
    y = n
    n = len(x)
    
    px = np.linspace(np.min(x),np.max(x),n)
    
    
    if method == 'half_gauss':
        popt, pcov = curve_fit(half_gauss,x,y)
        sigma = popt[0]
        py = half_gauss(px,sigma)
        
        spline = UnivariateSpline(px, py-np.max(py)/2,s=0) #Plane for HWHM
        root_x = float(spline.roots()[0])
        root_y = half_gauss(root_x,sigma)
     
        
    if method == 'gaussian':
        popt, pcov = curve_fit(gaussian,x,y)
        mu = popt[0]        
        sigma = popt[1]
        A = popt[2]
        py = gaussian(px,mu,sigma,A)
        
        spline = UnivariateSpline(px, py-np.max(py)/2,s=0) #Plane for HWHM
        root_x = float(spline.roots()[0])
        root_y = gaussian(root_x,mu,sigma)
        
    if method == 'gauss':
        popt, pcov = curve_fit(gauss,x,y)
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]
        py = gauss(px,a,b,c,d)
        
        spline = UnivariateSpline(px, py-np.max(py)/2,s=0) #Plane for HWHM
        if len(spline.roots()) > 1:
            print (spline.roots())
            root_x = float(spline.roots()[1])
        else:
            print (spline.roots())
            root_x = float(spline.roots())
        root_y = gauss(root_x,a,b,c,d)
        
    if method == 'bimodal':
        popt, pcov = curve_fit(bimodal,x,y)
        a = popt[0]
        b = popt[1]
        c = popt[2]
        d = popt[3]
        py = bimodal(px,*popt)
        
        spline = UnivariateSpline(px, py-np.max(py)/2,s=0) #Plane for HWHM
        if len(spline.roots()) > 1:
            print (spline.roots())
            root_x = float(spline.roots()[1])
        else:
            print (spline.roots())
            root_x = float(spline.roots())
        root_y = bimodal(root_x,*popt)
        

        
    plt.plot(px,py, linewidth=2,color='orange',label='gaussian fit')

    plt.scatter(root_x,root_y,label='HWHM={}'.format(round(root_x,2)),s=30,color='red')
    plt.plot([root_x,root_x],[0,root_y],linestyle='--', color='red')
    plt.plot([0,root_x],[root_y,root_y],linestyle='--', color='red')
    plt.legend(loc='best')
    
    if savefig == True:
        plt.savefig('{}/{}_noise_dist.png'.format(url,title))
    
    if closefig == True:
        plt.close()
        
    
    return root_x 

    
    
    
