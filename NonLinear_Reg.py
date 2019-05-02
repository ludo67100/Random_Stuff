# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 10:58:47 2019

This code computes linear regression on synaptic charges measured in uncaging experiments
Model (Exponential fit) is segregated between inputs depending on their location 
inspired by : http://apmonitor.com/che263/index.php/Main/PythonRegressionStatistics

@author: ludovic.spaeth
"""


import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import uncertainties as unc
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

#Pip install uncertainties if needed 
try :
    import uncertainties.unumpy as unp
    
except : 
    import pip
    pip.main(['install','uncertainties'])
    import uncertainties.unumpy as unp 
    
conditions = ['WT','CUFF','SHAM','ENR']
conditions = ['WT']



for cond in conditions : 
    
    savedir = 'U:/01_ANALYSIS/Linear Fit on Charges'
    savefig = False
    showfig = True
      
    #Import the data 
    url = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/PCA/00_1_MONTH_ALL.xlsx'
    condition = cond
    data = pd.read_excel(url,sheet_name=condition,index_col=0,header=0) 
    
    _coeff_ = []
    
    band_list = data.columns.values[1:9]
        
    for xband in band_list: 
        X = xband

          
        fig, axes = plt.subplots(3,3,figsize=[10,8])
        
        for yband,ax,idx in zip(data.columns.values[1:9],axes.ravel(),range(len(axes.ravel()))):
            Y = yband
        
            low,up=0.05,0.95 #limits for quantile
            x = data[[X]].values.ravel()

            y = data[[Y]].values.ravel()

            
            n = len(y)
            
            def f(x,a,b,c):
                #The linear or non linear correlation
                #Here the slope = a, the intercept = b
                return a*np.exp(b*x)+c
            
            #use curvefit to return the optimal a & b parameters
            popt, pcov = curve_fit(f,x,y)
            
            #Retrieve variance parameter values
            a = popt[0]
            b = popt[1]
            c = popt[2]
            print ('Optimal values')
            print ('a: ' + str(a))
            print ('b: ' + str(b))
            print ('c: ' + str(c))
            
            #Compute the r2 value
            # close to 1 = strong fit and variance relationship between x and y
            # close to 0 = variance in x does not explain variance in y 
            r2 = 1.0-(sum((y-f(x,a,b,c))**2)/((n-1.0)*np.var(y,ddof=1)))
            print ('R^2 ' + str(r2))
            
            #Calculate parameters confidence interval 
            a,b,c = unc.correlated_values(popt,pcov)
            print ('Uncertainty')
            print ('a :' + str(a))
            print ('b :' + str(b))
            print ('c :' + str(c))            
            
            #Plot the data
        
            ax.scatter(x,y,s=20,color='skyblue')

            
            #Calculate regression confidence interval
            px = np.linspace(0,max(x),100) #-Interval in which you expect the prediction 
            py = a*unp.exp(b*px)+c
            nom = unp.nominal_values(py)
            std = unp.std_devs(py)
            
            def prediction_bands(x,xd,yd,p,func,conf=0.95):
                #x=requested points
                #xd = x data
                #yd = ydata
                #p = parameters
                #func = function name
                
                alpha = 1.0-conf #Significance
                N = xd.size   # data sample size
                var_n = len(p) #Number of parameters
                
                #Quantile of student's distribution for p=(1-alpha/2)
                q = stats.t.ppf(1.0 - alpha / 2.0, N-var_n)
                
                #Std of an individual measurement 
                se = np.sqrt(1.0 / (N-var_n)*np.sum((yd-func(xd,*p))**2))
                
                #Auxiliary definitions
                sx = (x - xd.mean())**2
                sxd = np.sum((xd-xd.mean())**2)
                
                #Best fit model (predicted values)
                yp = func(x,*p)
                
                #Prediction band
                dy = q * se *np.sqrt(1.0+(1.0/N) + (sx/sxd))
                
                #Upper and Lower prediction band
                lower_band, upper_band = yp-dy,yp+dy
                return lower_band, upper_band 
            
            
            
            lpb,upb = prediction_bands(px,x,y,popt,f,conf=0.95)
            
            #plot the regression 
            ax.plot(px,nom,c='lightcoral',label='r2=%s'%round(r2,3),linewidth=2.0)
            
            #plot the uncertaincies lines, 95% confidence region
            ax.plot(px,nom+1.96*std,c='orange')
            ax.plot(px,nom-1.96*std,c='orange')
            
            #Prediction bands, 95% prediction bands
            ax.plot(px,lpb,color='gray',linestyle='--')  
            ax.plot(px,upb,color='gray',linestyle='--')
            
            ax.set_xlabel(str(X)+' charge')
            ax.set_ylabel(str(Y)+' charge')
            
            ax.legend(loc='best')
            #ax.set_title('Linear regression, r=%s, r^2=%s'%(round(np.sqrt(r2),4),round(r2,4)))
            if idx == 1:
                
                ax.set_title('%s linear fit (%s based)'%(condition,X))
           
            plt.tight_layout()
            
            if savefig == True:
                plt.savefig('%s/%s_linear_fit_%s_based.pdf'%(savedir,condition,X))
                    
            _coeff_.append(r2)
            
        if showfig == False:
            plt.close()
      
    #Corr matrix 
    coefficient = np.asarray(_coeff_).reshape((len(band_list),len(band_list)))
    
    plt.figure()
    plt.title('Coeffcients (r2) matrix for %s maps'%condition)
    plt.imshow(coefficient)
    for (j,i),value in np.ndenumerate(coefficient):
        plt.text(i,j,round(value,2),ha='center',va='center',color='white')
    plt.xticks(ticks=np.arange(0,len(band_list),1),labels=band_list,rotation=45)
    plt.yticks(ticks=np.arange(0,len(band_list),1),labels=band_list,rotation=45)
    plt.tight_layout()
    
    if savefig == True:
        plt.savefig('%s/%s_linear_corr_matrix.pdf'%(savedir,condition))
    
    cut_off = 0.55
    
    plt.figure()
    mask = cut_off <= coefficient
    
    plt.title('Coeffcients (r2) matrix for %s maps >= %s'%(condition,cut_off))
    plt.imshow(mask)
    plt.xticks(ticks=np.arange(0,len(band_list),1),labels=band_list,rotation=45)
    plt.yticks(ticks=np.arange(0,len(band_list),1),labels=band_list,rotation=45)
    plt.tight_layout()
    
    if savefig == True:
        plt.savefig('%s/%s_linear_cutoff_corr_matrix.pdf'%(savedir,condition))    


