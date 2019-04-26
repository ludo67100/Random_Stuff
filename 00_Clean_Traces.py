# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 13:13:29 2019

@author: ludovic.spaeth
"""

import numpy as np
import neo 
import os 
import pandas as pd 
from matplotlib import pyplot as plt 
import sys

grid = np.array([[1,2,3,4],[5,6,7,8]])

ch = 0

window = [0.4,0.6] #in seconds
baseline = [0.4,0.45] #in seconds
stim = [0.50,0.010] # start and duration in seconds

holdings = [-70,-50,0]
holdings = [-70]


savefig = False
savedata = False
savedir = 'U:/01_ANALYSIS/MF_BINDA'

bl_start = 0

url = 'U:/RAW DATA/MF_BINDA/Sorted/18-04-2019/Cell_1a'
name = '180419_1a'

RAW_SIGS = [] #For raw sigs as matrix and avg sigs as arrays

#-----------------------------------FUNCTIONs----------------------------------

#------------------------------------FIRST COLLECT THE DATA--------------------
for holding, idx in zip(holdings,range(len(holdings))):
    
    path = '{}/{}'.format(url,str(holding))
    file_list = sorted(os.listdir(path))
    
    print (file_list)
    
    #Basic info 
    r = neo.io.WinWcpIO('{}/{}'.format(path,file_list[0]))
    
    b=r.read_block()
    
    sampling_rate = float(b.segments[0].analogsignals[0].sampling_rate)
    time_vector = np.ravel(b.segments[0].analogsignals[0].times)
    
    win_begin = np.where(time_vector>=window[0])[0][0]
    win_compute_begin = np.where(time_vector>=0.5)[0][0]
    win_end = np.where(time_vector>=window[1])[0][0]   
    
    signals = np.zeros((grid.size,len(time_vector),len(file_list))) #The matrix for raw traces
    
    
    for file, idx in zip(file_list, range(len(file_list))):
        
        npath = '{}/{}'.format(path,file)
        
        reader = neo.io.WinWcpIO(npath)
        
        block = reader.read_block()
        
        for sweep in range(grid.size):
            
            sig = np.ravel(block.segments[sweep].analogsignals[ch].magnitude)
            
            bl_begin = np.where(time_vector>=baseline[0])[0][0]
            bl_end = np.where(time_vector>=baseline[1])[0][0]
            leak = np.mean(sig[bl_begin:bl_end])
            signal = sig - leak
            
#            trunc_time = time_vector[win_begin:win_end]
#            trunc_sig = signal[win_begin:win_end]
       
            if holding == 0 : #So inhibition
                
                plt.figure()
                plt.axvspan(stim[0],stim[0]+stim[1],color='skyblue',alpha=0.5) #For the stim
                plt.xlim(0.4,0.6)
                plt.plot(time_vector,sig)
                plt.pause(0.1)
                plt.show(block=False)
            
                decide = input('Keep this sweep ? : y/n : ')
                
                if decide=='break':
                    sys.exit()
                
                if decide == 'n':      
                    nans = np.zeros(len(signal))
                    nans[:] = np.nan
                    signals[sweep,:,idx] = nans
                    plt.close()
                    
                else:
                    signals[sweep,:,idx] = signal
                    plt.close()
            else:
                signals[sweep,:,idx] = signal

    #Append traces to matrices for later
    RAW_SIGS.append(signals)  
     
#------------------------------Now the serious shit----------------------------


AVG_SIGS = []
for site in range(grid.size):
    AVG_SIGS.append(np.nanmean(RAW_SIGS[0][site,:,:],axis=1))
plt.figure()
plt.plot(AVG_SIGS[0])







    