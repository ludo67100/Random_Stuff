# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:38:53 2019

@author: ludov
"""

from neo.rawio import WinWcpRawIO as win 

def remove_leak(signal,bl_start,bl_stop):
    
    import numpy as np
    
    leak = np.nanmean(signal[bl_start:bl_stop],axis=0)
    
    return signal-leak

    

file = 'C:/Users/ludov/Documents/Catwalk/Sham_Wt/Baseline2/SHAM/S1.wcp'

reader = win(file)
reader.parse_header() #HEADER NEEDS TO BE PARSED !!!!!!!!!!!!!!!

nb_sweeps = reader.header['nb_segment'][0]
sampling_rate = reader.get_signal_sampling_rate()

from matplotlib import pyplot as plt 
import numpy as np 

for sweep in range(nb_sweeps):    

    plt.figure()
    
    #Get the raw sigs, binary format 
    #i_start, i_stop = index to start, to stop
    raw_sigs = reader.get_analogsignal_chunk(block_index=0, seg_index=sweep, i_start=0, i_stop=-1,
    channel_indexes=[0,1])
    
    #Convert to float64 to rescale and get corret y values for both channels 
    left_sigs = reader.rescale_signal_raw_to_float(raw_sigs, dtype='float64')[:,0].ravel()
    right_sigs = reader.rescale_signal_raw_to_float(raw_sigs, dtype='float64')[:,1].ravel()
    
    #Check both signals for a given sweep have the same size... 
    assert len(left_sigs) == len(right_sigs),'left and right signals are not of equal length'
    
    time_vector = np.arange(0,len(left_sigs),1)*1./sampling_rate

    bl_start = np.where(time_vector>=0.0)[0][0]
    bl_stop =  np.where(time_vector>=0.1)[0][0]

    clean_left_signal = remove_leak(left_sigs,bl_start,bl_stop)
    clean_right_signal = remove_leak(right_sigs,bl_start,bl_stop)*-1

    sinusoid = clean_left_signal+clean_right_signal
    
    plt.plot(clean_left_signal,time_vector,color='skyblue',label='Left')
    plt.plot(clean_right_signal,time_vector,color='darkblue',label='right')
    
    plt.plot(sinusoid,time_vector, color='orange',label='sinusoid')
    
    plt.legend(loc='best')