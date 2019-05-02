# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 15:47:18 2019

@author: ludovic.spaeth
"""

#First the module we need
import neo 
from scipy.signal import butter, lfilter
from matplotlib import pyplot as plt
import numpy as np 
import peakutils

#File location
file = r'U:\RAW DATA\data for python course\PC_spiking.wcp'

#Create a reader
reader = neo.WinWcpIO(file)
#Call the block
bl = reader.read_block()

#Call the the time vector, the signal and the sampling rate
time = bl.segments[0].analogsignals[0].times
sampling_rate = float(bl.segments[0].analogsignals[0].sampling_rate)

sig=bl.segments[0].analogsignals[0].magnitude
#Remove the leak
sig = sig-np.mean(sig[0:1000])

#For faster computation, let's take only the first 5 seconds of the signal
limit = np.ravel(np.where(time>=5.0))[0]

sig = np.ravel(sig[0:limit]).astype('float')
time = np.ravel(time[0:limit])

#Get the peaks on the filtered data
threshold = np.std(sig[0:10000])*3
indexes = peakutils.indexes(sig,thres = threshold, min_dist =100)

#Plot the peaks on the raw data
plt.figure()
plt.title('Spike detection on the raw trace')
plt.xlabel('Time (s)')

plt.plot(time,sig.ravel(),color='skyblue')
for i in range(len(indexes)):
    plt.scatter(time[indexes[i]],sig[indexes[i]], color='gray')
    
    
peaks_x = peakutils.interpolate(time,sig,ind=indexes)
        
    

