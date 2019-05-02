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

sig = sig[0:limit]
time = time[0:limit]

#Define the filters
def butter_lowpass(cutoff, fs, order=10):
    nyq = 0.5 * fs #Nyquist frequencie is half the sampling rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=10):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# Filter requirements.
order = 1 #High value = high window for delayed inputs/outputs
fs = sampling_rate # sample rate, Hz
cutoff = 50 # desired cutoff frequency of the filter, Hz

# Filter the data, and plot both the original and filtered signals.
filtered_sig = butter_lowpass_filter(sig, cutoff, fs, order)

#Plot the filtered data
plt.figure()
plt.title('Filtered signal')
plt.xlabel('Time (s)')
plt.plot(time,filtered_sig)

#Get the peaks on the filtered data
threshold = np.std(filtered_sig[0:100])*5
indexes = peakutils.indexes(filtered_sig.ravel(),thres = threshold, min_dist = 100)

#Plot the peaks on the raw data
plt.figure()
plt.title('Spike detection on the raw trace')
plt.xlabel('Time (s)')

plt.plot(time,sig,color='skyblue')
for i in range(len(indexes)):
    plt.scatter(time[indexes[i]],sig[indexes[i]], color='gray')
    
    

        
    

