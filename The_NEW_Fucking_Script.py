# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:40:54 2019

@author: ludovic.spaeth
"""

import numpy as np 
from numpy import genfromtxt as gen
from matplotlib import pyplot as plt 
import glob
import os 
import pandas as pd 
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

#------------------------------------FUNCTIONS---------------------------------------------------------
#---------------------------------DO NOT MODIFY--------------------------------------------------------
def MAD(a,axis=None):
    '''
    Computes median absolute deviation of an array along given axis
    '''
    #Median along given axis but keep reduced axis so that result can still broadcast along a 
    
    med = np.nanmedian(a, axis=axis, keepdims=True)
    mad = np.median(np.abs(a-med),axis=axis) #MAD along the given axis
    
    return mad 

def stack_lines(_list):
    
    stacked_map = np.vstack((_list[0],
                            _list[1],
                            _list[2],
                            _list[3],
                            _list[4],
                            _list[5]))
    
    return stacked_map


groups = ['WT','CUFF_1_MONTH','SHAM_1_MONTH','ENR','CUFF_15_DAYS','SHAM_15_DAYS']


method = 'amplitude' #'amplitude' or 'charge'

ylim = 20

bin_for_median = 12 #In % of P1- : 12 c'est pas mal 

if method == 'charge':
    Z_score_limit = 3.09 #For Zscore limit
    vmin, vmax = -8,-0.7 #For maps plot to align every conditions

else:
    Z_score_limit = 2.09 #For Zscore limit
    vmin, vmax = -60,-10 #For maps plot to align every conditions



SAVING = True

for group in range(len(groups)):


    # FILES AND DIRECTORY ------------------------------------------------------------------------------------------------------------
    if method == 'amplitude':
        tag_str = 'Amp'
        #directory = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/AMP_ANALYSIS/00_MAPS/%s'%groups[group]
        directory = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/AMP_ANALYSIS/00_MAPS/%s'%groups[group]
    else:
        tag_str = 'Charge'
        directory = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/CHARGE_ANALYSIS/00_MAPS/%s'%groups[group]

    
    Bands = gen('U:/01_ANALYSIS/Map/AldolaseC/Control_96_sites/BANDS_96_SITES_CORRECTED_80_bins.csv', delimiter=',')
    
    zebrin_file = pd.read_excel(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\Mesures_ZII_HighRes_WT_Cuff_Sham_Enr.xlsx",sheet_name=groups[group],index_col=0,header=1)  #DataFrame containing zebrin band file 

    
    files = glob.glob(r'%s/*_%s_2D_OK.csv'%(directory,tag_str))  #Directory
    names = [os.path.basename(x) for x in files]  # Get maps name 
    
    for x in range(len(names)):
        names[x] = names[x].replace("_%s_2D_OK.csv"%tag_str,"")
        
    manips = names 
    
    if method == 'amplitude':
        savedir = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/AMP_ANALYSIS/00_MAPS/01_MEDIAN'
    else:
        savedir = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/CHARGE_ANALYSIS/00_MAPS/01_MEDIAN'
        
    
    condition = '%s'%groups[group]
    
    # LOADING MATRIX -----------------------------------------------------------------------------------------------------------------
    H = 6           #Map height in sites
    L = 64          #Map width in sites
    N = len(manips)
    
    all_maps, maps = plt.subplots(N,1,sharex=True,sharey=True,figsize=(8,9))
    all_maps.suptitle('{} zscores from {}'.format(groups[group],tag_str))
    
    _mat = np.zeros((H,L,N,3)) #[0] for map, [1] for position and [2] for Zscore
    
    #Positions_1D, Zscore_1D = [],[]
    
    for i in range(N):
        
        _mat[:,:,i,0]=gen(r'%s/%s_%s_2D_OK.csv'%(directory,manips[i],tag_str),delimiter=',')
        
        pos = gen(r'%s/%s_Positions_cp_centered_OK.csv'%(directory,manips[i]),delimiter=',')
        
        #Positions_1D.append(pos) #The 1D position array
        
        pos_2D = (pos,pos,pos,pos,pos,pos)
        
        _mat[:,:,i,1]=np.reshape(pos_2D,(6,64))
            
        _mat[:,:,i,2]=np.genfromtxt(r'%s/%s_%s_zscore_2D_OK.csv'%(directory,manips[i],tag_str),delimiter=',') 
        
        zscore_1D_for_plot = np.max(_mat[:,:,i,2], axis=0)
        
        
        maps[i].plot(pos,zscore_1D_for_plot,label='{}'.format(manips[i]))
        maps[i].legend(loc='best')
        
        
        #Zscore_1D.append(np.abs(np.max(_mat[:,:,-1,2],axis=0)))
        
#FOR 1D ANALYSIS----------------------------------------------------------------------------------------
#FOR 1D ANALYSIS----------------------------------------------------------------------------------------
#FOR 1D ANALYSIS-----------------------------------------------------------------------------------------
        #Create basis for concatenation at first loop
        if i == 0 :
            POSITIONS_1D = pos
            ZSCORES_1D = np.max(_mat[:,:,i,2], axis=0)
            AMPS_1D = np.min(_mat[:,:,i,0], axis=0)
        
        #Concatenate patterns for next loops
        else :
            POSITIONS_1D = np.concatenate((POSITIONS_1D,pos),axis=0)
            ZSCORES_1D = np.concatenate((ZSCORES_1D,np.max(_mat[:,:,i,2], axis=0)),axis=0)
            AMPS_1D = np.concatenate((AMPS_1D,np.min(_mat[:,:,i,0], axis=0)), axis=0)
            
    if SAVING == True:
        plt.savefig('{}/{}_{}_all_1D_zscores.pdf'.format(savedir,groups[group],method))
            
    #SORT AMPLS AND ZSCORE ACCORDING TO POSITIONS
    SORTED_1D_AMPS = [x for _, x in sorted(zip(POSITIONS_1D,AMPS_1D))]
    SORTED_1D_ZSCORES = [x for _, x in sorted(zip(POSITIONS_1D,ZSCORES_1D))]
            
    #Then plot
    fig, ax = plt.subplots(2,1, figsize=(7,8))
    plt.suptitle('%s raw sorting (max values per maps)'%groups[group])
    
    ax[0].plot(sorted(POSITIONS_1D),SORTED_1D_AMPS, label='%s'%tag_str)
    ax[0].set_ylabel('%s'%tag_str)
    
    ax[1].plot(sorted(POSITIONS_1D),np.ones(len(sorted(POSITIONS_1D)))*2.0, color='gray',linestyle='--')
    ax[1].plot(sorted(POSITIONS_1D),np.ones(len(sorted(POSITIONS_1D)))*3.0, color='black',linestyle='--')

    ax[1].plot(sorted(POSITIONS_1D),SORTED_1D_ZSCORES, label='Max Zscores')
    ax[1].set_ylabel('Zscore')
    ax[1].set_xlabel('Position (P1- norm)')
    
    raw_zebrin_values = zebrin_file.loc['MEAN (normalized)'].values[:8]
    raw_zebrin_stds = zebrin_file.loc['STD (normalized)'].values[:8]

    zebrin_color='green'
    
    #P2+ contra zebrin band
    ax[0].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)
    
    #P1+ zebrin band
    ax[0].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)
    
    #P2+ ipsi zebrin band
    ax[0].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)
    ax[1].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)
    
    #PC_cluster
    ax[0].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymin=0.9, alpha=0.5,color='red')
    ax[1].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymax=0.1, alpha=0.5,color='red')
   

    if SAVING == True:
        plt.savefig('{}/{}_{}_1D_raw_sorting.pdf'.format(savedir,groups[group],method))
    
    #BINNING FOR MEDIAN CALCUL
    step = bin_for_median #In % of P1- : 12 c'est pas mal 
    binning = np.arange(-300,300+step,step)
    
    _MEDS, _MADS, _POS, _COUNTS = [],[],[],[]
    
    for i in range(len(binning)):
        
        if i == len(binning)-1:
            break
        
        start, stop = binning[i],binning[i+1]
        _meds, _mads, _pos, _count = [],[],[],[]
        
        #print ('Bin %s to %s'%(start, stop))
        
        SORTED_POSITIONS = sorted(POSITIONS_1D)
        
        for j in range(len(SORTED_POSITIONS)):
            if start < SORTED_POSITIONS[j] <= stop:
                if np.isnan(SORTED_1D_ZSCORES[j])==False:
                    _meds.append(SORTED_1D_ZSCORES[j])
                    _pos.append(SORTED_POSITIONS[j])
                
        _MEDS.append(np.nanmedian(_meds))
        _COUNTS.append(np.count_nonzero(_meds))
        _POS.append(np.nanmedian(_pos))
        _MADS.append(MAD(_meds, axis=0))
        
#        print('--------Median values--------')
#        print(_meds)
#        print('------------MAD--------------')
#        print (MAD(_meds))
        
    #Then plot
    fig, ax = plt.subplots(2,1, sharex=True, figsize=(7,8))
    
    for h in range(len(_POS)):
        if _COUNTS[h] >= 5 : 
            color= 'orange'
        else:
            color='black'
        
        ax[0].bar(_POS[h],_COUNTS[h], color=color, width=10,alpha=0.8)  
    ax[0].set_xlabel('Distance (P1- norm)')
    ax[0].set_title('Count for median')
    ax[0].set_ylabel('Count')


    ax[1].set_title('%s Median Zscore (bin_size=%s)'%(groups[group], step))
    ax[1].set_ylim(0,ylim)
    ax[1].set_ylabel('Zscore')
  
    ax[1].fill_between(_POS,0,np.asarray(_MEDS)+np.abs(np.asarray(_MADS)),color='0.8')  
    ax[1].fill_between(_POS,0,_MEDS,color='0.3')  
    lim = np.ones(len(_POS))*Z_score_limit
    ax[1].fill_between(_POS,lim,_MEDS,where=_MEDS>=lim,color='crimson',interpolate=True)

    #P2+ contra zebrin band
    ax[1].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)
    
    #P1+ zebrin band
    ax[1].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)
   
    #P2+ ipsi zebrin band
    ax[1].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)
   
    #PC_cluster
    ax[1].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymax=0.13,ymin=0.1,color='white')
    

    if SAVING == True:
        plt.savefig('{}/{}_1D_median_{}_Zscore.pdf'.format(savedir,groups[group],method)) 

#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------
#FOR 2D ANALYSIS--------------same shit but line by line first----------------------------------------------        
#FOR 2D ANALYSIS--------------------------------------------------------------------------------------------                
  
    fig, ax = plt.subplots(H,1, figsize=(7,9),sharex=True,sharey=True)
    plt.suptitle('%s 2D maps line by line'%groups[group])
    
    _MEDIAN_ZSCORE_2D, _AVERAGE_AMP_2D, _COUNT_2D, _POSITIONS_2D, _SUM_2D = [], [],[],[],[]
    
    for j in range(H):
        for i in range(N):
            
            #Create basis for concatenation at first loop
            if i == 0 :
                POSITIONS_2D = _mat[j,:,i,1]
                ZSCORES_2D = _mat[j,:,i,2]
                AMPS_2D = _mat[j,:,i,0]
            
            #Concatenate patterns for next loops
            else :
                POSITIONS_2D = np.concatenate((POSITIONS_2D,_mat[j,:,i,1]),axis=0)
                ZSCORES_2D = np.concatenate((ZSCORES_2D,_mat[j,:,i,2]),axis=0)
                AMPS_2D = np.concatenate((AMPS_2D,_mat[j,:,i,0]), axis=0)
                
            #SORT AMPLS AND ZSCORE ACCORDING TO POSITIONS
        SORTED_2D_AMPS = [x for _, x in sorted(zip(POSITIONS_2D,AMPS_2D))]
        SORTED_2D_ZSCORES = [x for _, x in sorted(zip(POSITIONS_2D,ZSCORES_2D))]
        
        ax[j].plot(sorted(POSITIONS_2D),SORTED_2D_ZSCORES)
        ax[j].plot(sorted(POSITIONS_2D),np.ones(len(sorted(POSITIONS_2D)))*Z_score_limit,linestyle='--')
        label_line = j+1
        ax[j].set_ylabel('Zscore line %s'%label_line)
        
        if j == H-1:
            ax[j].set_xlabel('Distance (P1- norm)')
                
        zebrin_color='green'
        
        #P2+ contra zebrin band
        ax[j].axvspan(raw_zebrin_values[1], raw_zebrin_values[2], alpha=0.2,color=zebrin_color)        
        #P1+ zebrin band
        ax[j].axvspan(raw_zebrin_values[3], 0, alpha=0.2,color=zebrin_color)       
        #P2+ ipsi zebrin band
        ax[j].axvspan(raw_zebrin_values[4], raw_zebrin_values[6], alpha=0.2,color=zebrin_color)     
        #PC_cluster
        ax[j].axvspan(raw_zebrin_values[5]-raw_zebrin_stds[5], raw_zebrin_values[5]+raw_zebrin_stds[5],ymax=0.1, alpha=0.5,color='red')

        #BINNING FOR MEDIAN CALCUL
        step = bin_for_median #In % of P1- : 12 c'est pas mal 
        binning = np.arange(-300,300+step,step)
        
        _MEDS, _MADS, _POS, _COUNTS, _AMPS, _SUM = [],[],[],[],[],[]
        
        for i in range(len(binning)):
            
            if i == len(binning)-1:
                break
            
            start, stop = binning[i],binning[i+1]
            _meds, _mads, _pos, _count, _amps, _sum = [],[],[],[],[],[]
            
            #print ('Bin %s to %s'%(start, stop))
            
            SORTED_POSITIONS = sorted(POSITIONS_2D)
            
            for j in range(len(SORTED_POSITIONS)):
                if start < SORTED_POSITIONS[j] <= stop:
                    if np.isnan(SORTED_2D_ZSCORES[j])==False:
                        _meds.append(SORTED_2D_ZSCORES[j])
                        _pos.append(SORTED_POSITIONS[j])
                        _amps.append(SORTED_2D_AMPS[j])
                        _sum.append(SORTED_2D_AMPS[j])
                    
            _MEDS.append(np.nanmedian(_meds))
            _COUNTS.append(np.count_nonzero(_meds))
            _POS.append(np.nanmedian(_pos))
            _MADS.append(MAD(_meds, axis=0))
            _AMPS.append(np.nanmean(_amps,axis=0))
            _SUM.append(np.nansum(_sum,axis=0))
        
        _MEDIAN_ZSCORE_2D.append(np.asarray(_MEDS))
        _AVERAGE_AMP_2D.append(np.asarray(_AMPS))
        _COUNT_2D.append(np.asarray(_COUNTS))
        _POSITIONS_2D.append(np.asarray(_POS))
        _SUM_2D.append(np.asarray(_SUM))

    if SAVING == True:
        plt.savefig('{}/{}_{}_2D_raw_sorting.pdf'.format(savedir,groups[group],method))

    fig, ax = plt.subplots(2,1,figsize=(14,5))

    plt.suptitle('%s 2D maps'%groups[group])
    ax[0].set_title('Median Zscore')
    median_zscore_2d = ax[0].imshow(stack_lines(_MEDIAN_ZSCORE_2D),interpolation='kaiser', cmap='hot',vmin=0,aspect='auto')
    fig.colorbar(median_zscore_2d, ax=ax[0])
        
    ax[1].set_title('Average %s'%tag_str)
    ax[1].set_xticks(np.arange(0,len(_POSITIONS_2D[0]),1))
    ax[1].set_xticklabels(_POSITIONS_2D[0].astype(int),rotation=-90)
    mean_amplitude_2d = ax[1].imshow(stack_lines(_AVERAGE_AMP_2D),interpolation='kaiser', cmap= 'hot_r',vmax=vmax,vmin=vmin,aspect='auto')
    fig.colorbar(mean_amplitude_2d,ax=ax[1])   

    if SAVING == True:
        plt.savefig('{}/{}_2D_median_{}_Zscore.pdf'.format(savedir,groups[group],method))
    
    

