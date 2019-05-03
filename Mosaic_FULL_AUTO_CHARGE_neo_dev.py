# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:20:21 2019

@author: ludovic.spaeth
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:40:03 2019

Created on Tue Nov 27 16:14:33 2018

Script built to compute Synaptic maps from in vitro uncaging experiments, in a "automatic" way

What the script does : 
    - loads WinWcp (.wcp) files containing voltage clamp recordings
    - Gets info for spacing, orientation and so on from GC_PC_mappings_data_info.xlsx excel file
    - Gets Zebrin positions and info from Mesures_ZII_HighRes_WT_Cuff_Sham_Enr.xlsx excel file
    
    - Calculate average synaptic charge, average noise, and zscore from ephy data
    - Build synaptic map in consensus orientation 
    
    - Save figures and tables 

@author: ludovic.spaeth
"""

from neo.io.winwcpio import WinWcpIO as Win
from matplotlib import pyplot as plt
import numpy as np

from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.interpolate import UnivariateSpline

from scipy import trapz


import pandas as pd

#-----------------------DIRECTORY & INFO---------------------------------------
#--------------CHECK ALL THESE INFO BEFORE RUNNING----------------------------- 


DONE = [] #To append names of maps that are done 


FOLDERS = ['ENR 11-05-2017','ENR 16-05-2017','ENR 16-05-2017','ENR 16-05-2017',
           'ENR 17-05-2017','ENR 17-05-2017','ENR 18-05-2017','ENR 18-05-2017',
           'ENR 18-05-2017','ENR 19-05-2017','ENR 19-05-2017','ENR 21-03-2019',
           'ENR 21-03-2019','ENR 22-03-2019','ENR 22-03-2019','ENR 23-03-2019','ENR 23-03-2019']

NAMES = ['110517(1)','160517(1)','160517(2)','160517(3)',
         '170517(1)','170517(2)','180517(1)','180517(2)',
         '180517(3)','190517(1)','190517(2)','210319(2)',
         '210319(3)','220319(1)','220319(2)','230319(1)','230319(2)']


df  = pd.read_excel(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\Mesures_ZII_HighRes_WT_Cuff_Sham_Enr.xlsx",sheet_name='ENR',index_col=0,header=1)  #DataFrame containing zebrin band file 
    
filelist  = pd.read_excel(r"U:\01_ANALYSIS\01_BRAVE_NEW_WORLD\GC-PC_mappings_data_info.xlsx",sheet_name='ENR',index_col=0,header=0)  #DataFrame containing zebrin band file 


for NAME,FOLDER,IDX in zip(NAMES,FOLDERS,range(len(NAMES))):
    name = NAME
    manip = NAME.replace('WT_','')
    directory = 'U:/RAW DATA/Enrichissement/{}'.format(FOLDER)
    
    savedir = 'U:/01_ANALYSIS/01_BRAVE_NEW_WORLD/CHARGE_ANALYSIS/00_MAPS/ENR'
    
    SAVING = True
    CLOSE_FIG = True
     
    ch = int(filelist.loc[[manip],['Channel']].values.ravel()[0]) #Channel 0 for Im1 and 2 for Im2
    
    reshape = int(filelist.loc[[manip],['Ordre']].values.ravel()[0])
    global_reshape = np.array([int(x) for x in str(reshape)])

    #-------------------------GETS FILELIST----------------------------------------
    ZSOCRE_CUT = 3.0
    
    global_reshape = global_reshape-1
    
    scanspot_list = ['Scanspot 1','Scanspot 2','Scanspot 3','Scanspot 4']
    
    
    files = []
    
    for i in range(global_reshape.size):
        idx_list = np.ravel(filelist.loc[[manip],[scanspot_list[i]]].values)
        record_list  = []
        
        for idx in idx_list :
            print (type(idx))
            print (idx)
            if isinstance(idx,str) == True:
    
                record_list.append(idx)
                
        files.append(record_list)
        
    files = np.asarray(files)
    
    #----------------------ZEBRIN BANDS & ORIENTATION-------------------------------
    _BANDS_micron =  df.loc[['%s'%name],'P2- contra':'P2- ipsi']  #Selects values in xcl file from P2- contra to P2- ipsi
        
    _BANDS_norm =  df.loc[['%s norm_P1-'%name],'P2- contra':'P2- ipsi']  #Selects values in xcl file from P2- contra to P2- ipsi
    
    orientation= df.loc[['%s'%name],['Position']]   # Gets the piece of dataframe
    orientation= orientation.iloc[0,0]              # Gets the float value in the df 
    orientation= int(orientation)                   #Converts float in integer 
        
    position_left = df.loc[['%s'%name],['Pos Left']] # 32 if cell has been centered, may change in case of paired recordings 
    position_left = position_left.iloc[0,0]
    position_left = int(position_left)
        
    position_right = df.loc[['%s'%name],['Pos Right']] # 32 if cell has been centered, may change in case of paired recordings 
    position_right = position_right.iloc[0,0]
    position_right = int(position_right)
    
    
    P2mcontra = _BANDS_norm.iloc[0,0]
    P2pcontra = _BANDS_norm.iloc[0,1]
    P1mcontra = _BANDS_norm.iloc[0,2]
    P1p	= _BANDS_norm.iloc[0,3]
    BordgaucheP1m = 0.
    BorddroitP1m = 100.
    cp = _BANDS_norm.iloc[0,5]
    P2pipsi   = _BANDS_norm.iloc[0,6]
    P2mipsi   = _BANDS_norm.iloc[0,7]
            
    zebrin = np.array([P2pcontra,P1mcontra,P1p,BordgaucheP1m,cp,BorddroitP1m,P2pipsi,P2mipsi])
    
    #Position et bords de la carte ----------------------------------------------------------------------------------------------
    print ('Position de la cellule: ',cp,"%",'de P1-'	)						
    p1_minus = _BANDS_micron.iloc[0,4]			                               # Taille de P1- en microns
    print ('P1- mesure' ,p1_minus, ' microns')
    step = float(20)														# Taille de sites de photostim en microns
    print (step,' microns par site de photostim')
    step2 = float((step/p1_minus)*100)										# Calcul du step en % de P1-
    print ('soit ',step2," %",'de P1-')
    lb = position_left*step/p1_minus*100									# Calcul bord gauche de la map en % de P1-
    print ('Bord Gauche de la carte = ',lb,"%")
    rb = position_right*step/p1_minus*100								    # Calcul du bord droit de la map en % de P1-
    print ('Bord Droit de la carte = ',rb,"%")
    
    step_positions = (lb+rb)/63
    positions = np.linspace(-lb,rb,64)
    positions_cp_centered = np.linspace(-lb,rb,64)+cp
    
    
    grid = np.array([[1,17,5,21,9,25,13,29,2,18,6,22,10,26,14,30],
                     [50,35,54,39,96,42,91,46,79,36,53,40,85,45,56,48],
                     [31,3,19,7,23,11,27,15,32,4,20,8,24,12,28,16],
                     [59,66,89,60,75,58,84,71,88,62,95,77,82,57,74,87],
                     [33,73,37,70,41,65,44,80,34,68,38,93,43,94,47,61],
                     [64,49,78,52,76,69,86,63,92,51,72,55,83,67,90,81]])
                     
    sites = grid.size
    
    
    
    
    _GRIDS_OF_SIGNAL_AMPLITUDE = []
    _GRIDS_OF_NOISE_AMPLITUDE = []
    
    
    #-----------SCAN SPOT CALCUL---------------------------------------------------
    for scanspot in range(files.shape[0]):  
        
        records = files[scanspot]
    
        r = Win(filename='%s/%s.wcp'%(directory,records[0]))
        bl = r.read_block()
        
        sampling = float(bl.segments[0].analogsignals[ch].sampling_rate)
        points = len(bl.segments[0].analogsignals[ch])
        
        baseline_begin = 0
        baseline_end = int(0.045*sampling)
        
        if IDX <= 1 :
            window_begin = int(0.2*int(sampling))
            window_end = int(0.4*int(sampling))
            noise_begin = int(0.71*sampling)
            noise_end = int(0.91*sampling)  

        else:
            window_begin = int(0.5*int(sampling))
            window_end = int(0.7*int(sampling))
            noise_begin = int(0.21*sampling)
            noise_end = int(0.41*sampling)
            
        time = np.arange(0,points,1)/sampling*1000 #TIME IN MS
        
        
        _SWEEPS = np.zeros((96,len(records),int(window_end-window_begin))) #Remove +1 for paired recordings 
        _NOISE  = np.zeros((96,len(records),int(noise_end-noise_begin)))
    
               
        _GRID_SIGNAL_AMP,_GRID_NOISE_AMP = [],[]
        _GRID_SIGNAL_CHARGE,_GRID_NOISE_CHARGE = [],[]
        
    #---------------MEASURE WINDOW & LEAK REMOVING---------------------------------
        for record in range(len(records)):
            r = Win(filename='%s/%s.wcp'%(directory,records[record]))
            bl = r.read_block()
        
            for (row,col),site in np.ndenumerate(grid-1):
                sampling = float(bl.segments[site].analogsignals[ch].sampling_rate)
                leak = np.mean(bl.segments[site].analogsignals[ch].magnitude[baseline_begin:baseline_end])
                                    
                noise = np.squeeze(bl.segments[site].analogsignals[ch].magnitude[noise_begin:noise_end]-leak)
    
                signal = np.squeeze(bl.segments[site].analogsignals[ch].magnitude[window_begin:window_end]-leak)
           
                _SWEEPS[site,record,:]=signal
                _NOISE[site,record,:]=noise
    
    #-------------------------PLOT AND AMP MEASURES--------------------------------        
        f2,ax = plt.subplots(6,16,sharex=True,sharey=True,figsize=(12,12))
        plt.suptitle('%s scanspot %s'%(name,scanspot+1))
    
        time_ticks = [window_begin,(window_begin+window_end)/2,window_end]
        tick_labels = [0,100,200]
        
        for (row,col),site in np.ndenumerate(grid-1):
            signal_average = np.nanmean(_SWEEPS[site],axis=0)
            noise_average = np.nanmean(_NOISE[site],axis=0)
            
            _SIGNAL_CHARGE = trapz(signal_average,dx=(1./sampling))
            _NOISE_CHARGE = trapz(noise_average,dx=(1./sampling))
            
    
            ax[row,col].plot(noise_average,color='0.5')        
            ax[row,col].plot(signal_average,color='blue')
            ax[row,col].set_title(site+1)
    
            for index in range(len(signal_average)):
                if signal_average[index] == np.min(signal_average):
                    min_signal = index
                    signal_window = signal_average[min_signal-50:min_signal+50]
    
                    
            for index in range(len(noise_average)):            
                if noise_average[index] == np.min(noise_average):
                    min_noise = index
                    noise_window = noise_average[min_noise-50:min_noise+50]
                    
    
            _GRID_SIGNAL_AMP.append(_SIGNAL_CHARGE)
            _GRID_NOISE_AMP.append(_NOISE_CHARGE)
        
        if SAVING == True:
            plt.savefig(r'%s/%s_ScanSpot_%s.png'%(savedir,name,scanspot+1))
        else :
            print ('Fig is not saved')
        
        if CLOSE_FIG == True:
            plt.close()
    #----------------RESHAPE AND CO------------------------------------------------
        _GRID_SIGNAL_AMP = np.asarray(_GRID_SIGNAL_AMP)    
        _GRID_NOISE_AMP = np.asarray(_GRID_NOISE_AMP)
        
        _IMG_GRID_SIGNAL_AMP = np.reshape(_GRID_SIGNAL_AMP,(grid.shape[0],grid.shape[1]))    
        _IMG_GRID_NOISE_AMP = np.reshape(_GRID_NOISE_AMP,(grid.shape[0],grid.shape[1]))
        
        _GRIDS_OF_SIGNAL_AMPLITUDE.append(_IMG_GRID_SIGNAL_AMP)
        _GRIDS_OF_NOISE_AMPLITUDE.append(_IMG_GRID_NOISE_AMP)    
       
    #---------------------GLOBAL RESHAPE & CONCATENATE-----------------------------
    
    if global_reshape.size == 1 : #Only the first scanspot
        
        fake_map = np.empty((grid.shape[0],grid.shape[1]))
        fake_map.fill(np.nan)
        
        _FULL_SIGNAL_MAP = np.concatenate((fake_map,fake_map,_GRIDS_OF_SIGNAL_AMPLITUDE[0],fake_map),axis=1)
    
        _FULL_NOISE_MAP = np.concatenate((fake_map,fake_map,_GRIDS_OF_NOISE_AMPLITUDE[0],fake_map),axis=1)
    
        
      
    elif global_reshape.size == 4: #Full map
    
        _FULL_SIGNAL_MAP = np.concatenate((_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[0]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[1]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[2]],_GRIDS_OF_SIGNAL_AMPLITUDE[global_reshape[3]]),axis=1)
    
        _FULL_NOISE_MAP = np.concatenate((_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[0]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[1]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[2]],_GRIDS_OF_NOISE_AMPLITUDE[global_reshape[3]]),axis=1)
    
    else: #Adapt because of lacking scanspot
        fake_map = np.empty((grid.shape[0],grid.shape[1]))
        fake_map.fill(np.nan)
        
        _FULL_SIGNAL_MAP = np.array([])
        _FULL_NOISE_MAP = np.array([]) 
    
        for i in range(len(global_reshape)):
    
            _ind = global_reshape[i]   
            
            if _ind == global_reshape[0]:        
                _FULL_SIGNAL_MAP = _GRIDS_OF_SIGNAL_AMPLITUDE[_ind]
                _FULL_NOISE_MAP = _GRIDS_OF_NOISE_AMPLITUDE[_ind]
    
    
            if _ind != global_reshape[0]:                    
                _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,_GRIDS_OF_SIGNAL_AMPLITUDE[_ind]),axis=1)
                _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,_GRIDS_OF_NOISE_AMPLITUDE[_ind]),axis=1)
    
            
        if np.array_equal(global_reshape,np.array((1,0))) == True:   #IF SCANSPOT 3 & 4 ARE MISSING IN 2-1-3-4 CONFIG          
            _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,fake_map,fake_map),axis=1)            
            _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,fake_map,fake_map),axis=1)  
            
        elif np.array_equal(global_reshape,np.array((1,0,2))) == True:  #IF SCANSPOT 4 is MISSING IN 2-1-3-4 CONFIG          
            _FULL_SIGNAL_MAP = np.concatenate((_FULL_SIGNAL_MAP,fake_map),axis=1)            
            _FULL_NOISE_MAP = np.concatenate((_FULL_NOISE_MAP,fake_map),axis=1)  
            
        elif np.array_equal(global_reshape,np.array((0,1))) == True:  #IF SCANSPOT 3 & 4 are MISSING IN 4-3-1-2 CONFIG          
            _FULL_SIGNAL_MAP = np.concatenate((fake_map,fake_map,_FULL_SIGNAL_MAP),axis=1)            
            _FULL_NOISE_MAP = np.concatenate((fake_map,fake_map,_FULL_NOISE_MAP),axis=1)  
            
        elif np.array_equal(global_reshape,np.array((2,0,1))) == True:  #IF SCANSPOT 4 is MISSING IN 4-3-1-2 CONFIG          
            _FULL_SIGNAL_MAP = np.concatenate((fake_map,_FULL_SIGNAL_MAP),axis=1)            
            _FULL_NOISE_MAP = np.concatenate((fake_map,_FULL_NOISE_MAP),axis=1)  
            
        elif np.array_equal(global_reshape,np.array((0,1,2))) == True:  #IF SCANSPOT 4 is MISSING IN 4-1-2-3 CONFIG          
            _FULL_SIGNAL_MAP = np.concatenate((fake_map,_FULL_SIGNAL_MAP),axis=1)            
            _FULL_NOISE_MAP = np.concatenate((fake_map,_FULL_NOISE_MAP),axis=1)  
    
    
    #---------------------------MAP FLIP ACCORDING TO POSITION-----------------------
            
            
    if orientation== 1:
        print ("Map in position 1: nothing to change")
        
    if orientation== 2:
        #Flip et position----------------------------------------------------------
        _FULL_SIGNAL_MAP = np.fliplr(_FULL_SIGNAL_MAP)	            #Flips map left-Right axis
        _FULL_NOISE_MAP = np.fliplr(_FULL_NOISE_MAP)	            #Flips map left-Right axis
        grid = np.fliplr(grid)				                          #Flips grid on left-right axis
        print ("Map in position 2: Flip on medio-lateral axis")
        
    if orientation==3:
        _FULL_SIGNAL_MAP = np.flipud(_FULL_SIGNAL_MAP)		#Flips map on top-down axis
        _FULL_NOISE_MAP = np.flipud(_FULL_NOISE_MAP)		    #Flips map on top-down axis   
        grid = np.flipud(grid)	
        print ("Map in position 3: Flip on top-down axis")   
            
    if orientation== 4:
        #Flip et position----------------------------------------------------------
        _FULL_SIGNAL_MAP = np.fliplr(_FULL_SIGNAL_MAP)	
        _FULL_SIGNAL_MAP = np.flipud(_FULL_SIGNAL_MAP)
        _FULL_NOISE_MAP = np.fliplr(_FULL_NOISE_MAP)	
        _FULL_NOISE_MAP = np.flipud(_FULL_NOISE_MAP)    
        
        grid = np.fliplr(grid)	
        grid = np.flipud(grid)			#Flips on left-right axis and top bottom axis
        print ("Map in position 4: Flip on both axis")
    
    
    #----------------------NOISE HISTOGRAM & MEAN----------------------------------
    noise = np.ravel(_FULL_NOISE_MAP)
    noise = noise[np.logical_not(np.isnan(noise))]
    
    mean_noise = np.mean(noise)
    
    plt.figure()
    plt.suptitle('%s noise distribution'%name)
    
    (mu,sigma) = norm.fit(noise) #Best fit for data
    n,bins,patches, = plt.hist(noise, 100,normed=1,facecolor='gray',alpha=0.75) #Histogram
    fit = mlab.normpdf(bins,mu,sigma) #yaxis fit
    spline = UnivariateSpline(bins, fit-np.max(fit)/2,s=0) #Plane for HWHM
    plt.plot(bins,fit,'r',linewidth=2) 
    r1,r2 = spline.roots() #Finds the roots, r1 = HWHM
    
    if SAVING == True:
        plt.savefig(r'%s/%s_Charge_Map_Noise_Hist.png'%(savedir,name))
      
    if CLOSE_FIG == True:
        plt.close()
    
    
    #--------------------------------2D ZSCORE-------------------------------------
    sigma_2D = np.full((6,64),np.abs(r1))
    noise_2D = np.full((6,64),np.abs(mean_noise))
    Zscore_2D = (np.abs(_FULL_SIGNAL_MAP)-noise_2D)/sigma_2D
    
    #-----------------------------1D PATTERN & ZSCORE------------------------------
    Amp_pattern = np.max(np.abs(_FULL_SIGNAL_MAP),axis=0)
    
    mean_array = np.ones(Zscore_2D.shape[1])*np.abs(mean_noise)
    sigma_array = np.ones(Zscore_2D.shape[1])*np.abs(r1)
    Zscore_1D = (Amp_pattern-mean_array)/sigma_array
    
    #Zscore_1D = np.max(Zscore_2D,axis=0)
    #----------------------------------FIGURE PLOT---------------------------------
    
    fig1,ax = plt.subplots(5,1, sharex=False, sharey=False,figsize=(8.27,11.69))		
    					 
    complete_grid = np.concatenate((grid, grid, grid, grid), axis = 1)
    
    fig1.suptitle('%s' %name, fontsize=20)
    
    v = ax[4].imshow(np.abs(_FULL_SIGNAL_MAP),interpolation='nearest', cmap='hot',vmin = np.abs(r1))
    cbar1 = plt.colorbar(v,orientation='horizontal')
    ax[4].set_title('2D Max Charge',loc='right')
    
    mask = Zscore_2D>=ZSOCRE_CUT				#2D Zscore
    ax[3].imshow(mask, interpolation='none',cmap='hot')
    ax[3].set_title('2D Zscore',loc='right')
    
    zcolor = 'green'; zalpha = 0.5
    ax[2].axvspan(zebrin[0],zebrin[1],color=zcolor,alpha=zalpha)
    ax[2].axvspan(zebrin[2],zebrin[3],color=zcolor,alpha=zalpha)
    ax[2].axvspan(zebrin[2],zebrin[3],color=zcolor,alpha=zalpha)
    ax[2].axvspan(zebrin[5],zebrin[6],color=zcolor,alpha=zalpha)
    
    ax[2].axvspan(zebrin[4]-0.2,zebrin[4]+0.2,color='red',alpha=zalpha) #the cell
    
    
    ax[2].set_xlabel('Distance (P1- normalized - cp aligned)')	
    ax[2].set_title('Zebrin Bands',loc='right')
    ax[2].set_xlim(-lb, rb)
    
    
    x = np.arange(0,len(Amp_pattern),1)
    ax[0].fill_between(x,0,Amp_pattern,color='lightblue',label='Synaptic Charge',interpolate=True)
    ax[0].set_title('Charge pattern',loc='right')
    ax[0].set_xlim(0,len(Amp_pattern)-1)
    ax[0].set_ylabel('Absolute Max Charge (pC)')
    
    ax[1].fill_between(x,0,Zscore_1D,color='0.4',linewidth=2,label='Zscore')
    z_lim = np.ones(len(Zscore_1D))*ZSOCRE_CUT
    ax[1].fill_between(x,z_lim,Zscore_1D,where=Zscore_1D>=z_lim,color='lightcoral',linewidth=2,label='Zscore',interpolate=True)
    ax[1].set_title('1D Zscore',loc='right')
    ax[1].set_xlim(0,len(Zscore_1D)-1)
    ax[1].set_ylabel('Zscore')
    
    if SAVING == True:
        plt.savefig(r'%s/%s_Charge_Map_Zebrin.pdf'%(savedir,name))
    
    if CLOSE_FIG == True:
        plt.close()
    #-----------------------------SAVINGS------------------------------------------
    
    if SAVING == True :
        np.savetxt(r"%s/%s_Charge_zscore_max_OK.csv" %(savedir,name),Zscore_1D, delimiter=',')
        np.savetxt(r"%s/%s_Charge_zscore_2D_OK.csv" %(savedir,name),Zscore_2D, delimiter=',')
        
        np.savetxt(r"%s/%s_Charge_max_OK.csv" %(savedir,name),Amp_pattern, delimiter=',')
        np.savetxt(r"%s/%s_Charge_2D_OK.csv" %(savedir,name),_FULL_SIGNAL_MAP, delimiter=',')
        
        np.savetxt(r"%s/%s_Charge_Noisemap_OK.csv" %(savedir,name),_FULL_NOISE_MAP, delimiter=',')
        np.savetxt(r"%s/%s_Charge_Sigma_OK.csv" %(savedir,name),sigma_2D, delimiter=',')
        
        np.savetxt(r"%s/%s_Zebrin_OK.csv" %(savedir,name),zebrin, delimiter=',')
        
        np.savetxt(r"%s/%s_Positions_OK.csv" %(savedir,name),positions, delimiter=',')
        np.savetxt(r"%s/%s_Positions_cp_centered_OK.csv" %(savedir,name),positions_cp_centered, delimiter=',')
        
        np.savetxt(r"%s/%s_Files_List.csv" %(savedir,name),files,delimiter=',',fmt="%s")
        
    else : 
        print ("No data was saved")
        
    DONE.append(name)
    
print ('These maps have been computed successfully')
print (DONE)
