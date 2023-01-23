# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 17:56:19 2023

@author: klab
"""

import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams.update({'font.size': 7})
import random
import numpy as np

trialTypes = ['freeze_reward', 'freeze_error', 'freeze_reward_omit', 'freeze_error_give_reward']

prob_to_give_reward_in_correct_trials = 90
prob_to_give_reward_in_error_trials = 10

listOfCorrectTrials = []
listOfErrorTrials = []

correctTrialType = None
errorTrialType = None

for i in range(1000): 

    pick = random.choice([x+1 for x in range(100)])
    
    #Determine probability for correct trial
    
    if pick > prob_to_give_reward_in_correct_trials:
        
        correctTrialType = 'freeze_reward_omit'
        
    else:
        correctTrialType = 'freeze_reward'
        
    listOfCorrectTrials.append(correctTrialType)
        
    

    #Determine probability for error trial  
        
    if pick > prob_to_give_reward_in_error_trials:
        
        errorTrialType = 'freeze_error'
        
    else:
         
        errorTrialType = 'freeze_error_give_reward'
        
    listOfErrorTrials.append(errorTrialType)
    
fig, ax = plt.subplots()

ax.set_ylabel('% of trials')

ax.bar([0,1,2,3], 
       [len([x for x in listOfCorrectTrials if x == 'freeze_reward'])/len(listOfCorrectTrials)*100, 
        len([x for x in listOfCorrectTrials if x == 'freeze_reward_omit'])/len(listOfCorrectTrials)*100,
        len([x for x in listOfErrorTrials if x == 'freeze_error'])/len(listOfErrorTrials)*100,
        len([x for x in listOfErrorTrials if x == 'freeze_error_give_reward'])/len(listOfErrorTrials)*100], 
       width=0.5)

ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['freeze_reward','freeze_reward_omit','freeze_error','freeze_error_give_reward'])







