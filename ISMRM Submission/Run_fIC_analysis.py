# Python script to run fIC analysis functions

# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys

# Import fIC analysis functions
import fIC_analysis

# First define model type to analyse
ModelNum = 1

# ROI Name
ROIName = 'L1_b3000_NT'


# Find list of patients
fnames = glob.glob(f'VERDICT outputs/*/Model {ModelNum}/fIC.mat')

PatNums = [ os.path.split( os.path.split( os.path.split(fname)[0])[0])[1] for fname in fnames  ]


# Create list for 'good PatNums' which have a well defined ROI from Natasha
goodPatNums = []

# For each patient, save lesion ROI masks and extract fICs
for indx, PatNum in enumerate(PatNums):
     
     
    # # Check if ROI mask exists and extract fICs from lesion
    # if os.path.isfile(f'ROIs/{PatNum}/{ROIName}.npy'):
    #     goodPatNums.append(PatNum)
    #     # Extract fICs
    #     fIC_analysis.extractROIfICs(PatNum, ROIName = ROIName, ModelNum = ModelNum)
    #     continue
    # else:
    #     print(f'No existing ROI mask for {PatNum}')
    
    
    # If ROI mask doesn't exist, try and extract ROI 
    try:
        
        # Lesion mask
        fIC_analysis.saveROImask(PatNum, ROIName = ROIName)
        
        # Extract fICs
        fIC_analysis.extractROIfICs(PatNum, ROIName = ROIName, ModelNum = ModelNum)
        
        # Append to goodPatNums list
        goodPatNums.append(PatNum)
        
    # If no ROI RTstruct file, or error, remove patient number from list
    except:
        print(f'ROI error for {PatNum}')

   
PatNums = goodPatNums



# Calculate median fICs  
fIC_analysis.avgROIfICs(PatNums = PatNums, ROIName = 'L1_b3000_NT', ModelNum = ModelNum, avg_type = 'median')

# Apply ROC analysis
fpr, tpr, thresholds, roc_auc = fIC_analysis.fIC_ROC( ModelNum = ModelNum)


# == Save results

# Make folder
try:
    os.makedirs(f'ROC results/Model {ModelNum}')
except: 
    None
    
# Save ROC results
np.save(f'ROC results/Model {ModelNum}/fpr.npy', fpr)
np.save(f'ROC results/Model {ModelNum}/tpr.npy', tpr)
np.save(f'ROC results/Model {ModelNum}/thresholds.npy', thresholds)
np.save(f'ROC results/Model {ModelNum}/roc_auc.npy', roc_auc)

Youden_Indices = tpr-fpr

print(f'Thresholds: {thresholds}')
print(f' Best threshold: {thresholds[Youden_Indices == np.max(Youden_Indices)]}')

plt.figure()
plt.plot(fpr,tpr,marker = '*')
plt.show()