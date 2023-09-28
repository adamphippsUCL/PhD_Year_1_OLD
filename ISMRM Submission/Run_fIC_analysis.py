# Python script to run fIC analysis functions

# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Import fIC analysis functions
import fIC_analysis

# First define model type to analyse
ModelType = 'Original'

# Find list of patients
fnames = glob.glob(f'VERDICT outputs/*/{ModelType}/fIC.mat')
PatNums = [ os.path.split( os.path.split( os.path.split(fname)[0])[0])[1] for fname in fnames  ]


# # For each patient, save lesion ROI masks and extract fICs
# for PatNum in PatNums:
    
#     print(PatNum)
#     # Lesion mask
#     fIC_analysis.saveROImask(PatNum, ROIName = 'L1')
    
#     # Extract fICs
#     fIC_analysis.extractROIfICs(PatNum, ROIName = 'L1', ModelType = ModelType)
    
  
# Calculate median fICs  
fIC_analysis.avgROIfICs(ROIName = 'L1', ModelType = 'Original', avg_type = 'median')

# Apply ROC analysis
fpr, tpr, thresholds, roc_auc = fIC_analysis.fIC_ROC(ModelType = 'Original')

Youden_Indices = tpr-fpr

print(thresholds[Youden_Indices == np.max(Youden_Indices)])

print(roc_auc)
plt.figure()
plt.plot(fpr,tpr,marker = '*')
plt.show()