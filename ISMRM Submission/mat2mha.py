# Python script to convert .mat image volume files to mha so hey can be viewed in ITK SNAP

# Import relevant libraries
import SimpleITK as sitk
from scipy.io import loadmat
import numpy as np
import sys

# Define Patient ID, model type, and volume name
PatNum = 'BAR_038'
ModelNumbers = list(range(1,15))

VolumeName = 'vb0'

for ModelNum in ModelNumbers:
    
    try:
        # Load .mat file
        volume = loadmat(f'VERDICT outputs/{PatNum}/Model {ModelNum}/{VolumeName}.mat')[VolumeName]
        # remove infinities
        volume[volume == np.inf] = 0
        # Remove nan
        volume[np.isnan(volume)] = 0
        
        # Change image orientation
        volume = np.moveaxis( volume , -1, 0)  

        # Save as mha file
        sitk.WriteImage( sitk.GetImageFromArray(volume), f'VERDICT outputs/{PatNum}/Model {ModelNum}/{VolumeName}.mha' )
        
    except: None




