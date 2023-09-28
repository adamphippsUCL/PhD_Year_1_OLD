# Python script to convert .mat image volume files to mha so hey can be viewed in ITK SNAP

# Import relevant libraries
import SimpleITK as sitk
from scipy.io import loadmat
import numpy as np

# Define Patient ID, model type, and volume name
PatNum = 'BAR_003'
ModelType = 'No VASC Reduced Rs'
VolumeName = 'fIC'

# Load .mat file
volume = loadmat(f'VERDICT outputs/{PatNum}/{ModelType}/{VolumeName}.mat')[VolumeName]
# remove infinities
volume[volume == np.inf] = 0
# Change image orientation
volume = np.moveaxis( volume , -1, 0)  

# Save as mha file
sitk.WriteImage( sitk.GetImageFromArray(volume), f'VERDICT outputs/{PatNum}/{ModelType}/{VolumeName}.mha' )




