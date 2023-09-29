# Python script to convert .mat image volume files to mha so hey can be viewed in ITK SNAP

# Import relevant libraries
import SimpleITK as sitk
from scipy.io import loadmat
import numpy as np

# Define Patient ID, model type, and volume name
PatNum = 'BAR_009'
ModelTypes = ['Original', 'No VASC', 'No VASC Reduced Rs 1', 'No VASC Reduced Rs 2', 'No VASC Reduced Rs 3', 'No VASC Reduced Rs 4']
VolumeName = 'fIC'

for ModelType in ModelTypes:
    # Load .mat file
    volume = loadmat(f'VERDICT outputs/{PatNum}/{ModelType}/{VolumeName}.mat')[VolumeName]
    # remove infinities
    volume[volume == np.inf] = 0
    # Change image orientation
    volume = np.moveaxis( volume , -1, 0)  

    # Save as mha file
    sitk.WriteImage( sitk.GetImageFromArray(volume), f'VERDICT outputs/{PatNum}/{ModelType}/{VolumeName}.mha' )




