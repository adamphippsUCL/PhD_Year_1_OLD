# Python script to combine VERDICT images into array suitable for verdict_fit input

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
import sys

# Import imgtools functions
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM, EPI # type: ignore


# Define DICOM filenames
b90_fname = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\20230912 Simplified VERDICT model\Test DICOMs\I7'
b500_fname = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\20230912 Simplified VERDICT model\Test DICOMs\I6'
b1500_fname = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\20230912 Simplified VERDICT model\Test DICOMs\I5'
b2000_fname = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\20230912 Simplified VERDICT model\Test DICOMs\I4'
b3000_fname = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\20230912 Simplified VERDICT model\Test DICOMs\I3'

# List of filenames
bVals_fnames = [b90_fname, b500_fname, b1500_fname, b2000_fname, b3000_fname]


# Read each image as dcm object
NormalisedImages = []

for fname in bVals_fnames:
    
    # Make dcm object for b value
    bdcm = DICOM.MRIDICOM(fname)
    
    # Image array
    b_ImageArray = bdcm.constructImageArray()
    
    # b values and diffusion type of each slice
    b_bVals = bdcm.DICOM_df['Diffusion B Value'].values
    
    # Find b=0 and b!=0 images
    b_b0Image = b_ImageArray[b_bVals == 0]
    b_bImage = b_ImageArray[(b_bVals != 0)*(bdcm.DICOM_df['Diffusion Gradient Orientation'].isnull().values)]

    # Normalised image
    b_Normalised = b_bImage/b_b0Image
    
    NormalisedImages.append(b_Normalised.T)
    
    
Y = np.stack(NormalisedImages, axis = -1)

# Create dictionary
Ydata = dict(zip(['Y'], [Y]))

# Save as matlab file
savemat('20230912 Simplified VERDICT model/Test fIC/Y.mat', Ydata)


