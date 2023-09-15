# Python script to test loading of ROIs and measurement of mean fIC within ROIs

# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os

# Import DICOM from imgtools
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM # type: ignore

# ===================================================================

# Specify patient number
pat_num = 'INN_129'

# Define path to image DICOMs
Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{pat_num}\scans" 

# Find filename for b3000 image
b3000_DICOM_fname = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*')[0]

# Create dcm object 
b3000dcm = DICOM.MRIDICOM(b3000_DICOM_fname)
b3000_ImageArray = b3000dcm.constructImageArray()
b3000_bVals = b3000dcm.DICOM_df['Diffusion B Value'].values

# Define path to ROI DICOM
RTStruct_path = rf"D:\UCL PhD Imaging Data\INNOVATE ROIs\{pat_num}"

# Find RTstruct filename
RTStruct_fname = glob.glob(f'{RTStruct_path}/*')[0]



# === Create ROI mask

# Instantiate contours object
contours = DICOM.contours(RTStruct_fname)


# Define lesion structure number (hardcoded here but should be found automatically in future)
LesionStructNum = 2

# Make mask for lesion
LesionMask = contours.create_mask(Struct_Num = LesionStructNum, DICOM_fname = b3000_DICOM_fname)

# Remove duplicate spatial slices
LesionMask = LesionMask[b3000_bVals == 0]

# Save lesion mask
np.save(f'ROIs/{pat_num}/lesion.npy', LesionMask)



plt.figure()
plt.imshow(LesionMask[6], cmap = 'gray')
plt.figure()
plt.imshow(b3000_ImageArray[0::5][6], cmap = 'gray')
plt.show()


