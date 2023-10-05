# Python script to test loading of ROIs and measurement of mean fIC within ROIs

# Import relevant libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
import SimpleITK as sitk

# Import DICOM from imgtools
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM # type: ignore

# ===================================================================

# Specify patient number
pat_num = 'BAR_025'

try:
    # Define path to image DICOMs
    Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{pat_num}\scans" 

    # Find filename for b3000 image
    b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*')
    
    test_valid = b3000_DICOM_fnames[0]
    
except:
    # Define path to image DICOMs
    Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{pat_num}" 

    # Find filename for b3000 image
    b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}\*b3000_80\DICOM\*')
    
   

# Test if DICOM MF or SF
if len(b3000_DICOM_fnames) == 1:
    # MF
    b3000_DICOM_fname = b3000_DICOM_fnames[0]
    b3000dcm = DICOM.MRIDICOM(b3000_DICOM_fname)
    
elif len(b3000_DICOM_fnames) > 1:
    # SF
    multiframe = False
    b3000dcm = DICOM.MRIDICOM(DICOM_fnames = b3000_DICOM_fnames, multiframe = multiframe)
    
else:
    print('No')
    
    

b3000_ImageArray = b3000dcm.constructImageArray()
b3000_bVals = b3000dcm.DICOM_df['Diffusion B Value'].values

# Define path to ROI DICOM
RTStruct_path = rf"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE\INNOVATE ROIs\{pat_num}"

# Find RTstruct filename
RTStruct_fname = glob.glob(f'{RTStruct_path}/*')[0]



# === Create ROI mask

# Instantiate contours object
contours = DICOM.contours(RTStruct_fname)

LesionStructName = 'L1_b3000_NT'

# Define lesion structure number (hardcoded here but should be found automatically in future)
LesionStructNum = contours.Struct_Name_Num_dict[LesionStructName]

# Make mask for lesion
LesionMask = contours.create_mask(Struct_Num = LesionStructNum, DICOM_dcm = b3000dcm)



# Remove duplicate spatial slices
LesionMask = LesionMask[b3000_bVals == 0]

LesionSliceIndx = np.where( np.sum(LesionMask, axis = (1,2)) != 0)[0][0]


# Save lesion mask
np.save(f'ROIs/{pat_num}/{LesionStructName}.npy', LesionMask)

# Save mask an b0 from 3000 as mha
sitk.WriteImage(sitk.GetImageFromArray(LesionMask), 'LesionMask.mha')
sitk.WriteImage(sitk.GetImageFromArray(b3000_ImageArray[b3000_bVals == 0]), 'b0from3000.mha')


plt.figure()
plt.imshow(LesionMask[6], cmap = 'gray')
plt.figure()
plt.imshow(b3000_ImageArray[0::5][6], cmap = 'gray')
plt.show()


