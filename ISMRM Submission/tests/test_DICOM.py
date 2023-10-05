# Script to test SF DICOM function
import sys
import glob
import matplotlib.pyplot as plt

# Import DICOM from imgtools
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM # type: ignore



# Specify patient number
pat_num = 'BAR_012'


# Define path to image DICOMs
Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{pat_num}" 

# Find filename for b3000 image
b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*')


dcm = DICOM.MRIDICOM(DICOM_fnames = b3000_DICOM_fnames, multiframe = False )

dcm.df2excel()

array = dcm.constructImageArray()

print(dcm.DICOM_df['Diffusion B Value'].values)
print(array.shape)

# plt.figure()
# plt.imshow(array[4,:,:], cmap = 'gray')
# plt.show()


vc = dcm.constructVoxelCoordinates()
print(vc.shape)

