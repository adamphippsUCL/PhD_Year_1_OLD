# Python script to contain functions to analyse fIC values with lesion ROI

# Import relevnat libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from scipy.io import loadmat
import pandas as pd

# Import DICOM from imgtools
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM # type: ignore


# Function for saving ROI masks
def saveROImask(
                PatNum, 
                ROIName, 
                INNOVATE_path = r"D:\UCL PhD Imaging Data\INNOVATE", 
                INNOVATE_ROIs_path = r"D:\UCL PhD Imaging Data\INNOVATE ROIs",
                output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs"):
    
    '''
    Function to save mask of ROI in RTStruct file.
    
    ROIs drawn on b3000 VERDICT scan as this is used as the registration target
    
    '''
    
    # Define path to DICOMs for patient
    Img_DICOM_path = f'{INNOVATE_path}/{PatNum}\scans'
    
    # Find filename for b3000 image
    b3000_DICOM_fname = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*')[0]

    # Create dcm object 
    b3000dcm = DICOM.MRIDICOM(b3000_DICOM_fname)
    b3000_ImageArray = b3000dcm.constructImageArray()
    b3000_bVals = b3000dcm.DICOM_df['Diffusion B Value'].values

    # Define path to ROI DICOM
    RTStruct_path = f'{INNOVATE_ROIs_path}/{PatNum}'

    # Find RTstruct filename
    RTStruct_fname = glob.glob(f'{RTStruct_path}/*')[0]


    # === Create ROI mask

    # Instantiate contours object
    contours = DICOM.contours(RTStruct_fname)

    # Define lesion structure number (hardcoded here but should be found automatically in future)
    LesionStructNum = contours.Struct_Name_Num_dict[ROIName]
 
    # Make mask for lesion
    LesionMask = contours.create_mask(Struct_Num = LesionStructNum, DICOM_fname = b3000_DICOM_fname)

    # Remove duplicate spatial slices
    LesionMask = LesionMask[b3000_bVals == 0]

    # Save lesion mask
    np.save(f'{output_path}/{PatNum}/{ROIName}.npy', LesionMask)
                    
                                  
# Function for extracting fIC values from ROI
def extractROIfICs(
    PatNum,
    ROIName,
    ModelType,
    VERDICT_output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\VERDICT outputs",
    ROI_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs",
    output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\fIC results",
):
    
    '''
    Function to extract fIC values from within an ROI
    
    '''                
    
    # Load fIC array
    fIC = loadmat(f'{VERDICT_output_path}/{PatNum}/{ModelType}/fIC.mat')['fIC']
    # Permute array axes (account for MATLAB Python differences)
    fIC = np.moveaxis( fIC , -1, 0)   
    
    # Load ROI mask
    ROIMask = np.load(f'{ROI_path}/{PatNum}/{ROIName}.npy')
    # Make Bool
    ROIMask = (ROIMask != 0)
    
    # Extract fIC from ROI
    ROI_fIC = fIC[ROIMask]
    
    # Save as numpy array
    try:
        os.makedirs(f'{output_path}/{PatNum}/{ROIName}')
    except:
        None
        
    np.save(f'{output_path}/{PatNum}/{ROIName}/{ModelType}.npy', ROI_fIC)
    
       
def readBiopsyResults(
    biopsy_data_xlsx = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE\INNOVATE patient groups 2.0.xlsx'
):
    
    '''
    Function to read biopsy data excel sheet and save results as binary dataframe
    
    '''
    
    BiopsyDF = pd.read_excel(biopsy_data_xlsx)
    
    # Clinically significant patients
    csPats = list( (BiopsyDF['Clinically significant (Gleason >=7)'].values)[1:] )
    
    # Ones (binary 1 for cs)
    ones = list(np.ones(len(csPats)))
    
    # Non-clincially significant patients 
    ncsPats = list( (BiopsyDF['Clinically insignificant (Gleason <7)'].values)[1:] )
    
    # Zeros (binary 0 for ncs)
    zeros = list(np.zeros(len(ncsPats)))
    
    
    # Make dataframe
    Patients = csPats + ncsPats
    Results = ones + zeros
    
    BiopsyResultsDF = pd.DataFrame({'Patient ID': Patients, 'Biopsy Result': Results})
    
    return BiopsyResultsDF
         
       
# Function to calculate avg fIC over ROIs
def avgROIfICs(
    ROIName,
    ModelType,
    avg_type = 'mean',
    results_path = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\fIC results'
    
):
    
    # Extract list of fIC results filenames
    fIC_fnames = glob.glob(f'{results_path}/*/{ROIName}/{ModelType}.npy')
    
    # For each file, extract patient number and calculate average fIC
    PatNums = []
    avgfICs = []
    
    for fname in fIC_fnames:
        
        # Find patient number
        PatNum = os.path.split( 
                               os.path.split(
                                   os.path.split(fname)[0]
                               )[0]
                               )[1]

        # Append to list
        PatNums.append(PatNum)
        
        # Load fICs
        fICs = np.load(fname)
        
        # Calculate average
        if avg_type == 'mean':
            avg_fIC = np.mean(fICs)
        elif avg_type == 'median':
            avg_fIC = np.median(fICs)
        else:
            print('Incorrect input, default to mean')
            avg_fIC = np.mean(fICs)
            
        # Append to list
        avgfICs.append(avg_fIC)
        
        
        # Create dataframe
        fIC_DF = pd.DataFrame({'Patient ID': PatNums, 'Avg fIC': avgfICs})
        
        print(fIC_DF)