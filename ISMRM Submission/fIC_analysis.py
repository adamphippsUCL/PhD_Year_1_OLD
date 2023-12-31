# Python script to contain functions to analyse fIC values with lesion ROI

# Import relevnat libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from scipy.io import loadmat
import pandas as pd
import pickle
from sklearn.metrics import roc_curve, roc_auc_score
import SimpleITK as sitk
from skimage.metrics import mean_squared_error

# Import DICOM from imgtools
sys.path.insert(0, r'C:\Users\adam\OneDrive - University College London\UCL PhD\MRes Year\Project\MRes_Project\Prostate MRI Project')
from imgtools import DICOM # type: ignore


# Function for saving ROI masks
def saveROImask(
                PatNum, 
                ROIName, 
                INNOVATE_path = r"D:\UCL PhD Imaging Data\INNOVATE", 
                INNOVATE_ROIs_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE\INNOVATE ROIs",
                output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs"):
    
    '''
    Function to save mask of ROI in RTStruct file.
    
    ROIs drawn on b3000 VERDICT scan as this is used as the registration target
    
    '''
    
    try:
        # Define path to image DICOMs
        Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{PatNum}\scans" 

        # Find filename for b3000 image
        b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}/*b3000_80/DICOM/*.dcm')
        
        # Test to throw error if empty list
        test_valid = b3000_DICOM_fnames[0]
 
    except:
        
        # Define path to image DICOMs
        Img_DICOM_path = rf"D:\UCL PhD Imaging Data\INNOVATE\{PatNum}" 

        # Find filename for b3000 image
        b3000_DICOM_fnames = glob.glob(f'{Img_DICOM_path}\*b3000_80\DICOM\*.dcm')


    # Test if DICOM MF or SF
    if len(b3000_DICOM_fnames) == 1:
        # MF
        multiframe = True
        b3000_DICOM_fname = b3000_DICOM_fnames[0]
        b3000dcm = DICOM.MRIDICOM(b3000_DICOM_fname)
        
    elif len(b3000_DICOM_fnames) > 1:
        # SF
        multiframe = False
        b3000dcm = DICOM.MRIDICOM(DICOM_fnames = b3000_DICOM_fnames, multiframe = multiframe)
        
    else:
        print(f'No b3000 file for patient {PatNum}')
        sys.exit()
        
        
    # Create dcm object 
    b3000_ImageArray = b3000dcm.constructImageArray()
    b3000_bVals = b3000dcm.DICOM_df['Diffusion B Value'].values
    b0 = b3000_ImageArray[b3000_bVals == 0]
    
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
    LesionMask = contours.create_mask(Struct_Num = LesionStructNum, DICOM_dcm = b3000dcm)

    # Remove duplicate spatial slices
    LesionMask = LesionMask[b3000_bVals == 0]
    
    
    
    '''New code: accounting for z flips of fIC map'''
    
    # == First, load in b0 data from Matlab output
    
    # Load .mat file
    matb0 = loadmat(f'VERDICT outputs/{PatNum}/Model 1/vb0.mat')['vb0']
    
    # Permute matlab b0 volume
    matb0 = np.moveaxis( matb0 , -1, 0)  
    

    # == Calculate MSE for different relative orientations
    MSE0 = mean_squared_error(matb0, b0)
    MSE1 = mean_squared_error(matb0, b0[::-1, : , :])
    
    # If MSE0 > MSE1, python b3000 image and mask need to be flipped to match fIC
    if MSE0 > MSE1:
        LesionMask = LesionMask[::-1,:,:]
        b0 = b0[::-1,:,:]
    else:
        None
        
    
    # == Save lesion mask
    try:
        os.makedirs(f'{output_path}/{PatNum}')
    except:
        None
      
    # Save as npy  
    np.save(f'{output_path}/{PatNum}/{ROIName}.npy', LesionMask)
    
    # Save as mha
    sitk.WriteImage( sitk.GetImageFromArray(LesionMask), f'{output_path}/{PatNum}/{ROIName}.mha' )
    
    # == Save b=0 from b3000 image
    
    # Save as npy  
    np.save(f'{output_path}/{PatNum}/b0.npy', b0)
    
    # Save as mha
    sitk.WriteImage( sitk.GetImageFromArray(b0), f'{output_path}/{PatNum}/b0.mha' )     
 
 
saveROImask('BAR_038', 'L1_b3000_NT')

# Function for saving Nifti masks as npy arrays
def saveNiftimask(
    PatNum,
    Nifti_INNOVATE_ROIs_path = r"D:\UCL PhD Imaging Data\Nifti INNOVATE ROIs",
    output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs"):
    
    # Read Nifti image
    mask = sitk.GetArrayFromImage(
        sitk.ReadImage(f'{Nifti_INNOVATE_ROIs_path}/{PatNum}.nii.gz')
    )
    
    print(mask.shape)
    plt.figure()
    plt.imshow(mask[4])
    plt.show()
    
           
# Function for extracting fIC values from ROI
def extractROIfICs(
    PatNum,
    ROIName,
    ModelNum,
    VERDICT_output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\VERDICT outputs",
    ROI_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\ROIs",
    output_path = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\fIC results",
):
    
    '''
    Function to extract fIC values from within an ROI
    
    '''                
    
    # Load fIC array
    fIC = loadmat(f'{VERDICT_output_path}/{PatNum}/Model {ModelNum}/fIC.mat')['fIC']
    
    # remove infinities
    fIC[fIC == np.inf] = 0
    
    # Remove nan
    fIC[np.isnan(fIC)] = 0
    
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
        
    # Save as npy
    np.save(f'{output_path}/{PatNum}/{ROIName}/Model {ModelNum}.npy', ROI_fIC)
    
    # Save fIC image as mha
    sitk.WriteImage(sitk.GetImageFromArray(fIC), f'VERDICT outputs/{PatNum}/Model {ModelNum}/fIC.mha')
    

      
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
    PatNums, 
    ROIName,
    ModelNum,
    avg_type = 'mean',
    results_path = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\fIC results'
    
):
    
    # Extract list of fIC results filenames
    fIC_fnames = glob.glob(f'{results_path}/*/{ROIName}/Model {ModelNum}.npy')
    
    # For each file, extract patient number and calculate average fIC
    avgfICs = []
    thesePatNums = []
    
    for fname in fIC_fnames:
        
        # Find patient number
        PatNum = os.path.split( 
                               os.path.split(
                                   os.path.split(fname)[0]
                               )[0]
                               )[1]

        # Check if patient number in PatNums
        if PatNum not in PatNums:
            print(PatNum)
            continue
        else:
            thesePatNums.append(PatNum)
            
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
    fIC_DF = pd.DataFrame({'Patient ID': thesePatNums, 'Avg fIC': avgfICs})

    print(fIC_DF)
    # Create directory
    try:
        os.makedirs(f'{results_path}/Average fIC Dataframes/Model {ModelNum}')
    except:
        None
        
    # Save dataframe as pickle
    with open(f'{results_path}/Average fIC Dataframes/Model {ModelNum}/average_fIC_df.pickle', 'wb') as handle:
        pickle.dump(fIC_DF, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

        
    
def fIC_ROC(
    ModelNum,
    results_path = r'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission\fIC results'
    ):
    
    '''
    Python function to generate ROC curve for lesion classification from 
    a specified model type
    
    General methods:
    
    1. Read in average fIC and biopsy results dataframes
    2. Create corresponding arrays of average fIC and biopsy outcomes (significant or insignificant)
    3. Use sklearn to generate ROC curve
    
    '''
    
    # Read in average fIC dataframe
    with open(f'{results_path}/Average fIC Dataframes/Model {ModelNum}/average_fIC_df.pickle', 'rb') as handle:
        fIC_DF = pickle.load(handle)
        

    # Read in biopsy results dataframe
    BiopsyResults_DF = readBiopsyResults()

    # Construct list of common patients (in Biopsy DF and fIC DF)
    fIC_PatList = fIC_DF['Patient ID'].values
    Biopsy_PatList = BiopsyResults_DF['Patient ID'].values
    
    PatList = fIC_PatList[ np.array([Pat in Biopsy_PatList for Pat in fIC_PatList]) ]
    
    # Construct arrays for biopsy results and fIC
    BiopsyResults = []
    fICs = []
    
    for PatNum in PatList:
        # Extract fIC
        fIC_Bools = (fIC_DF['Patient ID'].values == PatNum)
        fICs.append(fIC_DF['Avg fIC'].values[fIC_Bools][0])
        # Extract biopsy result
        Biopsy_Bools = (BiopsyResults_DF['Patient ID'].values == PatNum)
        BiopsyResults.append(BiopsyResults_DF['Biopsy Result'].values[Biopsy_Bools][0])
        
    # Make arrays
    fICs = np.asarray(fICs)    
    BiopsyResults = np.asarray(BiopsyResults)
    print(fICs, BiopsyResults)
    
    # Create ROC curve
    fpr, tpr, thresholds = roc_curve(y_true = BiopsyResults, y_score = fICs)
    
    # Calculate roc_auc_score
    roc_auc = roc_auc_score(y_true = BiopsyResults, y_score = fICs)
    
    return fpr, tpr, thresholds, roc_auc



    
   
   
# # for PatNum in ['BAR_003', 'BAR_004', 'BAR_005', 'BAR_006', 'BAR_009', 'BAR_033', 'INN_019']:
# #     saveROImask(PatNum, 'L1')
# #     extractROIfICs(PatNum, ROIName = 'L1', ModelType = 'Original')
        
# avgROIfICs('L1', 'Original')
# print( fIC_ROC('Original') )