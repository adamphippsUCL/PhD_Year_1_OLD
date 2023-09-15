% MATLAB script to calculate fIC maps for patients from the INNOVATE
% clinical trial and then save the results.

% Author: Adam Phipps

%=====================================================================

% First, specify patient number 
pat_num = 'INN_129';

% Define DICOM folder path
DICOM_path = ['D:\UCL PhD Imaging Data\INNOVATE\' pat_num '\scans'];

% Define VERDICT model type
model_type = 'Reduced Rs';

% Create output folder
ISMRM_folder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission';

output_path = [ISMRM_folder '/VERDICT outputs/' pat_num '/' model_type];

if ~exist(output_path, 'dir')
   mkdir(output_path)
end

% Define options for model type
if strcmp(model_type, 'Original')
    ncompart = 2;
    excludebvals = [];
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam(DICOM_path, output_path);


elseif strcmp(model_type, 'No VASC')
    ncompart = 1;
    excludebvals = [90];
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam(DICOM_path, output_path, ncompart = ncompart, excludebvals = excludebvals);

elseif strcmp(model_type, 'Reduced Rs')
    ncompart = 2;
    excludebvals = [];
    Rs = 3:2:17;
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam(DICOM_path, output_path, Rs = Rs);

else
    disp('Incorrect model specification')
    exit()
end


% Save results
save([output_path '/scheme.mat'], 'scheme')
save([output_path '/Y.mat'], 'Y')   
save([output_path '/fIC.mat'], 'fIC')
save([output_path '/fEES.mat'], 'fEES')
save([output_path '/fVASC.mat'], 'fVASC')
save([output_path '/R.mat'], 'R')
