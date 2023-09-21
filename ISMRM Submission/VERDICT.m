
function VERDICT(pat_num, model_type, opts)

% MATLAB function to carry out VERDICT processing 

% INNOVATE patient number and VERDICT model type specified as inputs

% Model type can be: 'Original' or 'No VASC'

arguments
    pat_num % INNOVATE patient number
    model_type % Type of VERDICT model to fit
    opts.INNOVATE_path = 'D:\UCL PhD Imaging Data\INNOVATE\' % Path to INNOVATE imaging folder
    opts.parent_folder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission' % Parent folder to save outputs
end

% Define DICOM folder path
DICOM_path = [opts.INNOVATE_path pat_num '\scans'];

% Define output folder
output_path = [opts.parent_folder '/VERDICT outputs/' pat_num '/' model_type];

if ~exist(output_path, 'dir')
   mkdir(output_path)
end


% Define options for model type

% Original VERDICT model
if strcmp(model_type, 'Original')

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam(DICOM_path, output_path);

% No VASC compartment model
elseif strcmp(model_type, 'No VASC')
    ncompart = 1;
    excludebvals = [90];
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam(DICOM_path, output_path, ncompart = ncompart, excludebvals = excludebvals);

% Reduced Radii in fitting    
elseif strcmp(model_type, 'Reduced Rs')
    Rs = 5:2:15;
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

end
