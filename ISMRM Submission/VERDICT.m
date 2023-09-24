
function VERDICT(pat_num, model_type, opts)

% MATLAB function to carry out VERDICT processing 

% INNOVATE patient number and VERDICT model type specified as inputs

% Model type can be: 'Original' or 'No VASC'

arguments
    pat_num % INNOVATE patient number
    model_type % Type of VERDICT model to fit
    opts.INNOVATE_path = "D:\UCL PhD Imaging Data\INNOVATE\" % Path to INNOVATE imaging folder
    opts.parent_folder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\PhD_Year_1\ISMRM Submission" % Parent folder to save outputs
end

% Define DICOM folder path
DICOM_path = join( [opts.INNOVATE_path pat_num "\scans"], "");

% Define output folder
output_path = join([opts.parent_folder "/VERDICT outputs/" pat_num "/" model_type], "");

if ~exist(output_path, "dir")
   mkdir(output_path)
end


% Define options for model type

% Original VERDICT model
if strcmp(model_type, 'Original')
    
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path) ...
        );

% No VASC compartment model
elseif strcmp(model_type, 'No VASC')

    ncompart = 1;
    excludebvals = [90];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        excludebvals = excludebvals);

% Reduced Radii in fitting    
elseif strcmp(model_type, 'No VASC Reduced Rs')

    ncompart = 1;
    excludebvals = [90];
    Rs = linspace(7.5, 12.5, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        Rs = Rs);

else
    disp('Incorrect model specification')
    exit()
end


% Save results
save([convertStringsToChars(output_path) '/scheme.mat'], 'scheme')
save([convertStringsToChars(output_path) '/Y.mat'], 'Y')   
save([convertStringsToChars(output_path) '/fIC.mat'], 'fIC')
save([convertStringsToChars(output_path) '/fEES.mat'], 'fEES')
save([convertStringsToChars(output_path) '/fVASC.mat'], 'fVASC')
save([convertStringsToChars(output_path) '/R.mat'], 'R')

end
