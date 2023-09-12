function grappa_basic
% GRAPPA_BASIC Demo for reconstruction of GRAPPA acquired data.
% 
% Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
% See function grappa_detail.m for an example showing more of the 
% workings and functionality of the SIRF code.
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file from a GRAPPA MRI acquisition in the ISMRMRD format.
%    Example GRAPPA datasets:
%    a) 'meas_MID00108_FID57249_test_2D_2x.dat' is 
%       available from https://www.ccppetmr.ac.uk/downloads
%       This is in the manufacturer's raw data format and needs to be
%       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
%       This executable is installed on the Virtual Machine.
%
%    b) An ISMRMRD h5 file can be simulated using the MATLAB function
%    gen_us_data.m . This will reconstruct faster than the real data in a).
%
% Usage:
%  grappa_basic
%
%
% Adapted by David Atkinson (D.Atkinson@ucl.ac.uk) from original 
% code by Evgueni Ovtchinnikov
%
% See also GRAPPA_DETAIL GEN_US_DATA


% load the SIRF mutilities library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

% This demo uses gadgetron as the MR reconstruction engine. 
% Load the SIRF library and import the SIRF class.
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

import mGadgetron.*


try 
    % Get the filename of the input ISMRMRD h5 file
    disp('Select ISMRMRD H5 file')
    [fn,pn] = uigetfile('*','Select ISMRMRD H5 file') ;
    filein = fullfile(pn,fn) ;
    
    % Load this ISMRMRD h5 file
    input_data = AcquisitionData(filein);
    
    % Pre-process this input data using three preparation gadgets 
    % from gadgetron
    prep_gadgets = [{'NoiseAdjustGadget'}, ...
                    {'AsymmetricEchoAdjustROGadget'} ...
                    {'RemoveROOversamplingGadget'} ];
    
    preprocessed = input_data.process(prep_gadgets);
    
    
    % Perform reconstruction of the preprocessed data.
    % 1. set the reconstruction to be for Cartesian GRAPPA data.
    recon = GenericCartesianGRAPPAReconstruction();
    
    % 2. set the reconstruction input to be the data we just preprocessed.
    recon.set_input(preprocessed);
    
    % 3. run (i.e. 'process') the reconstruction.
    fprintf('---\n reconstructing...\n');
    recon.process();
    
    % Extract an image object from the reconstruction and convert this 
    % to a MATLAB array.
    image_obj = recon.get_output('image');
    idata = image_obj.as_array();  % returns a complex array
    
    sl = 5 ; % Number of the slice to be displayed.
    if size(idata,3) < sl
        sl = 1 ;
    end
    
    % Display the modulus and phase for this reconstructed slice.
    figure('Name',['idata, slice: ',num2str(sl)])
    subplot(1,2,1), imshow(abs(idata(:,:,sl)),[]), title('Abs')
    subplot(1,2,2), imshow(angle(idata(:,:,sl)),[-pi pi]), title('Phase')
    

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
