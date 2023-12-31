function verdict(dfolder, output_folder, opts)
% VERDICT A wrapper function for VERDICT processing. Intended either
% for MATLAB directly, or deployed as a Docker container.
%
% See the help in verdict_process for more details of the processing, 
% algorithms and options. The help here is mainly regarding the calling 
% syntax.
%
% MATLAB
% verdict(dfolder, output_folder, Name, Value, ...)
%
% COMMAND LINE (e.g. Docker)
% Name and Value are treated as character strings
% verdict 'inputfolder' 'outputfolder' solver SPAMS  addb0ToScheme true
%
% Examples
%  % This Example will add outputs to the folder dfolder, the same top-level 
%  folder as the input. This top-level could be the study
%  dfolder = '/User/myname/data/STUDYID-005' ;
%  verdict(dfolder, dfolder, outputDicom=true)
%
%  verdict(dfolder, tempdir, solver='SPAMS')    % MATLAB mode only
%  verdict(dfolder, tempdir, 'solver', 'SPAMS') 
%
% David Atkinson
%
% See also verdict_process bv2scheme verdict_fit sviewer getSeriesVERDICT 

% Note input arguments will be characters when called from command line in
% deployed mode hence the defaults below. verdict_process takes care of 
% this e.g. converting the string 'true' to the boolean true

arguments
    dfolder
    output_folder
    opts.outputDicom   = 'true' % output results in DICOM files
    opts.addPhilipsPrivate = 'true' % adds Private fields to Philips if available (for scanner upload)
    opts.outputMatFile = 'true'  % output variables in one .mat file
    opts.outputAInReport='false' % Output AMICO matrix A in report
    opts.register      = 'true'  % input data will have be registered
    opts.swapinvXNAT   = 'false' % swap in the data from a vXNAT.mat file
    opts.swapinvMAT     = 'false' % swap in MAT file for Guanda's project
    opts.usedirecdiff  = 'false' % use the 3 directional diffusion images
    opts.solver        = 'lsqnonnegTikonhov' % solver used
    opts.forceXNATscheme = 'false' % ignores scanner and uses XNAT scheme file
    opts.addb0ToScheme   = 'true'  % adds b=0 (signal of 1) to every series
    opts.quiet           = 'false' % suppress figures and report viewer
    opts.maskhwmm        = '48'    % mask half-width in mm used for registration
    opts.resultsFolderName = ''    % If empty will default to res-datetime
    opts.fICReconNumber  = '5'   % DICOM fIC recon number for saved fIC (specified in testing)
    opts.vBaseSeriesNumber  % DICOM Base Series number for saved fIC (specified in testing), default is from b=90
    % Final fIC Series Number will be 100*vBaseSeriesNumber + fICReconNumber
    opts.bvSortDirection = 'descend' % to correspond to XNAT pipeline
    opts.vADCbmax = '1600' % max b-value used in VERDICT ADC calculation
    opts.allowedSeriesNumbers = [] % SeriesNumbers in this set can be used
    opts.forcedSchemeName = ''     % Force scheme name for debugginh
    opts.ncompart = '2'            % Nuber of tissue compartment beyond sphere (1 for EES, 2 for EES and VASC)
end

% Set Release/Version tag here (will be output in report and .mat file)
opts.verdictVersion = '1.010 - R20230827' ;

makeDOMCompilable() % Needed for deployed Report Generator

verdict_process(dfolder, output_folder, opts) 

