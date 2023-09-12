function [dinfoout, dictJoined] = addPhilipsPrivate(dinfoout, locout, baseSeries, reconNum)
% ADDPHILIPSPRIVATE Adds private fields to enable DICOM import on Philips scanner
% Called from verdict_process before writeDicom
%
%  [dinfoout, dictJoined] = addPhilipsPrivate(dinfoout, locout, baseSeries, reconNum)
%
%  On input, dinfoout may contain private fields labelled 'Private...'
%  Remove these, add the required ones with their Philips' names (for
%  clarity and to ensure the VR and VM are known).
% 
% Philips use SeriesNumbers such as 601, 602 etc where the '6' is the scan
% number and the '1', '2' etc successive reconstructions. Here baseSeries
% would, for example, be 6 and reconNum 2.
%
% David Atkinson, University College London
%
% See Also verdict_process DICOMDICT writeDicom joindictprivate
%

fields = fieldnames(dinfoout) ;

for ifield = 1:length(fields)
    this_field = fields{ifield} ;
    if length(this_field) >= 7 && strcmp(this_field(1:7),'Private')
        % Private field in input - delete it.
        dinfoout = rmfield(dinfoout,this_field) ;
    end
end

DICOMManufacturer = 'Philips Medical Systems' ;

% Now switch to new dictionary, read in info, copy across desired fields.

[dictJoined, dictOrig] = joindictprivate ;
dicomdict("set",dictJoined) ;

for iloc = 1:length(locout)
    this_loc = locout(iloc) ;

    dinwp = dicominfo(dinfoout(this_loc).Filename) ; % dinfofull in with private

    dinfoout(this_loc).EchoNumbers = 1; 
    dinfoout(this_loc).NumberOfTemporalPositions = 1 ;
    dinfoout(this_loc).PresentationLUTShape = 'IDENTITY' ;
    dinfoout(this_loc).StudyDescription = '' ;
    dinfoout(this_loc).ImagePlaneOrientation = 'TRANSVERSAL' ;
    dinfoout(this_loc).ScanningSequence = 'SE' ;
    dinfoout(this_loc).SequenceVariant = 'SK' ;
    dinfoout(this_loc).VolumetricProperties = 'VOLUME' ;
    dinfoout(this_loc).ImagedNucleus = '1H' ;
    dinfoout(this_loc).Manufacturer = DICOMManufacturer ;
    dinfoout(this_loc).SpacingBetweenSlices = dinfoout(this_loc).slc2c ;


    % These also need to be added to mf2single_keep in writeDicom.m
    dinfoout = copyfield(dinfoout, this_loc, 'Private_2001_10xx_Creator', dinwp) ;
   
    dinfoout(this_loc).MRImagePhaseNumber = 1 ; % (2001,1008)
    
    dinfoout = copyfield(dinfoout, this_loc, 'ImagePlaneNumber', dinwp) ; % (2001,100a)
    %dinfoout = copyfield(dinfoout, this_loc, 'ImagePlaneOrientation', dinwp) ;

    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesNrOfEchoes', dinwp) ; % (2001,1014)
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesNrOfPhases', dinwp) ; % (2001,1017)
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesNrOfSlices', dinwp) ; % (2001,1018)
    
    dinfoout(this_loc).MRSeriesReconstructionNumber = reconNum ;  % (2001,101d)
    
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesScanningTechniqueDesc', dinwp) ; % (2001,1020)
    dinfoout = copyfield(dinfoout, this_loc, 'Stack', dinwp) ; % (2001,105f) This is a structure including offsets, angulations etc
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesNrOfStacks', dinwp) ; % (2001,1060)

    dinfoout(this_loc).MRSeriesAcquisitionNumber = baseSeries ;  % (2001,107B)
    
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesNrOfDynamicScans', dinwp) ; % (2001,1081)

    dinfoout = copyfield(dinfoout, this_loc, 'Private_2005_10xx_Creator', dinwp) ;

    dinfoout = copyfield(dinfoout, this_loc, 'MRImageTypeMR', dinwp) ; % (2005,1011)
    dinfoout = copyfield(dinfoout, this_loc, 'MRMeasurementScanResolution', dinwp) ; % (2005,101d)
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesScanDuration', dinwp) ; % (2005,1033) 
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesDataType', dinwp) ; % (2005,1035)
    dinfoout(this_loc).MRImageChemicalShiftNumber = 0 ;                % (2005,1040)
    dinfoout = copyfield(dinfoout, this_loc, 'MRImageScanningSequencePrivate', dinwp) ; % (2005,106e)
    dinfoout = copyfield(dinfoout, this_loc, 'MRSeriesGeometryCorrection', dinwp) ; % (2005,10a9)

end

dicomdict("set",dictOrig) ;

end

% - - - - - - - - - - 

function dout = copyfield(din, loc, field, dinwp)
dout = din ;

if isfield(dinwp,field)
    dout(loc).(field) = dinwp.(field) ;
else
    if ~isfield(din, field) % already there from dfparse
        warning([field,' is not a field'])
    end
end

end