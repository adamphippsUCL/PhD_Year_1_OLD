function explore_PhilipsPrivate
origf = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/ClemenceExampleMinimalSeries/DICOM/IM_0005' ;

outfolder = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/testdicoms' ;

d = joindictprivate;
dicomdict('set',d) 

dinfo = dicominfo(origf);
[X] = dicomread(origf) ;

dfn = fullfile(outfolder,'direct_copy.dcm') ;

dinfo = rmfield(dinfo,{'DiffusionGradientOrientation', ...
    'PositionReferenceIndicator', ...
    'TransmitCoilName', ...
    'TriggerTime', ...  
    'SpacingBetweenSlices', ... 
    'SliceLocation', ...  % checked
    'LowRRValue', ...
    'HighRRValue', ...
    'IntervalsAcquired', ...
    'IntervalsRejected', ...
    'HeartRate', ...
    'TriggerWindow', ...
    'RequestingService', ...
    'RequestedProcedureDescription', ...
    'RequestedContrastAgent', ...
    'PerformedStationAETitle', ...
    'PerformedStationName', ...
    'PerformedLocation', ...
    'AcquisitionDate', ...
    'AcquisitionTime', ...
    'PregnancyStatus', ...
    'RequestingPhysician', ...
    'ReferringPhysicianName', ...
    'CodeValue', ...
    'CodingSchemeDesignator', ...
    'CodeMeaning'}) ;

dinfo.StudyInstanceUID = dicomuid ;
dinfo.SeriesInstanceUID = dicomuid;
dinfo.FrameOfReferenceUID = dicomuid ;
dinfo.SOPInstanceUID = dicomuid ;
dinfo.InstanceCreatorUID = dicomuid ;
dinfo.ImplementationClassUID = dicomuid ;
dinfo.MediaStorageSOPInstanceUID = dicomuid ;
dinfo.ImageType= 'DERIVED\PRIMARY\DIFFUSION\FIC' ;
dinfo.AcquisitionDateTime = '20230201094950.55000' ;
dinfo.NumberOfPhaseEncodingSteps = 168 ;
dinfo.EchoTrainLength = 127;
dinfo.PercentSampling = 100.6 ;
dinfo.PercentPhaseFieldOfView= 100;
dinfo.PixelBandwidth= 1676;
 dinfo.Manufacturer= 'Philips Healthcare'; % FFS - has to be Philips Medical Systems !! ?? 
dinfo.ProtocolName = 'testy2' ;

status = dicomwrite(X, dfn, dinfo, ...
    'CreateMode','Copy', ...
    'MultiframeSingleFile', false, ...
    'UseMetadataBitDepths', false, ...  % OK if false
    'WritePrivate', true, ...
    'VR', 'explicit') 

disp(['Written: ',dfn])