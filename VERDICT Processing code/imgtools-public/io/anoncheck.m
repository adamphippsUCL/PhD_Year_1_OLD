function anoncheck(dirn)
% ANONCHECK Check the anonymisation of DICOM data
%
% anoncheck
% anoncheck(folder)
%
% Checks all DICOM files in folder and sub-folders. 
%
% FOR EACH FOLDER:
% Displays folder name.
%
% Checks for differences between values in specified fields.
%
% Checks PatientID is present and within the PatientName and Filename.
%   ie PatientName should include PatientID (which will be a studyID)
%   If PatientID is '1020102', a PatientName of 'ReImagine_1020102' is OK.
%
% Checks Deidentification has been performed (DicomCleaner run).
% 
% Checks Deidentification methods that were used.
% 
% Checks for presence of other 'dangerous' fields e.g. EthnicGroup.
% 
% NOTE most checks stop after warning from one file to prevent an overload of
% messages. This means later files in a folder might have unchecked issues.
%
% Dates that are 20000101 are not printed.
%
% Copyright 2021. David Atkinson, University College London
%


if nargin < 1 || ~exist(dirn,'dir') 
    dirn = pref_uigetdir('anoncheck','dirn','','Select folder') ;
end

[fpathA, this_dir, ext] = fileparts(dirn) ;
% Directory can have  a dot in the name which looks like an extension
% hence:
this_dir = [this_dir ext] ;


% Determine DICOM files in same folder to be checked
dlist = dir(dirn) ;
ifl = 0 ;
for id = 1:length(dlist)
    if ~dlist(id).isdir
        fn = fullfile(dlist(id).folder, dlist(id).name) ;
        if isdicom(fn)
            ifl = ifl + 1 ;
            flist{ifl} = fn ;
        end
    else
        % folder, call anoncheck recursively
        if ~strcmp(dlist(id).name,'.') && ~strcmp(dlist(id).name,'..')
             anoncheck(fullfile(dlist(id).folder, dlist(id).name) ) ;
        end
    end
end

% % If there were no DICOM files, return (recursive calls)
if ifl == 0
    return
end

% check each folder below calling anoncheck
% % recursively
% if ifl == 0
%     for id = 1:length(dlist)
%         if dlist(id).isdir && ~strcmp(dlist(id).name,'.') && ~strcmp(dlist(id).name,'..')
%             anoncheck(fullfile(dlist(id).folder, dlist(id).name) ) ;
%         end
%     end
%     return
% end

disp(' ')
disp(['- Will check ',num2str(length(flist)),' files in folder: ',this_dir])



% Fields for visual check that are unequivocal
fvis = {...
    'PatientName', ...
    'PatientID', ...
    'PatientBirthDate', ...
    'PatientAddress'} ;

% Aim is for the following list to be study specific. These are also displayed if
% present
study_vis = { ...
    'PatientAge', ...
    'PatientWeight', ...
    'PatientSex', ...
    'StudyDate', ...
    'SeriesDate', ...	
    'AcquisitionDate', ...
    'ContentDate', ...	
    'AccessionNumber', ...
    'StudyID', ...
    'PerformedProcedureStepStartDate', ...
    'IssueDateofImagingServiceRequest', ...
    'AcquisitionDateTime', ...
    'InstanceCreationDate', ...
    'InstitutionName', ...
    'InstitutionalDepartmentName' ...
    'InstitutionAddress', ...
    'ReferringPhysicianName', ...
    'ReferringPhysicianAddress', ...
    'ReferringPhysicianTelephoneNumbers', ...
    'ConsultingPhysicianName', ...
    'PerformingPhysicianName', ...
    'NameOfPhysiciansReadingStudy', ...
    } ;
    
 % the required DICOM deidentification methods (taken from DicomCleaner) 
 % Adapt for requirements
req_dmethods = { ...
    'Deidentified', ...
    'Descriptors removed except series', ...
    'Patient Characteristics removed', ...
    'Device identity removed', ...
    'Institution identity removed', ...
    'Dates modified', ...
    'Structured content removed', ...
    'Unsafe private removed', ...
    'UIDs remapped' } ;



% 'Dangerous' fields. Warning issued if present with a non-blank value
fnp = { ...
    'PatientBirthTime', ...
    'PatientBirthDateInAlternativeCalendar', ...
    'PatientInsurancePlanCodeSequence', ...
    'OtherPatientIDs', ...
    'OtherPatientNames', ...
    'OtherPatientIDsSequence', ...
    'PatientBirthName', ...
    'PatientMotherBirthName', ...
    'MilitaryRank', ...
    'MedicalAlerts', ...
    'Allergies', ...
    'AdmittingDiagnosesDescription', ...
    'CountryOfResidence', ...
    'RegionOfResidence', ...
    'PatientTelephoneNumbers', ...
    'PatientTelecomInformation', ...
    'EthnicGroup', ...
    'Occupation', ...	
    'SmokingStatus', ...
    'AdditionalPatientHistory', ...	
    'PregnancyStatus', ...
    'LastMenstrualDate', ...
    'PatientReligiousPreference', ...
    'PatientComments', ...
    'ImageComments', ...
    'OperatorsName'} ;

% Use first file in list as reference file
dinfo = dicominfo(flist{1}, 'UseDictionaryVR', true) ;


% Check all the files in flist
wdeident = false ; % Checks stopped after warning to prevent too many messages
widmatch = false ;
wfnp = false ;     % 'field not present'

dref = dinfo ; % dref has fields added if they are not in reference file
dref.cfviswarned = false ;

for ifile = 1:length(flist)
    d2 = dicominfo(flist{ifile}, 'UseDictionaryVR', true) ;
    % cfvis2 compares the relevant structures in two files and reports
    % differences 
    if dref.cfviswarned == false
        dref = cfvis2(fvis, dref, d2)    ;    % Name, DOB, ID, Address
        dref = cfvis2(study_vis, dref, d2)  ; % study fields
    end
    if ~widmatch
      widmatch = cidmatch(d2)   ;  % Check PatientID is present and 
    end                            % in the FamilyName and file name
    if ~wdeident
      wdeident = cdeident(req_dmethods, d2) ; % De-identification used
    end
    if ~wfnp
        wfnp = cfnp(fnp, d2)  ;              % Other dangerous fields
    end
end
cfvis(fvis, dref)       % Output values of fields in ref
cfvis(study_vis, dref)

end

function widmatch = cidmatch(d2)
% CIDMATCH Check PatientID is present and in the Family Name and File name

widmatch = false ;
if ~isfield(d2,'PatientID')
    ID = [] ;
else
    ID = d2.PatientID ;
end
if isempty(ID)
    fprintf('<strong>!! PatientID is empty or missing</strong>\n')
    widmatch = true ;
end

if ~isempty(ID)
    if ~contains(d2.PatientName.FamilyName, ID)
        fprintf('<strong>!! Patient ID is not in Patient Name</strong>\n')
        widmatch = true ;
    end
    
    if ~contains(d2.Filename, ID)
        fprintf('<strong>!! Patient ID is not in folder name</strong>\n')
        widmatch = true ;
    end
end


end

function wdeident = cdeident(req_dmethods, dinfo)
% CDEIDNT Check Deidentification is recorded in DICOM file
%
% DicomCleaner records its processing
%

wdeident = false ;
% Check that de-identification has been performed with required steps
if ~isfield(dinfo,'PatientIdentityRemoved') || ~strcmp(dinfo.PatientIdentityRemoved,'YES')
    fprintf('<strong>!! There is no record of Patient Identity Removed</strong>\n')
    wdeident = true ;
else
    val = dinfo.DeidentificationMethod ;
    aa = textscan(val,'%s','Delimiter','\') ;
    steps = aa{1} ;
    
    for imeth = 1: length(req_dmethods)
        locmeth = strcmp(steps,req_dmethods{imeth}) ;
        if max(locmeth) < 1
            fprintf('<strong>!! Required de-identification method NOT performed: %s</strong>\n', req_dmethods{imeth})
            wdeident = true ;
        end
    end
end

end

function wfnp = cfnp(fnp, dinfo)
% CFNP Warn if dangerous fields are present with non-blank values
%

wfnp = false ;

for ifield = 1:length(fnp)
    if isfield(dinfo,fnp{ifield})
        % skip if blank or 'Not Stated / Unk'
        val = dinfo.(fnp{ifield}) ;
        if ~isempty(val)
            if ischar(val) && ~strcmp(val,'Not Stated / Unk')
                fprintf('<strong>!! %s present with value: %s</strong>\n',fnp{ifield},val)
                wfnp = true ;
            end
        end
    end
end


end

function dref = cfvis2(fvis, dref, dfile)

for ifield = 1: length(fvis)
    fd = fvis{ifield} ;
    
    if isfield(dfile, fd) && isfield(dref, fd)
        % check they are equal
        if ~isequal(dref.(fd), dfile.(fd))
            disp(['!! Unequal values for ',fd,': '])
            disp(dref.(fd))
            disp(dfile.(fd))
            dref.cfviswarned = true ;
        end
        
    elseif isfield(dfile, fd) && ~isfield(dref, fd)
        % put into reference list
        dref.(fd) = dfile.(fd) ;
    elseif ~isfield(dfile,fd) && isfield(dref, fd)
        % do nothing
    else
        % do nothing
    end
    
end


end
function cfvis(fvis, dinfo)
% CFVIS Print field values for visual check
% Does not print Dates that are 20000101

for ifield = 1: length(fvis)
    if isfield(dinfo,fvis{ifield}) % PatientAddress not always present
        fv = dinfo.(fvis{ifield}) ; % field value
        
        if isstruct(fv)
            fields = fieldnames(fv) ;
            for isubf = 1: length(fields)
                val = fv.(fields{isubf}) ;
                if ~isempty(val)
                    fprintf('%33s : %s\n',[fvis{ifield},' ',fields{isubf}],val)
                end
            end
        elseif isnumeric(fv)
            fprintf('%33s : %s\n', fvis{ifield},num2str(fv))
        else
            if ~isempty(fv)
                fieldn = fvis{ifield} ;
                if length(fieldn)>4 && strcmp(fieldn(end-3:end),'Date') && ...
                        strcmp(fv,'20000101')
                    % don't print dates that are 20000101
                else
                    fprintf('%33s : %s\n',fvis{ifield},fv)
                end
            end
        end
    end
    
    
    
end
end

    
    
    
    

