function [dinfo_out, vSeriesNumbers, vBV, dinfoT2ax, T2axSeriesNumbers, vMAT, dscheme]  = getSeriesVERDICT(dataList, opts)
% getSeriesVERDICT Get VERDICT DICOM series, DICOM summary info and
% b-values (Classic or Enhanced). Also finds  T2 axial scan (Enhanced only).
%
% VERDICT scans identified as those Series with exactly 2 b-values.
%
% T2 axial are identified as those with ImageType ORIGINAL\PRIMARY\T2\NONE
% and approximately axial in orientation. Only tested for Enhanced DICOM
% T2s. Ignores a Luminal Water scan if present.
%
% [dinfo, vSeriesNumbers, vBV, dinfoT2ax, T2axSeriesNumbers, vMAT, dscheme]  =  ...
%                                getSeriesVERDICT(dataList, Name, Value)
% 
% dataList either a single file, a single folder, or, a cell array of 
% files and/or folders. DICOM files can be Classic MR or Enhanced MR
%
% Name, Value
% useSubfolders {true} | false  . Parses all subfolders below those in
% dataList
% nExpectedSeries {5} Number of VERDICT series expected
% opts.allowedSeriesNumbers {[]} Series Numbers that are allowed foe
% VERDICT
%
% dinfo structure with all VERDICT series included
% vSeriesNumbers  VERDICT scans' series numbers, typically [5 1]
% vBV             VERDICT scans' b-values, typically [5 2]
%  The above outputs are in the order of series number 
%
% dinfoT2ax  dinfo structure for T2ax scan(s)
% T2axSeriesNumbers  The Series number(s) of the T2 ax scans found.
% vMAT  A volume saved as either vXNAT.mat from the deniftify app, or vMAT.mat.
%       Expected format for vMAT is [row col slice nbv 2] where 
%       vb0 = vMAT(:,:,:,:,1) and vbv = vMAT(:,:,:,:,2)
% dscheme cell array where dscheme{SeriesNumber} = [small_delta big_delta]
%  if a Philips Data Storage DICOM file is in the folder tree.
%
% Copyright 2023  David Atkinson
%
% See also VERDICT_main dfparse

arguments
    dataList 
    opts.useSubfolders logical = true
    opts.nExpectedSeries = 5
    opts.allowedSeriesNumbers = []
    opts.excludebvals = []
end

vMAT = [] ; ivMAT = 0;

if ~iscell(dataList)
    dataList = {dataList} ;
end

ndataList = length(dataList) ;
ifn = 0 ;

for idata = 1:ndataList
    if isfolder(dataList{idata})
        if opts.useSubfolders
            dstruc = dir([dataList{idata},filesep,'**',filesep,'*']) ;
        else
            dstruc = dir(dataList{idata}) ; 
        end  

        for idstruc = 1:numel(dstruc)
            if dstruc(idstruc).isdir
                % ignore (if sub-folders were wanted, they will be in dstruc)
            elseif dstruc(idstruc).bytes>132
                ffn = fullfile(dstruc(idstruc).folder, dstruc(idstruc).name) ;
                if isdicom(ffn)
                    ifn=ifn+1;
                    fns{ifn}  = ffn ;
                end
                if strcmp(dstruc(idstruc).name, 'vXNAT.mat') || ...
                        strcmp(dstruc(idstruc).name, 'vMAT.mat')
                    vMATstruc = load(ffn) ;
                    if isfield(vMATstruc,'vXNAT')
                        vMAT = vMATstruc.vXNAT ;
                    else
                        vMAT = vMATstruc.vMAT ;
                    end
                    vMAT(isnan(vMAT)) = 0 ;
                    ivMAT = ivMAT + 1;
                end
            end
        end
    elseif exist(dataList{idata},'file')
        if isdicom(dataList{idata})
            ifn=ifn+1;
            fns{ifn}  = dataList{idata} ;
        end
    end
end

if ivMAT > 1
    warning(['Should be 0 or 1 vXNAT.mat or vMAT file, found: ',num2str(ivMAT)])
end

if ~exist('fns','var')
    warning('No files found in input.')
    disp(['Causes include: incorrect file or path names, network or external drive not available.'])
    
    dinfo_out=[]; vSeriesNumbers=[]; vBV=[]; dinfoT2ax=[]; T2axSeriesNumbers=[]; vMAT=[];
    return
end

% Now have a list of Classic or Enhanced DICOM files in fns
% Loop over the files to fill dinfo
dinfo = [] ;
dinfoT2ax = [] ;

UseDictionaryVR = [] ;

dscheme = {} ;

for ifn = 1:length(fns)
    [dinfo_this, UseDictionaryVR] = dfparse(fns{ifn}, false, UseDictionaryVR) ;
    if isfield(dinfo_this,'DiffusionBValue')
        if length(dinfo_this) > 1
            if isempty(opts.allowedSeriesNumbers) || ismember(dinfo_this(1).SeriesNumber , opts.allowedSeriesNumbers)
                % Enhanced DICOM
                nbv = length(unique([dinfo_this.DiffusionBValue])) ;
                if nbv == 2
                    dinfo = [dinfo dinfo_this] ;
                elseif nbv == 3
                    if contains(dinfo_this(1).ImageType,'DERIVED')
                        % skip
                    else
                        dinfo = [dinfo dinfo_this] ;
                    end
                end
            end
        else
            % Classic DICOM file
            if contains(dinfo_this.ImageType,'DERIVED\PRIMARY\DIFFUSION\ADC') || ...
                contains(dinfo_this.ImageType,'DERIVED\PRIMARY\DIFFUSION\EADC') || ...
                contains(dinfo_this.ImageType,'ORIGINAL\PRIMARY\M_FFE\M\FFE') || ...
                dinfo_this.RepetitionTime < 1500 || ...
                dinfo_this.EchoTime < 30
                % skip
            else
                dinfo = [dinfo dinfo_this] ;
            end
        end
    end  

    % Look for T2W axial also 
    if length(dinfo_this) < 40 % not an Enhanced DICOM LWI (has more frames)

        if isfield(dinfo_this(1),'ImageType') && contains(dinfo_this(1).ImageType,'ORIGINAL\PRIMARY\T2\NONE')
            IOP = dinfo_this(1).ImageOrientationPatient ;
            planeNorm = cross(IOP(1:3),IOP(4:6)) ;
            if abs(dot(planeNorm,[0 0 1])) > 0.8 % axial only
                dinfoT2ax = [dinfoT2ax  dinfo_this] ;
            end
        end
    end

    
    xxinfo = XXread(fns{ifn},'suppressNotXXWarning',true,'outFileInEditor',false) ;
    if isfield(xxinfo,'IF_DIFF_deltas')
        deltas = xxinfo.IF_DIFF_deltas{:} ;
        if ~isempty(deltas)
            deltas_cell = textscan(deltas, '%f/%f') ;
            sdelta = deltas_cell{2} ;
            BDELTA = deltas_cell{1} ;
            dscheme{xxinfo.SeriesNumber} = [sdelta BDELTA] ;
        end
    end
        
    

end

if ~isempty(dinfoT2ax)
    dinfoT2ax = data_process(dinfoT2ax) ;
    snT2ax = [dinfoT2ax.SeriesNumber] ;
    T2axSeriesNumbers = unique(snT2ax) ;
else
    T2axSeriesNumbers = [];
end

dinfo = process_duplicates(dinfo) ;

dinfo = data_process(dinfo) ;

% Parse to get Series Numbers and b-values for VERDICT scans 

sn = [dinfo.SeriesNumber] ;
bv = [dinfo.DiffusionBValue] ;

usn = unique(sn) ;
% disp([num2str(usn),' unique Series Numbers read in.'])

iSeries = 0 ;
vSeriesNumbers=[]; vBV=[] ; dinfo_out = [] ;

for iusn = 1:length(usn)
    if isempty(opts.allowedSeriesNumbers) || ismember(usn(iusn) , opts.allowedSeriesNumbers)
        loc_sn = sn==usn(iusn) ; % locations in dinfo of this series number
        ubv_this = unique(bv(loc_sn)) ; % unique b-values in this series

        % CODE ADDED TO EXCLUDE B VALUE
        if sum( ismember(ubv_this, opts.excludebvals))
            disp(['VERDICT series' usn(iusn) 'with b values ' num2str(ubv_this) ' removed'])
            continue
        end

        disp(['Series ',num2str(usn(iusn)),' contains b-values: ',num2str(ubv_this)])

        if length(ubv_this) == 2
            % assumed to be a VERDICT series
            iSeries = iSeries + 1 ;
            vSeriesNumbers(iSeries) = usn(iusn) ;
            vBV(iSeries,:) =  ubv_this ;
            dinfo_out = [dinfo_out dinfo(loc_sn)] ;
        elseif length(ubv_this) == 3
            iSeries = iSeries + 1 ;
            vSeriesNumbers(iSeries) = usn(iusn) ;
            vBV(iSeries,:) =  ubv_this([1 2])  ;
            dinfo_out = [dinfo_out dinfo(loc_sn)] ;

            iSeries = iSeries + 1 ;
            vSeriesNumbers(iSeries) = usn(iusn) ;
            vBV(iSeries,:) =  ubv_this([1 3]) ;
            
        end
    end
end

if isempty(vBV)
    warning('MATLAB:getSeriesVERDICT:noVERDICTFiles', ...
        'No VERDICT files found')
end

if length(vSeriesNumbers) ~= opts.nExpectedSeries 
    warning('MATLAB:getSeriesVERDICT:unexpectedNumberVERDICTFiles', ...
        [num2str(length(vSeriesNumbers)),' VERDICT series found, expected ', ...
         num2str(opts.nExpectedSeries)])
end

end

function dinfo = process_duplicates(dinfo)
% Process dinfo for duplicate scans
% 
% Needs modifying to put output in report

if isfield(dinfo,'StudyInstanceUID')
    uStudyUID = unique({dinfo.StudyInstanceUID}) ;

    nuStudy = length(uStudyUID) ;
    if nuStudy > 1
        warning(['Detected ',num2str(nuStudy),' studies.'])

        % In example seen, the whole study was duplicated with the same
        % Series appearing twice. Check number of files for each study,
        % keep the largest or first.
        %
        % Alternative example seen in XNAT with Enhanced DICOM with the
        % whole study saved as 2 studies, some VERDICT in one, and some in
        % another.

        for iStudy = 1:nuStudy
            cscans(iStudy) = sum(contains({dinfo.StudyInstanceUID},uStudyUID{iStudy})) ;
            disp(['Study ',num2str(iStudy),' has ',num2str(cscans(iStudy)),' entries'])
        end

        if sum(cscans) < 400
            % use all
            disp(['Using all studies'])
        else
            [mx,imx]=max(cscans) ;
            disp(['Using Study # ',num2str(imx(1)),' with ',num2str(mx),' images.'])

            tokeep = contains({dinfo.StudyInstanceUID},uStudyUID{imx(1)}) ;

            dinfo = dinfo(tokeep) ;
        end
    end
end

% Below on;y useful after filtering for VERDICT scans

% % Check for scan repeated due to forgetting "O images in DB". Unlike above, 
% % scans will have separate Series
% 
% [sn, id, isn] = unique([dinfo.SeriesNumber]) ;
% mostFrequent = mode(isn) ;
% % Keep Series that have the most common number of entries (assumes we are
% % discarding a lower number with incorrect numbers).
% 
% infoToKeep = find(isn==mostFrequent) ;
% if length(infoToKeep) < length(dinfo)
%     warning(['Keeping ',num2str(length(infoToKeep)),' images out of ',...
%         num2str(length(dinfo)),' read.'])
% 
%     dinfo = dinfo(infoToKeep) ;
% end





end
