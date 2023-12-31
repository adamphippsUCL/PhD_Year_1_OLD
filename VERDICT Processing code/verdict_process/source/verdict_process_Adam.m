function [scheme, Y,fIC, fEES, fVASC, R, rmse] = verdict_process_Adam(dfolder, output_folder, opts)
% VERDICT_PROCESS VERDICT processing 
%
% See wrapping function VERDICT for how to call this function
%
% Processing of VERDICT MRI data. Fits a ball, sphere and astrosticks
% tissue model to diffusion data from 5 scans. Outputs include fIC map.
% The acquisition scheme is currently hard-coded for the scanners listed in
% bv2scheme and any scheme file present is ignored. Results may differ 
% if the assumed scheme was not used for the actual measurement.
%
% verdict_process(dfolder, output_folder, Name, Value, ...)
%
% Loads data from dfolder (Enhanced or Classic DICOM files in any order).
% Other DICOMS may be present but there should be exactly one dataset for
% each of the 5 VERDICT scans (b=90, 500, 1500, 2000, 3000 each with a
% b=0).
% Optional registration (translations only, based on central region of
% central slice, b=0 to b=0 alignment with transformation applied to DW data).
% Fit of model to data. Model uses scheme (b-value, G, deltas)
% Default solver is lsqnonneg with Tikonhov regularisation. 
% Output of maps.
% Colormaps are experimental and not tested.
%
% Uses DICOM ManufacturerModelName to choose a scheme unless forceXNATscheme 
% is true. ("XNAT" here refers to reference processing used for Singh et al
% Radiology 2022)
%
% Outputs a PDF report. Logging previously used the Advanced Logger 
% for MATLAB  (available on FileExchange / GitHub) but removed to aid with
% deployment to Docker
%
%
% When called from the command line (e.g. in Docker), Name and Value must
% be character strings
%
% The divergent colourmap comes from colorcet:
% Reference:
% Peter Kovesi. Good Colour Maps: How to Design Them.
% arXiv:1509.03700 [cs.GR] 2015
% https://arxiv.org/abs/1509.03700
%
% Examples
%  % This Example will add DICOMs and other output to the folder dfolder
%  dfolder = '/User/myname/data/STUDYID-005' ;
%  verdict_process(dfolder, dfolder, outputDicom=true)
%
%  verdict_process(dfolder, tempdir, solver='SPAMS')    % MATLAB mode only
%  verdict_process(dfolder, tempdir, 'solver', 'SPAMS') 
%
% As a command line:
% verdict_process 'inputfolder' 'outputfolder' solver SPAMS  addb0ToScheme true
%
% Configurations
% --------------
% XNAT results closest (swapping in volumes from XNAT registrations)
% swapinvXNAT=true, register=false, usedirecdiff=true, solver='SPAMS', forceXNATscheme=true, addb0ToScheme=true 
% XNATmimic = {'register','false','swapinvXNAT','true','usedirecdiff','true','solver','SPAMS','forceXNATscheme','true','addb0ToScheme','true'}
%
% As above, but correct scanner (scheme) with b=0 added.
% swapinvXNAT=true; register=false; usedirecdiff=true; solver='SPAMS'; forceXNATscheme=false; addb0ToScheme=true ;
%
% XNAT registered input, correct scheme (with added b=0), lsqnoninTokonohv reg
% swapinvXNAT=true; register=false; usedirecdiff=true; solver='lsqnonnegTikonhov'; forceXNATscheme=false; addb0ToScheme=true ;
%
% MATLAB process (the default)
% swapinvXNAT=false, register=true, usedirecdiff=false, solver='lsqnonnegTikonhov', forceXNATscheme=false, addb0ToScheme=true 
%
% Swapping in a file vMAT.mat. Should be a file vMAT.mat containing the
% variable vMAT such that vb0 = vMAT(:,:,:,:,1)  and vbv = vMAT(:,:,:,:,2)
% corresponding to ascending b-values 90,500,1500,2000,3000
% swapinvMAT=true, usedirecdiff=false
%
%
% David Atkinson
%
% See also bv2scheme verdict_fit  sviewer  getSeriesVERDICT verdict

% input arguments will be characters when called from command line in
% deployed mode. Converted below to boolean or numerical if necessary

% The arguments block was moved from here to the wrapper function verdict to 
% allow verdict.m and pcode for verdict_process

opts = convertCharsToLogical(opts) ; 
if ischar(opts.maskhwmm)
    opts.maskhwmm = str2double(opts.maskhwmm) ;
end
if isfield(opts,'vBaseSeriesNumber') && ischar(opts.vBaseSeriesNumber)
    opts.vBaseSeriesNumber = str2double(opts.vBaseSeriesNumber) ;
end
if isfield(opts,'fICReconNumber') && ischar(opts.fICReconNumber)
    opts.fICReconNumber = str2double(opts.fICReconNumber) ;
end
if isfield(opts,'vADCbmax') && ischar(opts.vADCbmax)
    opts.vADCbmax = str2double(opts.vADCbmax) ;
end
if isfield(opts,'allowedSeriesNumbers') && ischar(opts.allowedSeriesNumbers)
    opts.allowedSeriesNumbers = str2double(opts.allowedSeriesNumbers) ;
end
if isfield(opts,'ncompart') && ischar(opts.ncompart)
    opts.ncompart = str2double(opts.ncompart) ;
end


% Set timestamp tstr for output folder naming 
tnow = datetime ;
tstr_nice = char(tnow) ;
tnow.Format='yyyy-MM-dd''T''HHmmss' ; % will form part of path name (no weird characters)
tstr = char(tnow) ;

% Create subfolder for results (name using time string)
if ~isfield(opts,'resultsFolderName') || isempty(opts.resultsFolderName)
    opts.resultsFolderName = ['res-',tstr];
end
resultsFolder = fullfile(output_folder, opts.resultsFolderName) ;


% % Set up Logging 
% logFileName = 'verdict.log' ;
% logFFN = fullfile(resultsFolder, logFileName) ;
% logObj = mlog.Logger("verdict_log") ;
% logObj.LogFolder = resultsFolder ;
% logObj.LogFile = logFFN ; 
% 
% logObj.info(['Timestamp: ',tstr_nice,' (',tstr,')'])
% % 
% % 
% % [status,msg,msgID] = mkdir(output_folder) ;
% % if status == 0
% %     warning(msgID, msg)
% %     return
% % end
% % 
% % [status,msg,msgID] = mkdir(resultsFolder) ;
% % if status == 0
% %     warning(msgID, msg)
% %     return
% % end
% % 
% % if opts.outputDicom == true
% %     [status,msg,msgID] = mkdir(fullfile(resultsFolder,'DICOM')) ;
% %     if status == 0
% %         warning(msgID, msg)
% %         return
% %     end
% % end

% % 
% % % Set up Reporting
% % rptFileName = 'rptverdict.pdf' ;
% % rptFFN = fullfile(resultsFolder, rptFileName) ;
% % 
% % import mlreportgen.report.*
% % import mlreportgen.dom.*
% % 
% % rpt    = Report(rptFFN,'pdf');
% % rpt.Layout.Landscape = true;
% % append(rpt,TableOfContents)
% % 
% % if isdeployed 
% %     append(rpt, Paragraph(['Code is running deployed. If Docker, the folder volume names in this report will have been mapped to a local folder ']))
% % end

verdictVersion = opts.verdictVersion ;
% % 
% % append(rpt, Paragraph(['This report was generated: ', tstr_nice]) )
% % append(rpt, Paragraph(' '))
% % append(rpt, Paragraph(['VERDICT version: ', verdictVersion ]))
% % append(rpt, Paragraph(' '))
% % append(rpt, Paragraph(['This report file name is: ', rptFFN]))
% % append(rpt, Paragraph(' '))
% % append(rpt, Paragraph(['The results folder is: ',resultsFolder]))

% % append(rpt,Chapter('Inputs'))

if opts.quiet
    figVisible = 'off' ; % figures do not display, but are still in report
else
    figVisible = 'on' ;
end
% % 
% % 
% % if opts.swapinvXNAT
% %     if opts.register == true || opts.usedirecdiff == false
% %         msg = 'Cannot register or use iso diff when swapinvXNAT is true';
% %         warning(msg);
% %         append(rpt, Paragraph(msg)) ;
% %         append(rpt, Paragraph(' '))
% %     end
% % end
% % 
% % if opts.swapinvMAT
% %     if opts.swapinvXNAT
% %         msg = 'Cannot have both vMAT and vXNAT swap in.';
% %         warning(msg);
% %         append(rpt, Paragraph(msg)) ;
% %         append(rpt, Paragraph(' '))
% %     end
% %     if opts.usedirecdiff == true
% %         msg = 'vMAT does not have directional diffusion, cannot usedirecdiff.';
% %         warning(msg);
% %         append(rpt, Paragraph(msg)) ;
% %         append(rpt, Paragraph(' '))
% %     end
% % end

if opts.Rs
    Rs = opts.Rs;
else
    Rs = linspace(0.01,15.1,18) ; % radii used in fit.
end
% VERDICT AMICO paper used 0.01:15.1 
% In online AMICO linspace(0.01,20.1,20);
% https://github.com/daducci/AMICO_matlab/blob/master/models/AMICO_VERDICTPROSTATE.m
% Exact choice seems to make little difference to fIC.

rmseThreshold = 0.05 ; % Root Mean Square Error above this coloured rmseRGB
rmseRGB = [0.8 0.8 0.8] ; % [0.8  0  0.8] is magenta-like

fw_pix = 1000 ; % width of larger figures in pixels


switch opts.usedirecdiff
    case true
        nbd = 3 ; % diffusion weighted images per file
    case false
        nbd = 1 ; % isotropic diffusion weighted image
    otherwise
        warning('Unknown value of usedirecdiff')
end

if opts.addb0ToScheme
    naddb0 = 1 ; % number of b=0 added per series 
else
    naddb0 = 0 ;
end

if opts.forceXNATscheme == true
    if opts.addb0ToScheme == false
        warning('forceXNATscheme hence adding b=0 scans to scheme')
    end
    naddb0 = 1 ; % number of b=0 added per series 
end


% % append(rpt, Paragraph('Input Options'))
% % append(rpt, Paragraph(' '))
% % append(rpt, Paragraph(['Number of diffusion images used per Series (nbd): ',num2str(nbd)]))
% % append(rpt, Paragraph(['Number of b=0 added per Series: ',num2str(naddb0)]))
% % append(rpt, Paragraph(' '))

% % Report parameters in a table
% optsout = opts;
% optsout = rmfield(optsout,'allowedSeriesNumbers') ;
% pTObj = MATLABTable(struct2table(optsout, 'AsArray',true)) ;
% pTObj.Style = [pTObj.Style { FontSize('8pt')} ] ;
% append(rpt, pTObj)
% 
% append(rpt,Paragraph(' '))
% append(rpt,Paragraph(['Input data folder: ',dfolder]))
% append(rpt,Paragraph(' '))


dicomdict("factory")

[dinfo, vSeriesNumbers, vBV, dinfoT2ax, T2axSeriesNumbers, vMAT, dscheme] = getSeriesVERDICT(dfolder, allowedSeriesNumbers = opts.allowedSeriesNumbers, ...
    excludebvals = opts.series_excludebvals) ;





if isempty(dinfo)
% %     append(rpt,Paragraph(' '))
% %     msg = 'NO INPUT FOUND. (Check file and path names / connection to networked or external drives)';
% %     append(rpt,Paragraph(msg))

    disp('Exiting as no input found.')
    return
end


% % % Report exam details
% % append(rpt,Paragraph(' '))
% % append(rpt,Paragraph('DICOM details'))

dfi = dicominfo(dinfo(1).Filename) ;
% % if isfield(dfi,'PatientID')
% %     append(rpt,Paragraph(['PatientID: ',dfi.PatientID]))
% % end
% % if isfield(dfi,'StudyDescription')
% %     append(rpt,Paragraph(['StudyDescription: ',dfi.StudyDescription]))
% % end
% % if isfield(dfi,'InstitutionName')
% %     append(rpt,Paragraph(['InstitutionName: ',dfi.InstitutionName]))
% % end
% % if isfield(dfi,'ProtocolName')
% %     append(rpt,Paragraph(['ProtocolName: ',dfi.ProtocolName]))
% % end
% % if isfield(dfi,'Manufacturer')
% %     append(rpt,Paragraph(['Manufacturer: ',dfi.Manufacturer]))
% % end
% % if isfield(dfi,'ManufacturerModelName')
% %     append(rpt,Paragraph(['ManufacturerModelName: ',dfi.ManufacturerModelName]))
% % end
% % if isfield(dfi,'SoftwareVersions')
% %     append(rpt,Paragraph(['SoftwareVersions: ',dfi.SoftwareVersions]))
% % end
% % append(rpt,Paragraph(' '))

% % % Check dinfo for reasonable parameter values
% % % Add to report and/or logging
% % pcheck = checkV(dinfo) ;
% % append(rpt, pcheck)

[sortedBV, indsortedBV] = sort(vBV(:,2), opts.bvSortDirection) ;

nSeries = length(vSeriesNumbers) ;

iptprefsOrig = iptgetpref ;
iptsetpref('ImshowBorder','tight')
iptsetpref('ImshowInitialMagnification','fit') 

if opts.register
    regTranslations = zeros([nSeries 2]) ;  % Registration transformations stored for report
    [optimizer,metric] = imregconfig('monomodal') ; % monomodal as b=0 to b=0
    fw = fw_pix ;    % figure width in pixels
    fh = fw/nSeries*3 ;
% %     hfreg_check = figure(Name='Pre and Post Registration', Position=[250 500 fw round(fh)],Visible=figVisible);
% %     treg = tiledlayout(3,nSeries,'TileSpacing','none','Padding','tight') ;
% %     axrs = [] ;
end

% Pre-allocate variables used for checking data
medianvb0 = zeros([1 nSeries]) ; medianvbv = zeros([1 nSeries]) ; plotbv = zeros([1 nSeries]) ;
TE = zeros([1 nSeries]) ; TR = zeros([1 nSeries]) ;

nDeltaFromDicom = 0 ;

for iSeries = 1: nSeries
    sn_this = vSeriesNumbers(indsortedBV(iSeries)) ;
    bv_this = vBV(indsortedBV(iSeries),2) ;

    % Get b=0 data
    [vb0, mb0, b0loc] = d2mat(dinfo,{'slice','bv','series'},'bv',0, ...
        'series',sn_this,'op','fp') ;



    TR(indsortedBV(iSeries)) = dinfo(b0loc(1)).RepetitionTime ;
    if isfield(dinfo,'EffectiveEchoTime')
        TE(indsortedBV(iSeries)) = dinfo(b0loc(1)).EffectiveEchoTime ;
    end

    % Get b>0 data
    if opts.usedirecdiff
        [vbv, mbv, bvloc] = d2mat(dinfo,{'slice','bdirec','ddty','series'}, ...
            'ddty',1,'series',sn_this,'op','fp') ;
        if size(vbv,4)~=3, warning('Expected 3 diffusion directions'), end
    else
        [vbv, mbv, bvloc] = d2mat(dinfo,{'slice','ddty','series','bv'}, ...
            'series',sn_this, 'bv', bv_this, ...
            'ddty', 2, 'op','fp') ;
    end

    if opts.swapinvXNAT
        % in the reference XNAT, b-value ordering is 3000 to 90
        switch bv_this
            case 90
                vb0 = vMAT(:,:,:,17) ;
                vbv = vMAT(:,:,:,18:20) ;
            case 500
                vb0 = vMAT(:,:,:,13) ;
                vbv = vMAT(:,:,:,14:16) ;
            case 1500
                vb0 = vMAT(:,:,:,9) ;
                vbv = vMAT(:,:,:,10:12) ;
            case 2000
                vb0 = vMAT(:,:,:,5) ;
                vbv = vMAT(:,:,:,6:8) ;
            case 3000
                vb0 = vMAT(:,:,:,1) ;
                vbv = vMAT(:,:,:,2:4) ;
            otherwise
                error(['Unrecognised bv_this: ',num2str(bv_this)])
        end
    end

    if opts.swapinvMAT
        % vMAT 
        switch bv_this
            case 90
                vb0 = vMAT(:,:,:,1,1) ;
                vbv = vMAT(:,:,:,1,2) ;
            case 500
                vb0 = vMAT(:,:,:,2,1) ;
                vbv = vMAT(:,:,:,2,2) ;
            case 1500
                vb0 = vMAT(:,:,:,3,1) ;
                vbv = vMAT(:,:,:,3,2) ;
            case 2000
                vb0 = vMAT(:,:,:,4,1) ;
                vbv = vMAT(:,:,:,4,2) ;
            case 3000
                vb0 = vMAT(:,:,:,5,1) ;
                vbv = vMAT(:,:,:,5,2) ;
            otherwise
                error(['Unrecognised bv_this: ',num2str(bv_this)])
        end
    end

    % compute medians to check signal
    vexcend = vb0(:,:,3:end-2) ;
    medianvb0(iSeries) = median(vexcend(:),'omitnan') ;
    vexcend = vbv(:,:,3:end-2,:) ;
    medianvbv(iSeries) = median(vexcend(:),'omitnan') ;
    plotbv(iSeries) = bv_this ;

    if iSeries == 1 

        %% ADAM SAVE b0 from b3000
    
        % Save vb0 so image flips can be dealt with!!!
        save([convertStringsToChars(output_folder) '/vb0.mat'], 'vb0')

        %%

        % first pass
        Y = zeros([size(vb0,[1 2 3]), (nbd + naddb0)*nSeries]) ;
        vb0tot = zeros(size(vb0));
        v2000 =  zeros(size(vb0)) ;
        v2000norm = zeros(size(v2000)) ;

        % these are output in a .mat file for testing/debugging
        vprereg  = zeros([size(vb0,[1 2 3]), (nbd + 1)*nSeries]) ;
        vpostreg = zeros([size(vb0,[1 2 3]), (nbd + 1)*nSeries]) ;
        

        % Find pixels in mask aroud image centre where prostate is assumed
        % to be located (within +/- maskhwmm) mm of centre.
        maskhw = floor(opts.maskhwmm/mb0.geom(1).PixelSpacing_HW(1)) ;
        [ny, nx, nz] = size(vb0,[1 2 3]) ;

        maskcentre = [ceil( (ny+1)/2 )  ceil((nx+1)/2) ] ;

        maskc = { max(1,maskcentre(1)-maskhw) : min(ny,maskcentre(1)+maskhw) , ...
            max(1,maskcentre(2)-maskhw) : min(nx,maskcentre(2)+maskhw) } ;

        reg_slice = ceil((nz+1)/2) ; % registration uses only one slice per series
% %         append(rpt, Paragraph(['Registration slice number: ',num2str(reg_slice)])) ;


        vfixed = vb0(maskc{:},reg_slice) ; % fixed image for registration
        vfixedUpper = prctile(vfixed(:),98) ;
        
        dfullinf = dicominfo(dinfo(1).Filename) ;
        if isfield(dfullinf,'ManufacturerModelName')
            scanner = dfullinf.ManufacturerModelName ;
% %             append(rpt, Paragraph(['scanner:', scanner])) ;
        else
            warning('Scanner (ManufacturerModelName) not known')
% %             append(rpt, Paragraph('Scanner (ManufacturerModelName) not known'))
        end

        if opts.forceXNATscheme
% %             append(rpt,Paragraph('! Forcing use of XNAT scheme file.'))
            scanner = 'XNAT' ;
        end

        if isfield(dfullinf,'BodyPartExamined')
            BodyPartExamined = dfullinf.BodyPartExamined ;
            if strcmp(BodyPartExamined,'KIDNEY' )
                scanner = [scanner,'Renal'] ;
% %                 append(rpt, Paragraph(['RENAL: scanner switched to:', scanner])) ;
            end
        else
            BodyPartExamined = 'UNKNOWN' ;
        end

        if ~isempty(opts.forcedSchemeName)
            scanner = opts.forcedSchemeName ;
% %             append(rpt, Paragraph(['FORCED SCHEME: scanner set to:', scanner])) ;
        end



% %         % Checks if a scheme file is present and reports only
% %         % Note code still uses bv2scheme hard-coded values
% %         [numFiles, fileNames] = findVerdictSchemeFiles(dfolder) ;
% %         if numFiles == 1
% %             [version, tableData] = readVerdictSchemeFile(fileNames{1}) ;
% %             [meetsRequirements, messg] = checkVerdictSchemeTable(tableData, scanner) ;
% %             if meetsRequirements
% %                 mRstr = ['Scheme file meets requirements. Version: ',num2str(version)] ;
% %                 append(rpt,Paragraph(mRstr))
% %             else
% %                 append(rpt,Paragraph(['Scheme file present but does not meet requirements: ', messg]))
% %             end
% %         elseif numFiles == 0
% %             append(rpt,Paragraph('No scheme file present.'))
% %         else
% %             append(rpt,Paragraph('More than one scheme file was found.'))
% %         end

        geomout = mb0.geom ; % used as reference and for writing DICOMs
        dinfoout = dinfo ; % for DICOM output
        locout = b0loc ;
        if ~isfield(opts,'vBaseSeriesNumber')
            if sn_this > 200
                vbsn = floor(sn_this/100) ;
            else
                vbsn = sn_this ;
            end

            opts.vBaseSeriesNumber = vbsn ;
        end
    end % iSeries == 1

% %     % Check data is consistent and report
% %     pcheck = checkV(geomout, mb0.geom) ;
% %     append(rpt, pcheck)
% %     pcheck = checkV(geomout, mbv.geom) ;
% %     append(rpt, pcheck)
% %     pcheck = checkV(vb0) ;
% %     append(rpt,pcheck)
% %     pcheck = checkV(vbv) ;
% %     append(rpt,pcheck)
% %     if ~strcmp(BodyPartExamined,'KIDNEY')
% %         pcheck = checkV(mb0.geom,'axial') ;
% %         append(rpt, pcheck)
% %     end
   

    if opts.register
        vmoving = vb0(maskc{:},reg_slice) ;
        ylinePos = ceil(size(vmoving,1)/2) ;
        xlinePos = ceil(size(vmoving,2)/2) ;
        nxlim = size(vmoving,2) ;
        nylim = size(vmoving,1) ;

        % different b=0's do not always have the same scale - 
        % normalising first seems more robust for registration

        vmoving_toreg = mat2gray(vmoving,[0 double(prctile(vmoving(:),98))]) ;
        vfixed_toreg  = mat2gray(vfixed, [0 double(prctile(vfixed(:),98))]) ;

        tform = imregtform(vmoving_toreg,vfixed_toreg,'translation',optimizer,metric);

        regTranslations(iSeries,:) = tform.Translation ; % store for report

        vb0_reg = zeros(size(vb0)) ;
        vbv_reg = zeros(size(vbv)) ;

        for islice =1:size(vb0,3)
            vb0_reg(:,:,islice) = imwarp(vb0(:,:,islice),tform,"OutputView",imref2d(size(vb0(:,:,1)))) ;
            for ibd = 1:nbd
                vbv_reg(:,:,islice,ibd) = imwarp(vbv(:,:,islice,ibd), tform,"OutputView",imref2d(size(vbv(:,:,1)))) ;
            end
        end

        % Figure to QA registration. Top row preserves scaling
        % Middle row is moving images
        % Bottom row is registered.
% %         axr = nexttile(treg,tilenum(treg,1,iSeries)) ;
% %         imshow(vmoving,[0 vfixedUpper],'Parent',axr) ;
% %         axrs = [axrs axr];
% % 
% %         axr = nexttile(treg,tilenum(treg,2,iSeries)) ;
% %         imshow(vmoving,[ ],'Parent',axr), hold on
% %         plot(axr,[1 nxlim],[ylinePos ylinePos],'Color','y')
% %         plot(axr,[xlinePos xlinePos],[1 nylim],'Color','y')
% %         axrs = [axrs axr];
% %         axr = nexttile(treg,tilenum(treg,3,iSeries)) ;
% %         imshow(vb0_reg(maskc{:},reg_slice),[ ],'Parent',axr), hold on
% %         plot(axr,[1 nxlim],[ylinePos ylinePos],'Color','y')
% %         plot(axr,[xlinePos xlinePos],[1 nylim],'Color','y')
% %         axrs = [axrs axr];
    else
        vb0_reg = vb0 ;
        vbv_reg = vbv ;
    end

    if abs(bv_this-2000) < 10
        v2000 = sum(vbv_reg,4) ;
        locb2000  = bvloc ;
        v2000norm = v2000 ./ vb0_reg ;
    end

    vb0tot = vb0tot + vb0_reg ;

    vnorm = vbv_reg./repmat(vb0_reg,[1 1 1 nbd]) ; % DW images, normalised by b=0
    bvinY_this = repmat(bv_this,[1 nbd]) ;


    if naddb0 == 1
        vnorm = cat(4,ones(size(vb0_reg)), vnorm ) ;
        bvinY_this = cat(2,0,bvinY_this) ;
    end

    bvinY(1, 1+(iSeries-1)*(nbd + naddb0) : iSeries*(nbd + naddb0)) = bvinY_this ;
    Y(:,:,:,1+(iSeries-1)*(nbd + naddb0) : iSeries*(nbd + naddb0)) = vnorm ; % Stack DW images in 4th dim

    % dscheme form getSeriesVERDICT is cell array with non-empty entries
    % only at SeriesNumbers where there is an XX file with IF_delta_Delta
    % (as SeriesNumbers typically jump in stesp of 100, this is largely
    % empty)
    if ~isempty(dscheme) && length(dscheme{sn_this}) > 1
        bvs = bv2scheme(bv_this, 'PDS', dscheme{sn_this}) ;
        nDeltaFromDicom  = nDeltaFromDicom  + 1 ;
    else
        bvs = bv2scheme(bv_this, scanner) ; 
    end
    scheme(1+(iSeries-1)*(nbd+naddb0)+1 : iSeries*(nbd+naddb0)) = bvs ;
    if naddb0 == 1
        scheme(1+(iSeries-1)*(nbd+naddb0)) = bv2scheme(0, scanner) ;
    end


    % for output/debugging
    vprereg(:,:,:,1+(iSeries-1)*(nbd + 1) : iSeries*(nbd + 1)) = cat(4,vb0,vbv) ;
    vpostreg(:,:,:,1+(iSeries-1)*(nbd + 1) : iSeries*(nbd + 1)) = cat(4,vb0_reg,vbv_reg) ;


end % iSeries

% % if opts.register
% %     linkaxes(axrs)
% % end
% % 
% % hfmedian = figure(Name='medians', Visible=figVisible) ;
% % plot(plotbv,medianvb0,'LineWidth',2,'DisplayName','Median B0'), hold on
% % plot(plotbv,medianvbv,'LineWidth',2,'DisplayName','Median vbv')
% % plot(plotbv,medianvbv./medianvb0*medianvb0(1),'LineWidth',2,'DisplayName','Normalised median vbv')
% % grid on
% % xlabel('b-value'), ylabel('median')
% % legend
% % 
% % figrpt = Figure(hfmedian);
% % figrpt.Snapshot.Caption = 'Medians excluding end slices' ;
% % figrpt.Snapshot.ScaleToFit = true ;
% % append(rpt,figrpt);

% Display inputs Y for QA. Also place in report.

Dref = 2e-9 ;  % A reference diffusivity for scaling the signals

% When "extra" b=0 are in Y, this leads to unnecessary blanks displayed so
% remove from display. These are added at the start of each 'b-value'
if naddb0 == 0
    row2Y4 = 1:size(Y,4) ;
else
    addedRows = 1:(nbd+naddb0):size(Y,4) ;
    row2Y4 = 1:size(Y,4) ;
    row2Y4(addedRows) = [] ;
end
nRow = length(row2Y4) ;

% % fw = fw_pix ;    % figure width in pixels
% % fh = fw/size(Y,3)*nRow ;
% % hfY = figure(Name='Y components',Position=[200 50 fw round(fh)],Visible=figVisible) ;
% % t=tiledlayout(nRow, size(Y,3),'TileSpacing','none','Padding','tight') ;
% % axs = [] ;
% % for iY3 = 1: size(Y,3)
% %     for iRow = 1: nRow
% %         iY4 = row2Y4(iRow) ;
% %         ax = nexttile(tilenum(t,iRow,iY3)) ;
% %         axs = [axs ax];
% %         upper = 1.2 * exp(-Dref*scheme(iY4).bval*1000) ;
% %         imshow(Y(:,:,iY3,iY4),[0 upper])
% %     end
% % end
% % linkaxes(axs)


% % % Report Generation
% % append(rpt, Paragraph(' '))
% % pobj = Paragraph(['Series Numbers and b-values in input.   Scanner: ',scanner, ...
% %     '. No. of delta/Delta read directly from DICOM: ',num2str(nDeltaFromDicom)]) ;
% % append(rpt, pobj)
% % 
% % tableStyle = { ...
% %     Width('100%'), ...
% %     Border('solid','black','1px'), ...
% %     ColSep('solid','black','1px'), ...
% %     RowSep('solid','black','1px') ...
% %     };
% % 
% % schemeReportOutput = scheme ;
% % if naddb0 ~= 0
% %     schemeReportOutput(addedRows) = [] ;
% % end
% % 
% % schemeTableObj = mlreportgen.dom.MATLABTable(struct2table(schemeReportOutput)) ;
% % schemeTableObj.Style = tableStyle ;
% % 
% % sT = table(vSeriesNumbers', vBV, TE', TR') ;
% % sT = sT(indsortedBV,:) ;  % put in b-value order
% % sT.Properties.VariableNames = {'Series', 'b-values', 'TE', 'TR'} ;
% % seriesTableObj = mlreportgen.dom.MATLABTable(sT) ;
% % seriesTableObj.Style = tableStyle ;
% % 
% % % lo_table is invisible - used for layout
% % lo_table = Table({seriesTableObj,' ', schemeTableObj}) ;
% % lo_table.entry(1,1).Style = {Width('3.2in')};
% % lo_table.entry(1,2).Style = {Width('.4in')};
% % lo_table.entry(1,3).Style = {Width('3.2in')};
% % lo_table.Style = {Width('100%'), ResizeToFitContents(false)};
% % 
% % append(rpt,lo_table) 
% % append(rpt,Paragraph(' '))

% % 
% % figrpt = mlreportgen.report.Figure(hfY);
% % figrpt.Snapshot.Caption = 'Entries in Y, the diffusion images normalized by their corresponding b=0' ;
% % figrpt.Snapshot.ScaleToFit = true ;
% % append(rpt,figrpt);

% % 
% % if opts.register
% %     hffixed = figure(Name='Reg Fixed',Visible=figVisible) ;
% %     imshow(vfixed,[0 vfixedUpper])
% %     figrpt = mlreportgen.report.Figure(hffixed);
% %     figrpt.Snapshot.Caption = 'The slice and the masked region treated as fixed for registration' ;
% %     append(rpt,figrpt);
% %     close(hffixed)
% % 
% %     append(rpt, Paragraph(' '))
% %     append(rpt, Paragraph('Translations found in registration (pixels)')) ;
% %     
% %     if any(abs(regTranslations) > 5, 'all')
% %         append(rpt,Paragraph('CHECK Registrations (large translation found)'))
% %     end
% % 
% %     regTableObj = mlreportgen.dom.MATLABTable(table(regTranslations)) ;
% %     append(rpt,regTableObj)
% % 
% %     figrpt = mlreportgen.report.Figure(hfreg_check);
% %     figrpt.Snapshot.Caption = 'Slice and mask region. Before registration, same windowing (top), scaled (middle) and after (bottom) registration' ;
% %     append(rpt,figrpt);
% % end

% % append(rpt,Chapter('Fitting'))


%% ADAM ADDITIONS

% == Mask for fitting

% Check if mask has been inputted
if exist('opts.mask', 'var')
    disp('Input mask');

else % Bounding box mask (Not sure why input isn't working)
    disp('no mask')
    opts.mask = zeros([size(vb0,[1 2 3])]) ;
    opts.mask(44:132,44:132,:) = 1;
end


% == Exclude b values from fitting!

for bval = opts.fitting_excludebvals
    
    disp(['b value ' num2str(bval) ' removed for fitting'])

    % Find scheme index for b value
    scheme_bools = ([scheme.bval] ~= bval);

    % Exclude b value from scheme
    scheme = scheme(scheme_bools);

    % Exclude b value image from Y
    Y = Y(:,:,:,scheme_bools);
    

end


%% Apply fitting

% ===== PERFORM FIT

[fIC, fEES, fVASC, R, rmse, A, tparams, vfopt] = verdict_fit(scheme,Y,Rs=Rs, ...
    solver=opts.solver, ncompart=opts.ncompart, mask = opts.mask) ;

% % if opts.outputAInReport
% %     append(rpt, Paragraph(' '))
% %     append(rpt, Paragraph(' Entries of AMICO system matrix A'))
% %     append(rpt, Paragraph(' '))
% %     ATableObj = mlreportgen.dom.MATLABTable(table(A)) ;
% %     ATableObj.Style = [ATableObj.Style {FontSize('8pt')}] ;
% %     append(rpt,ATableObj)
% % end
% % 
% % append(rpt, Paragraph(' '))
% % append(rpt, Paragraph(' Tissue Properties Used in Fit:'))
% % append(rpt, Paragraph(' '))
% % tissueObj = MATLABTable(struct2table(tparams)) ;
% % tissueObj.Style = [tissueObj.Style { FontSize('8pt')} ] ;
% % append(rpt, tissueObj)
% % 
% % append(rpt, Paragraph(' '))
% % pobj = Paragraph(['Solver used in verdict_fit was: ',vfopt.solver]);
% % append(rpt,pobj)
% % 
% % fsum = fIC + fEES + fVASC ;
% % 
% % m.geom = geomout ;
% % 
% % if anynan(R)
% %     append(rpt,Paragraph(' R contained NaNs')) ;
% %     R(isnan(R))=0 ;
% % end
% % 
% % nfIC = fIC./fsum ;
% % 
% % fQA = cat(4,fEES,fVASC, fsum, fIC, nfIC, rmse) ;
% % d4annotation = {'fEES','fVASC','fsum','fIC','nfIC','rmse'} ;
% % if ~opts.quiet
% %     sviewer(fQA,m,Name='fEES fVASC fsum, fIC, nfIC, rmse', CLim=[0.2 0.8], ...
% %     d4annotation=d4annotation, indexD4=4)
% % 
% %     sviewer(R,m)
% % end
% % 
% % % Calculate vADC, VERDICT ADC from b-values up to opts.vADCbmax
% % indb = find(bvinY <= opts.vADCbmax) ;
% % [vADC, vS0] = calcADC(Y(:,:,:,indb), bvinY(indb)) ;
% % 
% % 
% % % Builds a 5D RGB matrix out_rgb with fIC map, overlay, data.
% % %
% % % Colormap fIC using jet and with all fIC>1 clipped
% % % Set rmse> rmseThresh to magenta (not in jet color map)
% % 
% % nColor = 256 ;
% % cmap = jet(nColor) ;
% % out_rgb = zeros([size(fIC,[1 2 3]) 9 3]) ;
% % probe = zeros([size(fIC,[1 2 3]) 9]) ;
% % 
% % for islice = 1:size(fIC,3)
% %     % gray2ind clips fIC below 0 or above 1.
% %     ficrgb = ind2rgb(gray2ind(fIC(:,:,islice),nColor), cmap) ;
% %     % columnise to enable re-colouring of regions with poor fit rmse
% %     ficrgbcol = reshape(ficrgb,[size(ficrgb,1)*size(ficrgb,2) 3]) ;
% % 
% %     loc = find(rmse(:,:,islice) > rmseThreshold ) ;
% %     ficrgbcolnomask = ficrgbcol ;
% %     ficrgbcol(loc,:) = repmat(rmseRGB,[length(loc) 1])  ; % magenta
% %     ficrgb = reshape(ficrgbcol,[size(ficrgb,1) size(ficrgb,2) 3]) ;
% %     ficrgbnomask = reshape(ficrgbcolnomask,[size(ficrgb,1) size(ficrgb,2) 3]) ;
% %     out_rgb(:,:,islice,1,:) = reshape(ficrgb, [size(fIC,[1 2]) 1 1 3]) ;
% %     out_rgb(:,:,islice,2,:) = reshape(ficrgbnomask, [size(fIC,[1 2]) 1 1 3]) ;
% % end
% % d4annotation{1} = 'fIC masked' ;
% % d4annotation{2} = 'fIC' ;
% % probe(:,:,:,1) = fIC ;
% % probe(:,:,:,2) = fIC ;
% % 
% % imgb0 = mat2gray(vb0tot/nSeries, ...
% %     [0 double(prctile(vb0tot(:)/nSeries, 98))]) ;
% % 
% % imgb2000 = mat2gray(v2000, [0 double(max(v2000(:)))]) ;
% % imgb2000norm = mat2gray(v2000norm, [0 0.4]) ;
% % imgvADC = mat2gray(vADC*1000, [0 2]) ;
% % 
% % % fIC single colour overlay
% % overlay = fIC ;
% % overlay(overlay>1)=1;
% % overlay(overlay<0)=0 ;
% % 
% % % alphaOverlay = overlay ; % This should probably be more targeted at 0.42
% % xmid = 0.42 ;
% % grad = 2 ;
% % intc = xmid*(1-grad) ;
% % alphaOverlay = grad*overlay + intc ;
% % alphaOverlay(alphaOverlay<0)=0 ;
% % alphaOverlay(alphaOverlay>1)=1;
% % 
% % Rd = 1 ; G = 0.1 ; B = G ; % Slightly pink overlay colour
% % RBase = 0.9 ; GBase = RBase ; BBase = GBase ; % Base image 'colour'
% % 
% % overlayRGB = cat(4, Rd*overlay,G*ones(size(overlay)),B*ones(size(overlay))) ;
% %         
% % imgBaseRGB = cat(4, RBase*imgb0, GBase*imgb0, BBase*imgb0) ;
% % 
% % imgRGB = alphaOverlay.*overlayRGB + (1-alphaOverlay).*imgBaseRGB ;
% % 
% % 
% % imgb0RGB    = cat(5, imgb0, imgb0, imgb0) ;
% % imgb2000RGB = cat(5, imgb2000, imgb2000, imgb2000 ) ;
% % imgb2000normRGB = cat(5, imgb2000norm, imgb2000norm, imgb2000norm) ;
% % imgvADCRGB = cat(5, imgvADC, imgvADC, imgvADC) ;
% % 
% % % divergent
% % cmap_str = 'D09' ;
% % baseMap = colorcet(cmap_str) ;
% % % Skew colormap to put mid-white at fic=0.42 for divergent
% % whiteVal = 0.42 ;
% % 
% % cmv = linspace(0,1,size(baseMap,1)) ;
% % cmq = cmv/whiteVal*0.5 ;
% % 
% % extrapVal = baseMap(end,:) ;
% % 
% % redv   = interp1(cmv, baseMap(:,1), cmq,"linear",extrapVal(1)) ;
% % greenv = interp1(cmv, baseMap(:,2), cmq,"linear",extrapVal(2)) ;
% % bluev  = interp1(cmv, baseMap(:,3), cmq,"linear",extrapVal(3)) ;
% % 
% % divMap = cat(2, redv(:), greenv(:), bluev(:)) ;
% % 
% % for islice = 1:size(fIC,3)
% %     ficdivrgb = ind2rgb(gray2ind(fIC(:,:,islice),nColor), divMap) ;
% %     out_rgb(:,:,islice,7,:) =reshape(ficdivrgb,[size(fIC,[1 2]) 1 1 3]) ;
% % end
% % d4annotation{7}=cmap_str;
% % probe(:,:,:,7) = fIC;
% % 
% % % jet fIC overlay, constant alpha
% % calpha = 0.6 ;
% % out_rgb(:,:,:,3,:) = calpha*imgb0RGB + (1-calpha)*out_rgb(:,:,:,1,:) ;
% % d4annotation{3}='0.4*fIC and 0.6*b0 ';
% % probe(:,:,:,3) = fIC;
% % 
% % out_rgb(:,:,:,5, :) = imgb0RGB ; d4annotation{5} = 'b=0' ;
% % out_rgb(:,:,:,6, :) = imgb2000RGB ;  d4annotation{6} = 'b=2000' ;
% % out_rgb(:,:,:,4, :) = reshape(imgRGB, [size(imgRGB,[1 2 3]) 1 3]) ; d4annotation{4} = 'overlay' ;
% % out_rgb(:,:,:,8, :) = imgb2000normRGB; d4annotation{8} = 'b 2000 norm' ;
% % out_rgb(:,:,:,9, :) = imgvADCRGB; d4annotation{9}='vADC x10-3' ;
% % probe(:,:,:,4) = fIC ;
% % probe(:,:,:,5) = vb0tot/nSeries ;
% % probe(:,:,:,6) = v2000 ;
% % probe(:,:,:,8) = v2000norm ;
% % probe(:,:,:,9) = vADC*1000 ;
% % 
% % if ~opts.quiet
% %     sviewer(out_rgb, m, Name='fICmask fic fIC overlay b0 b2000', isrgb=true, ...
% %         d4annotation=d4annotation, indexD4=1, probe=probe)
% % 
% %     % sviewer of T2W if present
% %     for iT2 = 1:length(T2axSeriesNumbers)
% %         [vT2, mT2] = d2mat(dinfoT2ax,{'slice','series'}, ...
% %             'series',T2axSeriesNumbers(iT2), 'op','fp') ;
% %         sviewer(vT2, mT2, Name=['[',num2str(T2axSeriesNumbers(iT2)),'] T2W ax'], ...
% %             CLim=[0 prctile(vT2(:),98)], probe=vT2)
% %     end
% %     drawnow
% % end
% % 
% % % Produce montages of the above for the Report Generation
% % 
% % append(rpt,Chapter('Outputs'))
% % 
% % % Difficulties here setting colorbar position with timings and tiledlayout
% % hfcmap = figure(Name='The jet and divergent colormap',Visible=figVisible) ;
% % hfcmap.Units = 'normalized';
% % hfcmap.Position = [0.05 0.1 0.3 0.3];
% % 
% % t=tiledlayout(hfcmap,2,1);
% % ax=nexttile(t) ;
% % c = colorbar(ax);
% % c.Layout.Tile = 1 ;
% % c.Location = 'south';
% % c.FontSize = 12 ;
% % 
% % %c=colorbar(ax,'Location','south','FontSize',12,'Position',[0.05 0.6 0.9 0.3]) ;
% % colormap(ax,jet)
% % ax.Visible="off";
% % drawnow
% % 
% % ax2=nexttile(t) ;
% % c2 = colorbar(ax2) ;
% % c2.Layout.Tile=2;
% % c2.Location = 'south' ;
% % c2.FontSize = 12 ;
% % colormap(ax2,divMap)
% % ax2.Visible="off";
% % drawnow
% % 
% % cobj = mlreportgen.report.Figure(hfcmap);
% % cobj.Snapshot.Caption = ['The Jet and ',cmap_str,' divergent colormap'] ;
% % close(hfcmap)
% % append(rpt,cobj);
% % 
% % append(rpt,mlreportgen.dom.PageBreak())
% % 
% % rpt = montageReport(rpt,out_rgb, maskc, d4annotation, figVisible) ;
% % 
% % 
% % if opts.outputMatFile == true
% %     resultsFileMat = fullfile(resultsFolder,'outputs.mat') ;
% %     append(rpt,Paragraph(['Saving variables to MAT file: ',resultsFileMat]))
% % 
% %     save(resultsFileMat, 'fIC', 'fEES', 'fVASC', 'R', 'rmse', 'A', 'tparams', ...
% %         'vfopt', 'scheme', 'Y', 'm' , 'sortedBV', 'vprereg', 'vpostreg', ...
% %         'verdictVersion')
% % 
% %     if exist('vb0tot','var')
% %         save(resultsFileMat, 'vb0tot', '-append')
% %     end
% %     if exist('v2000','var')
% %         save(resultsFileMat, 'v2000', '-append')
% %     end
% % end
% % 
% % close(rpt)
% % 
% % if isdeployed
% %     disp('Deployed code: Folder path (volume) is within the deployed application')
% %     disp(' For Docker containers, this volume will be mapped to a local folder')
% % end
% % 
% % if ~opts.quiet && ~isdeployed
% %     rptview(rpt)
% % else
% %     disp(['Report ready in: ', rptFFN])
% % end
% % 
% % iptsetpref('ImshowBorder',iptprefsOrig.ImshowBorder)
% % iptsetpref('ImshowInitialMagnification',iptprefsOrig.ImshowInitialMagnification) 
% % 
% % if opts.outputDicom == true
% %     resultsFolderDICOM = fullfile(resultsFolder,'DICOM') ;
% % 
% %     disp('Starting DICOM outputs')
% % 
% %     dict2write = dicomdict("get") ;
% %     VR = 'implicit' ;
% %     writePrivate = false ;
% %     limitBitsStored = false ;
% %     useInputPrivatePhilips = false ;
% % 
% % 
% %     if opts.addPhilipsPrivate
% %         % Test original is Enhanced and Philips
% %         if strcmp(dfi.SOPClassUID,'1.2.840.10008.5.1.4.1.1.4.1') && ...
% %                 isfield(dfi,'Private_2001_10xx_Creator') && ...
% %                 strcmp(dfi.Private_2001_10xx_Creator,'Philips Imaging DD 001')
% %                 % input looks suitable
% %             disp('Adding some Philips private fields to output DICOMs.')
% % 
% %             [dinfoout, dict2write] = addPhilipsPrivate(dinfoout, locout, opts.vBaseSeriesNumber, opts.fICReconNumber ) ;
% %         
% %             VR = 'explicit' ;
% %             writePrivate = true ;
% %             limitBitsStored = true ; % For Philips (but not RGB)
% %             useInputPrivatePhilips = true ;
% % 
% %         else
% %             disp('NOT adding some Philips private fields to output DICOMs.')
% %             opts.addPhilipsPrivate = false ;
% %         end
% %     end
% %         
% %     wDargs = { ...
% %         'FrameOfReferenceUID','keep', ...
% %         'StudyInstanceUID','keep', ...
% %         'SeriesInstanceUID','new', ...
% %         'geom',geomout, ...
% %         'VR', VR, ...
% %         'writePrivate', writePrivate, ...
% %         'dict', dict2write, ...
% %         'useInputPrivatePhilips', useInputPrivatePhilips, ...
% %         'limitBitsStored', limitBitsStored, ...
% %         'rgbAsPalette', false} ;
% % 
% %     fICSeriesNumber = 100*opts.vBaseSeriesNumber + opts.fICReconNumber ;
% % 
% %     writeDicom(fIC, 'fICWithRescale','folder_name',resultsFolderDICOM, ...
% %         'header', {dinfoout, locout }, ...
% %         'quiet', true, ...
% %         'fnstem','fIC_withrescale', ...
% %         'SeriesDescription','fIC with RescaleSlope', ...
% %         'ProtocolName', 'fIC w rescale. Not for clinical or commercial use', ...
% %         'SeriesNumber', fICSeriesNumber, ...
% %         wDargs{:} ) ;
% % 
% %     newReconNumber  = opts.fICReconNumber +1 ; 
% %     newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %     if opts.addPhilipsPrivate
% %         [dinfoout(locout).MRSeriesReconstructionNumber] = deal(newReconNumber) ;
% %     end
% % 
% %     writeDicom(squeeze(256*out_rgb(:,:,:,4,:)), 'rgb', ... 
% %         'folder_name',resultsFolderDICOM, ...
% %         'header', {dinfoout, locout }, ...
% %         'fnstem','fIC_overlay', ...
% %         'quiet', true, ...
% %         'SeriesDescription','fIC overlay', ...
% %         'ProtocolName', 'fIC overlay. Not for clinical or commercial use', ...
% %         'SeriesNumber', newSeriesNumber, ...
% %         wDargs{:} ) ;
% % 
% %     newReconNumber  = newReconNumber +1 ;
% %     newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %     if opts.addPhilipsPrivate
% %         [dinfoout(locout).MRSeriesReconstructionNumber] = deal(newReconNumber) ;
% %     end
% % 
% %     writeDicom(squeeze(256*out_rgb(:,:,:,1,:)), 'rgb', ...
% %         'fnstem','fIC_jet_mask', ...
% %         'quiet', true, ...
% %         'SeriesDescription','fIC jet mask', ...
% %         'ProtocolName', 'fIC jet masked. Not for clinical or commercial use', ...
% %         'SeriesNumber', newSeriesNumber, ...
% %         'folder_name',resultsFolderDICOM, ...
% %         'header', {dinfoout, locout }, ...
% %         wDargs{:}) ;
% % 
% % 
% %     newReconNumber  = newReconNumber +1 ;
% %     newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %     if opts.addPhilipsPrivate
% %         [dinfoout(locout).MRSeriesReconstructionNumber] = deal(newReconNumber) ;
% %     end
% % 
% %     writeDicom(squeeze(256*out_rgb(:,:,:,2,:)), 'rgb', ...
% %         'fnstem','fIC_jet', ...
% %         'quiet', true, ...
% %         'SeriesDescription','fIC jet', ...
% %         'ProtocolName', 'fIC jet. Not for clinical or commercial use', ...
% %         'SeriesNumber', newSeriesNumber, ...
% %         'folder_name',resultsFolderDICOM, ...
% %         'header', {dinfoout, locout }, ...
% %         wDargs{:}) ;
% % 
% %     newReconNumber  = newReconNumber +1 ;
% %     newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %     if opts.addPhilipsPrivate
% %         [dinfoout(locout).MRSeriesReconstructionNumber] = deal(newReconNumber) ;
% %     end
% % 
% %     writeDicom(squeeze(256*out_rgb(:,:,:,7,:)), 'rgb', ...
% %         'fnstem','fIC_diverge', ...
% %         'quiet', true, ...
% %         'SeriesDescription','fIC diverge', ...
% %         'ProtocolName', 'fIC divergent. Not for clinical or commercial use', ...
% %         'SeriesNumber', newSeriesNumber, ...
% %         'folder_name',resultsFolderDICOM, ...
% %         'header', {dinfoout, locout }, ...
% %         wDargs{:}) ;
% % 
% % 
% % 
% %     if exist('imgb0','var')
% %         newReconNumber  = newReconNumber +1 ;
% %         newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %         if opts.addPhilipsPrivate
% %             [dinfoout(locout).MRSeriesReconstructionNumber] = deal(newReconNumber) ;
% %         end
% % 
% %         writeDicom(imgb0, 'grayWithRescale',...
% %             'quiet', true, ...
% %             'folder_name',resultsFolderDICOM, ...
% %             'header', {dinfoout, locout }, ...
% %             'fnstem','VERDICT_b0_aligned_summed', ...
% %             'ImageType',  'DERIVED\PRIMARY\DIFFUSION\ALIGNED' , ...
% %             'SeriesDescription','VERDICT b0 aligned summed', ...
% %             'ProtocolName', 'b0 aligned. Not for clinical or commercial use', ...
% %             'SeriesNumber', newSeriesNumber, ...
% %             wDargs{:}) ;
% %     end
% % 
% %     if exist('imgb2000','var')
% %         newReconNumber  = newReconNumber +1 ;
% %         newSeriesNumber = 100*opts.vBaseSeriesNumber + newReconNumber ;
% %         if opts.addPhilipsPrivate
% %             [dinfoout2000, dict2write] = addPhilipsPrivate(dinfoout, locb2000, opts.vBaseSeriesNumber, newReconNumber ) ;
% %         else
% %             dinfoout2000 = dinfoout ;
% %         end
% % 
% %         writeDicom(imgb2000, 'grayWithRescale','folder_name',resultsFolderDICOM, ...
% %             'header', {dinfoout2000, locb2000 }, ...
% %             'quiet', true, ...
% %             'fnstem','VERDICT_b2000_aligned', ...
% %             'ImageType',  'DERIVED\PRIMARY\DIFFUSION\ALIGNED' , ...
% %             'SeriesDescription','VERDICT b2000 aligned', ...
% %             'ProtocolName', 'b2000 aligned. Not for clinical or commercial use', ...
% %             'SeriesNumber', newSeriesNumber, ...
% %             wDargs{:}) ;
% %     end
% % 
% % end
% % 
% % disp('Finished.')
% % 
% % end
% % 
% % function rpt = montageReport(rpt,out_rgb, maskc, d4annotation, figVisible)
% % % montageReport Montage data for Report Generation
% % %
% % for ifig=1:size(out_rgb,4)
% %     hfmont = figure(Name='montage',Visible=figVisible) ;
% %     immont = squeeze(out_rgb(maskc{:},:,ifig,:)) ;
% %     immont = permute(immont,[1 2 4 3]) ;
% %     montage(immont)
% % 
% %     figrpt = mlreportgen.report.Figure(hfmont);
% %     figrpt.Snapshot.Caption = d4annotation{ifig} ;
% %     figrpt.Snapshot.ScaleToFit  = true ;
% %     append(rpt,figrpt);
% %     close(hfmont)
% % end
% % 
% % 
end

function pcheck = checkV(varargin)
% checkV Checks VERDICT data/parameters are as expected
%
% Intended to warn of situations where scans are mixed up, or re-planning has
% occurred. Checks slice locations, but does not check the PatientID.
% 
% pcheck = checkV(dinfo)  - checks certain tag values if present
% pcheck = checkV(v)      - checks for extreme data (no NaNs or entire slices zero)
% pcheck = checkV(geom1, geom2) - checks equality of structures (within tolerance)
% pcheck = checkV(geom, 'axial') - checks slices are close to non-angulated
%                                  axial
%
% Data should already have inividually been through data_process
%
% pcheck is a mlreportgen.dom.Paragraph or []

import mlreportgen.dom.*

if length(varargin) == 1 && isstruct(varargin{1})
    dinfo = varargin{1} ;

    % Check StudyDate, StudyTime, TR, TE, FlipAngle
    % Note StudyDate and StudyTime are often changed in anonymisation
    warned = checkDTag(dinfo, 'StudyDate', 'unique') ;
    warned = or(warned, checkDTag(dinfo, 'StudyTime', 'unique') );
    warned = or(warned, checkDTag(dinfo, 'RepetitionTime', 'geq', 2000) );
    warned = or(warned, checkDTag(dinfo, 'EffectiveEchoTime', 'geq', 30) );
    warned = or(warned, checkDTag(dinfo, 'FlipAngle', 'geq', 85) );

    if warned
        pcheck = Paragraph('Problem with one or more of StudyDate/Time, TR, TE, FA');
    else
        pcheck = [];
    end

elseif length(varargin) == 1 && isnumeric(varargin{1})
    % checkV(v)   
    v = varargin{1} ;
    warned = false ;
    for islice = 1:size(v,3)
        dslice = v(:,:,islice,:,:) ;
        sumSlice = sum(dslice(:)) ;
        if sumSlice <= 0 
            warned = true ;
            warning('Slice all 0')
        end
        if anynan(sumSlice)
            warned = true ;
            warning('NaN in data')
        end
    end
    if warned
        pcheck = Paragraph('! Data contains at least one slice all zero or NaNs present (maybe from registration)') ;
    else
        pcheck = [];
    end


elseif length(varargin) == 2 && isstruct(varargin{2})
    % checkV(geom1, geom2)
    geom1 = varargin{1} ;
    geom2 = varargin{2} ;
    n1 = length(geom1) ;

    warned = false ;
    if length(geom2) ~= n1
        warned = true ;
        warning('geoms not same size')
    end
    for ig = 1:n1
        if ~isequal([geom1(ig).Height geom1(ig).Width], ...
                    [geom2(ig).Height geom2(ig).Width])
            warned = true ;
            warning("Heights and widths must match")
        end
        if norm(geom1(ig).IPP - geom2(ig).IPP) > 1e-3*norm(geom1(ig).IPP)
            warned = true ;
            warning('IPPs not coincident')
        end
        if norm(geom1(ig).IOP(1:3) - geom2(ig).IOP(1:3) ) > 1e-3
            warned = true ;
            warning('IOP(1:3) not equal')
        end
        if norm(geom1(ig).IOP(4:6) - geom2(ig).IOP(4:6)) > 1e-3
            warned = true ;
            warning('IOP(4:6) not equal')
        end
        if ~isequal(geom1(ig).FrameOfReferenceUID, geom2(ig).FrameOfReferenceUID)
            warned = true ;
            warning('FrameOfReferenceUIDs not equal')
        end
    end

    if warned
        pcheck  = Paragraph('!! Geometries not consistent') ;
    else
        pcheck = [] ;
    end


elseif length(varargin) == 2 && ischar(varargin{2})
    % checkV(geom,'axial')
    if isequal(varargin{2},'axial')
        geom = varargin{1} ;
        n = length(geom) ;
        warned = false ;

        for ig = 1:n
            % 4e-2 threshold is about 2 degrees. Could use dgeom to get
            % angulations
            if norm(geom(ig).IOP(1:3) - [1; 0; 0]) > 4e-2
                warned = true ;
                warning(['IOP(1:3) not Left: ',num2str(geom(ig).IOP(1:3)')])
            end
            if norm(geom(ig).IOP(4:6) - [0; 1; 0]) > 4e-2 
                    warned = true ; 
                    warning(['IOP(4:6) not Posterior: ',num2str(geom(ig).IOP(4:6)')]) 
            end
        end

    else
        error('unknown calling pattern')
    end

    if warned
        pcheck  = Paragraph('!! Slice not true axial - check delta and DELTA.') ;
    else
        pcheck = [] ;
    end

else
    error('Unexpected parameters in call to checkV.')
end

end
    
function warned = checkDTag(dinfo, field, test, tval)
% checkDTag Check DICOM Tags are as expected
%
% Checks for either unique value, or, >=
%
% warned = checkDTag(dinfo, field, 'unique') (for character fields only)
% warned = checkDTag(dinfo, field, 'geq', testvalue)
%
%
% Examples
%  warned = checkDTag(dinfo,'RepetitionTime','geq',2000)
%  warned = checkDTag(dinfo,'StudyDate','unique')

warned = false ;

if isfield(dinfo,field)
     switch test
        case 'unique'
            vals = {dinfo.(field)} ;
            uvals = unique(vals) ;
            if length(uvals) ~= 1
                warned = true ;
                warning(['Found ',num2str(length(uvals)), ...
                    ' but expected exactly 1 value for field: ',field])
                uvals
            end
        case 'geq'
            vals = [dinfo.(field)] ;
            uvals = unique(vals) ;
            loc = find(uvals<tval) ;
            if ~isempty(loc)
                warned = true ;
                warning(['Found ',num2str(length(loc)),' values less than ', ...
                    num2str(tval),' for field: ',field])
            end
        otherwise
            error(['Unknown test: ',test])
    end
end
end


function opts = convertCharsToLogical(opts) 
% convertCharsToLogical Converts a field with char value to 
% logical if 'true' or 'false'
% 

fields = fieldnames(opts) ;

for ifield = 1:length(fields)
    fname = fields{ifield} ;
    if ischar(opts.(fname))
        switch opts.(fname)
            case {'true', 'True','TRUE'}
                opts.(fname) = true ;
            case {'false','False','FALSE'}
                opts.(fname) = false ;
        end
    end
end

end

