function try_VERDICT(dfolder, output_folder, opts)
% TRY_VERDICT VERDICT processing 
%
% Processing of VERDICT MRI data, fitting a ball, sphere and astrosticks
% model.
%
% Loads data from folder (Enhanced or Classic DICOM).
% Optional registration (translations only, based on central region of
% central slice, b=0 to b=0 with transformation applied to DW data).
% Fit of model to data. Model uses scheme (b-value, G, deltas)
% Default solver is lsqnonneg with Tikonhov regularisation. 
% Output of maps.
%
% Uses DICOM ManufacturerModelName to set the scheme unless forceXNATscheme 
% is true. ("XNAT" refers to reference processing used for Singh et al
% Radiology 2022)
%
% Outputs a PDF report and a log. Requires the Advanced Logger for MATLAB 
% (available on FileExchange / GitHub)
%
% try_VERDICT(dfolder, output_folder, Name, Value, ...)
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
%  try_VERDICT(dfolder, tempdir, solver='SPAMS')    % MATLAB mode only
%  try_VERDICT(dfolder, tempdir, 'solver', 'SPAMS') 
%
% As a command line:
% try_VERDICT 'inputfolder' 'outputfolder' solver SPAMS  addb0ToScheme true
%
% Configurations
% --------------
% XNAT closest
% swapinvXNAT=true, register=false, usedirecdiff=true, solver='SPAMS', forceXNATscheme=true, addb0ToScheme=true 
% XNATmimic = {'register','false','swapinvXNAT','true','usedirecdiff','true','solver','SPAMS','forceXNATscheme','true','addb0ToScheme','true'}
%
% XNAT but correct scanner (scheme) with b=0 added.
% swapinvXNAT=true; register=false; usedirecdiff=true; solver='SPAMS'; forceXNATscheme=false; addb0ToScheme=true ;
%
% XNAT registered input, correct scheme (with added b=0), lsqnoninTokonohv reg
% swapinvXNAT=true; register=false; usedirecdiff=true; solver='lsqnonnegTikonhov'; forceXNATscheme=false; addb0ToScheme=true ;
%
% MATLAB process (the default)
% swapinvXNAT=false, register=true, usedirecdiff=false, solver='lsqnonnegTikonhov', forceXNATscheme=false, addb0ToScheme=true 
%
%
% David Atkinson
%
% See also bv2scheme verdict_fit  sviewer  getSeriesVERDICT

% input arguments will be characters when called from command line in
% deployed mode.
arguments
    dfolder
    output_folder
    opts.outputDicom   = 'false' % output results in DICOM files
    opts.outputMatFile = 'true'  % output variables in one .mat file
    opts.register      = 'true'  % input data will have be registered
    opts.swapinvXNAT   = 'false' % swap in the data from a vXNAT.mat file
    opts.usedirecdiff  = 'false' % use the 3 directional diffusion images
    opts.solver        = 'lsqnonnegTikonhov' % solver used
    opts.forceXNATscheme = 'false' % ignores scanner and uses XNAT scheme file
    opts.addb0ToScheme   = 'true'  % adds b=0 (signal of 1) to every series
    opts.quiet           = 'false' % suppress figures and report viewer
    opts.maskhwmm        = '48'    % mask half-width in mm used for registration
end

opts = convertCharsToLogical(opts) ; 
if ischar(opts.maskhwmm)
    opts.maskhwmm = str2num(opts.maskhwmm) ;
end

% Set timestamp tstr for output folder naming and logging
tnow = datetime ;
tstr_nice = char(tnow) ;
tnow.Format='yyyy-MM-dd''T''HHmmss' ; % will form part of path name (no weird characters)
tstr = char(tnow) ;

% Create subfolder for results (name using time string)
resultsFolder = ['res-',tstr];
resultsFolder = fullfile(output_folder, resultsFolder) ;

% Set up Logging 
logFileName = 'verdict.log' ;
logFFN = fullfile(resultsFolder, logFileName) ;
logObj = mlog.Logger("verdict_log") ;
logObj.LogFolder = resultsFolder ;
logObj.LogFile = logFFN ; 

logObj.message(['Timestamp: ',tstr_nice,' (',tstr,')'])


[status,msg,msgID] = mkdir(output_folder) ;
if status == 0
    logObj.warning([msgID,' ', msg])
    return
end

[status,msg,msgID] = mkdir(resultsFolder) ;
if status == 0
    logObj.warning([msgID,' ', msg])
    return
end

[status,msg,msgID] = mkdir(fullfile(resultsFolder,'DICOM')) ;
if status == 0
    logObj.warning([msgID,' ', msg])
    return
end

% Set up Reporting
rptFileName = 'rptverdict.pdf' ;
rptFFN = fullfile(resultsFolder, rptFileName) ;

import mlreportgen.report.*
import mlreportgen.dom.*

rpt    = Report(rptFFN,'pdf');
rpt.Layout.Landscape = true;
append(rpt,TableOfContents)

append(rpt, Paragraph(['This report was generated: ', tstr_nice]) )
append(rpt, Paragraph(['This report is: ', rptFFN]))
append(rpt, Paragraph(['Log file name is: ',logFileName]))
append(rpt, Paragraph(['resultsFolder: ',resultsFolder]))

append(rpt,Chapter('Inputs'))

if opts.quiet
    figVisible = 'off' ; % figures do not display, but are still in report
else
    figVisible = 'on' ;
end


if opts.swapinvXNAT
    logObj.info('swapinvXNAT: true')
    if opts.register == true || opts.usedirecdiff == false
        logObj.warning(['Cannot register or use iso diff with swapinvXNAT']);
    end
else
    logObj.info('swapinvXNAT: false')
end

logObj.info(['register: ',num2str(opts.register)]) ;
logObj.info(['usedirecdiff: ',num2str(opts.usedirecdiff)])
logObj.info(['solver: ',opts.solver])


Rs = linspace(0.01,15.1,18) ; % radii used in fit. 
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
        logObj.warning('Unknown value of usedirecdiff')
end

if opts.addb0ToScheme
    naddb0 = 1 ; % number of b=0 added per series 
else
    naddb0 = 0 ;
end

if opts.forceXNATscheme == true
    if opts.addb0ToScheme == false
        logObj.warning('forceXNATscheme hence adding b=0 scans to scheme')
    end
    naddb0 = 1 ; % number of b=0 added per series 
end

logObj.info(['Number of b=0 added to scheme: ',num2str(naddb0)])
logObj.info(['Number of diffusion images per Series (nbd): ',num2str(nbd)])

% Report parameters in a table

pTObj = MATLABTable(struct2table(opts)) ;
pTObj.Style = [pTObj.Style { FontSize('8pt')} ] ;
append(rpt, pTObj)

logObj.info(['Calling getSeriesVERDICT with dfolder: ',dfolder])
[dinfo, vSeriesNumbers, vBV, dinfoT2ax, T2axSeriesNumbers, vXNAT] = getSeriesVERDICT(dfolder) ;

% Report exam details
dfi = dicominfo(dinfo(1).Filename) ;
if isfield(dfi,'PatientID')
    append(rpt,Paragraph(['PatientID: ',dfi.PatientID]))
end
if isfield(dfi,'StudyDescription')
    append(rpt,Paragraph(['StudyDescription: ',dfi.StudyDescription]))
end
if isfield(dfi,'InstitutionName')
    append(rpt,Paragraph(['InstitutionName: ',dfi.InstitutionName]))
end
if isfield(dfi,'ProtocolName')
    append(rpt,Paragraph(['ProtocolName: ',dfi.ProtocolName]))
end
if isfield(dfi,'Manufacturer')
    append(rpt,Paragraph(['Manufacturer: ',dfi.Manufacturer]))
end
if isfield(dfi,'ManufacturerModelName')
    append(rpt,Paragraph(['ManufacturerModelName: ',dfi.ManufacturerModelName]))
end
if isfield(dfi,'SoftwareVersions')
    append(rpt,Paragraph(['SoftwareVersions: ',dfi.SoftwareVersions]))
end


% Check dinfo for reasonable parameter values
% Add to report and/or logging
pcheck = checkV(dinfo) ;
append(rpt, pcheck)

[sortedBV, indsortedBV] = sort(vBV(:,2)) ;

nSeries = length(vSeriesNumbers) ;

if opts.register
    regTranslations = zeros([nSeries 2]) ;  % Registration transformations stored for report
    [optimizer,metric] = imregconfig('monomodal') ; % monomodal as b=0 to b=0
    fw = fw_pix ;    % figure width in pixels
    fh = fw/nSeries*3 ;
    hfreg_check = figure(Name='Pre and Post Registration', Position=[250 500 fw round(fh)],Visible=figVisible);
    treg = tiledlayout(3,nSeries,'TileSpacing','none','Padding','tight') ;
    axrs = [] ;
end

for iSeries = 1: nSeries
    sn_this = vSeriesNumbers(indsortedBV(iSeries)) ;
    bv_this = vBV(indsortedBV(iSeries),2) ;

    % Get b=0 data
    [vb0, mb0, b0loc] = d2mat(dinfo,{'slice','bv','series'},'bv',0, ...
        'series',sn_this,'op','fp') ;

    % Get b>0 data
    if opts.usedirecdiff
        [vbv, mbv, bvloc] = d2mat(dinfo,{'slice','bdirec','ddty','series'}, ...
            'ddty',1,'series',sn_this,'op','fp') ;
        if size(vbv,4)~=3, warning('Expected 3 diffusion directions'), end
    else
        [vbv, mbv, bvloc] = d2mat(dinfo,{'slice','ddty','series'},'series',sn_this, ...
            'ddty', 2, 'op','fp') ;
    end

    if opts.swapinvXNAT
        logObj.message('swapinvXNAT is true')
        % in the reference XNAT, b-value ordering is 3000 to 90
        switch bv_this
            case 90
                vb0 = vXNAT(:,:,:,17) ;
                vbv = vXNAT(:,:,:,18:20) ;
            case 500
                vb0 = vXNAT(:,:,:,13) ;
                vbv = vXNAT(:,:,:,14:16) ;
            case 1500
                vb0 = vXNAT(:,:,:,9) ;
                vbv = vXNAT(:,:,:,10:12) ;
            case 2000
                vb0 = vXNAT(:,:,:,5) ;
                vbv = vXNAT(:,:,:,6:8) ;
            case 3000
                vb0 = vXNAT(:,:,:,1) ;
                vbv = vXNAT(:,:,:,2:4) ;
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
        % first pass
        Y = zeros([size(vb0,[1 2 3]), (nbd + naddb0)*nSeries]) ;
        vb0tot = zeros(size(vb0));
        v2000 =  zeros(size(vb0)) ;

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
        logObj.info(['Registration slice number: ',num2str(reg_slice)]) ;


        vfixed = vb0(maskc{:},reg_slice) ; % fixed image for registration
        vfixedUpper = prctile(vfixed(:),98) ;
        
        dfullinf = dicominfo(dinfo(1).Filename) ;
        if isfield(dfullinf,'ManufacturerModelName')
            scanner = dfullinf.ManufacturerModelName ;
            logObj.info(['scanner:', scanner]) ;
        else
            logObj.warning('Scanner (ManufacturerModelName) not known')
        end

        if opts.forceXNATscheme
            logObj.message('Forcing XNAT scanner.');
            append(rpt,Paragraph('! Forcing use of XNAT scheme file.'))
            scanner = 'XNAT' ;
        end

        geomout = mb0.geom ; % used as reference and for writing DICOMs
        dinfoout = dinfo ; % for DICOM output
        locout = b0loc ;
    end % iSeries 1

    % Check data is consistent and report
    pcheck = checkV(geomout, mb0.geom) ;
    append(rpt, pcheck)
    pcheck = checkV(geomout, mbv.geom) ;
    append(rpt, pcheck)
    pcheck = checkV(vb0) ;
    append(rpt,pcheck)
    pcheck = checkV(vbv) ;
    append(rpt,pcheck)
    pcheck = checkV(mb0.geom,'axial') ;
    append(rpt, pcheck)
   

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
        axr = nexttile(treg,tilenum(treg,1,iSeries)) ;
        imshow(vmoving,[0 vfixedUpper],'Parent',axr) ;
        axrs = [axrs axr];

        axr = nexttile(treg,tilenum(treg,2,iSeries)) ;
        imshow(vmoving,[ ],'Parent',axr), hold on
        plot(axr,[1 nxlim],[ylinePos ylinePos],'Color','y')
        plot(axr,[xlinePos xlinePos],[1 nylim],'Color','y')
        axrs = [axrs axr];
        axr = nexttile(treg,tilenum(treg,3,iSeries)) ;
        imshow(vb0_reg(maskc{:},reg_slice),[ ],'Parent',axr), hold on
        plot(axr,[1 nxlim],[ylinePos ylinePos],'Color','y')
        plot(axr,[xlinePos xlinePos],[1 nylim],'Color','y')
        axrs = [axrs axr];
    else
        vb0_reg = vb0 ;
        vbv_reg = vbv ;
    end

    if abs(bv_this-2000) < 10
        v2000 = sum(vbv_reg,4) ;
        locb2000  = bvloc ;
    end

    vb0tot = vb0tot + vb0_reg ;

    vnorm = vbv_reg./repmat(vb0_reg,[1 1 1 nbd]) ; % DW images, normalised by b=0

    if naddb0 == 1
        vnorm = cat(4,ones(size(vb0_reg)), vnorm ) ;
    end

    Y(:,:,:,1+(iSeries-1)*(nbd + naddb0) : iSeries*(nbd + naddb0)) = vnorm ; % Stack DW images in 4th dim

    bvs = bv2scheme(bv_this, scanner) ; 
    scheme(1+(iSeries-1)*(nbd+naddb0)+1 : iSeries*(nbd+naddb0)) = bvs ;
    if naddb0 == 1
        scheme(1+(iSeries-1)*(nbd+naddb0)) = bv2scheme(0, scanner) ;
    end

    % for output/debugging
     vprereg(:,:,:,1+(iSeries-1)*(nbd + 1) : iSeries*(nbd + 1)) = cat(4,vb0,vbv) ;
    vpostreg(:,:,:,1+(iSeries-1)*(nbd + 1) : iSeries*(nbd + 1)) = cat(4,vb0_reg,vbv_reg) ;


end
if opts.register
    linkaxes(axrs)
end

hfmedian = figure(Name='medians', Visible=figVisible) ;
plot(plotbv,medianvb0,'DisplayName','Median B0'), hold on
plot(plotbv,medianvbv,'DisplayName','Median vbv')
plot(plotbv,medianvbv./medianvb0*medianvb0(1),'DisplayName','Normalised median vbv')
grid on
xlabel('b-value'), ylabel('median')
legend

figrpt = Figure(hfmedian);
figrpt.Snapshot.Caption = 'Medians excluding end slices' ;
figrpt.Snapshot.ScaleToFit = true ;
append(rpt,figrpt);

% Display inputs Y for QA. Also place in report.
Dref = 2e-9 ;  % A reference diffusivity for scaling the signals
fw = fw_pix ;    % figure width in pixels
fh = fw/size(Y,3)*size(Y,4) ;
hfY = figure(Name='Y components',Position=[200 50 fw round(fh)],Visible=figVisible) ;
t=tiledlayout(size(Y,4), size(Y,3),'TileSpacing','none','Padding','tight') ;
axs = [] ;
for iY3 = 1: size(Y,3)
    for iY4 = 1: size(Y,4)
        ax = nexttile(tilenum(t,iY4,iY3)) ;
        axs = [axs ax];
        upper = 1.2 * exp(-Dref*scheme(iY4).bval*1000) ;
        imshow(Y(:,:,iY3,iY4),[0 upper])
    end
end
linkaxes(axs)


% Report Generation
append(rpt, Paragraph(' '))
pobj = Paragraph(['Series Numbers in input.   Scheme for scanner: ',scanner]) ;
append(rpt, pobj)

tableStyle = { ...
    Width('100%'), ...
    Border('solid','black','1px'), ...
    ColSep('solid','black','1px'), ...
    RowSep('solid','black','1px') ...
    };

schemeTableObj = mlreportgen.dom.MATLABTable(struct2table(scheme)) ;
schemeTableObj.Style = tableStyle ;

sT = table(vSeriesNumbers', vBV) ;
sT.Properties.VariableNames = {'Series', 'b-values'} ;
seriesTableObj = mlreportgen.dom.MATLABTable(sT) ;
seriesTableObj.Style = tableStyle ;

% lo_table is invisible - used for layout
lo_table = Table({seriesTableObj,' ', schemeTableObj}) ;
lo_table.entry(1,1).Style = {Width('3.2in')};
lo_table.entry(1,2).Style = {Width('.4in')};
lo_table.entry(1,3).Style = {Width('3.2in')};
lo_table.Style = {Width('100%'), ResizeToFitContents(false)};

append(rpt,lo_table) 


figrpt = mlreportgen.report.Figure(hfY);
figrpt.Snapshot.Caption = 'Entries in Y, the diffusion images normalized by their corresponding b=0' ;
figrpt.Snapshot.ScaleToFit = true ;
append(rpt,figrpt);


if opts.register
    hffixed = figure(Name='Reg Fixed',Visible=figVisible) ;
    imshow(vfixed,[0 vfixedUpper])
    figrpt = mlreportgen.report.Figure(hffixed);
    figrpt.Snapshot.Caption = 'The slice and masked region treated as fixed for registration' ;
    append(rpt,figrpt);
    close(hffixed)

    append(rpt, Paragraph('Translations found in registration (pixels)')) ;
    
    if any(abs(regTranslations) > 5, 'all')
        append(rpt,Paragraph('CHECK Registrations (large translation found)'))
    end

    regTableObj = mlreportgen.dom.MATLABTable(table(regTranslations)) ;
    append(rpt,regTableObj)

    figrpt = mlreportgen.report.Figure(hfreg_check);
    figrpt.Snapshot.Caption = 'Slice and mask region. Before registration, same windowing (top), scaled (middle) and after (bottom) registration' ;
    append(rpt,figrpt);
end

append(rpt,Chapter('Fitting'))

% PERFORM FIT
logObj.message('Calling verdict_fit') ;

mask = ones([size(vb0,[1 2 3])]) ;

[fIC, fEES, fVASC, R, rmse, A, t, vfopt] = verdict_fit(scheme,Y,Rs=Rs, ...
    mask=mask, solver=opts.solver) ;

append(rpt, Paragraph(' '))
append(rpt, Paragraph(' Entries of AMICO system matrix A'))
append(rpt, Paragraph(' '))
ATableObj = mlreportgen.dom.MATLABTable(table(A)) ;
ATableObj.Style = [ATableObj.Style {FontSize('8pt')}] ;
append(rpt,ATableObj)

append(rpt, Paragraph(' '))
append(rpt, Paragraph(' Tissue Properties Used in Fit:'))
append(rpt, Paragraph(' '))
tissueObj = MATLABTable(struct2table(t)) ;
tissueObj.Style = [tissueObj.Style { FontSize('8pt')} ] ;
append(rpt, tissueObj)

append(rpt, Paragraph(' '))
pobj = Paragraph(['Solver used in verdict_fit was: ',vfopt.solver]);
append(rpt,pobj)

fsum = fIC + fEES + fVASC ;

m.geom = geomout ;

if anynan(R)
    logObj.warning('R contained NaNs') ;
    R(isnan(R))=0 ;
end

nfIC = fIC./fsum ;

fQA = cat(4,fEES,fVASC, fsum, fIC, nfIC, rmse) ;
d4annotation = {'fEES','fVASC','fsum','fIC','nfIC','rmse'} ;
if ~opts.quiet
    sviewer(fQA,m,Name='fEES fVASC fsum, fIC, nfIC, rmse', CLim=[0.2 0.8], ...
    d4annotation=d4annotation, indexD4=4)

    sviewer(R,m)
end



% Builds a 5D RGB matrix out_rgb with fIC map, overlay, data.
%
% Colormap fIC using jet and with all fIC>1 clipped
% Set rmse> rmseThresh to magenta (not in jet color map)

nColor = 256 ;
cmap = jet(nColor) ;
out_rgb = zeros([size(fIC,[1 2 3]) 7 3]) ;
for islice = 1:size(fIC,3)
    % gray2ind clips fIC below 0 or above 1.
    ficrgb = ind2rgb(gray2ind(fIC(:,:,islice),nColor), cmap) ;
    % columnise to enable re-colouring of regions with poor fit rmse
    ficrgbcol = reshape(ficrgb,[size(ficrgb,1)*size(ficrgb,2) 3]) ;

    loc = find(rmse(:,:,islice) > rmseThreshold ) ;
    ficrgbcolnomask = ficrgbcol ;
    ficrgbcol(loc,:) = repmat(rmseRGB,[length(loc) 1])  ; % magenta
    ficrgb = reshape(ficrgbcol,[size(ficrgb,1) size(ficrgb,2) 3]) ;
    ficrgbnomask = reshape(ficrgbcolnomask,[size(ficrgb,1) size(ficrgb,2) 3]) ;
    out_rgb(:,:,islice,1,:) = reshape(ficrgb, [size(fIC,[1 2]) 1 1 3]) ;
    out_rgb(:,:,islice,2,:) = reshape(ficrgbnomask, [size(fIC,[1 2]) 1 1 3]) ;
end
d4annotation{1} = 'fIC mask' ;
d4annotation{2} = 'fIC' ;

imgb0 = mat2gray(vb0tot/nSeries, ...
    [0 double(prctile(vb0tot(:)/nSeries, 98))]) ;

imgb2000 = mat2gray(v2000, [0 double(max(v2000(:)))]) ;


% fIC single colour overlay
overlay = fIC ;
overlay(overlay>1)=1;
overlay(overlay<0)=0 ;

% alphaOverlay = overlay ; % This should probably be more targeted at 0.42
xmid = 0.42 ;
grad = 2 ;
intc = xmid*(1-grad) ;
alphaOverlay = grad*overlay + intc ;
alphaOverlay(alphaOverlay<0)=0 ;
alphaOverlay(alphaOverlay>1)=1;

R = 1 ; G = 0.1 ; B = G ; % Slightly pink overlay colour
RBase = 0.9 ; GBase = RBase ; BBase = GBase ; % Base image 'colour'

overlayRGB = cat(4, R*overlay,G*ones(size(overlay)),B*ones(size(overlay))) ;
        
imgBaseRGB = cat(4, RBase*imgb0, GBase*imgb0, BBase*imgb0) ;

imgRGB = alphaOverlay.*overlayRGB + (1-alphaOverlay).*imgBaseRGB ;


imgb0RGB    = cat(5, imgb0, imgb0, imgb0) ;
imgb2000RGB = cat(5, imgb2000, imgb2000, imgb2000 ) ;

% divergent
cmap_str = 'D09' ;
baseMap = colorcet(cmap_str) ;
% Skew colormap to put mid-white at fic=0.42 for divergent
whiteVal = 0.42 ;

cmv = linspace(0,1,size(baseMap,1)) ;
cmq = cmv/whiteVal*0.5 ;

extrapVal = baseMap(end,:) ;

redv   = interp1(cmv, baseMap(:,1), cmq,"linear",extrapVal(1)) ;
greenv = interp1(cmv, baseMap(:,2), cmq,"linear",extrapVal(2)) ;
bluev  = interp1(cmv, baseMap(:,3), cmq,"linear",extrapVal(3)) ;

divMap = cat(2, redv(:), greenv(:), bluev(:)) ;

for islice = 1:size(fIC,3)
    ficdivrgb = ind2rgb(gray2ind(fIC(:,:,islice),nColor), divMap) ;
    out_rgb(:,:,islice,7,:) =reshape(ficdivrgb,[size(fIC,[1 2]) 1 1 3]) ;
end
d4annotation{7}=cmap_str;

% jet fIC overlay, constant alpha
calpha = 0.6 ;
out_rgb(:,:,:,3,:) = calpha*imgb0RGB + (1-calpha)*out_rgb(:,:,:,1,:) ;
d4annotation{3}='fIC ';

out_rgb(:,:,:,5, :) = imgb0RGB ; d4annotation{5} = 'b=0' ;
out_rgb(:,:,:,6, :) = imgb2000RGB ;  d4annotation{6} = 'b=2000' ;
out_rgb(:,:,:,4, :) = reshape(imgRGB, [size(imgRGB,[1 2 3]) 1 3]) ; d4annotation{4} = 'overlay' ;


if ~opts.quiet
    sviewer(out_rgb, m, Name='fICmask fic fIC overlay b0 b2000', isrgb=true, ...
        d4annotation=d4annotation, indexD4=1)

    % sviewer of T2W if present
    for iT2 = 1:length(T2axSeriesNumbers)
        [vT2, mT2] = d2mat(dinfoT2ax,{'slice','series'}, ...
            'series',T2axSeriesNumbers(iT2), 'op','fp') ;
        sviewer(vT2, mT2, Name=['[',num2str(T2axSeriesNumbers(iT2)),'] T2W ax'], ...
            CLim=[0 prctile(vT2(:),98)])
    end
end

% Produce montages of the above for the Report Generation

append(rpt,Chapter('Outputs'))

hfcmap = figure(Name='The jet colormap',Visible=figVisible) ;
hfcmap.Units = 'normalized';
hfcmap.Position = [0.05 0.1 0.3 0.1];
colormap jet

c=colorbar;
hax=gca;
hax.Visible="off";
c.Location = 'south' ;
c.FontSize=12;
c.Position = [0.05 0.1 0.9 0.7] ;
cobj = mlreportgen.report.Figure(hfcmap);
cobj.Snapshot.Caption = 'The Jet colormap' ;
close(hfcmap)
append(rpt,cobj);

append(rpt,mlreportgen.dom.PageBreak())

rpt = montageReport(rpt,out_rgb, maskc, d4annotation, figVisible) ;

close(rpt)

if ~opts.quiet
    rptview(rpt)
else
    logObj.message(['Report ready in: ', rptFFN])
end

if opts.outputMatFile == true
    resultsFileMat = fullfile(resultsFolder,'outputs.mat') ;
    logObj.info(['Saving variables to file: ',resultsFileMat])

    save(resultsFileMat, 'fIC', 'fEES', 'fVASC', 'R', 'rmse', 'A', 't', ...
        'vfopt', 'scheme', 'Y', 'm' , 'sortedBV', 'vprereg', 'vpostreg')

    if exist('vb0tot','var')
        save(resultsFileMat, 'vb0tot', '-append')
    end
    if exist('v2000','var')
        save(resultsFileMat, 'v2000', '-append')
    end
end

if opts.outputDicom == true
    resultsFolderDICOM = fullfile(resultsFolder,'DICOM') ;

    logObj.info(['Starting DICOM writing in folder: ',resultsFolderDICOM])

    writeDicom(fIC, 'fICx1000','folder_name',resultsFolderDICOM, ...
        'header', {dinfoout, locout }, ...
        'SeriesDescription','fIC x1000', ...
        'FrameOfReferenceUID','keep', ...
        'StudyInstanceUID','keep', ...
        'SeriesInstanceUID','new', 'geom',geomout) ;

    writeDicom(fIC, 'fICWithRescale','folder_name',resultsFolderDICOM, ...
        'header', {dinfoout, locout }, ...
        'quiet', true, ...
        'fnstem','fIC_withrescale', ...
        'SeriesDescription','fIC with RescaleSlope', ...
        'FrameOfReferenceUID','keep', ...
        'StudyInstanceUID','keep', ...
        'SeriesInstanceUID','new', 'geom',geomout) ;

    writeDicom(squeeze(out_rgb(:,:,:,1,:)), 'rgb', ...
        'fnstem','fIC_colour', ...
        'quiet', true, ...
        'SeriesDescription','fIC colour', ...
        'folder_name',resultsFolderDICOM, ...
        'header', {dinfoout, locout }, ...
        'FrameOfReferenceUID','keep', ...
        'StudyInstanceUID','keep', ...
        'SeriesInstanceUID','new', 'geom',geomout) ;

    if exist('imgb0','var')
        writeDicom(imgb0, 'grayWithRescale',...
            'quiet', true, ...
            'folder_name',resultsFolderDICOM, ...
            'header', {dinfoout, locout }, ...
            'fnstem','VERDICT_b0_aligned_summed', ...
            'ImageType',  'DERIVED\PRIMARY\DIFFUSION\ALIGNED' , ...
            'SeriesDescription','VERDICT b0 aligned summed', ...
            'FrameOfReferenceUID','keep', ...
            'StudyInstanceUID','keep', ...
            'SeriesInstanceUID','new', 'geom',geomout) ;
    end

    if exist('imgb2000','var')
        writeDicom(imgb2000, 'grayWithRescale','folder_name',resultsFolderDICOM, ...
            'header', {dinfoout, locb2000 }, ...
            'quiet', true, ...
            'fnstem','VERDICT_b2000_aligned', ...
            'ImageType',  'DERIVED\PRIMARY\DIFFUSION\ALIGNED' , ...
            'SeriesDescription','VERDICT b2000 aligned', ...
            'FrameOfReferenceUID','keep', ...
            'StudyInstanceUID','keep', ...
            'SeriesInstanceUID','new', 'geom',geomout) ;
    end

end

end

function rpt = montageReport(rpt,out_rgb, maskc, d4annotation, figVisible)
% montageReport Montage data for Report Generation
%
for ifig=1:size(out_rgb,4)
    hfmont = figure(Name='montage',Visible=figVisible) ;
    immont = squeeze(out_rgb(maskc{:},:,ifig,:)) ;
    immont = permute(immont,[1 2 4 3]) ;
    montage(immont)

    figrpt = mlreportgen.report.Figure(hfmont);
    figrpt.Snapshot.Caption = d4annotation{ifig} ;
    figrpt.Snapshot.ScaleToFit  = true ;
    append(rpt,figrpt);
    close(hfmont)
end


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
        pcheck  = Paragraph('!! Slice not true axial - implications for delta and DELTA.') ;
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

