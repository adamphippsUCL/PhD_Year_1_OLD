function writeDicom(data, data_type, varargin)
% writeDicom  Writes DICOM files 
%   Data from ADC, registration output, T1map, pCASL, resliced PETMR etc
% Note the structure of this function is likely to change.
%
%  writeDicom(data, data_type, param, value, ...)
%
% !! NOT FOR CLINICAL or COMMERCIAL USE !! (No checks on geometry, frame of ref, or patient)
%
%
% data:  [ny nx nslice nt]
%  or [ny nx 3] or [ny nx nz 3] if data_type is 'rgb'
%  or [ny nx nslice necho] for multi-echo
%
% data_type:
%   'ADCbody' expects that after multiplying by 1e6, ADC will have
%        values in range 0-4000.
%   'rsPET'  Rescaled PET
%   'positive' Rescales if too large, clips negative to zero.
%   'T1map' clips data outside of range 0 - 10000. Sets window 150 - 2000.
%   'rgb' is RGB, e.g. colour coded FA map.
%   'fatfraction' Converts input numbers in range 0-1 to integers in range
%   0 - 100.
%    'LWF'   uses a range 0 - 100
%    'LWF_2' uses a range 0 - 1000
%    'fICx1000'   uses a range 0 - 1000, does not provide a RescaleSlope
%    'fICWithRescale' uses a 0 - 1000 range, adds RescaleSlope and
%    RescaleIntercept (despite not strictly in classic MR?)
%    'grayWithRescale' Saves values between 0  and 1 using a Rescale of
%    1000. Clips outside this region
%    'RawRecon' Reconstruction from raw data.
%    'multi-echo'  data must then be [ny nx necho] or [ny nx necho] and TEs
%                  input
%
%  Param / value pairs
%    'geom'       REQUIRED (from  d2mat)
%    'fnstem', filename stem (default depends on data_type)
%    'folder_name', output folder name, default is to call uigetdir, will
%                   overwrite without warning
%    'SeriesDescription'   (default depends on data_type)
%    'SeriesNumber'  default is to add 20 to that in header
%    'FrameOfReferenceUID'  'keep', {'new'} or actual UID. 
%    'StudyInstanceUID'     {'keep'}, 'new' or actual UID.
%    'SeriesInstanceUID'     {'keep'}, 'new' or actual UID.
%    'anon' {'false'}, 'true' blanks only name, DOB and patientID
%    'header' cell array {dinfo, locs} or {'blankMR'}.  Determines
%         SOPClassUID. 
%         Note dinfo only used to provide original DICOM filenames , i.e. a
%         structure with the field Filename.
%         Data derived from multiframe MR will be 
%         converted to singleframe with some loss of information.
%         Siemens PETMRI SOPClassUID not supported by dicomwrite in Create
%         Mode
%            dinfo:  see above
%            locs: [nslice] or [nslice nt] locations in dinfo of slices
%    'quiet' {false} | true supress warning box
%    'burn_text' {''} text to burn into each image. For low res images, 
%                keep short e.g. 'Non-diagnostic'
%    'ImageComments' {''} text shown in Dicom viewer
%    'MRAquisitionType' {''}, if empty not set, otherwise always set
%    'ProtocolName' {''}
%    'TEs' required in seconds for data type 'multi-echo'
%    'PatientName'
%    'WindowWidth', 'WindowCenter' 
%    'ImageType'
%    'useInputPrivatePhilips' {false}
%    'limitBitsStored' {false}
%    'rgbAsPalette' {true}
%
%
% Requires text2im from https://uk.mathworks.com/matlabcentral/fileexchange/19896-convert-text-to-an-image
%
% Example use
% ===========
%  dinfo = datparse ;
%  [vol, mat, locs] = d2mat(dinfo,{'slice','bv','series'},'series',13,'op','fp') ;
%  [ADC,S0] = calcADC(vol,mat.bvVec) ;
%  writeDicom(ADC, 'ADCbody', 'header', {dinfo, locs(:,1)},'geom',mat.geom, 'FrameOfReferenceUID','keep' )
%
%  dinfo = datparse ;
%  [vol, mat, locs] = d2mat(dinfo,{'dyn'},'resize', [128 128],'op','fp') ;
%  [Ireg, globalres] = RDDR(vol) ;
%  writeDicom(Ireg,'positive', 'geom', mat.geom, 'header', {dinfo, locs} )
%
% T1map
%  dinfo = datparse ;  % MFA data
%  [volip, matp, locs] = d2mat(dinfo,{'slice','fa','wfio'},'wfio',3,'op','fp') ;
%  db1 = datparse ;  % B1 map
%  [vb1,matb1] = d2mat(db1,{'slice'},'op','fp') ;
%  vb1r = dreslice(vb1,matb1,matp) ;
%  [t1lipr, m0lipr] = MFAT1(volip(:,:,:,:),matp.faVec,4,vb1r(:,:,:)) ;
%  writeDicom(t1lipr, 'T1map', 'header', {dinfo, locs(:,1)}, 'geom', matp.geom, 'FrameOfReferenceUID','keep')
%
% pCASL
%   
% [cbf, cbf_mask, params, dcalib, mcalib, locs_calib] = pcasl_proc ;  
% writeDicom(cbf.*cbf_mask, 'pCASL', 'header', {dcalib,locs_calib}, 'geom',mcalib.geom,'FrameOfReferenceUID','keep')
%
% rsPET  Resliced PET from PET MRI scanner
%  Note FrameOfReferenceUID seems to be different between MR and PET?
%  Also, dicomwrite will not do full verification of PET data so have to
%  use CreateMode Copy
%
%  dpet = datparse ;  % select dynamic PET
%  ddsc = datparse ;  % select DSC
%  [vdsc,mdsc,ldsc] = d2mat(ddsc,{'slice','aqno'},'aqno',1,'op','dv') ;
%  [vdpet, mdpet,ldpet] = d2mat(dpet,{'slice','TinSeries'},'op','dv') ;
%  [volrs, matp_rs] = dreslice(vdpet,mdpet, mdsc) ;
%   dscfor = unique({ddsc.FrameOfReferenceUID}) ;
%  writeDicom(volrs,'rsPET','geom',matp_rs.geom,'FrameOfReferenceUID',dscfor{1}, 'header',{dpet, ldpet}) ;
%
% DTI 
%  
%  [D, op] = dmf2d('dinfo',dinfo,'slices','all') ;
%  cfa = d2cfa(D) ;
%  writeDicom(cfa,'rgb','folder_name',folder_name, 'FrameOfReferenceUID' , 'keep', ...
%     'header' ,{dinfo, op.loc1(:)} ,'geom', [op.mat.geom])
%
%
% SWI2D
%  dffe = datparse(dselect) ;  or  dffe = datparse ;
%  [vffe,mffe, locffe] = d2mat(dffe,{'slice','itype'},'op','dv') ;
%  cdata = vffe(:,:,:,1) .* exp(1i.* vffe(:,:,:,2)/1000) ;
%  [imp5, op] = swi2D(cdata,'fw',[0.5 0.5]) ;
%  writeDicom(abs(imp5),'positive','SeriesDescription','swi2D', ...
%     'header',{dffe locffe},'geom',mffe.geom,...
%     'burn_text','Non-Diagnostic','ImageComments',op.comment)
%
% Fat Fraction
% writeDicom(FF, 'fatfraction','geom', lFO.geom, ...
%  'FrameOfReferenceUID', 'keep', 'header', {FOinfo, lFO}, ...
%  'ImageComments','fat fraction','burn_text','FF non-diagnostic', 'SeriesInstanceUID','new')
%
% Luminal Water Fraction (use type LWF_2 for range 0 - 1000)
%   dinfo = dmfparse(multiTEfn) ;
%   [vme, mme, lme] = d2mat(dinfo,{'slice','echo'},'op','fp') ;
%   writeDicom(LWF, 'LWF','geom', mme.geom, ...
%     'FrameOfReferenceUID', 'keep', 'header', {dinfo, lme}, ...
%     'ImageComments','Luminal Water Fraction','burn_text',' LWF non-diagnostic ', ...
%     'SeriesInstanceUID','new')
%
%
% Creating a Multi-echo test data
%    writeDicom(img,'multi-echo','TEs',TEs, 'geom',geom, 'SeriesNumber',20, 'SeriesInstanceUID','new', ...
%      'StudyInstanceUID','new', 'MRAcquisitionType','2D', ...
%      'PatientName',pname)
%
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%
% Copyright 2020-2021. David Atkinson D.Atkinson@ucl.ac.uk
% University College London
%
% See also D2MAT calcADC DATPARSE text2im
%

% Set defaults (may be overwritten later)
isrgb = false ;
forceAddRescale = false ;
limitBitsStored = false ;

switch data_type
    case 'ADCbody'
        SeriesDescription = 'Calculated ADC' ;
        fnstem = 'ADC' ;
        ImageType = 'DERIVED\SECONDARY\ADC' ;  %?? Not sure of standard here
        
        scale = 1e6 ; % ADC scaling factor
        data = data * scale ;
        data_upper = 4000 ;  % ADC values over this are written as 4000.
        data_lower = 0 ;
        WindowWidth = 1999 ;
        WindowCenter = 1000 ;
    case 'pCASL'
        SeriesDescription = 'Calculated CBF' ;
        fnstem = 'ASL' ;
        ImageType = 'DERIVED\SECONDARY\CBF' ;
        
        data_upper = 100 ;  % CBF values
        data_lower = 0 ;
        WindowWidth = 49 ;
        WindowCenter = 50 ;
        
    case 'rsPET'
        SeriesDescription = 'Resliced PET' ;
        fnstem = 'rsPET' ;
        ImageType = 'DERIVED\SECONDARY\RSPET';
        data_max = max(data(:)) ;
        if data_max > intmax('uint16') 
            disp(['Data max exceeds uint16 - scaling'])
            scale = double(intmax('uint16')) / double(data_max) ;
        else
            scale = 1 ;
        end
        data = data * scale ;
        
        data_lower = 0 ;
        data_upper = max(data(:)) ;
        WindowWidth = data_upper ;
        WindowCenter = floor(data_upper /2) ;
        
    case 'positive'
        SeriesDescription = 'MATLAB Calculated Data' ;
        fnstem = 'res' ;
        ImageType = 'DERIVED\SECONDARY\CALCULATED' ;  %?? Not sure of standard here
        
        data_max = max(data(:)) ;
        if data_max > intmax('uint16') 
            disp(['Data max exceeds uint16 - scaling'])
            scale = double(intmax('uint16')) / double(data_max) ;
        else
            scale = 1 ;
        end
        data = data * scale ;
        data_upper = max(data(:)) ;
        data_lower = 0 ;
        WindowWidth = data_upper ;
        WindowCenter = round(data_upper/2) ;
    case 'multi-echo'
        ImageType = 'ORIGINAL\PRIMARY\M_SE\M\SE' ;
        SeriesDescription = 'MATLAB Generated Multi-echo' ;
        fnstem = 'ME' ;
        
        szdata = size(data) ;
        necho = szdata(end) ;
        if ndims(data) == 3
            data = reshape(data,[szdata(1) szdata(2) 1 necho]) ;
        end
        
        
        data_max = max(data(:)) ;
        if data_max > intmax('uint16') 
            disp(['Data max exceeds uint16 - scaling'])
            scale = double(intmax('uint16')) / double(data_max) ;
        else
            scale = 1 ;
        end
        data = data * scale ;
        data_upper = max(data(:)) ;
        data_lower = 0 ;
        WindowWidth = data_upper ;
        WindowCenter = round(data_upper/2) ;
    case 'T1map'
        SeriesDescription = 'MATLAB Calculated T1 map' ;
        fnstem = 'T1' ;
        ImageType = 'DERIVED\SECONDARY\T1' ;  %?? Not sure of standard here
        
        data_upper = 10000 ;  % T1 values over this are clipped.
        data_lower = 0 ;
        
        data_max = max(data(:)) ;
        if data_max > data_upper ;
            warning(['Data max exceeds ',num2str(data_upper), ...
                ' will clip values above.'])
        end
        
        WindowWidth = 1850 ;
        WindowCenter = 1075 ;
%     case 'newgeom'
%         SeriesDescription = 'MATLAB Calculated New Geometry' ;
%         fnstem = 'ng' ;
%         data_upper = max(data(:)) ;
%         data_lower = min(data(:)) ;

    case 'rgb'
        SeriesDescription = 'RGB' ;
        fnstem = 'RGB' ;
        ImageType = 'DERIVED\SECONDARY\RGB' ;
        
        data_upper = 256 ;
        data_lower = 0 ;
        WindowWidth = 256 ;
        WindowCenter = 128 ;
        isrgb = true ;
        
    case 'fatfraction'
        SeriesDescription = 'Fat Fraction';
        fnstem = 'ff' ;
        ImageType = 'DERIVED\SECONDARY\CALCULATED' ;  %?? Not sure of standard here
        
        scale = 100 ; % Convert to percent
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 100 ;  % Values over this are written as 100.
        data_lower = 0 ;
        WindowWidth = 100 ;
        WindowCenter = 50 ;
        
    case 'LWF'
        SeriesDescription = 'Luminal Water Fraction (%)';
        fnstem = 'lwf' ;
        ImageType = 'DERIVED\SECONDARY\CALCULATED' ;  %?? Not sure of standard here
        
        scale = 100 ; % Convert to percent
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 100 ;  % Values over this are written as 100.
        data_lower = 0 ;
        WindowWidth = 50 ;
        WindowCenter = 25 ;
        
     case 'LWF_2'
        SeriesDescription = 'Luminal Water Fraction (x1000)';
        fnstem = 'lwf' ;
        ImageType = 'DERIVED\SECONDARY\CALCULATED' ;  %?? Not sure of standard here
        
        scale = 1000 ; % Convert to percent
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 1000 ;  % Values over this are written as 1000.
        data_lower = 0 ;
        WindowWidth = 400 ;
        WindowCenter = 200 ;

     case 'fICx1000'
        SeriesDescription = 'fIC (x1000)';
        fnstem = 'fIC' ;
        ImageType = 'DERIVED\PRIMARY\DIFFUSION\FIC' ;  % Not sure if PRIMARY or SECONDARY
        
        scale = 1000 ; % 
        
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 1000 ;  % Values over this are written as 1000.
        data_lower = 0 ;
        WindowWidth = 1000 ;
        WindowCenter = 500 ;

        % RescaleSlope and RescaleIntercept appear to not structly be part
        % of the classic DICOM format for MR, hence dicomwrite does not
        % ouput them when in Create mode.

    case 'fICWithRescale'  % see also grayWithRescale
        % Technically non-standard
        SeriesDescription = 'fIC (with RescaleSlope)' ;
        fnstem = 'fIC' ;
        ImageType = 'DERIVED\PRIMARY\DIFFUSION\FIC' ;  % Not sure if PRIMARY or SECONDARY
        
        scale = 1000 ;
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 1000 ;  % Values over this are written as 1000.
        data_lower = 0 ;
        WindowWidth = 1 ;
        WindowCenter = 0.5 ;
        RescaleSlope = 1/scale ;  
        RescaleIntercept = 0 ;

        forceAddRescale = true ;
        

        % RealWorldValueMappingSequence = buildRealWorldValueMappingSequence(data_upper, ...
        %     data_lower, RescaleSlope, RescaleIntercept ) ;


    case 'grayWithRescale'
        % Technically non-standard
        SeriesDescription = 'grayscale (with RescaleSlope)' ;
        fnstem = 'gray' ;
        ImageType = 'DERIVED\PRIMARY\SECONDARY\GRAY' ;  % Not sure if PRIMARY or SECONDARY
        
        scale = 1000 ;
        data(isnan(data))= 0 ;
        data = data * scale ;
        data_upper = 1000 ;  % Values over this are written as 1000.
        data_lower = 0 ;
        WindowWidth = 1 ;
        WindowCenter = 0.5 ;
        RescaleSlope = 1/scale ;  
        RescaleIntercept = 0 ;

        forceAddRescale = true ;
        
    case 'RawRecon'
        SeriesDescription = 'RawRecon' ;
        fnstem = 'RR' ;
        ImageType = 'DERIVED\SECONDARY\RAWRECON' ;
        
        data_max = max(data(:)) ;
        if data_max > intmax('uint16') ;
            disp(['Data max exceeds uint16 - scaling'])
            scale = double(intmax('uint16')) / double(data_max) ;
        else
            scale = 1 ;
        end
        data = data * scale ;
        
        data_lower = 0 ;
        data_upper = max(data(:)) ;
        WindowWidth = data_upper ;
        WindowCenter = floor(data_upper /2) ;
        
    otherwise
        error(['Unknown data_type: ',data_type])
end


header = 'blankMR';
anon = false ;
UIfolder = true ;
quiet = false ; 
burn_text = '' ;
ImageComments = '' ;
MRAcquisitionType = '' ;
ProtocolName = '' ;
PatientName = '' ;
ContentQualification = 'RESEARCH' ;
writePrivate = false ;
VR = 'implicit' ;
dict = dicomdict("get") ;
useInputPrivatePhilips = false ;
rgbAsPalette = true ;

TEs = [] ;

for ipv = 1: 2 :length(varargin)
    switch varargin{ipv}
        case 'SeriesDescription'
            SeriesDescription = varargin{ipv+1} ;
        case 'fnstem'
            fnstem = varargin{ipv+1} ;
        case 'geom'
            geom = varargin{ipv+1} ;
        case 'FrameOfReferenceUID'
            FrameOfReferenceUID = varargin{ipv+1} ;
        case 'ImageType'
            ImageType = varargin{ipv+1} ;
        case 'SeriesNumber'
            SeriesNumber = varargin{ipv+1} ;
        case 'StudyInstanceUID'
            StudyInstanceUID = varargin{ipv+1} ;
        case 'SeriesInstanceUID'
            SeriesInstanceUID = varargin{ipv+1} ;
        case 'header'
            header = varargin{ipv+1} ;
        case 'anon'
            anon = varargin{ipv+1} ;
        case 'quiet'
            quiet = varargin{ipv+1} ;
        case 'folder_name'
            UIfolder = false ;
            folder_name = varargin{ipv+1} ;
        case 'burn_text'
            burn_text= varargin{ipv+1} ;
        case 'ImageComments'
            ImageComments= varargin{ipv+1} ;
        case 'MRAcquisitionType'
            MRAcquisitionType= varargin{ipv+1} ;
        case 'ProtocolName'
            ProtocolName = varargin{ipv+1} ;
        case 'TEs'
            TEs = varargin{ipv+1} ;
        case 'PatientName'
            PatientName = varargin{ipv+1} ;
        case 'useInputPrivatePhilips'
            useInputPrivatePhilips = varargin{ipv+1} ;
        case 'WindowWidth'
            WindowWidth = varargin{ipv+1} ;
        case 'WindowCenter'
            WindowCenter = varargin{ipv+1} ;
        case 'writePrivate'
            writePrivate = varargin{ipv+1} ;
        case 'VR'
            VR = varargin{ipv+1} ;
        case 'dict'
            dict = varargin{ipv+1} ;
        case 'limitBitsStored'
            limitBitsStored = varargin{ipv+1} ;
        case 'rgbAsPalette'
            rgbAsPalette = varargin{ipv+1} ;

        otherwise
            warning(['Unknown parameter: ',varargin{ipv}])
    end
end

if rgbAsPalette && ~isempty(burn_text)
    warning('writeDicom: RGB as PALETTE COLOR does not have burn text implemeted')
end

writeArgs = {'WritePrivate', writePrivate, 'VR', VR,...
    'MultiframeSingleFile', false, 'Dictionary',dict} ;

% data [ny nx 3] or [ny nx nz 3]

if isrgb
    nt =1 ; nslice = 1 ;
    if ndims(data)>4
        error(['RGB Data can only be [ny nx 3] or [ny nx nz 3]'])
    end
    
    % Im will be used as Indexed image if rgbAsPalette i.e.
    % Photometric Interpretation as Palette Color
    Im = mat2gray(data) ;
    
    if ndims(data)==4
        if size(data,4) ~= 3
            error(['RGB Data can only be [ny nx 3] or [ny nx nz 3]'])
        end
        
        nslice = size(data,3) ;  
        
        [Xm,map] = rgb2ind(squeeze(Im(:,:,1,:)), 256) ;
        for islice = 2:nslice
            [Xm(:,:,islice)]=rgb2ind(squeeze(Im(:,:,islice,:)), map) ;
        end
        
    else
        [Xm,map] =  rgb2ind(Im,256) ;
    end
else
    nt = size(data,4) ;
    nslice = size(data,3) ;
end


% At this point data should be at a scale where it can be converted to uint16
% Clip data to be within bounds
data( data < data_lower) = data_lower ;

data( data > data_upper) = data_upper ;

if ~quiet && ~isdeployed
    hw  = warndlg({'DICOM files not for clinical use.' ; ...
        'No check here for correct geometry.' ; ...
        'DICOM header information may be incorrect' ;...
        'Output will overwite existing files with same name'} ,...
        'Not for clinical or commercial use') ;
    uiwait(hw, 3)
end


if ispref('writeDicom','save_dir')
    defdir = getpref('writeDicom','save_dir');
else
    defdir = [] ;
end

if UIfolder
    [folder_name] = uigetdir(defdir,'Select folder for output') ;
    
    if folder_name == 0
        warning(['No folder specified'])
        return
    else
        setpref('writeDicom','save_dir',folder_name)
    end
else
    if ~exist(folder_name,'dir')
        warning(['Output folder: ',folder_name,', does not exist.'])
        return
    end
end

new_SeriesInstanceUID = dicomuid ; % new fake series UID.
new_StudyInstanceUID = dicomuid ; % new StudyInstanceUID if needed
new_FoRUID = dicomuid ; % FrameofReference if needed

% switch data_type
%     case 'newgeom'
%         amin = min(data(:)) ;
%         amax = max(data(:)) ;
%         
%         for islice=1:nslice
%             dinfo_out.FrameOfReferenceUID = FrameOfReferenceUID ;
%             dinfo_out.StudyInstanceUID = new_studyID;
%             dinfo_out.SeriesInstanceUID = new_dicomuid ;
%             dinfo_out.SeriesDescription = SeriesDescription ;
%             dinfo_out.Width = geom(islice).Width ;
%             dinfo_out.Height = geom(islice).Height ;
%             dinfo_out.PixelSpacing = [geom(islice).PixelSpacing_HW(1) geom(islice).PixelSpacing_HW(2)];
%             dinfo_out.ImageOrientationPatient = geom(islice).IOP ;
%             dinfo_out.ImagePositionPatient = geom(islice).IPP ;
%             dinfo_out.SliceThickness = geom(islice).SliceThickness ;
%             
%             fn_out = fullfile(folder_name,[fnstem,'_',num2str(islice,'%03u'),'.dcm']) ;
%             status = dicomwrite(mat2gray(data(:,:,islice),[amin amax]), fn_out, ...
%                 dinfo_out,'ObjectType','MR Image Storage')  ;
%         end
%         
%     otherwise

data = uint16(round(data)) ;

for it = 1:nt
    for islice = 1:nslice
        InstanceNumber = (it-1)*nslice + islice ;
        if iscell(header)
            dinfo = header{1};
            locs = header{2} ;
            if size(locs,1) ~= nslice && islice==1 && it==1
                disp(['#Slices in data and locs do not match: header info may be incorrect'])
            end

            basefn = dinfo(locs(islice,it)).Filename ;
            
            dinfull = dicominfo(basefn,'dictionary',dict) ;
            
            dinfo_out = dinfull ;

            if useInputPrivatePhilips
                % Also adds EchoNumbers (Seems to be needed) and
                % PresentationLUTShape (strictly should check Photometric
                % interpretation
                dthis = dinfo(locs(islice,it)) ;
                dinfo_out.EchoNumbers = dthis.EchoNumbers ;
                dinfo_out.NumberOfTemporalPositions = dthis.NumberOfTemporalPositions ;

                dinfo_out.Private_2001_10xx_Creator = dthis.Private_2001_10xx_Creator ;
                dinfo_out.MRImagePhaseNumber = dthis.MRImagePhaseNumber ;
                dinfo_out.ImagePlaneNumber = dthis.ImagePlaneNumber ;

                dinfo_out.ImagePlaneOrientation = dthis.ImagePlaneOrientation ;

                dinfo_out.MRSeriesNrOfEchoes = dthis.MRSeriesNrOfEchoes ;
                dinfo_out.MRSeriesNrOfPhases = dthis.MRSeriesNrOfPhases ;
                dinfo_out.MRSeriesNrOfSlices = dthis.MRSeriesNrOfSlices ;
                dinfo_out.MRSeriesReconstructionNumber = dthis.MRSeriesReconstructionNumber ;
                dinfo_out.MRSeriesScanningTechniqueDesc = dthis.MRSeriesScanningTechniqueDesc ;
                dinfo_out.Stack = dthis.Stack ;
                dinfo_out.MRSeriesNrOfStacks = dthis.MRSeriesNrOfStacks ;
                dinfo_out.MRSeriesAcquisitionNumber = dthis.MRSeriesAcquisitionNumber ;
                dinfo_out.MRSeriesNrOfDynamicScans = dthis.MRSeriesNrOfDynamicScans ;

                dinfo_out.Private_2005_10xx_Creator = dthis.Private_2005_10xx_Creator ;
                dinfo_out.MRImageTypeMR = dthis.MRImageTypeMR ;
                dinfo_out.MRSeriesDataType = dthis.MRSeriesDataType ;
                dinfo_out.MRImageChemicalShiftNumber = dthis.MRImageChemicalShiftNumber ;
                dinfo_out.MRImageScanningSequencePrivate = dthis.MRImageScanningSequencePrivate ;
                dinfo_out.MRSeriesScanDuration = dthis.MRSeriesScanDuration ;
                dinfo_out.MRSeriesGeometryCorrection = dthis.MRSeriesGeometryCorrection ;
                dinfo_out.MRMeasurementScanResolution = dthis.MRMeasurementScanResolution ;

                dinfo_out.PresentationLUTShape = dthis.PresentationLUTShape ;
                dinfo_out.StudyDescription = dthis.StudyDescription ;
                dinfo_out.ReceiveCoilName = dthis.ReceiveCoilName ;
                dinfo_out.TemporalPositionIdentifier = dthis.TemporalPositionIdentifier ;
                dinfo_out.ScanningSequence = dthis.ScanningSequence ;
                dinfo_out.SequenceVariant = dthis.SequenceVariant ;
                dinfo_out.VolumetricProperties = dthis.VolumetricProperties ;
                dinfo_out.ImagedNucleus = dthis.ImagedNucleus ;
                dinfo_out.ImagingFrequency = dthis.ImagingFrequency ;
                dinfo_out.EchoTime = dthis.EffectiveEchoTime ;
                dinfo_out.RepetitionTime = dthis.RepetitionTime ;
                dinfo_out.Manufacturer = dthis.Manufacturer ;
                dinfo_out.SpacingBetweenSlices = dthis.SpacingBetweenSlices ;

            end

        else
            switch header
                case 'blankMR'
                    % Blank dinfo
                    dinfo_out.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4' ;
                    dinfull = dinfo_out ;
            end
        end
        
        CreateMode = 'Create' ;
        switch dinfull.SOPClassUID
            % read from multi-frame MR, writing as singleframe
            case '1.2.840.10008.5.1.4.1.1.4.1'
                dinfo_out = mf2single_keep(dinfo_out) ;
                dinfo_out.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4' ;
            case '1.2.840.10008.5.1.4.1.1.128'  % PET
                CreateMode = 'Copy'  ;
                dinfo_out.RescaleSlope = 1;
                dinfo_out.RescaleIntercept = 0 ;
        end
        
        % FrameOfReferenceUID
        if ~exist('FrameOfReferenceUID','var')
            dinfo_out.FrameOfReferenceUID = new_FoRUID ;
        else
            switch FrameOfReferenceUID
                case 'keep' % copied from dinfull
                case 'new'
                    dinfo_out.FrameOfReferenceUID = new_FoRUID ;
                otherwise
                    dinfo_out.FrameOfReferenceUID = FrameOfReferenceUID ;
            end
        end
        
        % SeriesInstanceUID
        if ~exist('SeriesInstanceUID','var')
            dinfo_out.SeriesInstanceUID = new_SeriesInstanceUID ;
        else
            switch SeriesInstanceUID
                case 'keep' % copied from dinfull
                case 'new'
                    dinfo_out.SeriesInstanceUID = new_SeriesInstanceUID ;
                otherwise
                    dinfo_out.SeriesInstanceUID = SeriesInstanceUID ;
            end
        end
        
        if ~exist('StudyInstanceUID','var')
            %dinfo_out.StudyInstanceUID = new_StudyInstanceUID;
            % default is now to keep
        else
            switch StudyInstanceUID
                case 'keep'
                case 'new'
                    dinfo_out.StudyInstanceUID = new_StudyInstanceUID;
                otherwise
                    dinfo_out.StudyInstanceUID =  StudyInstanceUID ;
            end
        end
        
        
        if ~exist('SeriesNumber','var')
            if isfield(dinfull,'SeriesNumber') 
              dinfo_out.SeriesNumber = dinfull.SeriesNumber + 20;
            else
              dinfo_out.SeriesNumber = 20;
            end
        else
            dinfo_out.SeriesNumber = SeriesNumber ;
        end
        
        dinfo_out.Width = geom(islice).Width ;
        dinfo_out.Height = geom(islice).Height ;
        dinfo_out.PixelSpacing = [geom(islice).PixelSpacing_HW(1) geom(islice).PixelSpacing_HW(2)];
        dinfo_out.ImageOrientationPatient = geom(islice).IOP ;
        dinfo_out.ImagePositionPatient = geom(islice).IPP ;
        dinfo_out.SliceThickness = geom(islice).SliceThickness ;
        % May not correspond to originals dinfo_out.SliceLocation = geom(islice).D ;
        
        dinfo_out.ImageType = ImageType ;
        dinfo_out.SeriesDescription = SeriesDescription ;
        dinfo_out.WindowWidth = WindowWidth ;
        dinfo_out.WindowCenter = WindowCenter ;
        
        dinfo_out.ImageComments = ImageComments ;

        dinfo_out.InstanceNumber = InstanceNumber ;
        
        % Old code when not newgeom
        %dinfo_out.ImagePositionPatient = dinfo(locs(islice,1)).ImagePositionPatient ;
        %dinfo_out.ImageOrientationPatient = dinfo(locs(islice,1)).ImageOrientationPatient ;
        
        if exist('RescaleSlope','var') 
            dinfo_out.RescaleSlope = RescaleSlope ;
        end

        if exist('RescaleIntercept','var') 
            dinfo_out.RescaleIntercept = RescaleIntercept ;
        end

        if exist('RealWorldMappingValueSequence','var')
            dinfo_out.RealWorldMappingValueSequence = RealWorldMappingValueSequence ;
        end

        if ~isempty(ContentQualification)
            dinfo_out.ContentQualification = ContentQualification ;
        end

        if ~isempty(MRAcquisitionType)
            dinfo_out.MRAcquisitionType = MRAcquisitionType ;
        end
        
        if ~isempty(ProtocolName)
            dinfo_out.ProtocolName = ProtocolName ;
        end
        
        if ~isempty(PatientName)
            dinfo_out.PatientName = PatientName ;
        end
        
        if anon
            dinfo_out.PatientName = '';
            dinfo_out.PatientID = '' ;
            dinfo_out.PatientBirthDate = '' ;
        end
        
        if nt==1
            fn_out = fullfile(folder_name,[fnstem,'_',num2str(islice,'%03u'),'.dcm']) ;
        else
            fn_out = fullfile(folder_name,[fnstem,'_',num2str(islice,'%03u'),'_',num2str(it,'%03u'),'.dcm']) ;
        end
        if isrgb
            
            if rgbAsPalette
                if ndims(data)==4
                    % Xm should be indexed 
                    im2write = burn_in(squeeze(Xm(:,:,islice,:)), burn_text, map) ;
                    status = dicomwrite(im2write, map, fn_out, dinfo_out, ...
                        'CreateMode', CreateMode, writeArgs{:} ) ;
                else
                    im2write = burn_in(Xm, burn_text, map) ;
                    status = dicomwrite(im2write, map, fn_out, dinfo_out, ...
                        'CreateMode', CreateMode, writeArgs{:} ) ;
                end
            else % Output as DICOM RGB
                dinfo_out.ColorType = 'truecolor' ;
                dinfo_out.PhotometricInterpretation = 'RGB' ;
                % Needs to be uint8 for Philips scanner reload
                im2write = uint8(round(squeeze(data(:,:,islice,:)))) ;
                status = dicomwrite(im2write, fn_out, dinfo_out, ...
                    'CreateMode', CreateMode, writeArgs{:} ) ;

            end % end rgbAsPalette
        else 
            im2write = burn_in(data(:,:,islice,it), burn_text) ;
            switch data_type
                case 'multi-echo'
                    dinfo_this = dinfo_out ;
                    dinfo_this.EchoNumbers = it ;
                    dinfo_this.EchoTime = 1000*TEs(it) ;
                otherwise
                    dinfo_this = dinfo_out ;
            end
            if forceAddRescale
                tmpFile = tempname ;
                status = dicomwrite(im2write, tmpFile, dinfo_this, ...
                    'CreateMode', 'Create' , writeArgs{:}) ;
                
                dtemp = dicominfo(tmpFile) ;
                dtemp.RescaleSlope = dinfo_this.RescaleSlope ;
                dtemp.RescaleIntercept = dinfo_this.RescaleIntercept ;
                if isfield(dinfo_this,'PresentationLUTShape')
                    dtemp.PresentationLUTShape = dinfo_this.PresentationLUTShape ; % strictly this is not part of forceAddRescale
                end
                if isfield(dinfo_this,'VolumetricProperties')
                    dtemp.VolumetricProperties = dinfo_this.VolumetricProperties ; % strictly from multi-frame?
                end
                
                if limitBitsStored
                    dtemp.BitsStored = 12 ; 
                    dtemp.HighBit = 11 ;
                    useMetadataBitDepths = true ;
                else
                    useMetadataBitDepths = false ;
                end
                status = dicomwrite(im2write, fn_out, dtemp, ...
                    'CreateMode', 'Copy', writeArgs{:}, ...
                    'UseMetadataBitDepths', useMetadataBitDepths ) ;

                delete(tmpFile)

            else
                status = dicomwrite(im2write, fn_out, dinfo_this, ...
                    'CreateMode', CreateMode, writeArgs{:} ) ;
            end
        end
    end % islice
end % it


disp(['Written ',num2str(nslice*nt),' files to folder: ',folder_name])

end

function im2write = burn_in(im2D, text, map) 
im2write = im2D ;
if isempty(text)
    return
end
[ny, nx, nz] = size(im2D) ;

% text2im is from 
% https://uk.mathworks.com/matlabcentral/fileexchange/19896-convert-text-to-an-image
% 
imtext = text2im(text) ;
[nyt, nxt] = size(imtext) ;

ys = ny/nyt; xs = nx/nxt ;

scale = min([1 ys xs]) ;
imtexts = imresize(imtext,scale) ; 
imtexts = mat2gray(imtexts) ; % removes slight under or overshoot from interp

if nargin > 2 % map, 
    loclet = imtexts <= 0.5  ;
    locbak = imtexts > 0.5 ;
    
    mapsum = sum(map,2);
    [~,locmapb] = max(mapsum,[],1) ; % brightest color in map
    [~,locmapl] = min(mapsum,[],1) ; % darkest color in map
    
    if isa(im2D,'uint8') || isa(im2D,'uint16')
        locmapb = locmapb - 1 ;
        locmapl = locmapl - 1 ;
    end
    
    imtexts(loclet) = locmapl ;
    imtexts(locbak) = locmapb ;
    
else
    
    rng = double(max(im2D(:)) - min(im2D(:))) ;
    if rng == 0
        rng = 1;
    end
    imtexts = imtexts*rng + double(min(im2D(:))) ; % double
end

imtexts = cast(imtexts, 'like', im2D) ;

[ht, wt] = size(imtexts) ;

llr = ny; 
llc = max([1 floor(nx/2 - wt/2)]) ;
im2write(llr-ht+1:llr, llc:llc+wt-1) = imtexts ;


end

function dinfo_out = mf2single_keep(dinfo_in)

tokeep = { ...
    'AccessionNumber'
    'AcquisitionContrast'
    'AcquisitionDateTime'
    'AcquisitionDuration'
    'AcquisitionMatrix'
    'AcquisitionTime'
    'BodyPartExamined'
    'ColorType'
    'Columns'
    'DeidentificationMethod'
    'DeidentificationMethodCodeSequence' 
    'DeviceSerialNumber'
    'DiffusionBValue'
    'DiffusionGradientOrientation'
    'EchoNumbers'
    'EchoPulseSequence'
    'EchoTime'
    'EchoTrainLength'
    'FlipAngle'
    'Format'
    'FrameOfReferenceUID'
    'Height'
    'HighBit'
    'ImagePlaneOrientation'
    'ImageType'
    'ImagedNucleus'
    'ImagingFrequency'
    'InstanceNumber'
    'InstitutionalDepartmentName'
    'InstitutionName'
    'InversionTime'
    'LossyImageCompression'
    'MagneticFieldStrength'
    'Manufacturer'
    'ManufacturerModelName'
    'Modality'
    'MRAcquisitionType'
    'MultipleSpinEcho'
    'NumberOfAverages'
    'NumberOfPhaseEncodingSteps'
    'NumberOfTemporalPositions'
    'ParallelAcquisition'
    'ParallelAcquisitionTechnique'
    'ParallelReductionFactorInPlane'
    'PartialFourier'
    'PatientAge'
    'PatientBirthDate'
    'PatientID'
    'PatientIdentityRemoved'
    'PatientName'
    'PatientPosition'
    'PatientSex'
    'PatientWeight'
    'PercentPhaseFieldOfView'
    'PercentSampling'
    'PerformedStationAETitle'
    'PhotometricInterpretation'
    'PixelBandwidth'
    'PixelRepresentation'
    'PixelSpacing'
    'PresentationLUTShape'
    'ProtocolName'
    'PulseSequenceName'
    'ReceiveCoilName'
    'RepetitionTime'
    'RequestedContrastAgent'
    'ResonantNucleus'
    'Rows'
    'SamplesPerPixel'
    'ScanningSequence'
    'ScanOptions'
    'SeriesDate'
    'SeriesDescription'
    'SeriesNumber'
    'SeriesTime'
    'SequenceVariant'
    'SliceLocation'
    'SliceThickness'
    'SoftwareVersions'
    'SpacingBetweenSlices'
    'StationName'
    'StudyComments'
    'StudyDate'
    'StudyDescription'
    'StudyID'
    'StudyInstanceUID'
    'StudyTime'
    'TemporalPositionIdentifier'
    'TransmitCoilName'
    'VolumetricProperties'
    'Width' 
    'Private_2001_10xx_Creator'
    'MRImageChemicalShiftNumber'
    'MRImagePhaseNumber'
    'ImagePlaneNumber'
    'ImagePlaneOrientation'
    'MRSeriesNrOfEchoes'
    'MRSeriesNrOfPhases'
    'MRSeriesNrOfSlices'
    'MRSeriesReconstructionNumber'
    'MRSeriesScanningTechniqueDesc'
    'Stack'
    'MRSeriesNrOfStacks'
    'MRSeriesAcquisitionNumber'
    'MRSeriesNrOfDynamicScans'
    'Private_2005_10xx_Creator'
    'MRImageTypeMR'
    'MRSeriesDataType'
    'MRImageScanningSequencePrivate'
    'MRSeriesGeometryCorrection'
    'MRMeasurementScanResolution'} ;



for ifield = 1:length(tokeep)
    if isfield(dinfo_in, tokeep(ifield))
        dinfo_out.(tokeep{ifield}) = dinfo_in.(tokeep{ifield}) ;
    end
end
end


