function [dinfo ] = data_process(dinfo, THRESH, verbose)
% DATA_PROCESS Checks slices parallel, computes order and slice centre separations
%  Performs some other DICOM tidying up.
%
% [dinfo ] = data_process(dinfo)
% [dinfo ] = data_process(dinfo, THRESH)
% [dinfo ] = data_process(dinfo, THRESH, verbose)
% [dinfo ] = data_process(dinfo, [], verbose)
%
% THRESH defaults to 0.01. Previously 0.001 could fail for PET from mMR Biograph.
% Also, a higher value needed for PAR/REC due to rounding errrors?
%
% verbose {true} | false
%
% Code adapted from dicom2vol, called by datparse (SingleFrame) or 
% dmfparse (MultiFrame) or dfparse (from dselector).
%
% Provides an integer slice index that orders slices correctly in a
% right-handed coordinate system.
%
% itype enumerated ImageType (see code). Enumeration also listed in d2mat
%
% Some code repeated in geom_check
%
% Copyright, 2019-2021, David Atkinson.
% D.Atkinson@ucl.ac.uk
%
% See also DFPARSE D2MAT DSELECTOR DATPARSE DMFPARSE 

if nargin < 2 || isempty(THRESH)
    THRESH = 0.01 ; % Slice position differences less than THRESH are
    % considered to be the same slice
end

iopatol = 1e-4 ; % Absolute tolerance for differences in entries to IOP

if nargin < 3
    verbose = true ;
end

if ~isfield(dinfo,'SeriesNumber'), return, end

sn = [dinfo.SeriesNumber] ;
scans = unique(sn) ;
if verbose, disp(['The series numbers are: ',num2str(scans)]) , end

nseries = length(scans) ;
for iseries = 1:nseries
    fok = find(sn == scans(iseries)) ;
    if isempty(fok)
        disp(['No files for scan ',num2str(scans(iseries))])
        break
    end
    
    dinfo_s = dinfo(fok) ;
    isSurvey = false ;

    if isfield(dinfo_s,'ProtocolName') 
        if verbose
            disp(['Processing series: ',num2str(scans(iseries)),'. ', ...
                dinfo_s(1).ProtocolName])
        end

    end
    
    nslice = length(dinfo_s) ;
    iop = dinfo_s(1).ImageOrientationPatient;
    
    rdc = iop(1:3) ; % row direction cosine
    cdc = iop(4:6) ;
    sdc = cross(rdc,cdc) ; % slice direction cosine
    
    slicepos = zeros([nslice 1]) ;
    ipp = zeros([nslice 3]) ;
    
    for ifile = 1: nslice
        % Check IOPs are equal to within absolute tolerance of 1e-4
        iopadiff = abs(iop - dinfo_s(ifile).ImageOrientationPatient) ;
        if ~isempty(find(iopadiff> iopatol))
            % Slices are not parallel - issue warning unless survey scan

            if isfield(dinfo_s,'ProtocolName') && contains(dinfo_s(ifile).ProtocolName,'SURVEY')
            else
                warning([' Slices not parallel in: ', dinfo_s(ifile).Filename] )
            end

        end
        slicepos(ifile) = dot(sdc,dinfo_s(ifile).ImagePositionPatient) ;
    end
    
    [srtslicepos, idx] = sort(slicepos) ;
    % Assign integer slice index, check for same slice more than
    % once, e.g. in dynamics
    clear idx_sl
    idx_sl(1) = 1;
    firstpass = 1 ;
    slc2c = dinfo(fok(idx(1))).SliceThickness ; % overwritten if more than one slice
    
    for islice  = 2:nslice
        sdiff = srtslicepos(islice)-srtslicepos(islice-1) ;
        if sdiff < THRESH  % same slice pos
            idx_sl(islice) = idx_sl(islice-1) ;
        else
            if firstpass
                slsep = sdiff; % set slice separation
                slc2c = slsep ;
                firstpass = 0 ;
            else
                if abs(sdiff-slsep) > THRESH
                    slc2c = NaN;
                    if isfield(dinfo_s,'ProtocolName') && contains(dinfo_s(ifile).ProtocolName,'SURVEY')
                    else
                        warning([' Slice separations not equal, ',...
                            num2str(sdiff),'  slsep ',num2str(slsep)] )
                    end
                end
            end
            idx_sl(islice) = idx_sl(islice-1)+1;
        end
    end
    
    if verbose, disp([' ',num2str(idx_sl(end)),' unique slice positions out of ',...
        num2str(nslice),' total slices.']), end
    
    % now assign these slice indices to the relevant dinfo.
    for ifok = 1:length(fok)
        dinfo(fok(idx(ifok))).sl = idx_sl(ifok) ;
        dinfo(fok(idx(ifok))).slc2c = slc2c ;
    end
    
    % Diffusion data
    if isfield(dinfo_s,'DiffusionBValue')
        bv = [dinfo_s.DiffusionBValue] ;
        ubv = unique(bv) ;
        if ~isempty(ubv) && verbose
            disp([' Series contains b-values: ',num2str(ubv)])
        end
    elseif isfield(dinfo_s,'Private_0019_100c')  % Siemens single frame Imanova
        bv = [dinfo_s.Private_0019_100c] ;       % unclear of meaning from UCH PACs
        ubv = unique(bv) ;
        
        if ~isempty(ubv)
            if size(ubv,1) == 1 % For Siemens ADC files, this can be a column vector
                % (no idea why) exclude in this case.
                if verbose, disp([' Series contains b-values: ',num2str(ubv)]), end
                
                for iim = 1:nslice
                    dinfo(fok(iim)).DiffusionBValue = bv(iim) ;
                end
            end
        end
    end

    if isfield(dinfo_s,'DiffusionDirectionality')
        % use this
    elseif isfield(dinfo_s,'Private_0019_100d')
        ddty = {dinfo_s.Private_0019_100d} ;
        for iim = 1:nslice
            switch ddty{iim}
                case 'ISOTROPIC'
                    dinfo(fok(iim)).DiffusionDirectionality = 2 ;
                case 'DIRECTIONAL'
                    dinfo(fok(iim)).DiffusionDirectionality = 1 ;
                case 'NONE'
                    dinfo(fok(iim)).DiffusionDirectionality = 0 ;
                otherwise
                    warning(['Diff Directionality problems'])
                    dinfo(fok(iim)).DiffusionDirectionality = 0 ;
            end

        end
    end
    
    % Dixon and others
    if isfield(dinfo_s,'ImageType')
        it = {dinfo_s.ImageType};
        uit = unique(it) ;
        if ~isempty(uit)
            if verbose
                disp(['Unique Image Types present:'])
                for iuit = 1:length(uit)
                    disp(['  ',uit{iuit}])
                end
            end
            
            % rename Image Types for easier processing in d2mat
            loc_w = strcmp(it,'DERIVED\PRIMARY\W\W\DERIVED') ;
            loc_f = strcmp(it,'DERIVED\PRIMARY\F\F\DERIVED') ;
            loc_ip = strcmp(it,'DERIVED\PRIMARY\IP\IP\DERIVED') ;
            loc_op = strcmp(it,'DERIVED\PRIMARY\OP\OP\DERIVED') ; % see 22 for fat fraction
            
            loc_mffe = strcmp(it,'ORIGINAL\PRIMARY\M_FFE\M\FFE') ;
            loc_mir = strcmp(it,'ORIGINAL\PRIMARY\M_IR\M\IR') ;
            loc_t2map  = strcmp(it,'ORIGINAL\PRIMARY\T2 MAP\T2\UNSPECIFIED') ; %7
            loc_rir   =  strcmp(it,'ORIGINAL\PRIMARY\R_IR\R\IR') ; %8
            loc_pca = strcmp(it,'ORIGINAL\PRIMARY\VELOCITY MAP\P\PCA') ;
            loc_b1 = strcmp(it,'ORIGINAL\PRIMARY\M_B1\M\B1') ;  %10
            loc_se = strcmp(it,'ORIGINAL\PRIMARY\M_SE\M\SE') ;
            loc_pmap = strcmp(it,'ORIGINAL\PRIMARY\PHASE MAP\P\SE') ;
            loc_mpca = strcmp(it,'ORIGINAL\PRIMARY\M_PCA\M\PCA') ;
            loc_pffe = strcmp(it,'ORIGINAL\PRIMARY\PHASE MAP\P\FFE') ; %14
            loc_pb1  = strcmp(it,'ORIGINAL\PRIMARY\PHASE MAP\P\B1') ; % 15
            loc_pir  = strcmp(it,'ORIGINAL\PRIMARY\PHASE MAP\P\IR') ; % 16
            loc_rffe = strcmp(it,'ORIGINAL\PRIMARY\R_FFE\R\FFE') ; % 17
            loc_iir  = strcmp(it,'ORIGINAL\PRIMARY\I_IR\I\IR') ; % 18
            loc_iffe = strcmp(it,'ORIGINAL\PRIMARY\I_FFE\I\FFE') ; % 19
            loc_b0map = strcmp(it, 'ORIGINAL\PRIMARY\B0 MAP\B0\UNSPECIFIED') ; % 20
            loc_t2star = strcmp(it, 'ORIGINAL\PRIMARY\T2_STAR_UNSPECIF\T2_STAR\UNSPECIFIED') ; % 21
            loc_ff  = strcmp(it, 'DERIVED\PRIMARY\FF\FF\DERIVED') ; % 22 fat fraction
            loc_other  = strcmp(it, 'ORIGINAL\PRIMARY\OTHER') ; % 23 other (GE?)
            loc_petmr  = strcmp(it, 'ORIGINAL\PRIMARY\M\ND\NORM') ; % 24
            loc_proj   = strcmp(it, 'ORIGINAL\PRIMARY\PROJECTION IMAGE\M\FFE') ; % 25 MPR Philips
            loc_adc = strcmp(it, 'ORIGINAL\PRIMARY\ADC_UNSPECIFIED\ADC\UNSPECIFIED') ; % 26 ADC
            loc_axial = strcmp(it, 'ORIGINAL\SECONDARY\AXIAL' ); % 27 Axial CT
            loc_mfsplit = strcmp(it, 'ORIGINAL\PRIMARY\M\DIS2D\MFSPLIT') ; % 28 Siemens Mutiframe split(!
            loc_t1mod = strcmp(it, 'ORIGINAL\PRIMARY\T1\M') ; % 29 FFE modulus
            loc_t1phase = strcmp(it, 'ORIGINAL\PRIMARY\T1\P') ; % 30 FFE phase



            % Decided to enumerate to make d2mat easier. Needs refactoring.
            for jit = 1:length(it)
                if loc_w(jit)
                    dinfo(fok(jit)).wfio = 1 ;
                    dinfo(fok(jit)).itype = 1 ;
                elseif loc_f(jit)
                    dinfo(fok(jit)).wfio = 2 ;
                    dinfo(fok(jit)).itype = 2 ;
                elseif loc_ip(jit)
                    dinfo(fok(jit)).wfio = 3 ;
                    dinfo(fok(jit)).itype = 3 ;
                elseif loc_op(jit)
                    dinfo(fok(jit)).wfio = 4 ;
                    dinfo(fok(jit)).itype = 4 ;
                elseif loc_mffe(jit)
                    dinfo(fok(jit)).itype = 5 ;
                elseif loc_mir(jit)
                    dinfo(fok(jit)).itype = 6 ;
                elseif loc_t2map(jit)
                    dinfo(fok(jit)).itype = 7 ;
                elseif loc_rir(jit)
                    dinfo(fok(jit)).itype = 8 ;
                elseif loc_pca(jit)
                    dinfo(fok(jit)).itype = 9 ;
                elseif loc_b1(jit)
                    dinfo(fok(jit)).itype = 10 ;
                elseif loc_se(jit)
                    dinfo(fok(jit)).itype = 11 ;
                elseif loc_pmap(jit)
                    dinfo(fok(jit)).itype = 12 ;
                elseif loc_mpca(jit)
                    dinfo(fok(jit)).itype = 13 ;
                elseif loc_pffe(jit)
                    dinfo(fok(jit)).itype = 14 ;
                elseif loc_pb1(jit)
                    dinfo(fok(jit)).itype = 15 ;
                elseif loc_pir(jit)
                    dinfo(fok(jit)).itype = 16 ;
                elseif loc_rffe(jit)
                    dinfo(fok(jit)).itype = 17 ;
                elseif loc_iir(jit)
                    dinfo(fok(jit)).itype = 18 ;
                elseif loc_iffe(jit)
                    dinfo(fok(jit)).itype = 19 ;
                elseif loc_b0map(jit)
                    dinfo(fok(jit)).itype = 20 ;
                elseif loc_t2star(jit)
                    dinfo(fok(jit)).itype = 21 ;
                elseif loc_ff(jit)
                    dinfo(fok(jit)).itype = 22 ;
                elseif loc_other(jit)
                    dinfo(fok(jit)).itype = 23 ;
                elseif loc_petmr(jit)
                    dinfo(fok(jit)).itype = 24 ;
                elseif loc_proj(jit)
                    dinfo(fok(jit)).itype = 25 ;
                elseif loc_adc(jit)
                    dinfo(fok(jit)).itype = 26 ;
                elseif loc_axial(jit)
                    dinfo(fok(jit)).itype = 27 ;
                elseif loc_mfsplit(jit)
                    dinfo(fok(jit)).itype = 28 ;
                elseif loc_t1mod(jit)
                    dinfo(fok(jit)).itype = 29 ;
                elseif loc_t1phase(jit)
                    dinfo(fok(jit)).itype = 30 ;
                else
                    dinfo(fok(jit)).itype = -1 ;
                end
            end
            
            num_it = [dinfo.itype] ;
            unum_it = unique(num_it) ;
            if ~isempty(unum_it) && verbose
                disp(['Unique numerical itypes present: '])
                for iun = 1:length(unum_it)
                    disp(['   ',num2str(unum_it(iun))])
                end
            end
        end
    end
    
    % MRImageLabeType for ASL. See below for Single Frame DICOM
    if isfield(dinfo_s,'MRImageLabelType')
        it = {dinfo_s.MRImageLabelType};
        uit = unique(it) ;
        if ~isempty(uit)
            if verbose
                disp(['Unique Image Label (ASL) Types present:'])
                for iuit = 1:length(uit)
                    disp(['  ',uit{iuit}])
                end
            end
            
            % rename Image LabelTypes for easier processing in d2mat
            loc_lab = strcmp(it,'LABEL') ;
            loc_cont = strcmp(it,'CONTROL') ;
            
            for jit = 1:length(it)
                if loc_lab(jit)
                    dinfo(fok(jit)).ilabtype = 1 ;
                elseif loc_cont(jit)
                    dinfo(fok(jit)).ilabtype = 2 ;
                end
            end
        end
    end
    
    % MRImageLabeType for ASL  # MJS added Private_2005_1429 for single
    % frame images, the MRImageLabelType did not work for single frame
    if isfield(dinfo_s,'Private_2005_1429')
        it = {dinfo_s.Private_2005_1429};
        temp = it{:} ;
        if ~isempty(temp)
            uit = unique(it) ;
            if ~isempty(uit)
                if verbose
                    disp(['Unique Image Label Types present:'])
                    for iuit = 1:length(uit)
                        disp(['  ',uit{iuit}])
                    end
                end

                % rename Image LabelTypes for easier processing in d2mat
                loc_lab = strcmp(it,'LABEL') ;
                loc_cont = strcmp(it,'CONTROL') ;

                for jit = 1:length(it)
                    if loc_lab(jit)
                        dinfo(fok(jit)).ilabtype = 1 ;
                    elseif loc_cont(jit)
                        dinfo(fok(jit)).ilabtype = 2 ;
                    end
                end
            end
        end
    end
    
    % Repetition Time. Below is for single-frame.
    % For multi-frame, set in dmfparse.
    
    TRwarn = 1 ;
    if isfield(dinfo_s,'RepetitionTime')
        RTs = [dinfo_s.RepetitionTime] ;
        uRTs = unique(RTs) ;
        if verbose, disp(['RepetitionTimes are: ',num2str(uRTs(:)')]), end
        
        for itr = 1:length(RTs)
            if dinfo_s(itr).RepetitionTime == 0
                isPriv = isfield(dinfo_s,'Private_2005_1030') ;
                if TRwarn && isPriv
                    disp(['Setting TR to first value in Private Field'])
                    TRwarn = 0 ;
                end
                if isPriv
                    dinfo(fok(itr)).RepetitionTime = dinfo(fok(itr)).Private_2005_1030(1) ;
                end
            end
        end
    end
    
    %FlipAngle. Make more precise if Philips MRSeriesFlipAngle present
    %(otherwise can loose fractional parts)
    % Only for single-frame, multi is OK.
    if isfield(dinfo_s,'Private_2001_1023')
        for ifa = 1:nslice
            dinfo(fok(ifa)).FlipAngle = dinfo(fok(ifa)).Private_2001_1023 ;
        end
    end
    
    %RescaleSlope and Rescale Intercept. Can be missing in some PACs data.
    %Use RealWorkdValueMapping if present (does appear in some PACs). If 
    % the Private versions were extracted (currently code does not do this
    % for multi-frame), use Private value in preference. 
    % Note dmfparse currently creates a RescaleSlope and Intercept even 
    % if not in DICOM, but it seems they always are for Philips multi-frame.
    
    % Set Rescale from RealWorld if present
    % Set Rescale from Private if present.
    
    srwi = 'RealWorldValueIntercept';
    srws = 'RealWorldValueSlope' ;
    
    % dinfo = setf(dinfo_s, dinfo, fok, nslice, fpref, f2, fset, def)
    dinfo = setf(dinfo_s, dinfo, fok, nslice, srwi, 'RescaleIntercept', 'RescaleIntercept', 0, verbose) ;
    dinfo = setf(dinfo_s, dinfo, fok, nslice, 'Private_2005_1409', 'RescaleIntercept', 'RescaleIntercept', 0, verbose) ;
    
    dinfo = setf(dinfo_s, dinfo, fok, nslice, srws, 'RescaleSlope', 'RescaleSlope', 1, verbose) ;
    dinfo = setf(dinfo_s, dinfo, fok, nslice, 'Private_2005_140a', 'RescaleSlope', 'RescaleSlope', 1, verbose) ;
    
    % DA - set fok index to 1, not islice below
    if ~isfield(dinfo,'Private_2005_100e') || isempty(dinfo(fok(1)).Private_2005_100e) % ScaleSlope
        % Not in input, replace with 1/RescaleSlope (which is now set) for 'fp'
        % output from d2mat
        %
        % Additional test on empty is for Siemens data when there are
        % multiple series - first series wiil add the field, but it is
        % empty for 2nd unless this check
        for islice = 1:nslice
            dinfo(fok(islice)).Private_2005_100e = 1/dinfo(fok(islice)).RescaleSlope ;
        end
    end
    
    
    if isfield(dinfo,'AcquisitionTime') && isfield(dinfo,'SeriesTime')
        for islice =1:nslice
            dinfo(fok(islice)).TinSeries = round(1e6*(dinfo(fok(islice)).AcquisitionTime - ...
                dinfo(fok(islice)).SeriesTime)) ;
        end
        
        tis = [dinfo(fok(1:nslice)).TinSeries] ;
        dtis = diff(tis) ;
        loc_dtisnz = find(dtis > 0) ;
        if ~isempty(loc_dtisnz)
            dtisnz = dtis(loc_dtisnz) ;
            disp(['Mean dt for AcqTime: ',num2str(mean(dtisnz)/1e6),'s.'])
        end
    end
    
    if verbose, disp(' '), end
end % series

if isfield(dinfo,'Private_0019_100c')
    dinfo = rmfield(dinfo,'Private_0019_100c') ;
end

if isfield(dinfo,'Private_2001_1023')
    dinfo = rmfield(dinfo,'Private_2001_1023') ;
end

end



function dinfo = setf(dinfo_s, dinfo, fok, nslice, fpref, f2, fset, def, verbose)
% SETF Set field for RescaleIntercept and RescaleSlope
% Aim is to favour RealWorldValue over Rescale, then if Private present,
% use that.

% Need 2nd checks here because RescaleIntercept and Slope fields can get
% set when more that one series being read
isfpref = isfield(dinfo_s, fpref) && ~isempty([dinfo_s.(fpref)]) ;
isf2 = isfield(dinfo_s, f2) && ~isempty([dinfo_s.(f2)]) ;

if ~isfpref && ~isf2
    if verbose, disp([fset,' using default value: ',num2str(def)]), end
    for islice = 1:nslice
        dinfo(fok(islice)).(fset) = def ;
    end
end

if isfpref && ~isf2
    vals = zeros([1 nslice]) ;
    for islice = 1:nslice
        dinfo(fok(islice)).(fset) = dinfo(fok(islice)).(fpref) ;
        vals(islice) = dinfo(fok(islice)).(fpref) ;
    end
    if verbose, disp([ fset,' using ', fpref,', unique vals: ',num2str(unique(vals))]), end
end

if isf2 && ~isfpref
    vals = zeros([1 nslice]) ;
    for islice = 1:nslice
        dinfo(fok(islice)).(fset) = dinfo(fok(islice)).(f2) ;
        vals(islice) = dinfo(fok(islice)).(f2) ;
    end
    if verbose, disp([ fset,' using ',f2,', unique vals: ',num2str(unique(vals))]), end
end

if isf2 && isfpref
    disp([ fset,' both ',fpref,' and ',f2,' present'])
    warnf = 0 ;
    for islice = 1:nslice
        v1 = dinfo(fok(islice)).(fpref) ;
        v2 = dinfo(fok(islice)).(f2) ;
        if warnf ==0 && abs(v1-v2)/max([abs(v1) abs(v2) eps]) > 0.01
            warnf = 1;
            disp([fset,' being set to value in ',fpref])
            disp([f2,': ',num2str(v2)])
            disp([fpref,': ',num2str(v1)])
        end
        dinfo(fok(islice)).(fset) = v1 ;
    end
end
end


