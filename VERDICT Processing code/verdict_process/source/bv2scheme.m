function scheme = bv2scheme(bv,scanner, ddeltas)
% BV2SCHEME b-values to scheme
%
% scheme = bv2scheme(bv,scanner)
% scheme = bv2scheme(bv,scanner, ddeltas) when scanner='PDS'
%
% bv should be one of 90, 500, 1500, 2000, 3000 for Classic VERDICT
% scanner should be one of 'Ingenia', 'Achieva', 'SIGNA Premier', 'XNAT',
% 'Ingeniav2.0', 'IngeniaRenal', 'PDS'
%  'XNAT' is intended to be most similar to the prostate-XNAT processing 
%   up to at least January 2023.
%
%  'SIGNA Premier' data is from Cambridge (January 2023)
%  'PDS' requires ddeltas input (e.g. from Philips PDS DICOM XX file)
%
% scheme is a structure with fields bval, G, delta and DELTA.
%
% For all scanner values except 'XNAT':
%   scheme.bval is returned as the input bv
%   scheme.G    is computed using stejskal from delta, DETLA and bv
%
% For scanner 'XNAT':
%   scheme.G    is supplied here (taken from a scaner patch or seq dev mode)
%   scheme.bval is computed using stejskal from delta, DETLA and above G
%
%
% Scanners generally try hard to supply the requested b-value by taking
% into account the trapezoidal nature of gradients and possibly also the
% crushers around the 180 pulse. This means the typed-in b-value should be
% more accurate than one computed from G.
% The value of G from the scanner (diffusion grad strength) is not used for
% any scheme here, except 'XNAT' in order to match the prostate-XNAT processing.
%
%
% From Philips PPE
% From 1 sequence on simulator. Ingenia
% Seq            ACQ       TE/TR  Delta/delta  WFS
% b90_vx1.3   1.31/1.31/5 54/2024 26.9/3.5   57.5 
% b500_vx1.3  1.31/1.31/5 68/2220 33.9/10.5  57.5
% b1500_vx1.3 1.31/1.31/5 94/2976 46.9/23.5  57.5
% b2000_vx1.3 1.31/1.31/5 75/6699 37.4/14.0  57.5
% b3000_vx1.3 1.31/1.31/5 87/6291 43.4/20.0  57.4
%
% 
%
% See also stejskal verdict_fit

if bv == 0
    % used when a 'fake' b=0 is added to the scheme
    scheme.delta = 1 ;
    scheme.DELTA = 2 ;
    scheme.G = 0 ;
    scheme.bval = 0 ;

    disp('b-value 0. returned G and bval 0')
    return
end

switch scanner
    case 'PDS' % Philips XX file
        scheme.delta = ddeltas(1) ;
        scheme.DELTA = ddeltas(2) ;
    case 'Ingenia'
        switch bv
            case 90
                scheme.delta = 3.5 ;
                scheme.DELTA = 26.9 ;
            case 500
                scheme.delta = 10.5 ;
                scheme.DELTA = 33.9 ;
            case 1500
                scheme.delta = 23.5 ;
                scheme.DELTA = 46.9 ;
            case 2000
                scheme.delta = 14.0 ;
                scheme.DELTA = 37.4 ;
            case 3000
                scheme.delta = 20.0 ;
                scheme.DELTA = 43.4 ;
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end

    case 'Ingeniav2.0'
        switch bv
            case {90, 1000}
                scheme.delta = 11.5 ;
                scheme.DELTA = 28.4 ;
            case 500
                scheme.delta = 16.8 ;
                scheme.DELTA = 33.9 ;
            case 1500
                scheme.delta = 24.5 ;
                scheme.DELTA = 41.4 ;
            case 2000
                scheme.delta = 20.0 ;
                scheme.DELTA = 37.9 ;
            case 3000
                scheme.delta = 28.0 ;
                scheme.DELTA = 44.9 ;
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end


    case 'IngeniaRenal'
        switch bv 
            case 70
                scheme.delta = 4.8 ;
                scheme.DELTA = 27 ;
            case 90
                scheme.delta = 4.8 ;
                scheme.DELTA = 27 ;
            case 150
                scheme.delta = 4.8 ;
                scheme.DELTA = 27 ;
            case 500
                scheme.delta = 12 ;
                scheme.DELTA = 34 ;
            case 1000
                scheme.delta = 12 ;
                scheme.DELTA = 34 ;
            case 1500
                scheme.delta = 26.3 ;
                scheme.DELTA = 47 ;
            case 2000
                scheme.delta = 16.8 ;
                scheme.DELTA = 37.5 ;
            case 2200
                scheme.delta = 16.8 ;
                scheme.DELTA = 37.5 ;
            case 2500
                scheme.delta = 21.4 ;
                scheme.DELTA = 43.5 ;
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end
    case 'Achieva'
        % values taken from Bonet-Carne paper, amended with small change to
        % b2000
        switch bv
            case 90
                scheme.delta = 3.9 ;
                scheme.DELTA = 23.8 ;
            case 500
                scheme.delta = 11.4 ;
                scheme.DELTA = 31.3 ;
            case 1500
                scheme.delta = 23.9 ;
                scheme.DELTA = 43.8 ;
            case 2000
                 scheme.delta = 12.4 ; % 14.4 ;
                 scheme.DELTA = 32.3 ; % 34.3 ;
            case 3000
                scheme.delta = 18.9 ;
                scheme.DELTA = 38.8 ;
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end
    case 'SIGNA Premier'
        % B-values in GradAmps below correspond to Stejskal Tanner calc
        switch bv
            case 90
                scheme.delta = 4.716 ;
                scheme.DELTA = 23.808 ; % GradAmps 50.465, rmp 0.892
            case 500
                scheme.delta = 12.212 ;
                scheme.DELTA = 31.296 ; % GradAmps 41.496, ramp 0.892
            case 1500
                scheme.delta = 25.776 ;
                scheme.DELTA = 43.396 ; % GradAmps 30.113, ramp 0.892
            case 2000
                scheme.delta = 14.944 ;
                scheme.DELTA = 32.300 ; % GradAmps 67.701, ramp 0.892
            case 3000
                scheme.delta = 18.932 ;
                scheme.DELTA = 38.808 ; % GradAmps 60.06, ramp 0.892
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end
    case 'SIGNA PET/MR'
        switch bv
            case 90
                scheme.delta = 5.980 ;
                scheme.DELTA = 23.800 ; % GradAmps 39.755	
            case 500
                scheme.delta = 12.604 ;
                scheme.DELTA = 31.296 ; % GradAmps 40.150
            case 1500
                scheme.delta = 25.772;
                scheme.DELTA = 43.408 ; % GradAmps 30.111	
            case 2000
                scheme.delta = 25.832 ;
                scheme.DELTA = 34.432 ; % GradAmps 40.281
            case 3000
                scheme.delta = 27.164 ;
                scheme.DELTA = 43.828 ; % GradAmps 40.427
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end
    case 'XNAT' % The scheme on prostate-XNAT up to at least January 2023
        switch bv
            case 90
                scheme.delta = 3.9  ;
                scheme.DELTA = 23.8 ; % Grad 61.2
                scheme.G = 61.2 ;
            case 500
                scheme.delta = 11.4 ;
                scheme.DELTA = 31.3 ; % Grad 44.3
                scheme.G = 44.3 ;
            case 1500
                scheme.delta = 23.9 ;
                scheme.DELTA = 43.8 ; % Grad 32
                scheme.G = 32 ;
            case 2000
                scheme.delta = 14.4 ;
                scheme.DELTA = 34.3 ; % Grad 67.7
                scheme.G = 67.7 ;
            case 3000
                scheme.delta = 18.9 ;
                scheme.DELTA = 38.8 ; % Grad 60.0
                scheme.G = 60.0 ;
            otherwise
                warning('MATLAB:bv2scheme:unknownBValue',...
                    ['Unknown b-value for scanner: ', scanner])
        end
    otherwise
        error('MATLAB:bv2scheme:unknownManufacturerModel', ...
            ['ManufactuerModel (scanner) ',scanner,' not recognised'])
end

switch scanner
    case 'XNAT'
        scheme.bval = stejskal(scheme.delta, scheme.DELTA,G=scheme.G) ;
    otherwise
        scheme.bval = bv ;
        scheme.G = stejskal(scheme.delta, scheme.DELTA, bval=scheme.bval) ;
end

disp([scanner, ' bv ',num2str(bv),' Scheme: bval ',num2str(scheme.bval),' G ',num2str(scheme.G), ...
    ' delta ',num2str(scheme.delta),' DELTA ',num2str(scheme.DELTA)])

end