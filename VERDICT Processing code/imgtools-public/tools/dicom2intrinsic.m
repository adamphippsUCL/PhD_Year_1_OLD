function varargout = dicom2intrinsic(geom, varargin)
% DICOM2INTRINSIC Links DICOM and intrinsic coordinate spaces
% Converts DICOM parameters (IOP, IPP etc) to the affine matrix that 
% will convert an lph-coordinate to a MATLAB intrinsic spatial coordinate.
% Used directly, coordinates need to be homogeneous (1 at end) but function
% can perform transformations.
%
% A = dicom2intrinsic(geom)  
% A = dicom2intrinsic(geom, 'output', 'A')  
%
% S = dicom2intrinsic(geom, 'output', 'MATLAB_nifti_sform') 
%
% [L,P,H] = dicom2intrinsoc(geom, 'output', 'LPHcoords')
%
% [coords_lph] = dicom2intrinsic(geom, 'output', 'stackedLPH', ...
%                                                   'coords', coords)
%
% [coords_intrinsic] = dicom2intrinsic(geom, 'output', 'stackedIntrinsic', ...
%                                      'coords', coords_lph)
%
% geom is a structure and in most of the code its length is the number of
% slices and there is an IPP for every slice. 
%  
% To reduce the chance of mistakes, geom must be length 1 when coords are
% passed in. In this case, the intrinsic coords are with respect to the 
% slice representd by geom. 
% There are two likely scenarios;
%  a) geom refers to the slice under consideration and the third dimension
%     of the intrinsic coordinate will be 1.
%  b) geom refers to the first slice of a volume and the third dimension
%     of intrinsic coordinate is with respect to that first slice.
% Usage along the style of a) is preferred as this might be more robust
% against left-handed systems or non regular volumes.
%
% Example usage to get the affine matrix:
%   lph_coordH = [L P H 1] ;  % row vector
%   A = dicom2intrinsic(geom) ;
%   intrinsic_coordH = lph_coordH * A 
%
%   Snifti = dicom2intrinsic(geom, 'output', 'MATLAB_nifti_sform') ;
%    % Snifti is in form for direct allocation to sform.
%
% To get the L,P,H coordinates of every voxel.
%   [L,P,H] = dicom2intrinsoc(geom, 'output', 'LPHcoords')
%
% To convert intrinsic coordinates (stacked, not homogeneous) to LPH. 
%    [coords_lph] = dicom2intrinsic(geom, 'output', 'stackedLPH', ...
%                                        'coords', coords)
%
%  coords here can be [N x 3] or [N x 2] in which case 3rd dimension 
%   will be set to one (intrinsic z coordinate, not homogeneous 1).
% 
%
% To convert stacked LPH coords [N x 3] to intrinsic
%    [coords_intrinsic] = dicom2intrinsic(geom, 'output', 'stackedIntrinsic', ...
%                                      'coords', coords_lph)
%
%
%
% Copyright David Atkinson, 2020, University College London
% D.Atkinson@ucl.ac.uk
%
% See also VRESAMPLE

output = 'A';
for ipv = 1:2:length(varargin)
    switch varargin{ipv} 
        case 'output'
            output = varargin{ipv+1} ;
        case 'coords'
            coords = varargin{ipv+1} ;
        otherwise
            error('Unknown parameter name')
    end
end

% Builds to affine matrix A
[slsep, rdc, cdc, sdc] = geom_check(geom) ;
psh = geom(1).PixelSpacing_HW(1) ;
psw = geom(1).PixelSpacing_HW(2) ;

origin_pp = geom(1).IPP ;
origin_pp = origin_pp(:)' ;  % row vector

% set intrinsic first (top left centre) coordinate
switch output
    case 'MATLAB_nifti_sform'
        itlc = [ 0 0 0 ] ;
    otherwise
        itlc = [ 1 1 1] ;
end

Torig = [1 0 0 0 ; 
         0 1 0 0 ; 
         0 0 1 0 ; 
      -origin_pp 1] ;
  
Tintrorig = [1 0 0 0 ; 
             0 1 0 0 ; 
             0 0 1 0 ; 
             itlc  1 ]; % translation of intrinsic to origin 
    
Trim = [rdc(1) cdc(1) sdc(1) 0 ;
        rdc(2) cdc(2) sdc(2) 0 ;
        rdc(3) cdc(3) sdc(3) 0 ;
         0      0       0    1 ];
      
        
Tscale = [1/psw   0      0     0 ; 
          0     1/psh    0     0 ; 
          0       0   1/slsep  0 ; 
          0       0      0     1] ;
      
Tras2lph = diag([-1 -1 1 1]) ;
    
A = Torig *  Trim * Tscale * Tintrorig ;   

switch output
    case 'A'
        varargout{1} = A ;
    case 'MATLAB_nifti_sform'
        % first row-matrix form for RAS to intrinsic, then invert
        S = Tras2lph * A ;  % RAS to intrinsic
        S = inv(S) ;  % intrinsic to RAS
        % DO NOT TRANSPOSE as output is for affine3d in row-matrix form
        
        varargout{1} = S ;
    case 'LPHcoords'
        % (adapted from vresample)
        tformA = affine3d(A) ;
        nslice = length(geom) ;
        
        % The intrinsic coordinates
       [I1, I2, I3] = meshgrid(1:geom(1).Width, 1:geom(1).Height, 1:nslice) ;
       
       % Calculate the DICOM LPH coordinates 
       [L, P, H] = transformPointsInverse(tformA, I1, I2, I3) ;
        
       varargout{1} = L ; varargout{2} = P ; varargout{3} = H ; 
    case {'stackedLPH','stackedlph'}
        if length(geom) >1, warning('geom should refer to a single slice'),end
        tformA = affine3d(A) ;
        if size(coords,2) == 2 % Add z=0
            coords = cat(2, coords, ones([size(coords,1) 1])) ; 
        end
        [coords_lph] = transformPointsInverse(tformA, coords) ;
        varargout{1} = coords_lph ;
    case {'stackedIntrinsic'}
        if length(geom) >1, warning('geom should refer to a single slice'),end
        tformA = affine3d(A) ;
        [coords_intrinsic] = transformPointsForward(tformA, coords) ;
        varargout{1} = coords_intrinsic ;
    otherwise
        error(['Unknown output type: ',output])
end


% Check transform at some corner locations

thresh = 0.001 ; % threshold for vector difference

% top left centre
actv = [origin_pp 1] * A ;
if norm(actv(1:3) - itlc) > thresh
    actv
    error('TLC incorrect')
end

width = geom(1).Width ;
height = geom(1).Height ;

% bottom-right corner, this slice
brc = origin_pp + (rdc(:)' * psw * (width-1)) + (cdc(:)' * psh * (height-1)) ;
int_brc = [brc 1] * A ;

if norm(int_brc(1:3) - ([ width-1 height-1 0] + itlc) ) > thresh
    int_brc
    error('bottom right corner not correct')
end

% check corner opposite tlc, last slice
nslice = length(geom) ;

fc = brc + sdc(:)' * (nslice-1) * slsep ;

int_fc = [fc 1] * A ;

if norm(int_fc(1:3) - ([width-1 height-1 nslice-1]+ itlc) ) > thresh
    int_fc
    error(['far corner not correct'])
end




