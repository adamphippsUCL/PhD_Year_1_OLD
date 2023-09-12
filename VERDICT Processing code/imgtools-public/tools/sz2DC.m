function grid_DC = sz2DC(sz,dim)
% sz2DC Return the DC grid point based on a matrix size.
%
% grid_DC = sz2DC(sz)      - returns vector same length as sz
% grid_DC = sz2DC(sz,dim)  - return dimension dim only.
%
% Adapted from DCgp
% 
% Copyright 2020-2021. David Atkinson, University College London
%
% See also DC3D DCgp

numDims = length(sz) ;
if nargin > 1
    if dim > numDims
        grid_DC = 1 ;
        return
    end
end

grid_DC = zeros(1,numDims) ;

for idim = 1:numDims
  nx = sz(idim) ;
  if rem(nx,2) == 0
    grid_DC(idim) = nx/2 + 1 ;
  else
    grid_DC(idim) = (nx+1)/2 ;
  end
end

if nargin == 2
  grid_DC = grid_DC(dim) ;
end
