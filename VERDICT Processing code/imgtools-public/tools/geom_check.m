function [slsep, rdc, cdc, sdc] = geom_check(geom, varargin)
% GEOM_CHECK Checks slices are parallel, calculates slice centre separation
% or thickness if one slice. Checks slice direction
%
%  [slsep, rdc, cdc, sdc] = geom_check(geom, Name, value, ...)
%
%  geom is a structure with fields IOP, IPP, SliceThickness and possibly
%                                  MRAqType
%  Name value pairs:
%   'assertRH', {false} | true
%   'opstr',    {false} | true  outputs text
%                          'Centre' refers to image DC point of  mid slice
%
% slsep - slice centre separation (may differ from Slice thickness).
%
% rdc, cdc, sdc - row, column and slice direction cosines
%
% Copyright, 2019, David Atkinson 
% D.Atkinson@ucl.ac.uk
%
% See also DRESLICE  MPR_GEOM

THRESH = 0.001 ;
assertRH = false ; % error if not Right Handed
opstr = false ;

for ipv = 1:2:length(varargin)
    val = varargin{ipv+1} ;
    
    switch varargin{ipv}
        case 'assertRH'
            assertRH = val ;
        case 'opstr'
            opstr = val ;
        otherwise
            error(['Unknown param: ',varargin{ipv}])
    end
end
nsl = length(geom) ;

iop = geom(1).IOP ;
rdc = iop(1:3) ; % row direction cosine
cdc = iop(4:6) ;
sdc = cross(rdc,cdc) ; % slice direction cosine

slicepos = zeros([nsl 1]) ;

for isl = 1:nsl
    % Check IOPs are parallel
    
    iop_this = geom(isl).IOP ;
    if ~isequal(iop, iop_this)
      warning(' Slices not parallel ')
    end
   
    slicepos(isl) = dot(sdc, geom(isl).IPP) ;
end

[srtslicepos, idx] = sort(slicepos) ;
if nsl == 1
    slsep = geom(1).SliceThickness ;
else
    slsep = srtslicepos(2)-srtslicepos(1) ;
    
    % slices should be along positive sdc
    sv = geom(2).IPP - geom(1).IPP ;
    % check parallel
    if norm(sv./norm(sv) - sdc) > 0.01
        warning(['Slices not aligned with +ve slice direction'])
        disp(['sdc: ',num2str(sdc(:)'),' sv: ',num2str(sv(:)')])
        if assertRH
            error('Not Right Handed')
        end
    end
    if abs(norm(sv)-geom(1).SliceThickness) > 0.01
        warning(['Slice separation is not slice thickness in geom'])
        disp(['norm(sv): ',num2str(norm(sv)),' SliceTickness: ',...
            num2str(geom(1).SliceThickness)])
    end
    
end

for isl = 3:nsl
    sdiff = srtslicepos(isl)-srtslicepos(isl-1) ;
    if abs(sdiff-slsep) > THRESH
        warning([' Slice separations not equal, ',...
            num2str(sdiff),'  slsep ',num2str(slsep)] )
    end
end

if opstr
    midsl = ceil((nsl+1)/2) ;
    g = geom(midsl) ;
    
    oinfo = ori_info(g.IOP) ;
    
    % image centre point "DC"
    cp = g.IPP + ...
        (sz2DC(g.Height)-1) * g.PixelSpacing_HW(1) * g.IOP(4:6) + ...
        (sz2DC(g.Width)-1)  * g.PixelSpacing_HW(2) * g.IOP(1:3) ;
    
    %     height width      FOVmm        DC(L)  IPP(L)  IOP1L  IOP2L
    hdr  = '   Matrix        FOV mm               Centre     IPP     IOPr    IOPc';
    fmt1 = '[%4d x %4d]  %7.2f %7.2f mm  %8.2f   %8.2f  %6.3f  %6.3f' ;
    fmt2 = '                 |   ->%3s         %8.2f   %8.2f  %6.3f  %6.3f' ;
    fmt3 = '               %3s                 %8.2f   %8.2f  %6.3f  %6.3f' ;
    
    
    str1 = sprintf(fmt1, g.Height, g.Width, g.Height*g.PixelSpacing_HW(1), ...
           g.Width*g.PixelSpacing_HW(2), cp(1), g.IPP(1), g.IOP(1), g.IOP(4) ) ;
       
    str2 = sprintf(fmt2, oinfo.east_str, cp(2), g.IPP(2), g.IOP(2), g.IOP(5) ) ;
    str3 = sprintf(fmt3, oinfo.south_str, cp(3), g.IPP(3), g.IOP(3), g.IOP(6) ) ;
    
    disp(' ')
    disp(hdr)
    disp(str1)
    disp(str2)
    disp(str3)
end


end