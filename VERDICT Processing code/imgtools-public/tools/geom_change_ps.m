function [geom_out] = geom_change_ps(geom_q, PS_HW_new, varargin)
% GEOM_CHANGE_PS Change geom values to alter pixel size but keep FOV
%
%  geom_out = geom_change_ps(geom_q, PS_HW_new, Name, Value, ...) 
% 
%  geom_q is geom to be changed, PS_HW_new is the new (approx) PixelSpacing
%  geom_out contains the updated  parameters.
%  Note the actual PixelSpacing is chosen to maintain the original FOV.
%
% Use case is wanting high res data (e.g. T2) in plane of low res (e.g.
% DWI). Without modification, dreslice would return low-res slice of T2 in
% plane of DWI.
%
% User supplies high-res pixel spacing. Algorithm will maintain FOV and
% adjust actual pixel spacing to fit integer number of pixels (so that
% original DWI and be viewed side-by-side with T2 with same FOV).
% 
% Copyright, 2020, David Atkinson 
% D.Atkinson@ucl.ac.uk
%
% See also set_interpolant dgeomextract geom2sro plot_fov

test = false ;

geom_out = geom_q ;

PS_HW_q = geom_q(1).PixelSpacing_HW ;
Height = double(geom_q(1).Height) ;
Width  = double(geom_q(2).Width) ;

% Determine FOV (which will be held constant)
FOV_HW = [ PS_HW_q(1)*Height PS_HW_q(2)*Width ] ;

nHeight = round(FOV_HW(1)/PS_HW_new(1)) ;
nWidth  = round(FOV_HW(2)/PS_HW_new(2)) ;

PS_HW_n(1) = FOV_HW(1) / nHeight ;
PS_HW_n(2) = FOV_HW(2) / nWidth ;

geom_out(1).PixelSpacing_HW = PS_HW_n ;
geom_out(1).Width = nWidth ;
geom_out(1).Height = nHeight ;

% See dgeomextract for XData or geom2sro for spatial referencing
% New practice is NOT to keep updating these so strip out
if isfield(geom_out,'XData')
    geom_out = rmfield(geom_out,{'XData' ; 'YData' ; 'R2D'} ) ;
end

for islice = 1:length(geom_q)
    geom_out(islice).PixelSpacing_HW = PS_HW_n ;
    geom_out(islice).Width = geom_out(1).Width ;
    geom_out(islice).Height = geom_out(1).Height ;
    
    % Update IPP.
    % compute top left corner vertex (not pixel centre) for these calculations
    ipp = geom_out(islice).IPP;
    iop = geom_out(islice).IOP ;
    tlv = ipp - 0.5*(PS_HW_q(2)*iop(1:3) + PS_HW_q(1)*iop(4:6)) ;
    
    IPPnew = tlv + 0.5*(PS_HW_n(2)*iop(1:3) + PS_HW_n(1)*iop(4:6)) ;
    
    geom_out(islice).IPP = IPPnew ;
end

if test
    figure('Name','Test geom_change_ps')
    plot3(geom_q(1).IPP(1), geom_q(1).IPP(2), geom_q(1).IPP(3),'r+')
    hold on, grid on, xlabel('L'), ylabel('P'), zlabel('H'), axis square
    plot3(geom_out(1).IPP(1), geom_out(1).IPP(2), geom_out(1).IPP(3),'k+')
    
    plot_fov(geom_q,'r-')
    plot_fov(geom_out,'k--')
end

end

function plot_fov(geom, linespec)
% PLOT_FOV Plots outline of FOV from geom structure
%
% plot_fov(geom, linespec)
%
% See also geom_change_ps

% top left vertex (edge of corner of image)
tlv = geom(1).IPP -0.5* ( geom(1).PixelSpacing_HW(1)*geom(1).IOP(4:6) + ...
                          geom(1).PixelSpacing_HW(2)*geom(1).IOP(1:3) ) ;

% top right, bottom right and bottom left
trv = tlv + geom(1).PixelSpacing_HW(2)*geom(1).IOP(1:3)*(geom(1).Width+1) ;
brv = trv + geom(1).PixelSpacing_HW(1)*geom(1).IOP(4:6)*(geom(1).Height+1) ;
blv = tlv + geom(1).PixelSpacing_HW(1)*geom(1).IOP(4:6)*(geom(1).Height+1) ;

plot3([tlv(1) trv(1) brv(1) blv(1) tlv(1)], ...
      [tlv(2) trv(2) brv(2) blv(2) tlv(2)], ...
      [tlv(3) trv(3) brv(3) blv(3) tlv(3)], linespec)
end


                      
    
    
    