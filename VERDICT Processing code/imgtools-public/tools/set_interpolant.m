function interpolant = set_interpolant(geomV, geomQ, varargin)
% SET_INTERPOLANT Set separable interpolant for data reslicing / resampling
%
%  interpolant = set_interpolant(geomV, geomQ, Name, Value, ...)
%
% geomV describes the fixed volume image
% geomQ describes the locations queried, i.e. the geometry of the new image
%
% Name value pairs can be:
%  'plot_kernel' {true} | false
%  
% To change the pixel spacing of the output, use geom_change_ps first
%
% set_interpolant tries to choose a sensible kernel given the input data.
% It avoids aliasing by widening the kernel.
% The through plane kernel depends on whether the data is MR 3D or not.
%
% There may better kernel choices to preserve features. See articles by
% Thevanez & Unser for more details. 
%
% The MATLAB documentation for tformarray is sparse and testing is
% difficult!
%
% Copyright, 2020, David Atkinson 
% D.Atkinson@ucl.ac.uk
%
% See also geom_change_ps geom_check  tformarray makeresampler



% makeresampler takes as inputs {half_width, positive_half}
%  where positive_half is a vector of values
%
% This relates to the intrinsic space of the volume that is being
% interpolated into.

plot_kernel = true ;
if ~isfield(geomV, 'MRAqType'), geomV(1).MRAqType = ''; end
if ~isfield(geomQ, 'MRAqType'), geomQ(1).MRAqType = ''; end


for ipv = 1:2:length(varargin)
    switch varargin{ipv}
        case 'plot_kernel'
            plot_kernel = varargin{ipv+1} ;
        otherwise
            error('Unknown parameter')
    end
end


[slsepV, rdcV, cdcV, sdcV] = geom_check(geomV) ;
[slsepQ, rdcQ, cdcQ, sdcQ] = geom_check(geomQ) ;

ps_hwQ = geomQ(1).PixelSpacing_HW ;
ps_hwV = geomV(1).PixelSpacing_HW ;

% find directions in Q that give largest intervals along dimensions of V
[rdm, irdm] = mvec(rdcV, rdcQ, cdcQ, sdcQ, ps_hwQ, slsepQ) ;
[cdm, icdm] = mvec(cdcV, rdcQ, cdcQ, sdcQ, ps_hwQ, slsepQ) ;
[sdm, isdm] = mvec(sdcV, rdcQ, cdcQ, sdcQ, ps_hwQ, slsepQ) ;

% rcs.  Note rdm means ALONG a row
scale = [cdm/ps_hwV(1)  rdm/ps_hwV(2)  sdm/slsepV] ;
rscale = max([1 1 1], scale) ;

% cubic half kernel from peak to 2nd zero crossing
kv = cubic(linspace(0,2,200)) ; 


% Possible through slice scenarios:
%   1) True 3D input data: should use cubic or similar
%   2) 2D MS or M2D with no gaps: linear
%   3) Slices with gaps: linear
%   4) Single slice: nearest neighbour with width of slice thickness

% Through Slice
if length(geomV) == 1 % Single slice: nearest neighbour (box)
    khw = 1 / 2 ; % half_width is half of intrinsic 
    kph = [1 1] ; % box
else
    switch geomV(1).MRAqType
        case '3D' % cubic
            khw = 2*rscale(3) ;
            kph = kv/rscale(3) ;
        case '2D' % linear
            khw = rscale(3) ;     % This is slsep_in*rscale(3) in real units.
            kph = linspace(1,0,100)./rscale(3) ;   
        otherwise % Non MR, hence use linear
            khw = rscale(3) ;     % This is slsep_in*rscale(3) in real units.
            kph = linspace(1,0,100)./rscale(3) ;   
    end
end

% In a test where this was manually altered, the second entry appeared to
% affect the image vertically, suggesting the second entry should have the
% first element of rscale. See tformarray for "help".
interpolant = {    {2*rscale(2) kv/rscale(2)}, ...
                   {2*rscale(1) kv/rscale(1)} , ...
                   {khw         kph         }} ;

% For large angulations, e.g. TRA to SAG
switch geomQ(1).MRAqType
    case '2D'
        % Q is 2D and if the through-slice direction is dominating interpolation
        % in V (probably large angulation differences) - kernel
        % probably should be a slice profile but linear used here.
        if irdm == 3 
            interpolant{1} = {rscale(2) linspace(1,0,100)./rscale(2)};
        end
        if icdm == 3 
            interpolant{2} = {rscale(1) linspace(1,0,100)./rscale(1)};
        end
end

if plot_kernel
    figure('Name','kernels')
    y = interpolant{1}{2} ; x = linspace(0,interpolant{1}{1}, length(y)) ;
    plot(x,y,'LineWidth',2,'DisplayName','row'), hold on, grid on
    
    y = interpolant{2}{2} ; x = linspace(0,interpolant{2}{1}, length(y)) ;
    plot(x,y,'LineWidth',2,'DisplayName','col')
    
    y = interpolant{3}{2} ; x = linspace(0,interpolant{3}{1}, length(y)) ;
    plot(x,y,'LineWidth',2,'DisplayName','slice')
    
    legend
    
    disp(['The dim in Q that has largest sampling along row in V is: ',...
        num2str(irdm), ' with projection ',num2str(rdm),'mm.'])
    disp(['  largest sampling along col in V is: ',...
        num2str(icdm), ' with projection ',num2str(cdm),'mm.'])
    disp(['  largest sampling along slice in V is: ',...
        num2str(isdm), ' with projection ',num2str(sdm),'mm.'])
    
end

end
%---------------------------------
function [dm, idm] = mvec(vecV, rdcQ, cdcQ, sdcQ, ps_hwQ, slsepQ)
[dm, idm] = max( abs([ dot(vecV, rdcQ * ps_hwQ(2)) ...
                       dot(vecV, cdcQ * ps_hwQ(1)) ...
                       dot(vecV, sdcQ * slsepQ)     ])) ;
end

%---------------------------------
function f = cubic(x)
% From MATLAB imresize function:
% See Keys, "Cubic Convolution Interpolation for Digital Image
% Processing," IEEE Transactions on Acoustics, Speech, and Signal
% Processing, Vol. ASSP-29, No. 6, December 1981, p. 1155.

absx = abs(x);
absx2 = absx.^2;
absx3 = absx.^3;

f = (1.5*absx3 - 2.5*absx2 + 1) .* (absx <= 1) + ...
                (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) .* ...
                ((1 < absx) & (absx <= 2));
end