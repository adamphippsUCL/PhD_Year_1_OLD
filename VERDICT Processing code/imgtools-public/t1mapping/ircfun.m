function cf = ircfun(x, times, opt, Mt)
% IRFN Inverstion Recovery Cost Function 
%
%  cf = ircfun(x,times, opt, Mt)
%
% opt can be 
%  'lsqnl_mag2'  lsqnonlin magnitude, 2 parameters 
%                   M0(1-2*exp(-TI/T1))
%   x = [ T1 M0 ]
%
%  'lsqnl_mag3' Barral eqn 15   r_a + r_b.exp(-TI/T1)
%   x = [ T1  r_a  r_b/(-2r_a)] 
%
% For magnitude data 'itype' is 6. (real is 8)
%
% See Barral  'A Robust Methodology for In Vivo T1 Mapping' 
%     MRM 64:1057–1067 (2010)
%
%
% Copyright 2020-2021. David Atkinson, University College London
% D.Atkinson@ucl.ac.uk
%
% See also LLCFN IRT1 IRT1_CHECK

if nargin < 4
    Mt = zeros([1 length(times)]) ;
end

switch opt
    case 'lsqnl_mag2'  % magnitude 2 parameters
        cf = abs(x(2) * (1 - 2.*exp(-times(:)./x(1)))) - Mt(:)  ; 
    case { 'lsqnl_mag3', 'RD-NLS-PR3' }
        cf = abs(x(2) *(1 - 2.* x(3)*exp(-times(:)./x(1)))) - Mt(:) ;
    otherwise
        error(['Unknown opt: ',opt])
end

% Old code for itype 8 (real), but does not correctly handle pahse
   % cf = (x(2) * (1 - 2.*exp(-times./x(1))))  ; 
