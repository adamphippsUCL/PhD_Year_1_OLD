function fitt1 = IRT1(Mt, itimes, opt)
% IRT1 Inversion Recovery T1 calculation.
% lsqnonlin fit to magnitude or complex data Mt at inversion times itimes
% fit = IRT1(Mt, itimes, opt)
%
% opt is a string with the optimisation method. 
%    ( Note MATLAB lsqnonlin may not be Levenberg-Marquadt alg.)
%   'lsqnl_mag2' lsqnonlin magnitude data, 2 unknowns: T1, M0
%   'lsqnl_mag3' lsqnonlin magnitude data, 3 unknowns: T1, M0, inv
%      Barral has: r_a + r_b.exp(-TI/T1), here we use
%           abs( x(2) *(1 - 2.* x(3)*exp(-times./x(1))) ) 
%      so r_b = -2*x(2)*x(3) 
%
%   'RD-NLS-PR3' Reduced dimension, NLS, polarity restored, 3 unknowns
%      grid search
%
% fitt1 is a structure with fields:
%   T1
%   M0
%   RB (for opt 'RD-NLS-PR3' or 'lsqnl_mag3')
%   tau (for polarity restored, signal from times 1:tau are reversed in sign)
%
% See Barral et al
% Magnetic Resonance in Medicine 64:1057–1067 (2010)
%
% 
% Example (MultiFrame)
%  dinfoIR = datparse(dselect) ;  % select  Inversion Recovery.
%  ITYPE = 6 ;  % 6 for magnitude, 8 for real (don't use real) 
%  [vir, mir] = d2mat(dinfoIR,{'InversionTime','itype'},'itype',ITYPE,'op','fp') ; 
%  opt = 'lsqnl_mag3' ; % or 'lsqnl_mag2' or 'RD-NLS-PR3'
%  fitt1 = IRT1(vir, mir.tiVec, opt) ;
%
% Example (SingleFrame)
%  dinfo = datparse ;
%  ITYPE = 6 ;  % 6 for magnitude, 8 for real (don't use real)
%  [vir, mir] = d2mat(dinfo,{'series','itype'},'series', [ ], 'itype',ITYPE,'op','fp') ; 
%  opt = 'lsqnl_mag3' ; % or 'lsqnl_mag2' or 'RD-NLS-PR3' 
%  fitt1 = IRT1(vir, mir.tiVec_indata, opt) ;
%
%  
%
% Copyright 2020-2022. David Atkinson  University College London
% D.Atkinson@ucl.ac.uk
%
% See also LLT1 T1DISPLAY D2MAT DATPARSE GRIDS 


thresh = 15 ; % threshold for T1 calc

itimes = itimes(:) ; % column vector
Mt = double(Mt) ;
nt = length(itimes) ;
nd = ndims(Mt) ;
if size(Mt,nd) ~= nt
    error(['Number of times and last dimension of Mt must agree'])
end

sz_Mt = size(Mt) ;

np = prod(sz_Mt(1:end-1)) ;
T1 = zeros([np 1]) ;
M0 = zeros([np 1]) ;
RB = zeros([np 1]) ;  % only needed for 3 param fit
tau = zeros([np 1]) ;  % only needed for polarity reversed fits

Mt = reshape(Mt,[np nt]) ;
Mts = sum(abs(Mt),2) ;
loc = find(Mts > max(Mts(:))/thresh) ;
disp([num2str(length(loc)),' of ',num2str(np),' pixels and ',num2str(nt), ...
    ' time points.'])

options = optimset('Display','off') ;
nloc = length(loc) ;


if exist('parpool','file')
    T1_par = zeros([nloc 1]) ;
    M0_par = zeros([nloc 1]) ;
    RB_par = zeros([nloc 1]) ;
    tau_par = zeros([nloc 1]) ;
    T1_est = zeros([nloc 1]) ;
    M0_est = zeros([nloc 1]) ;
    RB_est = zeros([nloc 1]) ;
    
    % Set starting estimates (needs to be outside parfor loop)
    for il = 1:nloc
        ip = loc(il) ;
        [mint, lmin] = min(abs(Mt(ip,:))) ;
        % Basic eqn is M(t) = M0( 1 - 2*exp(t/T1) )
        T1_est(il) = itimes(lmin)/log(2) ; % T1 estimated at zero crossing
        M0_est(il) = max(abs(Mt(ip,:))) ; % M0 was estimated from final time, here added 
        % estimate to be from max in case recovery times do not go
        % far enough.

        switch opt
            case 'lsqnl_mag2'
                % do nothing more
            case 'lsqnl_mag3'
                RB_est(il) = 1 ;
            case 'RD-NLS-PR3'
            otherwise
                error(['Unknown opt: ',opt])
        end
    end
    
    switch opt
        case 'lsqnl_mag2'
            parfor il = 1:nloc
                ip = loc(il) ;
                f = @(x)ircfun(x,itimes,opt, Mt(ip,:)) ;
                [x,resnorm,residual,exitflag] = lsqnonlin(f,[T1_est(il) M0_est(il)],[0 0],[],options) ;
                T1_par(il) = x(1) ;  M0_par(il) = x(2) ;
            end
            T1(loc) = T1_par ;
            M0(loc) = M0_par ;
        case 'lsqnl_mag3'
            parfor il = 1:nloc
                ip = loc(il) ;
                f = @(x)ircfun(x,itimes,opt, Mt(ip,:)) ;
                [x,resnorm,residual,exitflag] = lsqnonlin(f,...
                    [T1_est(il) M0_est(il) RB_est(il)],[0 M0_est(il)/2 0], [6000 2*M0_est(il) 1], options) ;
                T1_par(il) = x(1) ;  M0_par(il) = x(2) ; RB_par(il) = x(3) ;
            end
            T1(loc) = T1_par ;
            M0(loc) = M0_par ;
            RB(loc) = RB_par ;
        case 'RD-NLS-PR3' 
            parfor il = 1:nloc
                ip = loc(il) ;
                [T1g, rag, rbg, taug] = grids(Mt(ip,:), itimes) ;
                T1_par(il) = T1g ; M0_par(il) = rag ; RB_par(il) = rbg ;
                tau_par(il) = taug ;
            end
            T1(loc) = T1_par ;
            M0(loc) = M0_par ;
            RB(loc) = RB_par ; 
            tau(loc) = tau_par ;
    end
    
else
    warning(['Old code no longer checked'])
    hw = waitbar(0,'Calculating T1s') ;
    upint = floor(length(loc)/100) ;
    for il = 1:length(loc)
        if rem(il,upint)==0
            waitbar(il/length(loc),hw) ;
        end
        ip = loc(il) ;
        % starting estimate
        
        % Eqn is M(t) = M0( 1 - 2*exp(t/T1) )
        [mint, lmin] = min(abs(Mt(ip,:))) ;
        x0(1) = itimes(lmin)/log(2) ; % T1 estimated at zero crossing
        x0(2) = Mt(ip,nt) ; % M0 estimated from final time
        
        f = @(x)ircfun(x,itimes, 'lsqnl_mag2', Mt(ip,:)) ;
        
        [x,resnorm,residual,exitflag] = lsqnonlin(f,x0,[0 0],[],options) ;
        T1(ip) = x(1) ;
        M0(ip) = x(2) ;
    end
    close(hw) ;
end %parpool

T1 = reshape(T1,[sz_Mt(1:end-1)]) ;
M0 = reshape(M0,[sz_Mt(1:end-1)]) ;
RB = reshape(RB,[sz_Mt(1:end-1)]) ;
tau = reshape(tau,[sz_Mt(1:end-1)]) ;

fitt1.T1 = T1;
fitt1.M0 = M0 ;
fitt1.RB = RB ;
fitt1.tau = tau ;
fitt1.opt = opt ;


% ! Now replaced by separate function !!
% function cf = ircfun(x,times, Mt, itype)
% % IRCFN Inverstion Recovery cost function for lsqnonlin
% % See also LLCFN
% 
% if itype == 6
%    cf = abs(x(2) * (1 - 2.*exp(-times./x(1))))  - Mt(:) ; 
% elseif itype == 8
%    cf = (x(2) * (1 - 2.*exp(-times./x(1))))  - Mt(:) ; 
% else
%     error
% end
