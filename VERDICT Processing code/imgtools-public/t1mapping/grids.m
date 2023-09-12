function [T1, ra, rb, tau] = grids(Mt, itimes)
% GRIDS Grid search for T1 fitting
% [T1, ra, rb, tau] = grids(Mt, itimes)
%
% Implements RD-NLS-PR3 from Barral et al, MRM 64:1057–1067 (2010)
% forfit to magnitude data.
%
% Copyright 2020-2021. David Atkinson, University College London
% D.Atkinson@ucl.ac.uk
%
% See also IRT1
%
T1s = [100:10:9000] ;

itimes = itimes(:)' ;
[mint, lmin] = min(abs(Mt)) ;
        
tau1 = lmin-1 ; if tau1== 0 ; tau1 = 1 ; end
tau2 = lmin  ;
cf1_max = realmin ;
cf2_max = realmin ;

for it1 =1:length(T1s)
    T1_curr = T1s(it1) ;
    
    [gamma1, psi1, ra1, rb1, Jtau1] = rd_params(T1_curr, tau1, Mt, itimes) ;
    cf1 = gamma1*gamma1/psi1 ;
    if cf1>cf1_max
        T1_est1 = T1_curr ;
        ra_est1 = ra1 ;
        rb_est1 = rb1 ;
        tau_est1 = tau1 ;
        Jtau_est1 = Jtau1 ;
        cf1_max = cf1 ;
    end
    
    [gamma2, psi2, ra2, rb2, Jtau2] = rd_params(T1_curr, tau2, Mt, itimes) ;
    cf2 = gamma2*gamma2/psi2 ;
    if cf2>cf2_max
        T1_est2 = T1_curr ;
        ra_est2 = ra2 ;
        rb_est2 = rb2 ;
        tau_est2 = tau2 ;
        Jtau_est2 = Jtau2 ;
        cf2_max = cf2 ;
    end
    
    
end

if Jtau_est1 < Jtau_est2
    T1 = T1_est1;
    ra = ra_est1;
    rb = rb_est1;
    tau = tau_est1 ;
else
    T1 = T1_est2;
    ra = ra_est2;
    rb = rb_est2;
    tau = tau_est2 ;
end


% Barral eqn 24
% gamma

% single T1_curr
function [gamma, psi, ra, rb, Jtau] = rd_params(T1_curr, tau, Mt, itimes)
PR = ones(size(Mt)) ;
PR(1:tau) = -1 ;
np = length(itimes) ;

sum_exp = sum( exp(-itimes/T1_curr)) ;

% eqn 
gamma = sum( exp(-itimes/T1_curr).* PR.*Mt ) - ...
                       1/np * sum_exp * (sum( PR.*Mt)) ;

% eqn 14
psi = sum( exp(-2*itimes/T1_curr)) - 1/np * sum_exp*sum_exp ;

% eqns 22, 23
rb = gamma / psi ;
    
ra = 1/np * ( sum(PR.*Mt ) - rb * sum_exp) ;

% Barral eqn 19
Jtau = sum( (PR .* Mt - (ra+rb*exp(-itimes/T1_curr))).^2 ) ;
