function explore_VERDICT_signal
% Signals from sphere, ball and astrosticks
%
%
% See also sphereGPD  astrosticks stejskal ball

delta = 23.9 ; % ms
DELTA = 43.8 ; % ms
G = 32 ;     % mT/m 
dIC = 2 ;    % 1e-9 m^2/s diffusivity for sphere (intra-cellular)
dVASC = 8 ;  % pseduo-diffusion
dEES = 2 ;

bvals = [0:1500] ; % for astrosticks in [s/mm^2]? e.g. 0, 150, 500, 1000


Rs = [1:30] ; % radii in um

sig = zeros([1 length(Rs)]) ;
for ir = 1:length(Rs)
    r_this = Rs(ir) ;
    sig(ir) = sphereGPD(delta, DELTA, G, r_this, dIC) ;
end

figure
plot(Rs,sig)
axis([0 max(Rs) 0 1]);
xlabel('radius')

% AstroSticks and Ball
sig_as = zeros([1 length(bvals)]) ;
sig_ball = zeros([1 length(bvals)]) ;
for ibval = 1:length(bvals)
   sig_as(ibval) = astrosticks(bvals(ibval), dVASC) ;
   sig_ball(ibval) = ball(bvals(ibval), dEES) ; 
end

figure
plot(bvals, sig_as,'DisplayName','astrostick')
axis([0 max(bvals) 0 1]);
xlabel('b-value')

hold on
plot(bvals,sig_ball,'DisplayName','ball')
legend

% Attempt to replicate Fig 1 in taxonony paper (for 9.4T system)
% Note the curves here show similar patterns, but the absolute values
% differ. Not clear why - possibly the diffusivity is not just the parallel
% diffusivity in the paper??

Gs = [10 : 10: 400];
d=0.8 ;
r=4;

sig_Grstar = zeros([1 length(Gs)]); % red star sphere, DELTA 20, delta 3, 
sig_Gbsq = zeros([1 length(Gs)]) ; % light blue square 30,30
sig_Gysq = zeros([1 length(Gs)]) ; % astro yellow square 30,30

for iG = 1:length(Gs)
    G_this = Gs(iG) ;
    sig_Grstar(iG) = sphereGPD(3,20,G_this,r,d) ;
    sig_Gbsq(iG)   = sphereGPD(30,30,G_this,r,d) ;

    bval = stejskal(3, 30, G_this) ;
    sig_Gysq(iG) = astrosticks(bval,d);
end

figure
plot(Gs,sig_Grstar,'r*')
hold on
plot(Gs,sig_Gbsq,'bs')
plot(Gs,sig_Gysq,'s')
axis([0 max(Gs) 0 1])
xlabel('G')


