% Matlab script to find Delta which maximises dS/dR at R = 10
% 
% === Code structure
% 
% 1. For a range of R, find signal as function of Delta
% 2. Evaluate dS/dR at R = 10 for each Delta.
% 3. Plot dS/dR as function of Delta to choose range of Delta's
%
% ----------------------------------------------------------------


% === 1. find signal as function of Delta for a range of R

% Define diffusion coefficient
d = 2; % 10^-3 mm/s

% Define b Value
b = 500; % s/mm^2

% Define range of gradient timings
delta = 15;  % ms
Dmin = delta;
dD = 5;
Dmax = 80;
Deltas = Dmin:dD:Dmax;  % ms

% Define range of radii
Rmin = 5;
Rmax = 18;
dR = 0.1;
Rs = Rmin:dR:Rmax;

% Array to be filled with signals
signals = zeros( [length(Rs), length(Deltas)] );

% Simulate signals
Rindx = 0;

for R = Rs

    Rindx = Rindx + 1;

    Dindx = 0;
   
    for Delta = Deltas
    
        Dindx = Dindx  + 1;
    
        % Calculate gradient strength
        G = stejskal(delta, Delta, 'bval', b);
    
        signalSphere = sphereGPD(delta, Delta, G, R, d) ;

        signals(Rindx, Dindx) = signalSphere;
      
    end

end


% Calculate dS/dR

[dSdD, dSdR] = gradient(signals, dD, -dR);

% Define radius
R = 10;
Rindx = round( (R-Rmin)/dR );
plot(Deltas, dSdR(Rindx,:))
hold on







