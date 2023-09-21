% MATLAB function to find b value which maxmises radii sensitivity

% Define diffusion coefficient
d = 2; % 10^-3 mm/s


% Define range of radii
Rmin = 4;
Rmax = 18;
dR = 0.1;
Rs = Rmin:dR:Rmax;

% Define range of gradient timings
delta = 10;  % ms
Dmin = delta;
dD = 1; % ms
Dmax = 100; % ms
Deltas = Dmin:dD:Dmax; % ms


% Range of b values to consider
bs = 500:50:1500;

% To be filled with radii sensitivities
Sensitivities = zeros(length(bs),1);

bindx = 0;
% Iterate over b values
for b = bs
    bindx = bindx + 1;

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


    % Find maximum dSdR at R=10
    R = 10;
    Rindx = ceil( (R-Rmin)/dR );
    [dSdR_max, maxindx] = max( dSdR(Rindx,:)) ;

    % Append to radii sensitivities
    ADC = 0.001;
    Sensitivities(bindx) = dSdR_max*exp(-b*ADC);


end

figure;
plot(bs, Sensitivities)
xlabel('b value')
ylabel('Radii sensitivity')
grid on
title(['Radii sensitivity, R = ' num2str(R) ', delta = ' num2str(delta) ', Delta = ' num2str(Deltas(maxindx)) ', ADC = ' num2str(ADC)] )

