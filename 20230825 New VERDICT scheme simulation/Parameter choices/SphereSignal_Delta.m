% MATLAB script to asssess sphere signal dependence on delta and Delta


% Define diffusion coefficient
d = 2; % 10^-3 mm/s

% Define b Value
b = 3000; % s/mm^2


%% Varying Delta

% Define range of gradient timings
delta = 1;  % ms
Dmin = delta;
dD = 1; % ms
Dmax = 60; % ms
Deltas = Dmin:dD:Dmax; % ms


% Define radii
Rs = 5:1:17;

fig = figure;

for R = Rs

    % Find signal as function of Delta
    signals = zeros(length(Deltas), 1);
    
    Dindx = 0;
    for Delta = Deltas
        Dindx = Dindx + 1;
        G = stejskal(delta, Delta, 'bval', b);
        signalSphere = sphereGPD(delta, Delta, G, R, d) ;
        signals(Dindx) = signalSphere;
    
    end
    
    plot(Deltas, signals, DisplayName = ['R = ' num2str(R) ' um'])
    hold on
end

legend(Location = 'southeast');
xlabel('Delta')
ylabel('Signal')
ylim([0,1])
title(['Sphere signal dependence on Delta, delta = ' num2str(delta) ' us, b = ' num2str(b)])
saveas(fig, 'Figures/SphereSignal_BigDelta.png')