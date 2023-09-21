% MATLAB script to asssess sphere signal dependence on delta and Delta


% Define diffusion coefficient
d = 2; % 10^-3 mm/s

% Define b Value
b = 1000; % s/mm^2


%% Varying delta

% Define range of gradient timings
Delta = 25;  % ms
dmin = 1;
dd = 1; % ms
dmax = Delta; % ms
deltas = dmin:dd:dmax; % ms


% Define radii
Rs = 5:1:17;

fig = figure;

for R = Rs

    % Find signal as function of Delta
    signals = zeros(length(deltas), 1);
    
    dindx = 0;
    for delta = deltas
        dindx = dindx + 1;
        G = stejskal(delta, Delta, 'bval', b);
        signalSphere = sphereGPD(delta, Delta, G, R, d) ;
        signals(dindx) = signalSphere;
    
    end
    
    plot(deltas, signals, DisplayName = ['R = ' num2str(R) ' um'])
    hold on
end

legend(Location = 'southeast');
xlabel('delta')
ylabel('Signal')
ylim([0,1])
title(['Sphere signal dependence on delta, Delta = ' num2str(Delta) ' us, b = ' num2str(b)])
saveas(fig, 'Figures/SphereSignal_smalldelta.png')