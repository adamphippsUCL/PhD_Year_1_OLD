% MATLAB script to analyse imapct of delta and Delta on sphere signal.
% 
% Want to compare magnitude of maximised dS/dR for delta and Delta


% Define diffusion coefficient
d = 2; % 10^-3 mm/s

% Define b Value
b = 1000; % s/mm^2

% Define range of radii
Rmin = 4;
Rmax = 18;
dR = 0.1;
Rs = Rmin:dR:Rmax;

%% Varying Delta

% Define range of gradient timings
delta = 10;  % ms
Dmin = delta;
dD = 1; % ms
Dmax = 100; % ms
Deltas = Dmin:dD:Dmax; % ms

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

% Plot dS/dR for a range of radii
fig = figure(1);
ax1 = axes;

% fig2 = figure(2);
% ax2 = axes;

for R = 6:2:15
    Rindx = ceil( (R-Rmin)/dR );
    plot(ax1, Deltas, dSdR(Rindx,:), DisplayName = ['R = ' num2str(R) ' um']);
    hold on

    % Maximum
    [M, indx] = max( dSdR(Rindx,:) );
%     scatter(ax1, Deltas(indx), M, 'c', 'black', 'HandleVisibility','off');
    xline(ax1, Deltas(indx), 'HandleVisibility','off');

end

legend(ax1)
xlabel(ax1, 'Delta')
ylabel(ax1, 'dSdR')
ylim(ax1, [0,0.11])
xlim(ax1, [0,Dmax])
grid on
title(ax1, ['dSdR, delta = ' num2str(delta) ', b = ' num2str(b)] )
saveas(fig, ['Figures/dSdR_delta' num2str(delta) '_b_' num2str(b) '.png'])

%% Dependence on b value

% For R = 10, find maximum dS/dR

R = 10;
Rindx = ceil( (R-Rmin)/dR );

dSdR_max = max( dSdR(Rindx,:)) ;

% %% Varying delta
% 
% 
% % Define range of gradient timings
% Delta = 40;  % ms
% dmin = 0;
% dd = 1; % ms
% dmax = Delta; % ms
% deltas = dmin:dd:dmax; % ms
% 
% 
% % Array to be filled with signals
% signals = zeros( [length(Rs), length(deltas)] );
% 
% % Simulate signals
% Rindx = 0;
% for R = Rs
%     Rindx = Rindx + 1;
%     dindx = 0;
%     for delta = deltas
%         dindx = dindx  + 1;
%         % Calculate gradient strength
%         G = stejskal(delta, Delta, 'bval', b);
%         signalSphere = sphereGPD(delta, Delta, G, R, d) ;
%         signals(Rindx, dindx) = signalSphere;
%     end
% end
% 
% % Calculate dS/dR
% [dSdd, dSdR] = gradient(signals, dd, -dR);
% 
% % Plot dS/dR for a range of radii
% figure;
% for R = 6:2:15
%     Rindx = ceil( (R-Rmin)/dR );
%     plot(deltas, dSdR(Rindx,:), DisplayName = ['R = ' num2str(R) ' um'])
%     hold on
% end
% legend
% xlabel('delta')
% ylabel('dSdR')
% 
