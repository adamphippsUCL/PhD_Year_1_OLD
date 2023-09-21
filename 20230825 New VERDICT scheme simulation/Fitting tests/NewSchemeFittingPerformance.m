% MATLAB script to test fitting performance of new VERDICT scheme

% Define New VERDICT scheme
V0 = [1,2,0];
V1 = [10, 20, 500];
V2 = [15, 25, 1200];
V3 = [20, 30, 1000];
V4 = [30, 45, 800];
% V5 = [25, 35, 2000];

Vs = [
    V0; V1; ...
    V0; V2; ...
    V0; V3;...
    V0; V4 ...
    ];%; V5];

% % Define old scheme
% V1 = [3.9, 23.8, 90];
% V2 = [11.4, 31.3, 500];
% V3 = [23.9, 43.8, 1500];
% V4 = [12.4, 32.3, 2000];
% V5 = [18.9, 38.8, 3000];
% 
% Vs = [V2; V3; V4; V5];


%% Define fitting grid

fICs = linspace(0,1,21);
Rs = linspace(7,13,7);
Nreps = 100;


% Test on specific fIC and R
fIC = 0.4;
fEES = 1-fIC;
fVASC = 0;
Ravg = 10;
Rstd = 0.001;

tissue_params = [fIC, fEES, fVASC, Ravg, Rstd];

% Array to be filled in with 
new_fIC_fits = zeros(Nreps,1);

for indx = 1:Nreps
    

    % Simulate signal
    [signals, Vscheme, params] = verdict_simulate_Adam(Vs, tissue_params);
   
    % Define SNRs
    bs = [Vscheme.bval];
    SNRs = 100*(exp(-bs*0.0008));

    % Add noise to signal
    signalsNoisy = abs(addnoise_Adam(signals, NoiseSigma = 0.02));

    signalsNoisy(1:2:end) = 1;

    scatter([Vscheme.bval], signalsNoisy)
    hold on

    % Fit to signal
    Y = zeros([1,1,length(signalsNoisy)]);
    Y(1,1,:) = signalsNoisy;
    
    [fIC_fit, fEES_fit, fVASC_fit, R_fit, rmse, A, t, opt] = verdict_fit(Vscheme, Y, ncompart = 1);

    % Append to fIC fits
    new_fIC_fits(indx) = fIC_fit;

end

ylim([0,1])

new_error = sqrt( mean((new_fIC_fits-fIC).^2) )
new_var = std(new_fIC_fits)
