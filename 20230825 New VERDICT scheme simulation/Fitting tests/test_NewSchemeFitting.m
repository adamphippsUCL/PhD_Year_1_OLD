% Matlab script to test fitting ability of proposed new VERDICT scheme


% Define New VERDICT scheme
V1 = [10, 20, 500];
V2 = [15, 25, 1000];
V3 = [20, 30, 1250];
V4 = [30, 45, 750];
V5 = [25, 35, 2000];

Vs = [V1; V2; V3; V4];%; V5];


% % Define old scheme
% V1 = [3.9, 23.8, 90];
% V2 = [11.4, 31.3, 500];
% V3 = [23.9, 43.8, 1500];
% V4 = [12.4, 32.3, 2000];
% V5 = [18.9, 38.8, 3000];
% 
% Vs = [V2; V3; V4; V5];




% Define parameter values
fIC = 0.3 ;
fEES = 1-fIC;
fVASC = 0;
Ravg = 10;
Rstd = 0.1;

tissue_params = [fIC, fEES, fVASC, Ravg, Rstd];


% Simulate signal
[signals, Vscheme, params] = verdict_simulate_Adam(Vs, tissue_params);


% Add noise to signal
signalsNoisy = abs(addnoise_Adam(signals, 20));
% signalsNoisy = abs(addnoise(signals,'snr',20,'signal',1) );

scatter([Vscheme.bval], signalsNoisy)
ylim([0,1])
hold on

Y = zeros([1,1,length(signalsNoisy)]);
Y(1,1,:) = signalsNoisy;

[fIC_fit, fEES_fit, fVASC_fit, R_fit, rmse, A, t, opt] = verdict_fit(Vscheme, Y, ncompart = 1);



