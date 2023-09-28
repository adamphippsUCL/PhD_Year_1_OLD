% MATLAB to investigate how fitted volume fraction error depends on 
% radius and noise level

% Outline:

% Randomise fRs
% Simulate signal
% Apply verdict fitting
% Investigate how fitted fRs differ from defined fRs
% Investigate how this depends on noise level


%% Define VERDICT scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
    V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    V0; V5
    ];

%% Define tissue parameters

% Volume fractions
fIC = 0.4;
fVASC = 0.1;
fEES = 1-fIC-fVASC;

% Spheres
tissueRmin = 0.1;
tissueRmax = 15;
tissuenR = 100; % Number of radii

% Noise level
NoiseSigma = 0.05;




%% Apply VERDICT fitting

% Define fitting
ncompart = 1;
fitRs = linspace(6,12,4);

for indx = 1:8


    %% Simulate signal
    
    [signals, Vscheme, Rs, fRs] = SimulateSignal(Vs, fIC, fVASC, tissueRmin, tissueRmax, tissuenR );
    
    % Add noise
    SignalsNoisy = abs( addnoise_Adam(signals, NoiseSigma = NoiseSigma) );
    SignalsNoisy(1:2:end) = 1;

    
    SignalsNoisy = SignalsNoisy((5-2*ncompart):end);
    Vscheme = Vscheme((5-2*ncompart):end);
    
    Y = zeros([1,1,length(SignalsNoisy)]);
    Y(1,1,:) = SignalsNoisy;
    
    [fIC_fit, fEES_fit, fVASC_fit, R, rmse, A, t, opt, x] = verdict_fit(Vscheme, Y, ncompart = ncompart, Rs = fitRs, solver = 'lsqlin_well_determined');
    
    fitfRs = x(1:length(fitRs));
    
    fig = figure;
    plot(Rs, fRs, '*')
    hold on
    plot(fitRs, fitfRs, '-p')
    legend(["Simulated fRs", "Fitted fRs"])
    saveas(fig, ['Figures/fR fitting/' num2str(indx) '.png'])

end