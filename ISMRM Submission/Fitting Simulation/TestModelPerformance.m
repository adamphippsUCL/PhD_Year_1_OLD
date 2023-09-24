% MATLAB script to test model performance

% For each number of model parameters, simulate signal and fitting for
% multiple random sets of radii

%% Define VERDICT scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];

Vs = [...
    V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    V0; V5];


%% Define tissue parameters

% Volume fractions
fIC = 0.4;
fVASC = 0.1;
fEES = 1-fIC-fVASC;

% Spheres
tissueRmin = 6;
tissueRmax = 14;
tissuenR = 100; % Number of radii

% Noise level
NoiseSigma = 0.02;


%% Define iterations 

% Number of simulation iterations for each number of model parameters
Niter = 1;


%% Define fitting parameters

% Radii range
fitRmin = 8;
fitRmax = 12;

% N compartments
ncompart = 1;

% Scheme
Vs = Vs( (5-2*ncompart):end , : );

%% Iterate over number of model parameters

% Range of model parameters
Nparams = 4:1:19;

% Empty array for results (
results = zeros(length(Nparams), Niter, 3);


for iterindx = 1:Niter

    % Simulate signal
    [signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, tissueRmin, tissueRmax, tissuenR );
    signals(1:2:end) = 1;

    Npindx = 0;
    for Nparam = Nparams
        Npindx = Npindx + 1;

        % Test fitting
        [rmse, bias, variance] = TestFitting(signals, Vscheme, fIC, fVASC, NoiseSigma, ncompart, fitRmin, fitRmax, (Nparam-ncompart), 200);

        % Append to results
        results(Npindx, iterindx, :) = [rmse, bias, variance];

    end
end

figure;
for iterindx = 1:Niter
    plot(Nparams, results(:,iterindx,1), '-*')
    hold on
end
