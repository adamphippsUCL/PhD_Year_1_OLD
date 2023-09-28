% MATLAB script to compare results of lsqnonneg and lslqlin fitting

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
tissueRmin = 0.5;
tissueRmax = 15;
tissuenR = 100; % Number of radii

% Noise level
NoiseSigma = 0.05;


%% Iterations
Niter = 50;

% Results for original fitting method
Results0 = zeros(Niter, 3);
% Results for constrained fitting
Results1 = zeros(Niter, 3);



for iterindx = 1:Niter


    [signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, tissueRmin, tissueRmax, tissuenR );

    % Needed if b=0 scans included in scheme!!!
    signals(1:2:end) = 1;
    
    % lsqnonneg fitting
    ncompart = 1;
    [rmse0, bias0, variance0, errors0, param_sums0] = TestFitting( ...
        signals( (5-2*ncompart):end ), ...
        Vscheme( (5-2*ncompart):end ), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        15.1, ...
        17, ...
        200, ...
        'lsqnonnegTikonhov' ...
        );

    Results0(iterindx,:) = [rmse0, bias0, variance0];
    

    % lsqlin fitting
    ncompart = 1;
    [rmse1, bias1, variance1, errors1, param_sums1] = TestFitting( ...
        signals( (5-2*ncompart):end ), ...
        Vscheme( (5-2*ncompart):end ), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        15.1, ...
        17, ...
        200, ...
        'lsqlin_underdetermined' ...
        );

    % Append results
    Results1(iterindx,:) = [rmse1, bias1, variance1];


end

fig = figure;
scatter( zeros(Niter,1), Results0(:,1), '*')
hold on
scatter( ones(Niter,1), Results1(:,1), '*')
xlim([-0.5,1.5])
ylabel('rmse')
xticks([0 1])
xticklabels({'lsqnoneg','lsqlin constrained'})
saveas(fig, 'Figures/Constrained Fitting/rmse_results.png')

fig = figure;
scatter( zeros(Niter,1), Results0(:,2), '*')
hold on
scatter( ones(Niter,1), Results1(:,2), '*')
xlim([-0.5,1.5])
ylabel('bias')
xticks([0 1])
xticklabels({'lsqnoneg','lsqlin constrained'})
saveas(fig, 'Figures/Constrained Fitting/bias_results.png')

fig = figure;
scatter( zeros(Niter,1), Results0(:,3), '*')
hold on
scatter( ones(Niter,1), Results1(:,3), '*')
xlim([-0.5,1.5])
ylabel('variance')
xticks([0 1])
xticklabels({'lsqnoneg','lsqlin constrained'})
saveas(fig, 'Figures/Constrained Fitting/variance_results.png')

fig = figure;
plot(param_sums0)
hold on
plot(param_sums1)
saveas(fig, 'Figures/Constrained Fitting/parameter_sums.png')

