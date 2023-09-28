% Trying to recreate results I saw yesterday

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

% % Without b=0
% Vs = [...
%     V1;
%     V2;
%     V3;
%     V4;
%     V5
%     ];


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

% Define number of iterations
Niter = 20;

Results0 = zeros(Niter, 3);
Results1 = zeros(Niter, 3);
Results2 = zeros(Niter, 3);

for iterindx = 1:Niter


    [signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, tissueRmin, tissueRmax, tissuenR );

    % Needed if b=0 scans included in scheme!!!
    signals(1:2:end) = 1;
    
    % Original VERDICT
    ncompart = 2;
    [rmse0, bias0, variance0] = TestFitting( ...
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

    Results0(iterindx,:) = [rmse0, bias0, variance0];
    
    % Original VERDICT Rs
    ncompart = 1;
    [rmse1, bias1, variance1] = TestFitting( ...
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
%     
    % Reduced Rs
    ncompart = 1;
    [rmse2, bias2, variance2] = TestFitting( ...
        signals( (5-2*ncompart):end), ...
        Vscheme( (5-2*ncompart):end), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        6, ...
        12, ...
        4,  ...
        200,...
        'lsqlin_well_determined');

    % Append results
    Results2(iterindx,:) = [rmse2, bias2, variance2];

end

fig = figure;
scatter( zeros(Niter,1), Results0(:,1), '*')
hold on
scatter( ones(Niter,1), Results1(:,1), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,1), '*')
xlim([-0.5,2.5])
ylabel('rmse')
xticks([0 1 2])
xticklabels({'Original','No VASC','No VASC reduced Rs'})
saveas(fig, 'Figures/rmse_results.png')

fig = figure;
scatter( zeros(Niter,1), Results0(:,2), '*')
hold on
scatter( ones(Niter,1), Results1(:,2), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,2), '*')
xlim([-0.5,2.5])
ylabel('bias')
xticks([0 1 2])
xticklabels({'Original','No VASC','No VASC reduced Rs'})
saveas(fig, 'Figures/bias_results.png')

fig = figure;
scatter( zeros(Niter,1), Results0(:,3), '*')
hold on
scatter( ones(Niter,1), Results1(:,3), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,3), '*')
xlim([-0.5,2.5])
ylabel('Variance')
xticks([0 1 2])
xticklabels({'Original','No VASC','No VASC reduced Rs'})
saveas(fig, 'Figures/variance_results.png')