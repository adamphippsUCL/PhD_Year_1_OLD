% MATLAB script to simulate VERDICT fitting for original and reduced models

%% Define VERDICT scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];


% % New scheme?
% V0 = [1,2,0];
% V2 = [10,20,500];
% V3 = [15,25,1000];
% V4 = [20,30,1500];
% V5 = [30,45,750];
% V1 = [18.9, 38.8, 3000];

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
fIC = 0.5;
fVASC = 0.15;
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

% Original model
Results0 = zeros(Niter, 3);

% No VASC, original Rs
Results1 = zeros(Niter, 3);

% No VASC, Rs = [0.1,5.1,10.1,15.1]
Results2 = zeros(Niter, 3);

% No VASC, Rs = [3, 6, 9, 12]
Results3 = zeros(Niter, 3);

% No VASC, Rs = [6,9,12,15]
Results4 = zeros(Niter, 3);

% No VASC, Rs = [0.1,3.1,6.1,9.1]
Results5 = zeros(Niter, 3);

% No VASC, No 3000, Rs = [4,8,12]
Results6 = zeros(Niter, 3);

% No VASC, No 3000, Rs = [4,8,12]
Results7 = zeros(Niter, 3);

% No VASC, No 3000, Rs = [4,8,12]
Results8 = zeros(Niter, 3);

for iterindx = 1:Niter

    iterindx

    [signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, tissueRmin, tissueRmax, tissuenR );
 
    
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
        'lsqnonnegTikonhov' ...
        );

    Results0(iterindx,:) = [rmse0, bias0, variance0];


    % No VASC, VERDICT Rs
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
        'lsqnonnegTikonhov' ...
        );

   % Append results
    Results1(iterindx,:) = [rmse1, bias1, variance1];    


    % No VASC, Rs = [0.1,5.1,10.1,15.1]
    ncompart = 1;
    [rmse2, bias2, variance2] = TestFitting( ...
        signals( (5-2*ncompart):end), ...
        Vscheme( (5-2*ncompart):end), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        15.1, ...
        4,  ...
        200,...
        'lsqnonneg');

   % Append results
    Results2(iterindx,:) = [rmse2, bias2, variance2];  


    % No VASC, Rs = [3,6,9,12]
    ncompart = 1;
    [rmse3, bias3, variance3] = TestFitting( ...
        signals( (5-2*ncompart):end), ...
        Vscheme( (5-2*ncompart):end), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        3, ...
        12, ...
        4,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results3(iterindx,:) = [rmse3, bias3, variance3];


     % No VASC, Rs = [6,9,12,15]
    ncompart = 1;
    [rmse4, bias4, variance4] = TestFitting( ...
        signals( (5-2*ncompart):end), ...
        Vscheme( (5-2*ncompart):end), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        6, ...
        15, ...
        4,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results4(iterindx,:) = [rmse4, bias4, variance4];   


    % No VASC, Rs = [0.1,3.1,6.1,9.1]
    ncompart = 1;
    [rmse5, bias5, variance5] = TestFitting( ...
        signals( (5-2*ncompart):end), ...
        Vscheme( (5-2*ncompart):end), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        9.1, ...
        4,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results5(iterindx,:) = [rmse5, bias5, variance5];   


    % No VASC, No 3000, Rs = [4,8,12], 
    ncompart = 1;
    [rmse6, bias6, variance6] = TestFitting( ...
        signals( (5-2*ncompart):end-2), ...
        Vscheme( (5-2*ncompart):end-2), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        4, ...
        12, ...
        3,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results6(iterindx,:) = [rmse6, bias6, variance6]; 



    % No VASC, No 3000, Rs = [0.1, 4.1 ,8.1], 
    ncompart = 1;
    [rmse7, bias7, variance7] = TestFitting( ...
        signals( (5-2*ncompart):end-2), ...
        Vscheme( (5-2*ncompart):end-2), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        8, ...
        3,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results7(iterindx,:) = [rmse7, bias7, variance7]; 
    

    % No VASC, No 3000, Rs = [0.1, 4.1 ,8.1], 
    ncompart = 1;
    [rmse8, bias8, variance8] = TestFitting( ...
        signals( (5-2*ncompart):end-2), ...
        Vscheme( (5-2*ncompart):end-2), ...
        fIC, ...
        fVASC, ...
        NoiseSigma, ...
        ncompart, ...
        0.1, ...
        8, ...
        3,  ...
        200,...
        'lsqnonneg');

    % Append results
    Results8(iterindx,:) = [rmse8, bias8, variance8]; 
end


fig = figure;
scatter( zeros(Niter,1), Results0(:,1), '*')
hold on
scatter( ones(Niter,1), Results1(:,1), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,1), '*')
hold on
scatter( 3*ones(Niter,1), Results3(:,1), '*')
hold on
scatter( 4*ones(Niter,1), Results4(:,1), '*')
hold on
scatter( 5*ones(Niter,1), Results5(:,1), '*')
hold on
scatter( 6*ones(Niter,1), Results6(:,1), '*')
xlim([-0.5,6.5])
ylabel('rmse')
xticks([0 1 2 3 4 5 6])
xticklabels({ ...
    'Original', ...
    'No VASC', ...
    'No VASC [0.1,5.1,10.1,15.1]', ...
    'No VASC, Rs = [3, 6, 9, 12]',...
    'No VASC, Rs = [6,9,12,15]',...
    'No VASC, Rs = [0.1,3.1,6.1,9.1]',...
    'No VASC, Rs = [4,8,12], No 3000',...

    })
grid on
saveas(fig, 'Figures/rmse_results.png')


fig = figure;
scatter( zeros(Niter,1), Results0(:,2), '*')
hold on
scatter( ones(Niter,1), Results1(:,2), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,2), '*')
hold on
scatter( 3*ones(Niter,1), Results3(:,2), '*')
hold on
scatter( 4*ones(Niter,1), Results4(:,2), '*')
hold on
scatter( 5*ones(Niter,1), Results5(:,2), '*')
hold on
scatter( 6*ones(Niter,1), Results6(:,2), '*')
xlim([-0.5,6.5])
ylabel('bias')
xticks([0 1 2 3 4 5 6])
xticklabels({ ...
    'Original', ...
    'No VASC', ...
    'No VASC [0.1,5.1,10.1,15.1]', ...
    'No VASC, Rs = [3, 6, 9, 12]',...
    'No VASC, Rs = [6,9,12,15]',...
    'No VASC, Rs = [0.1,3.1,6.1,9.1]',...
    'No VASC, Rs = [4,8,12], No 3000',...
    })
grid on
saveas(fig, 'Figures/bias_results.png')

fig = figure;
scatter( zeros(Niter,1), Results0(:,3), '*')
hold on
scatter( ones(Niter,1), Results1(:,3), '*')
hold on
scatter( 2*ones(Niter,1), Results2(:,3), '*')
hold on
scatter( 3*ones(Niter,1), Results3(:,3), '*')
hold on
scatter( 4*ones(Niter,1), Results4(:,3), '*')
hold on
scatter( 5*ones(Niter,1), Results5(:,3), '*')
hold on
scatter( 6*ones(Niter,1), Results6(:,3), '*')
xlim([-0.5,6.5])
ylabel('variance')
xticks([0 1 2 3 4 5 6])
xticklabels({ ...
    'Original', ...
    'No VASC', ...
    'No VASC [0.1,5.1,10.1,15.1]', ...
    'No VASC, Rs = [3, 6, 9, 12]',...
    'No VASC, Rs = [6,9,12,15]',...
    'No VASC, Rs = [0.1,3.1,6.1,9.1]',...
    'No VASC, Rs = [4,8,12], No 3000',...
    })
grid on
saveas(fig, 'Figures/variance_results.png')