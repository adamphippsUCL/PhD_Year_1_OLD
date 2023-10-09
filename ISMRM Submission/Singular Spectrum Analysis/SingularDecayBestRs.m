% MATLAB script to find best range of Rs which have slowest decaying
% singular values

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
% V4 = [20,30,1200];
% V5 = [30,40,750];
% V1 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
    V0; V5;...
    V0; V4;...
    V0; V3;...
    V0; V2;...
    V0; V1
    ];

% Number of compartments
ncompart = 1;

% Edit Vs
Vs = Vs(1:6+2*ncompart,:);

% Make scheme
Vscheme = MakeScheme(Vs);


%% Grid search over Rs

Rmin = 0.01;
Rmax = 40.1;

rmins = linspace(Rmin, Rmax, 100);
rmaxs = linspace(Rmin, Rmax, 100);

[Rmins, Rmaxs ] = ndgrid(rmins, rmaxs);

smins = zeros(size(Rmins));
CNs = ones(size(Rmins));

% Loop over grid
for minindx = 1:length(rmins)
    for maxindx = 1:length(rmaxs)

        rmin = Rmins(minindx, maxindx);
        rmax = Rmaxs(minindx, maxindx);

        % Check allowed
        if rmax > rmin
            continue

        else
            Rs = linspace(rmin, rmax, 4);
            A = GenMatrixA(Vscheme, ncompart, Rs);

            [U,S,V] = svd(A);
            s = diag(S);
            smin = min(s);
            smins(minindx, maxindx) = smin;
            
            CN = cond(A);
            CNs(minindx, maxindx) = CN;
        end





    end

end

fig = figure;
smins(1,1) = 0;
smins(1,2) = 0.01;
pcolor(Rmins, Rmaxs, smins)
colorbar()
title('Smallest singular value')
xlabel('Rmax')
ylabel('Rmin')
saveas(fig, 'Figures/smallest_singular_value.png')

logCNs = log(CNs);
logCNs(logCNs == 0) = NaN;
logCNs(logCNs > 10) = 10;

fig = figure;
pcolor(Rmins, Rmaxs, logCNs)
xlabel('Rmax')
ylabel('Rmin')
colorbar()
title('Condition number')
saveas(fig, 'Figures/condition_number.png')