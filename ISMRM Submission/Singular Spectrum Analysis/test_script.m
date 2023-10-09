% Script to test matrix A generation

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
    V0; V5;...
    V0; V4;...
    V0; V3;...
    V0; V2;...
%     V0; V1
    ];

Vscheme = MakeScheme(Vs);


%% Define fitting parameters


% Rs = linspace(0.1,15.1,17);
Rs = linspace(0.1,15.1,4);
% Rs = linspace(3,12,4);
% Rs = linspace(6,15,4);
% Rs = linspace(0.1,9.1,4);

ncompart = 1;


%% Generate matrix A

A = GenMatrixA(Vscheme, ncompart, Rs);

%% Generate singular values and condition number

[U,S,V] = svd(A);
s = diag(S);

CN = cond(A);

plot(log(s))
hold on

% %% Plot V
% x = [""];
% % Generate bar chart name
% indx = 0;
% for R = Rs
%     indx = indx + 1;
%     x(indx) = num2str(R);
% end
% x(indx+1) = "fEES";
% if ncompart == 2
%     x(indx+2) = "fVASC";
% end
% y = V(:,end)
% bar( y)