% Matlab script to test the verdict fitting

% Need to define scheme and combine VERDICT images

% Combined VERDICT image
load('Y.mat');


% Define VERDICT scheme
V1 = [3.2, 28.3, 90];
V2 = [9.7, 34.8, 500];
V3 = [22.7, 47.8, 1500];
V4 = [13.5, 38.5, 2000];
V5 = [18.9, 38.8, 3000];

Vs = [V1; V2; V3; V4; V5];

% Make VERDICT scheme
nscheme = length(Vs);
for i = 1:nscheme
    delta = Vs(i,1);
    Delta = Vs(i,2);
    bval = Vs(i,3);
    scheme(i) = fill_scheme(delta,Delta, bval);
end


% Apply VERDICT fitting
[fIC_fit_new, fEES_fit, fVASC_fit, R_fit, rmse, A, t, opt] = verdict_fit(scheme, Y);





function scheme = fill_scheme(delta, DELTA, bval)
% FILL_SCHEME Fills scheme structure and checks parameters
%
% See also stejskal

scheme.delta = delta;
scheme.DELTA = DELTA ;

scheme.bval = bval ;

scheme.G = stejskal(delta,DELTA,bval=bval);

    
end
