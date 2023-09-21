% MATLAB function to simulate VERDICT signal given fIC, fVASC, fEES, Ravg,
% Rstd

function [signals, Vscheme, params] = verdict_simulate_Adam(scheme, tissue_params, Rs, fRs)

% This function outputs the simulated signals for each scan in a specified
% VERDICT scheme.

% === Inputs:

% scheme = specification of VERDICT scheme parameters

% e.g. scheme = [V1, V2, V3, V4, V5] where Vi = [delta, Delta, G]

% tissue_params = specified tissue parameters

% e.g. tissue_params = [fIC, fEES, fVASC, Ravg, Rstd]


% === Outputs

% signals = signal magnitudes for each VERDICT scan

% -------------------------------------------------


% Deconstruct inputs

% VERDICT scheme
nscheme = length(scheme);
for i = 1:nscheme
    delta = scheme(i,1);
    Delta = scheme(i,2);
    bval = scheme(i,3);
    Vscheme(i) = fill_scheme(delta, Delta, bval);
end

% Tissue parameters
params.fIC = tissue_params(1);
params.fEES = tissue_params(2);
params.fVASC = tissue_params(3);
params.Rs = Rs;
params.fRs = fRs;
params.dEES = 2 ;
params.dIC = 2 ;
params.dVASC = 8 ;


% Radii used in fitting
nR = length(params.Rs);

% == First calculate volume fractions for different radii

% % Radii used in fitting (in micrometers)
% Rmin = 0.1;
% Rmax = 20;
% dR = 0.1;
% Rs = Rmin:dR:Rmax; 
% nR = length(Rs);
% 
% % Caclualte volume fractions (normal distribution around Ravg with std
% % Rstd, sum of fRs is fIC)
% fRs = gaussmf(Rs, [params.Rstd, params.Ravg]);
% fRs_sum = sum(fRs);
% fRs = (params.fIC/fRs_sum)*fRs;



% === Simulating signals

% Empty arrays for signals
sIC = zeros([nscheme nR]) ;
sEES = zeros([nscheme 1]) ;
sVASC = zeros([nscheme 1]) ;
stot = zeros([nscheme 1]) ;


for ischeme = 1:nscheme

    % SphereGPD signal for different radii
    for ir = 1:nR
        sIC(ischeme,ir) = sphereGPD(Vscheme(ischeme).delta, Vscheme(ischeme).DELTA, ...
            Vscheme(ischeme).G, params.Rs(ir), params.dIC);
    end

    % EES signal
    sEES(ischeme) = ball(Vscheme(ischeme).bval, params.dEES) ;
    
    % VASC signal
    sVASC(ischeme) = astrosticks(Vscheme(ischeme).bval, params.dVASC) ;

    % Total signal
    stot(ischeme) = sum(params.fRs.*sIC(ischeme,:)) + params.fEES*sEES(ischeme) + ...
        params.fVASC*sVASC(ischeme);

    
end

% scatter([Vscheme.bval], stot,'k')
% xlabel('bval')
% hold on

signals = stot;

end


function scheme = fill_scheme(delta, DELTA, bval)
% FILL_SCHEME Fills scheme structure and checks parameters
%
% See also stejskal

scheme.delta = delta;
scheme.DELTA = DELTA ;

scheme.bval = bval ;

scheme.G = stejskal(delta,DELTA,bval=bval);

    
end


