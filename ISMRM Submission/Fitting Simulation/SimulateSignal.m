% MATLAB script to define VERDICT scheme and simulate signal

function [signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, Rmin, Rmax, nR)

arguments
    Vs % VERDICT scheme parameters
    fIC % Intracellular volume fraction
    fVASC % Vascular volume fraction
    Rmin % Minimum radius
    Rmax % maximum radius
    nR % Number of different radii to simulate

end

%% Define tissue parameters
fEES = 1-fIC-fVASC;

% Random radii
Rs = Rmin + (Rmax - Rmin)*rand(1,nR);

% Random volume fractions
fRs = zeros(1,length(Rs));
for indx =1:length(Rs)
    fRs(indx) = (1-sum(fRs))*rand(1,1);
end
fRs = fIC*fRs;

% fR1 = rand(1,1);
% fR2 = (1-fR1)*rand(1,1);
% fR3 = 1-fR1-fR2;
% fRs = fIC*[fR1,fR2,fR3];

tissue_params = [fIC, fEES, fVASC];


%% Simulate signal over scheme

[signals, Vscheme] = verdict_simulate_Adam(Vs, tissue_params, Rs, fRs);


end


