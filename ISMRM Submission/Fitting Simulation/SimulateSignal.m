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

% % % Linspace radii
% Rs = linspace(Rmin, Rmax, nR);

% Random volume fractions
% fRs = zeros(1,length(Rs));
% for indx =1:length(Rs)
%     fRs(indx) = (1-sum(fRs))*rand(1,1);
% end

fRs = rand(1,length(Rs));
% fRs = ones(1, length(Rs));
 
% % Normal distributions adaptation!
% fRs = normpdf(Rs,8,2);

fRs = fIC*fRs/sum(fRs);


tissue_params = [fIC, fEES, fVASC];


%% Simulate signal over scheme

[signals, Vscheme] = verdict_simulate_Adam(Vs, tissue_params, Rs, fRs);


end


