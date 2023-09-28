% MATLAB function to test fitting
function [rmse, bias, variance, errors, param_sums] = TestFitting(signals, Vscheme, fIC, fVASC, NoiseSigma, ncompart, fitRmin, fitRmax, fitnR, Nrep, solver)

arguments
    signals % simulated signals
    Vscheme % VERDICT scheme
    fIC % IC volume fraction
    fVASC % VASC volume fraction
    NoiseSigma % Noise level (as fraction of b=0 signal intensity)
    ncompart % Number of compartments beyond IC (=1 for No VASC)
    fitRmin % Min radius in fitting
    fitRmax % Max radius in fitting
    fitnR % Number of radii used in fitting
    Nrep % Number of fitting repetitions
    solver % Type of solver used in optimisation
end

% Radii
fitRs = linspace(fitRmin,fitRmax,fitnR);

disp(length(signals))
%% Apply fitting multiple times

% Blank array to be filled with fitting results
fIC_fits = zeros(Nrep,1);
param_sums = zeros(Nrep,1);

for indx = 1:Nrep

    % Add noise to signals
    SignalsNoisy = abs( addnoise_Adam(signals, NoiseSigma = NoiseSigma) );
    SignalsNoisy(1:2:end) = 1;
    
    Y = zeros([1,1,length(SignalsNoisy)]);
    Y(1,1,:) = SignalsNoisy;

    % Apply fitting
    [fIC_fit, fEES_fit, fVASC_fit, R, rmse] = verdict_fit(Vscheme, Y, ncompart = ncompart, Rs = fitRs, solver = solver);
    fIC_fits(indx) = fIC_fit;
    param_sums(indx) = fIC_fit + fEES_fit + fVASC_fit;

end


%% Analyse fitting performance
% if ncompart == 1
%     fIC_this = fIC/(1-fVASC);
%     disp('yeeeet')
% else
%     fIC_this = fIC;
% end
fIC_this = fIC;
bias =  mean( (fIC_fits - fIC_this)  );
rmse = sqrt( mean( (fIC_fits - fIC_this).^2 ));
variance = var(fIC_fits);
errors = fIC_fits - fIC_this;

end