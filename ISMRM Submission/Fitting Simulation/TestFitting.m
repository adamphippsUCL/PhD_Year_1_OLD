% MATLAB function to test fitting
function [rmse, bias, variance] = TestFitting(signals, Vscheme, fIC, fVASC, NoiseSigma, ncompart, fitRmin, fitRmax, fitnR, Nrep)

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
end

% Radii
fitRs = linspace(fitRmin,fitRmax,fitnR);

disp(length(signals))
%% Apply fitting multiple times

% Blank array to be filled with fitting results
fIC_fits = zeros(Nrep,1);

for indx = 1:Nrep

    % Add noise to signals
    SignalsNoisy = abs( addnoise_Adam(signals, NoiseSigma = NoiseSigma) );
    SignalsNoisy(1:2:end) = 1;
    
    Y = zeros([1,1,length(SignalsNoisy)]);
    Y(1,1,:) = SignalsNoisy;

    % Apply fitting
    [fIC_fit, fEES_fit, fVASC, R, rmse] = verdict_fit(Vscheme, Y, ncompart = ncompart, Rs = fitRs);
    fIC_fits(indx) = fIC_fit;

end


%% Analyse fitting performance
bias =  mean( (fIC_fits - fIC)  );
rmse = sqrt( mean( (fIC_fits - fIC).^2 ));
variance = var(fIC_fits);

end