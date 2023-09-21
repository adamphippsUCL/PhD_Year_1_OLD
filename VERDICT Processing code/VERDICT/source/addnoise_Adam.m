% Function to add complex noise to signal

function signalsNoisy = addnoise_Adam(signals, opts)

% ---------------------------

% INPUTS

% signals = vector of signals from a voxel (over VERDICT scheme)

% SNRs = float/vector of SNR values for each VERDICT scan

% OUTPUTS

% signalsNoisy = vector of signals with added noise (complex noise)

% ---------------------------

arguments
    signals % vector of signals
    opts.SNRs = [] % float/vector of SNR values for each VERDICT scan
    opts.NoiseSigma = [] % Standard deviation of noise

end

signalsNoisy = zeros(size(signals));


% If single SNR defined
if length(opts.SNRs) == 1

    % Add noise to each signal
    for indx = 1:length(signals)
        
        % Calculate noise standard deviation
        NoiseSigma = signals(indx)/opts.SNRs;

        % Generate complex noise
        NoiseComplex = NoiseSigma.*randn + 1i*NoiseSigma.*randn;

        % Construct noisy signal
        signalsNoisy(indx) = signals(indx) + NoiseComplex;

    end


% If SNR defined for each signal
elseif length(opts.SNRs) == length(signals)
    
    % Add noise to each signal
    for indx = 1:length(signals)
        
        % Calculate noise standard deviation
        NoiseSigmas = signals(indx)/opts.SNRs(indx);

        % Generate complex noise
        NoiseComplex = NoiseSigmas.*randn + 1i*NoiseSigmas.*randn;

        % Construct noisy signal
        signalsNoisy(indx) = signals(indx) + NoiseComplex;

    end


elseif opts.NoiseSigma

    for indx = 1:length(signals)

        % Generate complex noise
        NoiseComplex = opts.NoiseSigma.*randn + 1i*opts.NoiseSigma.*randn;

        % Construct noisy signal
        signalsNoisy(indx) = signals(indx) + NoiseComplex;

    end

else
    signalsNoisy = signals;
end


end