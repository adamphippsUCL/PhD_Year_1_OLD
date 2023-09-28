% MATLAB script to carry out grid search over range of possible delta,
% Delta, and b. Find maximum dSdR over this grid


% Define diffusion coefficient
d = 2; % 10^-3 mm/s

% Define ADC
ADC = 0.0005;

% Define range of radii
Rmin = 5;
Rmax = 15;
dR = 0.1;
Rs = Rmin:dR:Rmax;


% Radius to evaluate dSdR at
Rcell = 10;


% == Configure parameter grid

% b values
bMin = 500;
bMax = 2000;
db = 50;
bs = bMin:db:bMax;

% Deltas
DeltaMin = 10;
DeltaMax = 60;
dD = 2;
Deltas = DeltaMin:dD:DeltaMax;

% deltas
deltaMin = 5;
deltaMax = 50;
dd = 2;

deltas = deltaMin:dd:deltaMax;

% Define grid of paramter values
[Ds, ds, Bs] = ndgrid(Deltas, deltas, bs);

% Empty array to be filled in
Signal_array = zeros(size(Ds));
dSdR_array = zeros(size(Ds));

% Loop over 
for Dindx = 1:length(Deltas)
    Delta = Deltas(Dindx);

    for dindx  = 1:length(deltas)
        delta = deltas(dindx);

        % Check allowed delta
        if delta > Delta
           delta = 0;
        end

        for bindx = 1:length(bs)
            b = bs(bindx);


            % Calculate gradient strength
            G = stejskal(delta, Delta, 'bval', b);
  
            % Empty signal array
            signals = zeros(length(Rs),1);

            Rindx = 0;
            for R = Rs
                Rindx = Rindx + 1;

                signalSphere = sphereGPD(delta, Delta, G, R, d) ;
                signals(Rindx) = signalSphere;
            end

            % Calculate gradient
            dSdR = gradient(signals, -dR);

            % Append to arrays
            Signal_array(Dindx, dindx, bindx) = signals( ceil( (Rcell-Rmin)/dR ) );
            Sensitivity_array(Dindx, dindx, bindx) = dSdR( ceil( (Rcell-Rmin)/dR ) )*exp(-b*ADC);
           
        end
    end
end


[M,I] = max(Sensitivity_array(:));
[I1, I2, I3] = ind2sub( size(Sensitivity_array), I);

DeltaBest = Deltas(I1)
deltaBest = deltas(I2)
bBest = bs(I3)

