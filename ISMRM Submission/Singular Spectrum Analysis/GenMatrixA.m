function A = GenMatrixA(scheme, ncompart, Rs)

% Function to generate Matrix A for a given verdict scheme and given
% fitting parameter combination.

% Nummber of radii
nr = length(Rs);

% Number of scheme elements
nscheme = length(scheme);

% Tissue parameters
t.dEES = 2;
t.dIC = 2;
t.dVASC = 8;
t.Rs = Rs;


A = zeros([nscheme nr+ncompart]) ; 
sIC  = zeros([nscheme nr]) ;
sEES = zeros([1 nscheme]) ;
sVASC = zeros([1 nscheme]) ;

for ischeme = 1:nscheme
    for ir = 1:nr
        sIC(ischeme,ir) = sphereGPD(scheme(ischeme).delta, scheme(ischeme).DELTA, ...
           scheme(ischeme).G, t.Rs(ir), t.dIC);
    end

    sEES(ischeme)  = ball(scheme(ischeme).bval, t.dEES) ;
    sVASC(ischeme) = astrosticks(scheme(ischeme).bval, t.dVASC) ;

    if ncompart == 2
        A(ischeme,:)   = [sIC(ischeme,:) sEES(ischeme) sVASC(ischeme) ] ;
    else
        A(ischeme,:)   = [sIC(ischeme,:) sEES(ischeme)] ;
    end
end

end