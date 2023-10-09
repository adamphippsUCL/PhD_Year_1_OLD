function Vscheme = MakeScheme(Vs)


    % VERDICT scheme
    nscheme = length(Vs);
    for i = 1:nscheme
        delta = Vs(i,1);
        Delta = Vs(i,2);
        bval = Vs(i,3);
        Vscheme(i) = fill_scheme(delta, Delta, bval);
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


end