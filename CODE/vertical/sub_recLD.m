%%% Biomass recruiting to size-class (g m-3 d-1)
function rec = sub_recLD(X,bio,ENVR)
    %X: growth rate out of Med size
    %bio: Med biomass
    
    %total growth out of Med size
    trec = X .* bio;

    %if shallow shelf, distribute in water column 
    if ENVR.H <= 200
        nz = length(X);
        frac = ones(size(X)) .* (1.0/nz);
        trec = sum(trec,'omitnan');
        rec = trec .* frac;
    else
        rec = trec;
    end

end