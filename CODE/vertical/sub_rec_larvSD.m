%%% Biomass recruiting to size-class (g m-3 d-1)
function rec = sub_rec_larvSD(X,bio,RE,Zm)
    %X: repro rate
    %bio: adult biomass
    %RE=rfrac: repro efficiency ~ sex ratio * egg survival
    %Zm: medium zoop biomass
    
    %total reproductive output
    trec = RE .* X .* bio;
    trec = sum(trec,'omitnan');
    %distribute in water column in proportion to MZ bio
    Zm(Zm<1e-5) = 0.0;
    zfrac = Zm ./ sum(Zm,'omitnan');
    
    rec = trec .* zfrac;

end