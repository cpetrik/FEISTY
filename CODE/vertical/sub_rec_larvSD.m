%%% Biomass recruiting to size-class (g m-3 d-1)
function rec = sub_rec_larvSD(X,bio,RE,ENVR,param)
    %X: repro rate
    %bio: adult biomass
    %RE=rfrac: repro efficiency ~ sex ratio * egg survival
    %Zm: medium zoop biomass

    Zm = ENVR.Zm;

    %Ld_Repro_zi = FEISTY%Ld_Repro(i, j, nk) / (nk * dzt(k))
    
    if ENVR.H <= 200
        rec = RE .* X .* bio;
    else
        %total reproductive output
        trec = RE .* X .* bio;
        
        % account for grid cell thickness (m) - Remy did not do this
        %trec = sum((trec.*param.dz),'omitnan');
        
        trec = sum((trec),'omitnan');

        %distribure evenly (how Remy did it)
        nz = length(X);
        zfrac = ones(size(X)) ./ (nz .* param.dz);

        %distribute in water column in proportion to MZ bio
        % Zm(Zm<1e-5) = 0.0;
        % zfrac = Zm ./ sum(Zm,'omitnan');

        rec = trec .* zfrac;
    end

end