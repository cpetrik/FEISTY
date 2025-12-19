%%% Biomass recruiting to size-class (g m-3 d-1)
function rec = sub_recMD(X,bio,td)
    % X = biomass specific maturation rate of smaller size class (gamma)
    % bio = biomass of smaller size class
    % td = pel/dem

    %growth out of size class at each depth
    drec = X .* bio;
    %total recruitment
    trec = sum(drec,'omitnan');
    %put in bottom layer - need to figure out how to do this if bottom !=75
    rec = trec .* td;

end
