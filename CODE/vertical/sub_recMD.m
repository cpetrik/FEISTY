%%% Biomass recruiting to size-class (g m-3 d-1)
function rec = sub_recMD(X,bio,td,param)
    % X = biomass specific maturation rate of smaller size class (gamma)
    % bio = biomass of smaller size class
    % td = pel/dem

    %fish(MD)%dBdt_fish(i,j,k) = sum(fish(SD)%Fout(i,j,1:nk) * FEISTY%Sd_B(i,j,1:nk) * dzt(1:nk))

    %growth out of size class at each depth
    drec = X .* bio;
    %total recruitment
    trec = sum((drec.* param.dz),'omitnan');
    %put in bottom layer - need to figure out how to do this if bottom !=75
    rec = trec .* td;

end
