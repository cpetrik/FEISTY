%%% Offline coupling
function [fish, ENVR] = sub_off_zm_fishvec(fish,ENVR)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    dZ = ENVR.dZm;
    
    %total biomass consumption of MZ by each fn type (g/m2)
    con = squeeze(fish.con(:,:,1)) .* fish.bio;
    %total biomass consumption of MZ by all (g/m2)
    tcon = sum(con,2);
    %original biomass-specific consumption rate of MZ of each type (g/g/m2)
    out = fish.con;
    %where consumption exceeds zoop mort
    id = (tcon > dZ);
    %frac of total MZ consumption by each type
    frac = ones(size(con));
    frac(id,:) = con(id,:) ./ tcon(id);
    %reduced biomass-specific consumption rate of MZ of each type (g/g/m2)
    %         = frac con * lost to HP (g/m2) / fish biomass (g)
    out(id,:,1) = (frac(id,:) .* dZ(id)) ./ fish.bio(id,:);
    fish.con = out;
    %fraction of HPloss consumed without adjusting
    ENVR.fZm = tcon ./ (dZ+eps);
    ENVR.fZm(id) = ones;
end
