%%% Offline coupling
function [fish, ENVR] = sub_off_zl_fishvec(fish,ENVR)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    
    % HPloss
    dZ = ENVR.dZl;
    %total biomass consumption of LZ by each fn type (g/m2)
    con = squeeze(fish.con(:,:,2)) .* fish.bio;
    %total biomass consumption of LZ by all (g/m2)
    tcon = sum(con,2);
    %original biomass-specific consumption rate of MZ of each type (g/g/m2)
    out = fish.con;
    %where consumption exceeds zoop mort
    id = (tcon > dZ);
    %frac of total LZ consumption by each type
    frac = ones(size(con));
    frac(id,:) = con(id,:) ./ tcon(id);
    %reduced biomass-specific consumption rate of LZ of each type (g/g/m2)
    %         = frac con * lost to HP (g/m2) / fish biomass (g)
    out(id,:,2) = (frac(id,:) .* dZ(id)) ./ fish.bio(id,:);
    fish.con = out;
    %fraction of HPloss consumed without adjusting
    ENVR.fZl = tcon ./ (dZ+eps);
    ENVR.fZl(id) = ones;
end
