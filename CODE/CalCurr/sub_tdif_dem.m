%%% Fraction of time spent in pelagic (for demersal)
function tdif = sub_tdif_dem(Z,param,bio1,bio2,bio3,bio4)
    % bio1, bio2: pelagic prey
    % bio3, bio4: demersal prey

    % use preferences in calculation
    biop = param.LD_phi_MF * bio1 + param.LD_phi_MP * bio2;
    biod = param.LD_phi_MD * bio3 + param.LD_phi_BE * bio4;
    
    tdif = zeros(size(Z));
    id = (Z < param.PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
    
end
