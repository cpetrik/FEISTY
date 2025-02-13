%%% Fraction of time spent in pelagic (for demersal)
function tdif = sub_tdif_dem(Z,param,bio1,bio2,bio3,bio4)
    % bio1: medium/adult forage prey biomass
    % bio2: medium/juvenile large pelagic prey biomass
    % bio3: medium/juvenile demersal prey biomass
    % bio4: benthic invertebrate prey biomass
    
    % use preferences (e.g. param.LD_phi_MF) in calculation
    biop = param.LD_phi_MF * bio1 + param.LD_phi_MP * bio2;
    biod = param.LD_phi_MD * bio3 + param.LD_phi_BE * bio4;
    
    tdif = zeros(size(Z));
    id = (Z < param.PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
    
end
