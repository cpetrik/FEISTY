%%% Fraction of time spent in pelagic (for demersal)
function tdif = sub_tdif_dem(Z,bio1,bio2,bio3,bio4)
    % bio1, bio2: pelagic prey
    % bio3, bio4: demersal prey

    global PI_be_cutoff LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE

    % use preferences in calculation
    biop = LD_phi_MF*bio1 + LD_phi_MP*bio2;
    biod = LD_phi_MD*bio3 + LD_phi_BE*bio4;
    
    tdif = zeros(size(Z));
    id = (Z < PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
    
end
