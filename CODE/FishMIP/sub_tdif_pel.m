%%% Fraction of time spent in pelagic (for piscivore)
function tdif = sub_tdif_pel(Z,bio1,bio2,biod)
    % bio1, bio2: pelagic prey
    % biod: demersal prey
    biop = bio1+bio2;
    tdif = ones(size(Z));
    id = (Z < PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
end
