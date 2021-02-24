%%% Fraction of time spent in pelagic (for piscivore)
function tdif = sub_tdif_pel_fishvec(Z,param,bio)
    % pelagic prey param.ixM(1:2)
    % demersal prey param.ixM(3)

    % use preferences in calculation
    biop = param.theta(param.ixL(2),param.ixM(1)) * bio(:,param.ixM(1)) + ...
        param.theta(param.ixL(2),param.ixM(2)) * bio(:,param.ixM(2));
    biod = param.theta(param.ixL(2),param.ixM(3)) * bio(:,param.ixM(3));
    
    tdif = ones(size(Z));
    id = (Z < param.PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
end
