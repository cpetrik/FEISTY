%%% Fraction of time spent in pelagic (for demersal)
function tdif = sub_tdif_dem_fishvec(Z,param,bio,bent)
    % pelagic prey param.ixM(1:2)
    % demersal prey param.ixM(3)

    % use preferences in calculation
    biop = param.theta(param.ixL(2),param.ixM(1)) * bio(:,param.ixM(1)) + ...
        param.theta(param.ixL(2),param.ixM(2)) * bio(:,param.ixM(2));
    biod = param.theta(param.ixL(2),param.ixM(3)) * bio(:,param.ixM(3)) + ...
        param.theta(param.ixL(2),3) * bent;
    
    tdif = zeros(size(Z));
    id = (Z < param.PI_be_cutoff);
    tdif(id,1) = biop(id,1) ./ (biop(id,1) + biod(id,1));
    
end
