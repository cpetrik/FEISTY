%%% BIOMASS MADE FROM REPRODUCTION
function [gamma, nu, rep, egg] = sub_rep(NX,gamma,nu,K,egg)
%nu: energy for growth or spawning
%K: proportion allocated to growth
%S: fraction = fraction of pop spawning at that time
%S: 0 or 1 = indicates if spawning season (1=yes)
%egg: energy stored for later repro
S = 1.0;

% NOTE: Still never going to accumulate biomass as muscle tissue
% If it is spawning season, it gets spawned
% If it is not spawning season, it gets stored as repro energy "egg"

    if K<1.0
        rho = zeros(NX,1);
        id = (nu > 0.0);
        rho(id,1) = (1.0-K) .* nu(id,1);  %energy available for from eating
        
        %rep = eggs from energy now + eggs from stored energy
        %those with energy spawn
        rep(id,1) = S .* (rho(id,1)+egg(id,1));         %fraction of pop reproducing now
        egg(id,1) = (1.0-S) .* (rho(id,1)+egg(id,1));   %rest gets stored for later
        
        %those without energy do not repro
        rep(~id,1) = zeros;
        egg(~id,1) = zeros;
        
        %add what would be growth to next size up as repro
        rep = rep + gamma;
        gamma = zeros(NX,1);
        
        %nu is now split into used for repro (nu) and stored (egg)
        nu(id,1) = rep(id,1);
    else
        rep = zeros(NX,1);
        egg = zeros(NX,1);
    end


end
